// Nurbs_API.cpp : Defines the basic NURBS API
//
#include "pch.h"
#include "NurbsLibSt.h"
#include <tchar.h>

#ifdef USE_IMSL
#include "IMSL_LINK.h" // IMSL
#endif

// Compute curve parametric values through averaging
void CompCurveParamsAverag(int p, int m, const CDoubleVector& _uj, CDoubleVector& _u)
{
	int n = m-p-1; // m = n + p + 1;

	ASSERT(p >= 1);
	ASSERT(m > p);
	_u.ReSize(m+1);
	double *u = _u.GetData();
	const double *uj = _uj.GetData();
	double div_p = 1.0 / p;
	int i,j;
	for(j = 0; j <= p; j++)	u[j] = 0.0;
	for(j = m-p; j <= m; j++)	u[j] = 1.0;

	for(j = 1; j <= n-p; j++)
	{
		double total = 0;
		for(i=j; i<=j+p-1; i++)
			total += uj[i];
		u[j+p] = div_p * total;		
	}
}


BOOL ComputeKnotVectorChordLength(int n, int p, const CPointEx3DVector& P, double* U)
{
	ASSERT(U != NULL);
	int nPoints = n+1;

	// find the total chord length
	double sum = 0.0;
	for(int i = 0; i < nPoints-1; i++)
	{
		sum += distance(P[i], P[i+1]);
	}

	if(sum == 0)
	{
		ASSERT(FALSE);
		return FALSE;
	}

	int mu = n + p + 1;

	// Compute an equally-spaced knot vectors in the intreval [0...1]
	
	double curSum = 0.0;

	for(int i = 0; i <= mu; i++)
	{
		if(i >= 0 && i <= p)
			U[i] = 0.0;
		else
		if(i > p && i <= (mu - p - 1))
		{
			curSum += distance(P[i-1], P[i]);
			U[i] = curSum / sum;
		}
		else
		if(i > (mu - p - 1) )
			U[i] =  1.0;
	}

	return TRUE;
}

BOOL ComputeChordLengthParams(int n, const CPointEx3DVector& P, 
							  CDoubleVector& U)
{
	U.ReSize(n+1);

	// Compute an chord-length knots in the intreval [0...1]

	U[0] = 0;
	double sum = 0;
	for(int i = 1; i <= n; i++)
	{
		sum += distance(P[i-1], P[i]);
		U[i] = sum;
	}

	if(sum == 0)
	{
		ASSERT(FALSE);
		return FALSE;
	}

	for(int i = 1; i <= n; i++)
	{
		U[i] /= sum;
		ASSERT(U[i]<=1);
	}
	U[n] = 1; // make sure last knot = 1 (avoid round off errors)

	return TRUE;
}

BOOL ComputeChordLengthParams(int n, const CPointExVector& P, 
							  CDoubleVector& U)
{
	U.ReSize(n+1);

	// Compute an chord-length knots in the intreval [0...1]

	U[0] = 0;
	double sum = 0;
	for(int i = 1; i <= n; i++)
	{
		double sumi = distance(P[i-1], P[i]);
		ASSERT(sumi>0);
		sum += sumi;
		U[i] = sum;
	}

	if(sum == 0)
	{
		ASSERT(FALSE);
		return FALSE;
	}

	for(int i = 1; i <= n; i++)
	{
		U[i] /= sum;
		ASSERT(U[i]<=1);
		ASSERT(U[i-1]<U[i]);
	}
	U[n] = 1; // make sure last knot = 1 (avoid round off errors)


	return TRUE;
}

int FindSpan(int n, int p, double u, const double* U, int gspan)
{
	ASSERT(U != NULL);
	int m = n + p + 1;
	if(fabs(u - U[0]) < 1e-9)
		u = U[0];
	else
	if(fabs(u - U[m]) < 1e-9)
		u = U[m];
	ASSERT( u >= U[0] );
//	ASSERT( u >= U[p] );
	ASSERT( u <= U[m] ); // check u e [u0,...,um]
//	ASSERT( u <= U[n+1] ); // check u e [u0,...,um]

	if( u>=U[n+1])
		return (n);
	if( u <= U[p+1])
		return (p);

	int low = p;
	int high = n+1;
	int mid = (gspan>=low && gspan<=high) ? gspan : (low+high) / 2;
//	int mid = (low+high) / 2;
	while( (u < U[mid] || u >= U[mid+1]))
	{
		if( u < U[mid] ) high = mid;
		else low = mid;
		mid = (low+high) / 2;
		ASSERT(low<high);
	}

	return mid;
}

void BasisFuns(int i, double u, int p, const double* U, double* N)
{
	ASSERT(U != NULL);
	ASSERT(N != NULL);

	if(p <= 3)
	{

		// This implementation pre-supposes that degree p <= 3
		// If p > 3, then left and right must be allocated dynamically
		// and not statically as follows.
		// Since, the degree of the curve is always 3, we choosed 
		// the second type of implementation to increase speed performance.
		ASSERT(p <= 3 );

		double left[4];
		double right[4];

		double saved, temp, denom;
	
		left[0] = right[0] = 0.0;

		N[0] = 1.0;
		for(int j = 1; j <= p; j++)
		{
			left[j] = u - U[i+1-j];
			right[j] = U[i+j] - u;
			saved = 0.0;
			for(int r=0; r<j; r++)
			{
				denom = (right[r+1] + left[j-r]);
				temp = N[r] / denom;
				N[r] = saved + right[r+1] * temp;
				saved = left[j-r] * temp;
			}
			N[j] = saved;
		}
	} // if p<=3

	else

	{

		// This is a general implementation for arbitrary p

		static CDoubleVector left;
		static CDoubleVector right;
		left.ReSize(p + 1);
		right.ReSize(p + 1);

		double saved, temp, denom;

		left[0] = right[0] = 0.0;

		N[0] = 1.0;
		for(int j = 1; j <= p; j++)
		{
			left[j] = u - U[i+1-j];
			right[j] = U[i+j] - u;
			saved = 0.0;
			for(int r=0; r<j; r++)
			{
				denom = (right[r+1] + left[j-r]);
				temp = N[r] / denom;
				N[r] = saved + right[r+1] * temp;
				saved = left[j-r] * temp;
			}
			N[j] = saved;
		}
	} // p>3
}

void OneBasisFun(int p, int m, double* U, int i, double u, double& Nip)
{
	ASSERT(U != NULL);
	ASSERT(u >= 0 && u <= 1);
	ASSERT(m > 0);

	if(p<=3)
	{
		// This implementation pre-supposes that degree p <= 3
		// If p > 3, then N must be allocated dynamically
		// and not statically as follows.
		// Since, the degree of the curve is always 3, we choosed 
		// the second type of implementation to increase speed performance.
		ASSERT(p <= 3 );

		double Uleft = 0;
		double Uright = 0;
		double N[4];
		double saved = 0, temp = 0;

		if((i == 0 && u == U[0]) ||
			(i == m-p-1 && u == U[m])) // special cases
		{
			Nip = 1.; return;
		}

		if(u < U[i] || u >= U[i+p+1]) // local property
		{
			Nip =0.; return;
		}

		for(int j = 0; j <= p; j++) // Initialize zeroth-degree functs
			if(u >= U[i+j] && u < U[i+j+1]) N[j] = 1.0;
			else N[j] = 0.0;

		for(int k = 1; k <= p; k++) // Compute triangular table
		{
			if(N[0] == 0) saved = 0.;
			else saved = ((u-U[i])*N[0])/(U[i+k]-U[i]);
			for(int j = 0; j < p-k+1; j++)
			{
				Uleft = U[i+j+1];
				Uright = U[i+j+k+1];
				if(N[j+1] == 0)
				{
					N[j] = saved; saved = 0.;
				}
				else
				{
					temp = N[j+1]/(Uright-Uleft);
					N[j] = saved + (Uright-u)*temp;
					saved = (u-Uleft)*temp;
				}
			}
		}
		Nip = N[0];
	}
	else
	{
		double Uleft = 0;
		double Uright = 0;
		static CDoubleVector N; N.ReSize(p+1);
		double saved = 0, temp = 0;

		if((i == 0 && u == U[0]) ||
			(i == m-p-1 && u == U[m])) // special cases
		{
			Nip = 1.; return;
		}

		if(u < U[i] || u >= U[i+p+1]) // local property
		{
			Nip =0.; return;
		}

		for(int j = 0; j <= p; j++) // Initialize zeroth-degree functs
			if(u >= U[i+j] && u < U[i+j+1]) N[j] = 1.0;
			else N[j] = 0.0;

		for(int k = 1; k <= p; k++) // Compute triangular table
		{
			if(N[0] == 0) saved = 0.;
			else saved = ((u-U[i])*N[0])/(U[i+k]-U[i]);
			for(int j = 0; j < p-k+1; j++)
			{
				Uleft = U[i+j+1];
				Uright = U[i+j+k+1];
				if(N[j+1] == 0)
				{
					N[j] = saved; saved = 0.;
				}
				else
				{
					temp = N[j+1]/(Uright-Uleft);
					N[j] = saved + (Uright-u)*temp;
					saved = (u-Uleft)*temp;
				}
			}
		}
		Nip = N[0];
	}
}

void DersBasisFuns(int i, double u, int p, int n, const double* U, ddouble ders)
{
	// This algorithm is based on the algorithm A2.3 of
	// 'The NURBS Book', pp.72.

	ASSERT(U != NULL);
	ASSERT(ders != NULL);
	ASSERT(n<=p);

	// This implementation pre-supposes that degree p <= 3
	// If p > 3, then left and right must be allocated dynamically
	// and not statically as follows.
	// We choosed the second type of implementation due to speed performance.
	ASSERT(p <= 3 );

	// for arbitrary p
//	double** ndu = alloc_array(p + 1., p + 1.);
//	double** a = alloc_array(2., p + 1.);

	double ndu[4][4];
	double a[2][4];

	double left[4];
	double right[4];
	double saved, temp;
	memset(ders, 0, sizeof(double)*4*4);

	ndu[0][0] = 1.0;
	left[0] = right[0] = 0.0;

	int j, r, k;
	for(j = 1; j <= p; j++)
	{
		left[j] = u - U[i+1-j];
		right[j] = U[i+j] - u;
		saved = 0.0;
		for(r=0; r<j; r++)
		{
			ndu[j][r] = right[r+1] + left[j-r];
			temp = ndu[r][j-1] / ndu[j][r];

			ndu[r][j] = saved + right[r+1] * temp;
			saved = left[j-r] * temp;
		}
		ndu[j][j] = saved;
	}

	for(j = 0; j <= p; j++) // Load the basis functions
		ders[0][j] = ndu[j][p];

	// This section computes the derivates.
	// 'The NURBS Book', Eq.(2.9), pp.61

	double d;
	int s1,s2, rk, pk, j1, j2;
	for(r = 0; r <= p; r++)
	{
		s1 = 0; s2 = 1; // Alternate rows in array a
		a[0][0] = 1.0;

		// Loop to compute kth derivative
		for(k = 1; k<=n; k++)
		{
			d = 0.0;
			rk = r - k;
			pk = p - k;
			if(r >= k)
			{
				a[s2][0] = a[s1][0] / ndu[pk+1][rk];
				d = a[s2][0] * ndu[rk][pk];
			}

			if(rk >= -1)
				j1 = 1;
			else
				j1 = -rk;

			if(r-1 <= pk)
				j2 = k-1;
			else
				j2 = p-r;

			for(j = j1; j <= j2; j++)
			{
				a[s2][j] = (a[s1][j]-a[s1][j-1]) / ndu[pk+1][rk+j];
				d += a[s2][j] * ndu[rk+j][pk];
			}

			if(r <= pk)
			{
				a[s2][k] = -a[s1][k-1] / ndu[pk+1][r];
				d += a[s2][k] * ndu[r][pk];
			}

			ders[k][r] = d;
			j = s1; s1 = s2; s2 = j; // Switch rows
			ASSERT(s1 >= 0 && s1 <=1);
			ASSERT(s2 >= 0 && s2 <=1);
		}// k loop
	}// r loop

	// Multiply through by the correct factors
	// Eq.(2.9)
	r = p;
	for(k = 1; k <= n; k++)
	{
		for(j = 0; j<= p; j++)
			ders[k][j] *= r;
		r *= (p-k);
	}

//	disp_array(ndu, p+1.);
//	disp_array(a, 2.);
}

void DersBasisFunsGen(int i, double u, int p, int n, const double* U, CDoubleMatrix& ders)
{
	// This algorithm is based on the algorithm A2.3 of
	// 'The NURBS Book', pp.72.
	ASSERT(U != NULL);
	ASSERT(ders.GetData() != NULL);
	ASSERT(n<=p);


	// This algorithm is based on the algorithm A2.3 of
	// 'The NURBS Book', pp.72.
	// This is a general implementation for arbitrary p

	ASSERT(U != NULL);
	ASSERT(n<=p);

	// for arbitrary p
	static CDoubleMatrix ndu;
	static CDoubleMatrix a;
	static CDoubleVector left;
	static CDoubleVector right;
	ndu.ReSize(p+1, p+1);
	a.ReSize(2, p+1);
	left.ReSize(p + 1);
	right.ReSize(p + 1);

	double saved, temp;

	ndu(0,0) = 1.0;
	left[0] = right[0] = 0.0;
	ders.SetVal(0);

	int j, r, k;
	for(j = 1; j <= p; j++)
	{
		left[j] = u - U[i+1-j];
		right[j] = U[i+j] - u;
		saved = 0.0;
		for(r=0; r<j; r++)
		{
			ndu(j,r) = right[r+1] + left[j-r];
			temp = ndu(r,j-1) / ndu(j,r);
				ndu(r,j) = saved + right[r+1] * temp;
			saved = left[j-r] * temp;
		}
		ndu(j,j) = saved;
	}

	for(j = 0; j <= p; j++) // Load the basis functions
		ders(0,j) = ndu(j,p);

	// This section computes the derivates.
	// 'The NURBS Book', Eq.(2.9), pp.61

	double d;
	int s1,s2, rk, pk, j1, j2;
	for(r = 0; r <= p; r++)
	{
		s1 = 0; s2 = 1; // Alternate rows in array a
		a(0,0) = 1.0;

		// Loop to compute kth derivative
		for(k = 1; k<=n; k++)
		{
			d = 0.0;
			rk = r - k;
			pk = p - k;
			if(r >= k)
			{
				a(s2,0) = a(s1,0) / ndu(pk+1,rk);
				d = a(s2,0) * ndu(rk,pk);
			}

			if(rk >= -1)
				j1 = 1;
			else
				j1 = -rk;

			if(r-1 <= pk)
				j2 = k-1;
			else
				j2 = p-r;

			for(j = j1; j <= j2; j++)
			{
				a(s2,j) = (a(s1,j)-a(s1,j-1)) / ndu(pk+1,rk+j);
				d += a(s2,j) * ndu(rk+j,pk);
			}

			if(r <= pk)
			{
				a(s2,k) = -a(s1,k-1) / ndu(pk+1,r);
				d += a(s2,k) * ndu(r,pk);
			}

			ders(k,r) = d;
			j = s1; s1 = s2; s2 = j; // Switch rows
			ASSERT(s1 >= 0 && s1 <=1);
			ASSERT(s2 >= 0 && s2 <=1);
		}// k loop
	}// r loop

	// Multiply through by the correct factors
	// Eq.(2.9)
	r = p;
	for(k = 1; k <= n; k++)
	{
		for(j = 0; j<= p; j++)
			ders(k,j) *= r;
		r *= (p-k);
	}
}


CPointEx CurvePoint(int n, int p, const double* U, const CPointEx* P, double u, double* N)
{
	static int gspan = -1;
	int span = FindSpan(n, p, u, U, gspan); gspan =span;
	BasisFuns(span, u, p, U, N);
	CPointEx C;
	int k;
	for(int i = 0; i <= p; i++)
	{
		k = span-p+i;
		C.x += (N[i] * P[k].x);
		C.y += (N[i] * P[k].y);
	}

	return C;
}

CPointEx CurvePoint(int n, int p, const double* U, const CPointExVector& P, double u)
{
	CPointEx C(0,0);
	if(p <= 3)
	{
		double N[4];
		static int gspan = -1;
		int span = FindSpan(n, p, u, U, gspan); gspan =span;
		BasisFuns(span, u, p, U, N);
		int k;
		for(int i = 0; i <= p; i++)
		{
			k = span-p+i;
			C.x += (N[i] * P[k].x);
			C.y += (N[i] * P[k].y);
		}
	}
	else
	{
		static CDoubleVector N;
		N.ReSize(p+1);
		static int gspan = -1;
		int span = FindSpan(n, p, u, U, gspan); gspan =span;
		BasisFuns(span, u, p, U, N.GetData());
		int k;
		for(int i = 0; i <= p; i++)
		{
			k = span-p+i;
			C.x += (N[i] * P[k].x);
			C.y += (N[i] * P[k].y);
		}
	}

	return C;
}

CPointEx3D CurvePoint(int n, int p, const double* U, const CPointEx3D* P, double u, double* N)
{
	static int gspan = -1;
	int span = FindSpan(n, p, u, U, gspan); gspan =span;
	BasisFuns(span, u, p, U, N);
	CPointEx3D C;
	int k;
	for(int i = 0; i <= p; i++)
	{
		k = span-p+i;
		C.x += (N[i] * P[k].x);
		C.y += (N[i] * P[k].y);
		C.z += (N[i] * P[k].z);
	}

	return C;
}

CPointEx* CurveDerivs(int n, int p, const double* U, const CPointEx* P, double u, int d)
{
	int du = min(d,p); // ensure that we do not exceed splines degree
	ASSERT(du > 0 && du <= p);

	CPointEx* CK = NULL;
	CK = new CPointEx[du+1];
	ASSERT(CK != NULL);
	int k;

	for(k = p+1; k <= d; k++)
	{
		CK[k].x = CK[k].y = 0.0;
	}

	static int gspan = -1;
	int span = FindSpan(n, p, u, U, gspan); gspan =span;
	ddouble Derivs;
	DersBasisFuns(span, u, p, du, U, Derivs );

	for(k = 0; k <= du; k++)
	{
		CK[k].x = 0.0;
		CK[k].y = 0.0;
		for(int j = 0; j <= p; j++)
		{
			int i = span-p+j;
			CK[k].x += (Derivs[k][j]* P[i].x);
			CK[k].y += (Derivs[k][j]* P[i].y);
		}
	}

	return CK;
}

CPointEx3D* CurveDerivs(int n, int p, const double* U, const CPointEx3D* P, double u, int d)
{
	int du = min(d,p); // ensure that we do not exceed splines degree
	ASSERT(du > 0 && du <= p);

	CPointEx3D* CK = NULL;
	CK = new CPointEx3D[du+1];
	ASSERT(CK != NULL);
	int k;

	for(k = p+1; k <= d; k++)
	{
		CK[k].x = CK[k].y = CK[k].z = 0.0;
	}

	static int gspan = -1;
	int span = FindSpan(n, p, u, U, gspan); gspan =span;
	ddouble Derivs;
	DersBasisFuns(span, u, p, du, U, Derivs );

	for(k = 0; k <= du; k++)
	{
		CK[k].x = 0.0;
		CK[k].y = 0.0;
		CK[k].z = 0.0;
		for(int j = 0; j <= p; j++)
		{
			int i = span-p+j;
			CK[k].x += (Derivs[k][j]* P[i].x);
			CK[k].y += (Derivs[k][j]* P[i].y);
			CK[k].z += (Derivs[k][j]* P[i].z);
		}
	}

	return CK;
}

void CurveDerivs(int n, int p, const double* U, const CPointEx3DVector& P, double u, int d, CPointEx3DVector& CK)
{
	int du = min(d,p); // ensure that we do not exceed splines degree
	ASSERT(du > 0 && du <= p);

	int k;

	for(k = p+1; k <= d; k++)
	{
		CK[k].x = CK[k].y = CK[k].z = 0.0;
	}

	if(p <= 3)
	{
		static int gspan = -1;
		int span = FindSpan(n, p, u, U, gspan); gspan =span;
		ddouble Derivs;
		DersBasisFuns(span, u, p, du, U, Derivs );

		for(k = 0; k <= du; k++)
		{
			CK[k].x = 0.0;
			CK[k].y = 0.0;
			CK[k].z = 0.0;
			for(int j = 0; j <= p; j++)
			{
				int i = span-p+j;
				CK[k].x += (Derivs[k][j]* P[i].x);
				CK[k].y += (Derivs[k][j]* P[i].y);
				CK[k].z += (Derivs[k][j]* P[i].z);
			}
		}
	}
	else
	{
		static int gspan = -1;
		int span = FindSpan(n, p, u, U, gspan); gspan =span;
		static CDoubleMatrix Derivs;
		Derivs.ReSize(p+1, p+1);
		DersBasisFunsGen(span, u, p, du, U, Derivs );

		for(k = 0; k <= du; k++)
		{
			CK[k].x = 0.0;
			CK[k].y = 0.0;
			CK[k].z = 0.0;
			for(int j = 0; j <= p; j++)
			{
				int i = span-p+j;
				CK[k].x += (Derivs(k,j)* P[i].x);
				CK[k].y += (Derivs(k,j)* P[i].y);
				CK[k].z += (Derivs(k,j)* P[i].z);
			}
		}
	}
}

void CurveDerivs(int n, int p, const double* U, const CPointEx3D* P, double u, int d, CPointEx3D* CK)
{
	int du = min(d,p); // ensure that we do not exceed splines degree
	ASSERT(du > 0 && du <= p);

	int k;

	for(k = p+1; k <= d; k++)
	{
		CK[k].x = CK[k].y = CK[k].z = 0.0;
	}

	if(p <= 3)
	{
		static int gspan = -1;
		int span = FindSpan(n, p, u, U, gspan); gspan =span;
		ddouble Derivs;
		DersBasisFuns(span, u, p, du, U, Derivs );

		for(k = 0; k <= du; k++)
		{
			CK[k].x = 0.0;
			CK[k].y = 0.0;
			CK[k].z = 0.0;
			for(int j = 0; j <= p; j++)
			{
				int i = span-p+j;
				CK[k].x += (Derivs[k][j]* P[i].x);
				CK[k].y += (Derivs[k][j]* P[i].y);
				CK[k].z += (Derivs[k][j]* P[i].z);
			}
		}
	}
	else
	{
		static int gspan = -1;
		int span = FindSpan(n, p, u, U, gspan); gspan =span;
		static CDoubleMatrix Derivs;
		Derivs.ReSize(p+1, p+1);
		DersBasisFunsGen(span, u, p, du, U, Derivs );

		for(k = 0; k <= du; k++)
		{
			CK[k].x = 0.0;
			CK[k].y = 0.0;
			CK[k].z = 0.0;
			for(int j = 0; j <= p; j++)
			{
				int i = span-p+j;
				CK[k].x += (Derivs(k,j)* P[i].x);
				CK[k].y += (Derivs(k,j)* P[i].y);
				CK[k].z += (Derivs(k,j)* P[i].z);
			}
		}
	}
}

void CurveDerivs(int n, int p, const double* U, const CPointExVector& P, double u, int d, CPointExVector& CK)
{
	int du = min(d,p); // ensure that we do not exceed splines degree
	ASSERT(du > 0 && du <= p);

	int k;

	for(k = p+1; k <= d; k++)
	{
		CK[k].x = CK[k].y = 0.0;
	}

	if(p <= 3)
	{
		static int gspan = -1;
		int span = FindSpan(n, p, u, U, gspan); gspan =span;
		ddouble Derivs;
		DersBasisFuns(span, u, p, du, U, Derivs );

		for(k = 0; k <= du; k++)
		{
			CK[k].x = 0.0;
			CK[k].y = 0.0;
			for(int j = 0; j <= p; j++)
			{
				int i = span-p+j;
				CK[k].x += (Derivs[k][j]* P[i].x);
				CK[k].y += (Derivs[k][j]* P[i].y);
			}
		}
	}
	else
	{
		static int gspan = -1;
		int span = FindSpan(n, p, u, U, gspan); gspan =span;
		static CDoubleMatrix Derivs;
		Derivs.ReSize(p+1, p+1);
		DersBasisFunsGen(span, u, p, du, U, Derivs );

		for(k = 0; k <= du; k++)
		{
			CK[k].x = 0.0;
			CK[k].y = 0.0;
			for(int j = 0; j <= p; j++)
			{
				int i = span-p+j;
				CK[k].x += (Derivs(k,j)* P[i].x);
				CK[k].y += (Derivs(k,j)* P[i].y);
			}
		}
	}
}

void CurveDerivs(int n, int p, const double* U, const CPointEx* P, double u, int d, CPointEx* CK)
{
	int du = min(d,p); // ensure that we do not exceed splines degree
	ASSERT(du > 0 && du <= p);

	int k;

	for(k = p+1; k <= d; k++)
	{
		CK[k].x = CK[k].y = 0.0;
	}

	if(p <= 3)
	{
		static int gspan = -1;
		int span = FindSpan(n, p, u, U, gspan); gspan =span;
		ddouble Derivs;
		DersBasisFuns(span, u, p, du, U, Derivs );

		for(k = 0; k <= du; k++)
		{
			CK[k].x = 0.0;
			CK[k].y = 0.0;
			for(int j = 0; j <= p; j++)
			{
				int i = span-p+j;
				CK[k].x += (Derivs[k][j]* P[i].x);
				CK[k].y += (Derivs[k][j]* P[i].y);
			}
		}
	}
	else
	{
		static int gspan = -1;
		int span = FindSpan(n, p, u, U, gspan); gspan =span;
		static CDoubleMatrix Derivs;
		Derivs.ReSize(p+1, p+1);
		DersBasisFunsGen(span, u, p, du, U, Derivs );

		for(k = 0; k <= du; k++)
		{
			CK[k].x = 0.0;
			CK[k].y = 0.0;
			for(int j = 0; j <= p; j++)
			{
				int i = span-p+j;
				CK[k].x += (Derivs(k,j)* P[i].x);
				CK[k].y += (Derivs(k,j)* P[i].y);
			}
		}
	}
}

CPointEx CurveFirstDeriv(int n, int p, const double* U, const CPointEx* P, double u)
{
	if(p > 3)
		return CurveFirstDerivGen(n, p, U, P, u);

	int du = min(1,p); // ensure that we do not exceed splines degree
	ASSERT(du > 0 && du <= p);

	CPointEx CK(0,0);
	static int gspan = -1;
	int span = FindSpan(n, p, u, U, gspan); gspan =span;
	ddouble Derivs;
	DersBasisFuns(span, u, p, du, U, Derivs );

	for(int j = 0; j <= p; j++)
	{
		int i = span-p+j;
		CK.x += (Derivs[1][j]* P[i].x);
		CK.y += (Derivs[1][j]* P[i].y);
	}

	return CK;
}

CPointEx CurveFirstDerivGen(int n, int p, const double* U, const CPointEx* P, double u)
{
	int du = min(1,p); // ensure that we do not exceed splines degree
	ASSERT(du > 0 && du <= p);

	CPointEx CK(0,0);

	int k = 1;
	static int gspan = -1;
	int span = FindSpan(n, p, u, U, gspan); gspan =span;
	static CDoubleMatrix Derivs;
	Derivs.ReSize(p+1, p+1);
	DersBasisFunsGen(span, u, p, du, U, Derivs );

	for(int j = 0; j <= p; j++)
	{
		int i = span-p+j;
		CK.x += (Derivs(1,j)* P[i].x);
		CK.y += (Derivs(1,j)* P[i].y);
	}

	return CK;
}

CPointEx3D CurveFirstDeriv(int n, int p, const double* U, const CPointEx3D* P, double u)
{
	int du = min(1,p); // ensure that we do not exceed splines degree
	ASSERT(du > 0 && du <= p);

	CPointEx3D CK(0,0,0);
	static int gspan = -1;
	int span = FindSpan(n, p, u, U, gspan); gspan =span;
	ddouble Derivs;
	DersBasisFuns(span, u, p, du, U, Derivs );

	for(int j = 0; j <= p; j++)
	{
		int i = span-p+j;
		CK.x += (Derivs[1][j]* P[i].x);
		CK.y += (Derivs[1][j]* P[i].y);
		CK.z += (Derivs[1][j]* P[i].z);
	}

	return CK;
}

CPointEx CurveNormal(int n, int p, const double* U, const CPointEx* P, double u)
{
	ASSERT(P != NULL);
	ASSERT(U != NULL);

	CPointEx Normal;
	CPointEx3D N, C1, C2;
	double d;

	// This formula ensures that the normal vector will always rest
	// at one side of the curve (left side always).

	CPointEx deriv;
	if(p <= 3)
		deriv = CurveFirstDeriv(n, p, U, P, u);// use only the first deriv
	else
		deriv = CurveFirstDerivGen(n, p, U, P, u);// use only the first deriv

	C1.x = deriv.x; C1.y = deriv.y; C1.z = 0.0;  
	C2 = CPointEx3D(0,0,1);
	N = cross(C1,C2);
	d = N.norm();
	if(d == 0)
		N = CPointEx3D(0,0,0);
	else
		N *= (1.0 / d);

	// N is a vector in the plane
	Normal.x = N.x;
	Normal.y = N.y;

	return Normal;
}


// This functions finds a point in the curve which is the
// closest to the given point. It is used in the GetParam
// function to give a good intitial value for performing 
// the necessary optimization.
double GetClosestPoint(const CPointEx& point, int nDataPoints,
						int n, int p, const double* U, const CPointEx* P, 
						double* N)
{
	double umin, distmin, dist;

	CPointEx C;

	double u = 0;
	double step = 1.0 / nDataPoints; // compute an equally spaced parameter's step
	C = CurvePoint(n, p, U, P, u, N);
	distmin = SQR_DISTANCE_EX(C, point);
	umin = u;
	u += step;
	for (int i = 1; i < nDataPoints; i++ ) // evaluate points on curve
	{
		C = CurvePoint(n, p, U, P, u, N);

		dist = SQR_DISTANCE_EX(C, point);
		if(dist < distmin)
		{
			distmin = dist;
			umin = u;
		}

		u += step;
	}
	
	return umin;
}

// This functions finds a point in the curve which is the
// closest to the given point. It is used in the GetParam
// function to give a good intitial value for performing 
// the necessary optimization.
double GetClosestPoint(const CPointEx3D& point, int nDataPoints,
						int n, int p, const double* U, const CPointEx3D* P, 
						double* N)
{
	double umin, distmin, dist;

	CPointEx3D C;

	double u = 0;
	double step = 1.0 / nDataPoints; // compute an equally spaced parameter's step
	C = CurvePoint(n, p, U, P, u, N);
	distmin = sqr_distance(C, point);
	umin = u;
	u += step;
	for (int i = 1; i < nDataPoints; i++ ) // evaluate points on curve
	{
		C = CurvePoint(n, p, U, P, u, N);

		dist = sqr_distance(C, point);
		if(dist < distmin)
		{
			distmin = dist;
			umin = u;
		}

		u += step;
	}
	
	return umin;
}

BOOL CurveParam(int n, int p, const double* U, 
				 const CPointEx* P, double* N, 
				 const CPointEx& point, int nDataPoints,
				 double& u, double ug)
{
	if(p>3)
		return CurveParamGen(n, p, U, P, point, nDataPoints, u, ug);
	ASSERT(P != NULL);
	ASSERT(U != NULL);

	const double e1 = 1e-6;
	const double e2 = 1e-6;

	BOOL IsClosed = (P[0] == P[n]);

	double ui, ui1;

	double num, denom;

	CPointEx derivs[3];
	CPointEx C;

	// Get a good start value to achive reliable convergence.
	// In fact this is the most time consuming operation.
	ui = (ug>=0 && ug<=1)?ug:GetClosestPoint(point, nDataPoints, n, p, U, P, N);
	ASSERT(0 <= ui && ui <= 1);

	C = CurvePoint(n, p, U, P, ui, N);
	CPointEx P_to_C = C-point;
	double dist_P_to_C = P_to_C.norm();
	if(dist_P_to_C<=e1)
	{
		u = ui;
		return TRUE;
	}
	CurveDerivs(n, p, U, P, ui, 2, derivs);
	num = derivs[1] * P_to_C;
	denom = derivs[2] * P_to_C + derivs[1].sqr_norm();

	ASSERT(denom != 0);
	ui1 = ui - num / denom;

	if(!IsClosed)
	{
		if(ui1 < 0) ui1 = 0;
		else
		if(ui1 > 1) ui1 = 1;
	}
	else
	{
		if(ui1 < 0) ui1 += 1;
		else
		if(ui1 > 1) ui1 -= 1;
	}

	if(ui1 < 0 || ui1 > 1)
	{
		return FALSE; // something went wrong; abort operation
	}

	ASSERT(0 <= ui1 && ui1 <= 1);
	
	u = ui1; // assume ui1 is the solution

	double crit1, crit2, crit3;
	CPointEx temp;
	
	temp = derivs[1] * (ui1-ui);
	crit1 = temp.norm();
	if(crit1 <= e1)
	{
		return TRUE;
	}
	
	crit2 = dist_P_to_C;
	if(crit2 <= e1)
	{
		return TRUE;
	}

	crit3 = fabs(num) / ( derivs[1].norm() * crit2 );
	if(crit3 <= e2)
	{
		return TRUE;
	}

	int nCount = 0; // a counter to avoid infinite loop
	// In general 2 or 3 iterations will suffice
	while( !(crit1 <= e1) && !(crit2 <= e1) && !(crit3 <= e2) && nCount++ < 100)
	{
		ui = ui1;

		CurveDerivs(n, p, U, P, ui, 2, derivs);
		C = CurvePoint(n, p, U, P, ui, N);
		P_to_C = C-point;
		num = derivs[1] * P_to_C;
		denom = derivs[2] * P_to_C + derivs[1].sqr_norm();

		if(fabs(denom) < MACHINE_EPSILON)
			return FALSE;

		ASSERT(denom != 0);
		ui1 = ui - num / denom;

		if(!IsClosed)
		{
			if(ui1 < 0) ui1 = 0;
			else
			if(ui1 > 1) ui1 = 1;
		}
		else
		{
			if(ui1 < 0) ui1 += 1;
			else
			if(ui1 > 1) ui1 -= 1;
		}

		if(ui1 < 0 || ui1 > 1)
			return FALSE;

		u = ui1; // assume ui1 is the solution
		
		temp = derivs[1] * (ui1-ui);
		crit1 = temp.norm();
		if(crit1 <= e1)
			return TRUE;
	
		temp = P_to_C;
		crit2 = temp.norm();
		if(crit2 <= e1)
			return TRUE;

		crit3 = fabs(num) / ( derivs[1].norm() * crit2 );
		if(crit3 <= e2)
			return TRUE;
	}

	return FALSE;
}

BOOL CurveParam(int n, int p, const double* U, const CPointEx* P, 
				const CPointEx& point, int nDataPoints, double& u, double ug)
{
	if(p>3)
		return CurveParamGen(n, p, U, P, point, nDataPoints, u, ug);

	ASSERT(P != NULL);
	ASSERT(U != NULL);

	const double e1 = 1e-6;
	const double e2 = 1e-6;

	BOOL IsClosed = (P[0] == P[n]);

	double ui, ui1;

	double num, denom;

	CPointEx derivs[3];
	CPointEx C;

	double N[4];
//	int span = FindSpan(n, p, u, U);
//	BasisFuns(span, u, p, U, N);

	// Get a good start value to achive reliable convergence.
	// In fact this is the most time consuming operation.
	ui = (ug>=0 && ug<=1)?ug:GetClosestPoint(point, nDataPoints, n, p, U, P, N);
	ASSERT(0 <= ui && ui <= 1);

	C = CurvePoint(n, p, U, P, ui, N);
	CPointEx P_to_C = C-point;
	double dist_P_to_C = P_to_C.norm();
	if(dist_P_to_C<=e1)
	{
		u = ui;
		return TRUE;
	}
	CurveDerivs(n, p, U, P, ui, 2, derivs);
	num = derivs[1] * P_to_C;
	denom = derivs[2] * P_to_C + derivs[1].sqr_norm();

	ASSERT(denom != 0);
	ui1 = ui - num / denom;

	if(!IsClosed)
	{
		if(ui1 < 0) ui1 = 0;
		else
		if(ui1 > 1) ui1 = 1;
	}
	else
	{
		if(ui1 < 0) ui1 += 1;
		else
		if(ui1 > 1) ui1 -= 1;
	}

	if(ui1 < 0 || ui1 > 1)
	{
		return FALSE; // something went wrong; abort operation
	}

	ASSERT(0 <= ui1 && ui1 <= 1);
	
	u = ui1; // assume ui1 is the solution

	double crit1, crit2, crit3;
	CPointEx temp;
	
	temp = derivs[1] * (ui1-ui);
	crit1 = temp.norm();
	if(crit1 <= e1)
	{
		return TRUE;
	}
	
	crit2 = dist_P_to_C;
	if(crit2 <= e1)
	{
		return TRUE;
	}

	crit3 = fabs(num) / ( derivs[1].norm() * crit2 );
	if(crit3 <= e2)
	{
		return TRUE;
	}

	int nCount = 0; // a counter to avoid infinite loop
	// In general 2 or 3 iterations will suffice
	while( !(crit1 <= e1) && !(crit2 <= e1) && !(crit3 <= e2) && nCount++ < 100)
	{
		ui = ui1;

		CurveDerivs(n, p, U, P, ui, 2, derivs);
		C = CurvePoint(n, p, U, P, ui, N);
		P_to_C = C-point;
		num = derivs[1] * P_to_C;
		denom = derivs[2] * P_to_C + derivs[1].sqr_norm();

		if(fabs(denom) < MACHINE_EPSILON)
		{
			return FALSE;
		}
		ASSERT(denom != 0);
		ui1 = ui - num / denom;

		if(!IsClosed)
		{
			if(ui1 < 0) ui1 = 0;
			else
			if(ui1 > 1) ui1 = 1;
		}
		else
		{
			if(ui1 < 0) ui1 += 1;
			else
			if(ui1 > 1) ui1 -= 1;
		}

		if(ui1 < 0 || ui1 > 1)
			return FALSE;

		u = ui1; // assume ui1 is the solution
		
		temp = derivs[1] * (ui1-ui);
		crit1 = temp.norm();
		if(crit1 <= e1)
			return TRUE;
	
		temp = P_to_C;
		crit2 = temp.norm();
		if(crit2 <= e1)
			return TRUE;

		crit3 = fabs(num) / ( derivs[1].norm() * crit2 );
		if(crit3 <= e2)
			return TRUE;
	}

	return FALSE;
}

BOOL CurveParamGen(int n, int p, const double* U, const CPointEx* P, 
				const CPointEx& point, int nDataPoints, double& u, double ug)
{
	ASSERT(P != NULL);
	ASSERT(U != NULL);

	const double e1 = 1e-6;
	const double e2 = 1e-6;

	BOOL IsClosed = (P[0] == P[n]);

	double ui, ui1;

	double num, denom;

	CPointEx derivs[3];
	CPointEx C;

	static CDoubleVector N;
	N.ReSize(p+1);
//	int span = FindSpan(n, p, u, U);
//	BasisFuns(span, u, p, U, N.GetData());

	// Get a good start value to achive reliable convergence.
	// In fact this is the most time consuming operation.
	ui = (ug>=0 && ug<=1)?ug:GetClosestPoint(point, nDataPoints, n, p, U, P, N);
	ASSERT(0 <= ui && ui <= 1);

	C = CurvePoint(n, p, U, P, ui, N.GetData());
	CPointEx P_to_C = C-point;
	double dist_P_to_C = P_to_C.norm();
	if(dist_P_to_C<=e1)
	{
		u = ui;
		return TRUE;
	}
	CurveDerivs(n, p, U, P, ui, 2, derivs);
	num = derivs[1] * P_to_C;
	denom = derivs[2] * P_to_C + derivs[1].sqr_norm();

	ASSERT(denom != 0);
	ui1 = ui - num / denom;

	if(!IsClosed)
	{
		if(ui1 < 0) ui1 = 0;
		else
		if(ui1 > 1) ui1 = 1;
	}
	else
	{
		if(ui1 < 0) ui1 += 1;
		else
		if(ui1 > 1) ui1 -= 1;
	}

	if(ui1 < 0 || ui1 > 1)
	{
		return FALSE; // something went wrong; abort operation
	}

	ASSERT(0 <= ui1 && ui1 <= 1);
	
	u = ui1; // assume ui1 is the solution

	double crit1, crit2, crit3;
	CPointEx temp;
	
	temp = derivs[1] * (ui1-ui);
	crit1 = temp.norm();
	if(crit1 <= e1)
	{
		return TRUE;
	}
	
	crit2 = dist_P_to_C;
	if(crit2 <= e1)
	{
		return TRUE;
	}

	crit3 = fabs(num) / ( derivs[1].norm() * crit2 );
	if(crit3 <= e2)
	{
		return TRUE;
	}

	int nCount = 0; // a counter to avoid infinite loop
	// In general 2 or 3 iterations will suffice
	while( !(crit1 <= e1) && !(crit2 <= e1) && !(crit3 <= e2) && nCount++ < 100)
	{
		ui = ui1;

		CurveDerivs(n, p, U, P, ui, 2, derivs);
		C = CurvePoint(n, p, U, P, ui, N.GetData());
		P_to_C = C-point;
		num = derivs[1] * P_to_C;
		denom = derivs[2] * P_to_C + derivs[1].sqr_norm();

		if(fabs(denom) < MACHINE_EPSILON)
			return FALSE;

		ASSERT(denom != 0);
		ui1 = ui - num / denom;

		if(!IsClosed)
		{
			if(ui1 < 0) ui1 = 0;
			else
			if(ui1 > 1) ui1 = 1;
		}
		else
		{
			if(ui1 < 0) ui1 += 1;
			else
			if(ui1 > 1) ui1 -= 1;
		}

		if(ui1 < 0 || ui1 > 1)
		{
			return FALSE;
		}

		u = ui1; // assume ui1 is the solution
		
		temp = derivs[1] * (ui1-ui);
		crit1 = temp.norm();
		if(crit1 <= e1)
			return TRUE;
	
		temp = P_to_C;
		crit2 = temp.norm();
		if(crit2 <= e1)
			return TRUE;

		crit3 = fabs(num) / ( derivs[1].norm() * crit2 );
		if(crit3 <= e2)
			return TRUE;
	}

	return FALSE;
}

BOOL CurveParam(int n, int p, const double* U, 
				 const CPointEx3D* P, double* N, 
				 const CPointEx3D& point, int nDataPoints,
				 double& u, double ug)
{
	ASSERT(P != NULL);
	ASSERT(U != NULL);

	const double e1 = 1e-6;
	const double e2 = 1e-6;

	BOOL IsClosed = (P[0] == P[n]);

	double ui, ui1;

	double num, denom;

	CPointEx3D derivs[3];
	CPointEx3D C;

	// Get a good start value to achive reliable convergence.
	// In fact this is the most time consuming operation.
	ui = (ug>=0 && ug<=1)?ug:GetClosestPoint(point, nDataPoints, n, p, U, P, N);
	ASSERT(0 <= ui && ui <= 1);

	C = CurvePoint(n, p, U, P, ui, N);
	CPointEx3D P_to_C = C-point;
	double dist_P_to_C = P_to_C.norm();
	if(dist_P_to_C<=e1)
	{
		u = ui;
		return TRUE;
	}
	CurveDerivs(n, p, U, P, ui, 2, derivs);
	num = derivs[1] * P_to_C;
	denom = derivs[2] * P_to_C + derivs[1].sqr_norm();

	ASSERT(denom != 0);
	ui1 = ui - num / denom;

	if(!IsClosed)
	{
		if(ui1 < 0) ui1 = 0;
		else
		if(ui1 > 1) ui1 = 1;
	}
	else
	{
		if(ui1 < 0) ui1 += 1;
		else
		if(ui1 > 1) ui1 -= 1;
	}

	if(ui1 < 0 || ui1 > 1)
	{
		return FALSE; // something went wrong; abort operation
	}

	ASSERT(0 <= ui1 && ui1 <= 1);
	
	u = ui1; // assume ui1 is the solution

	double crit1, crit2, crit3;
	CPointEx3D temp;
	
	temp = derivs[1] * (ui1-ui);
	crit1 = temp.norm();
	if(crit1 <= e1)
	{
		return TRUE;
	}
	
	crit2 = dist_P_to_C;
	if(crit2 <= e1)
	{
		return TRUE;
	}

	crit3 = fabs(num) / ( derivs[1].norm() * crit2 );
	if(crit3 <= e2)
	{
		return TRUE;
	}

	int nCount = 0; // a counter to avoid infinite loop
	// In general 2 or 3 iterations will suffice
	while( !(crit1 <= e1) && !(crit2 <= e1) && !(crit3 <= e2) && nCount++ < 100)
	{
		ui = ui1;

		CurveDerivs(n, p, U, P, ui, 2, derivs);
		C = CurvePoint(n, p, U, P, ui, N);
		P_to_C = C-point;
		num = derivs[1] * P_to_C;
		denom = derivs[2] * P_to_C + derivs[1].sqr_norm();

		if(fabs(denom) < 1e-16)
			return FALSE;

		ASSERT(denom != 0);
		ui1 = ui - num / denom;

		if(!IsClosed)
		{
			if(ui1 < 0) ui1 = 0;
			else
			if(ui1 > 1) ui1 = 1;
		}
		else
		{
			if(ui1 < 0) ui1 += 1;
			else
			if(ui1 > 1) ui1 -= 1;
		}

		if(ui1 < 0 || ui1 > 1)
			return FALSE;

		u = ui1; // assume ui1 is the solution
		
		temp = derivs[1] * (ui1-ui);
		crit1 = temp.norm();
		if(crit1 <= e1)
			return TRUE;
	
		temp = P_to_C;
		crit2 = temp.norm();
		if(crit2 <= e1)
			return TRUE;

		crit3 = fabs(num) / ( derivs[1].norm() * crit2 );
		if(crit3 <= e2)
			return TRUE;
	}

	return FALSE;
}

void FindCtrlPointsK_K1(const CPointEx& Q0, const CPointEx& Q3, 
						const CPointEx& T0, const CPointEx& T3,
						CPointEx& P1, CPointEx& P2)
{
	// This function computes two points P[1] and P[2] in order
	// a cubic B-spline to pass through P[0] and P[3].
	// See, 'The NURBS Book', pp395-397, for the following formulae.

	CPointEx T03(T0.x + T3.x, T0.y + T3.y);

	double a = 16 - SQR_NORM_EX(T03);
	double b = 12 * ( Q3 - Q0) * (T03);
	double c = -36 * SQR_NORM_EX( (Q3 - Q0) );
	double two_a = 2 * a;

	double D = b*b - 2 * two_a * c;
	ASSERT(D >= 0);
	ASSERT( two_a != 0);
	double sqrt_D = sqrt(D);

	double x1 = ( -b + sqrt_D ) / two_a;
	double x2 = ( -b - sqrt_D ) / two_a;

	ASSERT( !(x1<0 && x2<0) );
	ASSERT( !(x1>0 && x2>0) );
	
	double x = (x1 >= 0) ? x1 : x2;
	ASSERT(x >= 0);

	double x_div_3 = x / 3.0;
	P1.x = Q0.x + T0.x * x_div_3;
	P1.y = Q0.y + T0.y * x_div_3;
	P2.x = Q3.x - T3.x * x_div_3;
	P2.y = Q3.y - T3.y * x_div_3;
}

void FindCtrlPointsK_K1(const CPointEx3D& Q0, const CPointEx3D& Q3, 
						const CPointEx3D& T0, const CPointEx3D& T3,
						CPointEx3D& P1, CPointEx3D& P2)
{
	// This function computes two points P[1] and P[2] in order
	// a cubic B-spline to pass through P[0] and P[3].
	// See, 'The NURBS Book', pp395-397, for the following formulae.

	CPointEx3D T03(T0.x + T3.x, T0.y + T3.y, T0.z + T3.z);
	CPointEx3D Q30(Q3.x - Q0.x, Q3.y - Q0.y, Q3.z - Q0.z);

	double a = 16 - T03.sqr_norm();
	double b = 12 * (Q30 * T03);
	double c = -36 * sqr_norm( Q30 );
	double two_a = 2*a;

	double D = b*b - 2 * two_a * c;
	ASSERT(D >= 0);
	ASSERT(a != 0);

	double sqrt_D = sqrt(D);
	double x1 = ( -b + sqrt_D ) / two_a;
	double x2 = ( -b - sqrt_D ) / two_a;

	ASSERT( !(x1<0 && x2<0) );
	ASSERT( !(x1>0 && x2>0) );
	
	double x = (x1 >= 0) ? x1 : x2;
	ASSERT(x >= 0);

	double x_div_3 = x / 3.0;
	P1.x = Q0.x + T0.x * x_div_3;
	P1.y = Q0.y + T0.y * x_div_3;
	P1.z = Q0.z + T0.z * x_div_3;
	P2.x = Q3.x - T3.x * x_div_3;
	P2.y = Q3.y - T3.y * x_div_3;
	P2.z = Q3.z - T3.z * x_div_3;
}

void FindCtrlPoints(int nPoints, const CPointEx* Q, const CPointEx* T,
					 int nCtrlPoints, int m,
					 LPPOINTEX* P, CDoubleVector& U)
{
	// We construct a curve with equal speed at each Qk and at the 
	// midpoint of each Bezier segment. See, Eq.(9.52), Eq.(9.53) and Eq.(9.54)
	// See, 'The NURBS Book', pp398, for the following formulae.


	ASSERT(nPoints > 2);
	ASSERT(Q != NULL);
	ASSERT(T != NULL);
	ASSERT(U.GetData() != NULL);
	ASSERT(P != NULL);

	int n = nPoints - 1;
	(*P)[0] = Q[0];
	(*P)[nCtrlPoints-1] = Q[n];

	CPointEx PP1;
	CPointEx PP2;

	double u = 0;

	int k, two_k;
	for(k = 0; k < n; k++)
	{
		FindCtrlPointsK_K1(Q[k], Q[k+1], T[k], T[k+1], PP1, PP2 );
		two_k = 2*k;
		(*P)[two_k+1] = PP1;
		(*P)[two_k+2] = PP2;

		u = u + 3 * NORM_EX((PP1 - Q[k]));

		U[two_k + 4] = u;
		U[two_k + 5] = u;//(*U)[two_k + 4];
	}

	ASSERT(u > 0);
	for(k=0; k<n; k++)
	{
		two_k = 2*k;
		U[two_k + 4] = U[two_k + 4] / u;
		U[two_k + 5] = U[two_k + 4];
	}

	for(k=0; k<=3; k++)
	{
		U[k] = 0.0;
		U[m-k] = 1.0;
	}
}

void FindCtrlPoints(int nPoints, const CPointEx3D* Q, const CPointEx3D* T,
					 int nCtrlPoints, int m,
					 LPPOINTEX3D* P, CDoubleVector& U)
{
	ASSERT(nPoints > 2);
	ASSERT(Q != NULL);
	ASSERT(T != NULL);
	ASSERT(U.GetData() != NULL);
	ASSERT(P != NULL);

	int n = nPoints - 1;
	(*P)[0] = Q[0];
	(*P)[nCtrlPoints-1] = Q[n];

	CPointEx3D PP1;
	CPointEx3D PP2;

	double u = 0;

	int k, two_k;
	for(k = 0; k < n; k++)
	{
		FindCtrlPointsK_K1(Q[k], Q[k+1], T[k], T[k+1], PP1, PP2 );
		two_k = 2*k;
		(*P)[two_k+1] = PP1;
		(*P)[two_k+2] = PP2;

		u = u + 3 * (PP1 - Q[k]).norm();

		U[two_k + 4] = u;
		U[two_k + 5] = u; //(*U)[two_k + 4];
	}

	ASSERT(u > 0);
	for(k=0; k<n; k++)
	{
		two_k = 2*k;
		U[two_k + 4] = U[two_k + 4] / u;
		U[two_k + 5] = U[two_k + 4];
	}

	for(k=0; k<=3; k++)
	{
		U[k] = 0.0;
		U[m-k] = 1.0;
	}
}

void FindTans(int nPoints, const CPointEx* Q, LPPOINTEX* T)
{
	int n = nPoints - 1;
	ASSERT(n >= 2);

	ASSERT(Q != NULL);
	ASSERT(T != NULL);

	static CPointExVector q;
	q.ReSize(n+4);

	int k;
	
	for(k = 1; k <= n; k++)
	{
		q[k + 1].x = Q[k].x - Q[k-1].x;
		q[k + 1].y = Q[k].y - Q[k-1].y;
	}

	q[0 +1] = 2 * q[1 +1] - q[2 +1];
	q[-1 +1] = 2 * q[0 +1] - q[1 +1];

	q[n+1 +1] = 2 * q[n +1] - q[n-1 +1];
	q[n+2 +1] = 2 * q[n+1 +1] - q[n +1];

	double ak;
	CPointEx Vk;
	double qq1, qq2;
	double sum, norm;

	for(k=0; k<=n; k++)
	{
		qq1 = fabs( CROSS_EX(q[k-1 +1], q[k +1]));
		qq2 = fabs( CROSS_EX(q[k+1 +1], q[k+2 +1]));
		sum = qq1 + qq2;
		if(sum < 1e-6)
		{
			ak = 1; // this produces a corner.
		}
		else
		{
			ak = qq1 / sum; // formula (9.31)
		}

		Vk = (1-ak) * q[k +1] + ak * q[k+1 +1];

		norm = Vk.norm();
		if(norm > 0)
		{
			(*T)[k] = Vk * (1.0 / norm);
		}
	}

	// Special case; closed polygon
	// Set the last and first tangent = (T_0 + T_n) / 2
	if(Q[0] == Q[n])
	{
		(*T)[n] += (*T)[0];
		(*T)[n] *= 0.5;
		norm = (*T)[n].norm();
		if(norm > 0)
		{
			(*T)[n] /= norm;
		}
		(*T)[0] = (*T)[n];
	}
}


void FindTans(int nPoints, const CPointEx3D* Q, LPPOINTEX3D* T)
{
	int n = nPoints - 1;
	ASSERT(n >= 2);

	ASSERT(Q != NULL);
	ASSERT(T != NULL);

	static CPointEx3DVector q;
	q.ReSize(n+4);

	int k;
	
	for(k = 1; k <= n; k++)
	{
		q[k + 1].x = Q[k].x - Q[k-1].x;
		q[k + 1].y = Q[k].y - Q[k-1].y;
		q[k + 1].z = Q[k].z - Q[k-1].z;
	}

	q[0 +1] = 2 * q[1 +1] - q[2 +1];
	q[-1 +1] = 2 * q[0 +1] - q[1 +1];

	q[n+1 +1] = 2 * q[n +1] - q[n-1 +1];
	q[n+2 +1] = 2 * q[n+1 +1] - q[n +1];

	double ak;
	CPointEx3D Vk;
	double qq1, qq2;
	double sum, norm;

	for(k=0; k<=n; k++)
	{
		qq1 = ( cross(q[k-1 +1], q[k +1])).norm();
		qq2 = ( cross(q[k+1 +1], q[k+2 +1])).norm();
		sum = qq1 + qq2;
		if(sum < 1e-6)
		{
			ak = 1; // this produces a corner.
		}
		else
		{
			ak = qq1 / sum; // formula (9.31)
		}

		Vk = (1-ak) * q[k +1] + ak * q[k+1 +1];

		norm = Vk.norm();
		if(norm > 0)
		{
			(*T)[k].x = Vk.x / norm;
			(*T)[k].y = Vk.y / norm;
			(*T)[k].z = Vk.z / norm;
		}
	}

	// Special case; closed polygon
	// Set the last and first tangent = (T_0 + T_n) / 2
	if(Q[0] == Q[n])
	{
		(*T)[n] += (*T)[0];
		(*T)[n] *= 0.5;
		norm = (*T)[n].norm();
		if(norm > 0)
		{
			(*T)[n] /= norm;
		}
		(*T)[0] = (*T)[n];
	}
}

void RefineKnotVectCurve(int n, int p, const double* U, const CPointEx* P, const double* X, int r, double* Ubar, CPointEx* Q)
{
	// This is a generalized knot insertion algorithm, which allows
	// to insert the elements of vector X into the knot vector.
	// See, Algorithm (A5.4) of 'The NURBS Book', pp.162-167

	ASSERT(U != NULL);
	ASSERT(Ubar != NULL);
	ASSERT(P != NULL);
	ASSERT(Q != NULL);
	ASSERT(X != NULL);

	ASSERT(X[0] >=0 && X[0] <= 1);
	ASSERT(X[r] >=0 && X[r] <= 1);

	int m = n + p + 1;
	int a = FindSpan(n, p, X[0], U);
	int b = FindSpan(n, p, X[r], U);

	b = b + 1;
	int j;
	for(j = 0; j <= a - p; j++)	Q[j] = P[j];
	for(j = b - 1; j <= n; j++)	Q[j + r + 1] = P[j];

	for(j = 0; j <= a; j++)	Ubar[j] = U[j];
	for(j = b + p; j <= m; j++)	Ubar[j + r + 1] = U[j];
	
	int i = b + p - 1;
	int k = b + p + r;
	for(j = r; j >= 0; j--)
	{
		while( X[j] <= U[i] && i > a)
		{
			Q[k - p - 1] = P[i - p - 1];
			Ubar[k] = U[i];
			k = k - 1;
			i = i - 1;
		}

		Q[k - p - 1] = Q[k - p];
		for(int l = 1; l <= p; l++)
		{
			int ind = k - p + l;
			double alfa = Ubar[k+l] - X[j];

			if(fabs(alfa) == 0.0)
				Q[ind - 1] = Q[ind];
			else
			{
				alfa = alfa / (Ubar[k+l] - U[i-p+l]);
				Q[ind - 1] = alfa * Q[ind - 1] + ( 1.0 - alfa ) * Q[ind];
			}
		}
		Ubar[k] = X[j];
		k = k - 1;
	}
}

void RefineKnotVectCurve(int n, int p, const double* U, const CPointEx3D* P, const double* X, int r, double* Ubar, CPointEx3D* Q)
{
	// This is a generalized knot insertion algorithm, which allows
	// to insert the elements of vector X into the knot vector.
	// See, Algorithm (A5.4) of 'The NURBS Book', pp.162-167

	ASSERT(U != NULL);
	ASSERT(Ubar != NULL);
	ASSERT(P != NULL);
	ASSERT(Q != NULL);
	ASSERT(X != NULL);

	ASSERT(X[0] >=0 && X[0] <= 1);
	ASSERT(X[r] >=0 && X[r] <= 1);

	int m = n + p + 1;
	int a = FindSpan(n, p, X[0], U);
	int b = FindSpan(n, p, X[r], U);

	b = b + 1;
	int j;
	for(j = 0; j <= a - p; j++)	Q[j] = P[j];
	for(j = b - 1; j <= n; j++)	Q[j + r + 1] = P[j];

	for(j = 0; j <= a; j++)	Ubar[j] = U[j];
	for(j = b + p; j <= m; j++)	Ubar[j + r + 1] = U[j];
	
	int i = b + p - 1;
	int k = b + p + r;
	for(j = r; j >= 0; j--)
	{
		while( X[j] <= U[i] && i > a)
		{
			Q[k - p - 1] = P[i - p - 1];
			Ubar[k] = U[i];
			k = k - 1;
			i = i - 1;
		}

		Q[k - p - 1] = Q[k - p];
		for(int l = 1; l <= p; l++)
		{
			int ind = k - p + l;
			double alfa = Ubar[k+l] - X[j];

			if(fabs(alfa) == 0.0)
				Q[ind - 1] = Q[ind];
			else
			{
				alfa = alfa / (Ubar[k+l] - U[i-p+l]);
				Q[ind - 1] = alfa * Q[ind - 1] + ( 1.0 - alfa ) * Q[ind];
			}
		}
		Ubar[k] = X[j];
		k = k - 1;
	}
}

void RemoveCurveKnot(int n, int p, double* U, CPointEx* P, 
					  double u, int r, int s, int num, int& t, double error)
{
	// This is algorithm A5.8 of "The NURBS Book", pp.179-188.
	// Using this algorithm we remove any unnecassery knot
	// and control points of the curve.
	// This algorithm is perticularly usefull after
	// the curve-curve intersections (or any recursive subdivision)
	// algorithm in order to remove the inserted knots.

	int m = n + p + 1;
	int ord = p + 1;
	int fout = (2*r - s - p) / 2; // first ctrl point out
	int last = r - s;
	int first = r - p;
	int off, i, j, ii, jj, remflag, k;
	double alfi, alfj;

//	static CPointEx temp[10]; // this implies that 2p+1 < 10 (p<=4)
	static CPointExVector temp; temp.ReSize(2*p+1);

	for(t = 0; t < num; t++)
	{
		// This loop is Eq.(5.28)
		off = first - 1;
		temp[0] = P[off];
		temp[last+1-off] = P[last+1];
		i = first; j = last;
		ii = 1; jj = last - off;
		remflag = 0;

		while(j-i > t)
		{
			// Compute new ctrl points for one removal step
			alfi = (u - U[i]) / (U[i+ord+t] - U[i]);
			alfj = (u - U[j-t]) / (U[j+ord] - U[j-t]);
			temp[ii] = (P[i] - (1.0-alfi)*temp[ii-1] ) / alfi;
			temp[jj] = (P[j] - alfj*temp[jj+1]) / (1.0-alfj);
			i = i+1; ii = ii + 1;
			j = j - 1; jj = jj - 1;
		}// end while-loop

		if(j-i < t) // check if knot removable
		{
			if(DISTANCE_EX(temp[ii-1], temp[jj+1]) <= error)
				remflag = 1;
		}
		else
		{
			alfi = (u - U[i]) / (U[i+ord+t] - U[i]);
			if(DISTANCE_EX(P[i], alfi*temp[ii+t+1] + (1.0-alfi)*temp[ii-1]) <= error)
				remflag = 1;
		}

		if(remflag == 0) // cannot remove any more knots
			break; // get out of for-loop
		else
		{
			// successful removal. Seve new ctrl pts.
			i = first; j = last;
			while(j-i > t)
			{
				P[i] = temp[i-off];
				P[j] = temp[j-off];
				i += 1; j -= 1;
			}
		}
		first = first - 1; last = last + 1;
	}// end of for-loop

	if(t == 0)
		return;

	for(k = r+1; k<=m; k++) U[k-t] = U[k]; // shift knots left
	j = fout; i= j; // Pj through Pi will be overwritten
	for(k = 1; k < t; k++)
		if( (k % 2) == 1) // k modulo 2
			i = i + 1;
		else
			j = j - 1;
	for(k = i+1; k <= n; k++) // Shift
	{
		P[j] = P[k];
		j = j + 1;
	}
}

void RemoveCurveKnot(int n, int p, double* U, CPointEx3D* P, 
					  double u, int r, int s, int num, int& t, double error)
{
	// This is algorithm A5.8 of "The NURBS Book", pp.179-188.
	// Using this algorithm we remove any unnecassery knot
	// and control points of the curve.
	// This algorithm is perticularly usefull after
	// the curve-curve intersections (or any recursive subdivision)
	// algorithm in order to remove the inserted knots.

	int m = n + p + 1;
	int ord = p + 1;
	int fout = (2*r - s - p) / 2; // first ctrl point out
	int last = r - s;
	int first = r - p;
	int off, i, j, ii, jj, remflag, k;
	double alfi, alfj;

//	CPointEx3D temp[10]; // this implies that 2p+1 < 10 (p<=4)
	static CPointEx3DVector temp; temp.ReSize(2*p+1);

	for(t = 0; t < num; t++)
	{
		// This loop is Eq.(5.28)
		off = first - 1;
		temp[0] = P[off];
		temp[last+1-off] = P[last+1];
		i = first; j = last;
		ii = 1; jj = last - off;
		remflag = 0;

		while(j-i > t)
		{
			// Compute new ctrl points for one removal step
			alfi = (u - U[i]) / (U[i+ord+t] - U[i]);
			alfj = (u - U[j-t]) / (U[j+ord] - U[j-t]);
			temp[ii] = (P[i] - (1.0-alfi)*temp[ii-1] ) / alfi;
			temp[jj] = (P[j] - alfj*temp[jj+1]) / (1.0-alfj);
			i = i+1; ii = ii + 1;
			j = j - 1; jj = jj - 1;
		}// end while-loop

		if(j-i < t) // check if knot removable
		{
			if(distance(temp[ii-1], temp[jj+1]) <= error)
				remflag = 1;
		}
		else
		{
			alfi = (u - U[i]) / (U[i+ord+t] - U[i]);
			if(distance(P[i], alfi*temp[ii+t+1] + (1.0-alfi)*temp[ii-1]) <= error)
				remflag = 1;
		}

		if(remflag == 0) // cannot remove any more knots
			break; // get out of for-loop
		else
		{
			// successful removal. Seve new ctrl pts.
			i = first; j = last;
			while(j-i > t)
			{
				P[i] = temp[i-off];
				P[j] = temp[j-off];
				i += 1; j -= 1;
			}
		}
		first = first - 1; last = last + 1;
	}// end of for-loop

	if(t == 0)
		return;

	for(k = r+1; k<=m; k++) U[k-t] = U[k]; // shift knots left
	j = fout; i= j; // Pj through Pi will be overwritten
	for(k = 1; k < t; k++)
		if( (k % 2) == 1) // k modulo 2
			i = i + 1;
		else
			j = j - 1;
	for(k = i+1; k <= n; k++) // Shift
	{
		P[j] = P[k];
		j = j + 1;
	}
}

CPointEx3D CurvePoint(int n, int p, const double* U, const CPointEx3DVector& P, double u)
{
	CPointEx3D C(0,0,0);
	if(p <= 3)
	{
		double N[4];
		static int gspan = -1;
		int span = FindSpan(n, p, u, U, gspan); gspan =span;
		BasisFuns(span, u, p, U, N);
		int k;
		for(int i = 0; i <= p; i++)
		{
			k = span-p+i;
			C.x += (N[i] * P[k].x);
			C.y += (N[i] * P[k].y);
			C.z += (N[i] * P[k].z);
		}
	}
	else
	{
		static CDoubleVector N;
		N.ReSize(p+1);
		static int gspan = -1;
		int span = FindSpan(n, p, u, U, gspan); gspan =span;
		BasisFuns(span, u, p, U, N.GetData());
		int k;
		for(int i = 0; i <= p; i++)
		{
			k = span-p+i;
			C.x += (N[i] * P[k].x);
			C.y += (N[i] * P[k].y);
			C.z += (N[i] * P[k].z);
		}
	}

	return C;
}

/////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////
/* SPECIAL FUNCTIONS ///////////////////////////////////////////////////////*/
/* Used to Evaluate a B-spline basis and their derivatives at a given control point */

double NURBS(int i, int p, double u, double* U)
{
	if(u < U[i] || u >= U[i+p+1])
		return 0.0;

	if (p==0)
	{
		if ( (U[i]<=u) && (u<U[i+1]) )
			return (1.0);
		else 
			return (0.0);
	}
	
	double N1,N2;
	double DU1,DU2;
	
	DU1 = U[i+p] - U[i];
	DU2 = U[i+p+1] - U[i+1];
	
	if(DU1!=0) 
		N1 = NURBS(i, p-1, u, U) / DU1;
	else 
		N1 = 0.0;
	
	if(DU2!=0) 
		N2 = NURBS(i+1,p-1,u, U) / DU2;
	else 
		N2 = 0.0;
	
	return ( (u - U[i]) * N1 + (U[i+p+1] - u) * N2 );
}

double NURBS_K(int k, int n, int p, const double* U, double u, double* N)
{
	static int gspan = -1;
	int span = FindSpan(n, p, u, U, gspan); gspan = span;
//	int span = FindSpan(n, p, u, U);
	BasisFuns(span, u, p, U, N);
	int ki;
	for(int i = 0; i <= p; i++)
	{
		ki = span-p+i;
		if(ki == k)
			return N[i];
	}
/*	}
	else
	{
		static int gspan = -1;
		int span = FindSpan(n, p, u, U, gspan); gspan = span;
//		int span = FindSpan(n, p, u, U);
		BasisFuns(span, u, p, U, N);
		int ki;
		for(int i = 0; i <= p; i++)
		{
			ki = span-p+i;
			if(ki == k)
				return N[i];
		}
	}*/

	return 0.0;
}

void NURBS_K_DERIV(int K, int n, int p, const double* U, double u, int d, 
				   const CPointExVector& P, 
				   double& DERIV, CPointExVector& CK)
{
	int du = min(d,p); // ensure that we do not exceed splines degree
	ASSERT(du > 0 && du <= p);

	DERIV = 0;
	int k;
	for(k = p+1; k <= d; k++)
	{
		CK[k].x = CK[k].y = 0.0;
	}

	if(p <= 3)
	{
		static int gspan = -1;
		int span = FindSpan(n, p, u, U, gspan); gspan = span;
//		int span = FindSpan(n, p, u, U);
		static ddouble Derivs;
		DersBasisFuns(span, u, p, du, U, Derivs );

		for(k = 0; k <= du; k++)
		{
			CK[k].x = 0.0;
			CK[k].y = 0.0;
			for(int j = 0; j <= p; j++)
			{
				int i = span-p+j;
				if(K == i)
				{
					DERIV = Derivs[k][j];
				}
				CK[k].x += (Derivs[k][j]* P[i].x);
				CK[k].y += (Derivs[k][j]* P[i].y);
			}
		}
	}
	else
	{
		static int gspan = -1;
		int span = FindSpan(n, p, u, U, gspan); gspan = span;
//		int span = FindSpan(n, p, u, U);
		static CDoubleMatrix Derivs;
		Derivs.ReSize(p+1, p+1);
		DersBasisFunsGen(span, u, p, du, U, Derivs );

		for(k = 0; k <= du; k++)
		{
			CK[k].x = 0.0;
			CK[k].y = 0.0;
			for(int j = 0; j <= p; j++)
			{
				int i = span-p+j;
				if(K == i)
				{
					DERIV = Derivs(k,j);
				}
				CK[k].x += (Derivs(k,j)* P[i].x);
				CK[k].y += (Derivs(k,j)* P[i].y);
			}
		}
	}
}

void NURBS_K_DERIVS(int K, int n, int p, const double* U, double u, int d, 
				   const CPointExVector& P, 
				   CDoubleVector& DERIV, CPointExVector& CK)
{
	int du = min(d,p); // ensure that we do not exceed splines degree
	ASSERT(du > 0 && du <= p);

	int k;
	for(k = p+1; k <= d; k++)
	{
		CK[k].x = CK[k].y = 0.0;
	}

	if(p <= 3)
	{
		static int gspan = -1;
		int span = FindSpan(n, p, u, U, gspan); gspan = span;
//		int span = FindSpan(n, p, u, U);
		static ddouble Derivs;
		DersBasisFuns(span, u, p, du, U, Derivs );

		for(k = 0; k <= du; k++)
		{
			CK[k].x = 0.0;
			CK[k].y = 0.0;
			for(int j = 0; j <= p; j++)
			{
				int i = span-p+j;
				if(K == i)
				{
					DERIV[k] = Derivs[k][j];
				}
				CK[k].x += (Derivs[k][j]* P[i].x);
				CK[k].y += (Derivs[k][j]* P[i].y);
			}
		}
	}
	else
	{
		static int gspan = -1;
		int span = FindSpan(n, p, u, U, gspan); gspan = span;
//		int span = FindSpan(n, p, u, U);
		static CDoubleMatrix Derivs;
		Derivs.ReSize(p+1, p+1);
		DersBasisFunsGen(span, u, p, du, U, Derivs );

		for(k = 0; k <= du; k++)
		{
			CK[k].x = 0.0;
			CK[k].y = 0.0;
			for(int j = 0; j <= p; j++)
			{
				int i = span-p+j;
				if(K == i)
				{
					DERIV[k] = Derivs(k,j);
				}
				CK[k].x += (Derivs(k,j)* P[i].x);
				CK[k].y += (Derivs(k,j)* P[i].y);
			}
		}
	}
}

/////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////
// DISCRETE COMPUTATIONS
void DiscreteCurvature(const CPointExVector& p, CDoubleVector& curv)
{
	ASSERT(p.GetSize()==curv.GetSize());
	ASSERT(p.GetSize()>=3);
	int N = p.GetSize();
	curv[0]=0; curv[N-1] = 0; // assign zero to first and last pnts
	for(int i=1; i<N-1;i++)
		curv[i] = (p[i-1]-2*p[i]+p[i+1]).norm();
}

void DiscreteCurvature(const CPointExVector& p, CDoubleVector& curv, CPointExVector& curvVec)
{
	ASSERT(p.GetSize()==curv.GetSize());
	ASSERT(p.GetSize()==curvVec.GetSize());
	ASSERT(p.GetSize()>=3);
	int N = p.GetSize();
	curv[0]=0; curv[N-1] = 0; // assign zero to first and last pnts
	for(int i=1; i<N-1;i++)
	{
		curvVec[i] = (p[i-1]-2*p[i]+p[i+1]);
		curv[i] = curvVec[i].norm();
	}
}

void DiscreteCurvature(const CPointEx3DVector& p, CDoubleVector& curv, CPointEx3DVector& curvVec)
{
	ASSERT(p.GetSize()==curv.GetSize());
	ASSERT(p.GetSize()==curvVec.GetSize());
	ASSERT(p.GetSize()>=3);
	int N = p.GetSize();
	curv[0]=0; curv[N-1] = 0; // assign zero to first and last pnts
	for(int i=1; i<N-1;i++)
	{
		curvVec[i] = (p[i-1]-2*p[i]+p[i+1]);
		curv[i] = curvVec[i].norm();
	}
}

void DiscreteCurvature(const CPointEx3DVector& p, CDoubleVector& curv)
{
	ASSERT(p.GetSize()==curv.GetSize());
	ASSERT(p.GetSize()>=3);
	int N = p.GetSize();
	curv[0]=0; curv[N-1] = 0; // assign zero to first and last pnts
	for(int i=1; i<N-1;i++)
		curv[i] = (p[i-1]-2*p[i]+p[i+1]).norm();
}

void SqrDiscreteCurvature(const CPointExVector& p, CDoubleVector& curv)
{
	ASSERT(p.GetSize()==curv.GetSize());
	ASSERT(p.GetSize()>=3);
	int N = p.GetSize();
	curv[0]=0; curv[N-1] = 0; // assign zero to first and last pnts
	for(int i=1; i<N-1;i++)
		curv[i] = (p[i-1]-2*p[i]+p[i+1]).sqr_norm();
}

void SqrDiscreteCurvature(const CPointEx3DVector& p, CDoubleVector& curv)
{
	ASSERT(p.GetSize()==curv.GetSize());
	ASSERT(p.GetSize()>=3);
	int N = p.GetSize();
	curv[0]=0; curv[N-1] = 0; // assign zero to first and last pnts
	for(int i=1; i<N-1;i++)
		curv[i] = (p[i-1]-2*p[i]+p[i+1]).sqr_norm();
}

void SqrDiscreteCurvature(const CPointExVector& p, CDoubleVector& curv, CPointExVector& curvVec)
{
	ASSERT(p.GetSize()==curv.GetSize());
	ASSERT(p.GetSize()==curvVec.GetSize());
	ASSERT(p.GetSize()>=3);
	int N = p.GetSize();
	curv[0]=0; curv[N-1] = 0; // assign zero to first and last pnts
	for(int i=1; i<N-1;i++)
	{
		curvVec[i] = (p[i-1]-2*p[i]+p[i+1]);
		curv[i] = curvVec[i].sqr_norm();
	}
}

void SqrDiscreteCurvature(const CPointEx3DVector& p, CDoubleVector& curv, CPointEx3DVector& curvVec)
{
	ASSERT(p.GetSize()==curv.GetSize());
	ASSERT(p.GetSize()==curvVec.GetSize());
	ASSERT(p.GetSize()>=3);
	int N = p.GetSize();
	curv[0]=0; curv[N-1] = 0; // assign zero to first and last pnts
	for(int i=1; i<N-1;i++)
	{
		curvVec[i] = (p[i-1]-2*p[i]+p[i+1]);
		curv[i] = curvVec[i].sqr_norm();
	}
}

#ifdef USE_IMSL
void SolveLeastSquaresIter(CDoubleMatrix& A, CDoubleVector& b, CDoubleVector& x, CDoubleVector& RES, double TOL)
{
	// set error handling parameters first
	long lzero = 0, lone = 1, lthree1 = -3, fd = 3;
	UMACH(&lthree1, &fd);
	ERSET(&lzero, &lone, &lzero);

	long NRA = A.GetRows(),
		NCA = A.GetCols(),
		LDA = NRA;
	long KBASIS = 0;
	RES.SetSize(NRA);
	DLSBRR(&NRA, &NCA, A.GetData(), &LDA, b.GetData(), &TOL, x.GetData(), RES.GetData(), &KBASIS);
	//	CALL LSBRR (NRA, NCA, A, LDA, B, TOL, X, RES, KBASIS)
}
#endif


// Global Curve Interpolation [using LU Decomposition (should be faster)]
// Input : n = num of ctrl points is n+1
//		 : p = curve degree 
//		 : Q = vector of size m+1 with the fit data
//		 : Qparams = if not NULL will be used as the fit-data parameters,
//						otherwise will be computed by the function used chord-length parameterization
// Output: P = the computed control points of the approximating curve
//		 : knot = the computed knot vector of the approximating curve
//		 : Qparams = the fit-data parameters
//		 : maxSqrResidual = the maximum sqr residual of vector B-Ax. WARNING: only solution with MKL supports that
// p = degree
// n = number of control points is n+1
// m = number of fit data is m+1 (m>n)
BOOL LeastSquaresCurvApprox(int p, int n, const CPointExVector& Q, CPointExVector& P, CDoubleVector& knot, 
							CDoubleVector& U, double& maxSqrResidual)
{
	int m = Q.GetSize()-1;
	P.ReSize(n+1);
	P[0] = Q[0];
	P[n] = Q[m];

	// compute data parameters based on chord length in the interval [0,1]
	if(U.GetData()==NULL)
		ComputeChordLengthParams(m,Q,U);
	ASSERT(U.GetSize()==m+1);
	ASSERT(U[0]==0);
	ASSERT(U[m] == 1);

	// *** Knot vector
	// compute knot dimensions
	int mu = n + p + 1;
	knot.ReSize(mu+1); // allocate knots

	// Compute Knot vector by averaging
//	CompCurveParamsAverag(p, mu, U, knot);

	// Compute knot (p+1) "left" and "right" elements
	for(int i = 0; i <= p; i++)
	{
		knot[i] = 0;
		knot[n+i+1] =  1;
	}

	// construct knots acording to (9.69) of "The NURBS Book", page 412.
	// this will result in a positive definite and well-conditioned matrix
	double d = (m+1.0)/(n-p+1.0);
	for(int j=1; j<=n-p; j++)
	{
		int i = (int)floor(j*d); ASSERT(i<=j*d);
		double a = j*d-i; ASSERT(a>=0 && a<=1);
		ASSERT(U[i-1]<U[i]);
		knot[p+j] = (1.0-a)*U[i-1] + a*U[i];
		ASSERT(knot[p+j]>0 && knot[p+j]<1);
		ASSERT(knot[p+j-1]<knot[p+j]);
	}

	// construct system (9.65) of "The NURBS Book", page 411.
	static CDoubleMatrix N; N.ReSize(m-1,n-1);
#ifndef _DEBUG
	#pragma omp parallel default(none) firstprivate(p,m,n) shared(N,knot,U)
#endif
	for(int r=1; r<=m-1; r++)
	{
		int bNotZeroRow = 0;
		for(int c=1; c<=n-1; c++)
		{
			double basis;
			OneBasisFun(p,m,knot.GetData(),c,U[r], basis);
			N(r-1,c-1) = basis;
			bNotZeroRow += (basis>0);
		}
		ASSERT(bNotZeroRow>0);
	}

	// Eq.(9.63)
	static CDoubleMatrix R; R.ReSize(m-1,2); // right-hand columnwise vectors
	for(int k=1; k<=m-1; k++)
	{
		double basis0 = 0, basisn = 0;
		OneBasisFun(p,m,knot.GetData(),0,U[k], basis0);
		OneBasisFun(p,m,knot.GetData(),n,U[k], basisn);
		CPointEx Rk = Q[k] - basis0*Q[0] - basisn*Q[m];
		R(k-1,0) = Rk.x;
		R(k-1,1) = Rk.y;
	}

// ------------------------------------------
// New First method: use EIGEN Least Squares Solver
	CDoubleMatrix x;
	if (!LIN_SYS_EIGEN::SolveLeastSquares(N, R, x)) // x will contain the solution
		return false;

	// copy solution and compute residual
	for (int r = 1; r <= n - 1; r++)
	{
		P[r] = CPointEx(x(r - 1, 0), x(r - 1, 1));
	}

	// compute max residual
	maxSqrResidual = -1;
#ifdef _DEBUG
	double resRx = 0, resRy = 0;
	double resRxMax = -1, resRyMax = -1;
	double resRxMin = 1e6, resRyMin = 1e6;
#endif
	for (int k = 1; k <= m - 1; k++)
	{
		const CPointEx pk = CurvePoint(n, p, knot.GetData(), P, U[k]);
		const double diffx = Q[k].x - pk.x;
		const double diffy = Q[k].y - pk.y;

#ifdef _DEBUG
		resRx += diffx * diffx;
		resRy += diffy * diffy;
		resRxMax = max(resRxMax, diffx * diffx);
		resRxMin = min(resRxMin, diffx * diffx);
		resRyMax = max(resRyMax, diffy * diffy);
		resRyMin = min(resRyMin, diffy * diffy);
#endif

		maxSqrResidual = max(maxSqrResidual, diffx * diffx);
		maxSqrResidual = max(maxSqrResidual, diffy * diffy);
	}

#ifdef _DEBUG
	_tprintf(_T("\nEIGEN::LeastSquaresCurvApprox:\n\tSquare Sum Residual X-dim = %.4f\n\tSquare Sum Residual Y-dim = %.4f"
		"\n\tMax SqrResidual X-dim = %.4f\n\tMin SqrResidual X-dim = %.4f"
		"\n\tMax SqrResidual Y-dim = %.4f\n\tMin SqrResidual Y-dim = %.4f"
		"\n\tMax SqrResidual = %.4f\n"),
		resRx, resRy, resRxMax, resRxMin, resRyMax, resRyMin, maxSqrResidual);
#endif


	return true;
// ------------------------------------------
}


// *****************************************************************
// ** CPlane
// ** ~~~~~~~~~~~
// ** This is plane in the Euclidean space
// **
// *****************************************************************

	// default constructor. need to use Init(..) afterwards to initialize the plane
	// and use the functions of the interface
CPlane::CPlane() {
	m_a = CPointEx3D(0, 0, 0); m_b = CPointEx3D(0, 0, 0);
	m_c = CPointEx3D(0, 0, 0); m_r = CPointEx3D(0, 0, 0);
	m_A = m_B = m_C = m_D = 0.0;
	m_u_interval[0] = -1; m_u_interval[1] = 1;
	m_v_interval[0] = -1; m_v_interval[1] = 1;
}

CPlane::CPlane(const CPlane& plane) { Init(plane); }

// construct plane from three district points
CPlane::CPlane(const CPointEx3D& p0, const CPointEx3D& p1, const CPointEx3D& p2) {
	ASSERT(p0 != p1); ASSERT(p1 != p2); ASSERT(p0 != p2);
	VERIFY(Init(p0, p1, p2));
}

void CPlane::Destroy() {
	m_a = CPointEx3D(0, 0, 0); m_b = CPointEx3D(0, 0, 0);
	m_c = CPointEx3D(0, 0, 0); m_r = CPointEx3D(0, 0, 0);
	m_A = m_B = m_C = m_D = 0.0;
	m_u_interval[0] = -1; m_u_interval[1] = 1;
	m_v_interval[0] = -1; m_v_interval[1] = 1;
}

// init plane after serialization
BOOL CPlane::Init() {
	// See, paragraph 6.7.6 of Zeid's book, pp.320-321.
	ASSERT(m_r.sqr_norm() == 1);
	ASSERT(m_b.sqr_norm() == 1);
	ASSERT(m_c.sqr_norm() == 1);

	VERIFY(InitAnalytic());

	ASSERT(IsValid());
	return TRUE;
}

// init plane from three district points
BOOL CPlane::Init(const CPointEx3D& p0, const CPointEx3D& p1, const CPointEx3D& p2) {
	ASSERT(p0 != p1); ASSERT(p1 != p2); ASSERT(p0 != p2);
	CTriangle plane = { p0, p1, p2 };

	// See, paragraph 6.7.6 of Zeid's book, pp.320-321.
	m_r = GetPolyNormal(plane); // the plane normal ||r||=1
	ASSERT(fabs(m_r.sqr_norm() - 1) < 1e-6);

	// plane's directions
	m_a = plane[0];
	m_b = plane[1] - plane[0]; m_b.normalize();
	m_c = plane[2] - plane[0]; m_c.normalize();

	m_u_interval[0] = -1; m_u_interval[1] = 1;
	m_v_interval[0] = -1; m_v_interval[1] = 1;

	VERIFY(InitAnalytic());

	ASSERT(IsValid());
	return TRUE;
}

BOOL CPlane::Init(const CPointEx3D& p0, const CPointEx3D& p1, const CPointEx3D& p2,
	double u_interval[2], double v_interval[2]) {
	if (!Init(p0, p1, p2))
		return FALSE;
	m_u_interval[0] = u_interval[0]; m_u_interval[1] = u_interval[1];
	m_v_interval[0] = v_interval[0]; m_v_interval[1] = v_interval[1];
	return TRUE;
}


// init plane from three district points and its normal
BOOL CPlane::Init(const CPointEx3D& p0, const CPointEx3D& p1, const CPointEx3D& p2, const CPointEx3D& n) {
	ASSERT(p0 != p1); ASSERT(p1 != p2); ASSERT(p0 != p2);

	// See, paragraph 6.7.6 of Zeid's book, pp.320-321.
	m_r = n;
	ASSERT(fabs(m_r.sqr_norm() - 1) < 1e-6);

	// plane's directions
	m_a = p0;
	m_b = p1 - p0; m_b.normalize();
	m_c = p2 - p0; m_c.normalize();

	ASSERT(fabs(m_b * m_c) < 1e-6);
	ASSERT(fabs(m_b * m_r) < 1e-6);
	ASSERT(fabs(m_c * m_r) < 1e-6);

	m_u_interval[0] = -1; m_u_interval[1] = 1;
	m_v_interval[0] = -1; m_v_interval[1] = 1;

	VERIFY(InitAnalytic());

	ASSERT(IsValid());
	return TRUE;
}

BOOL CPlane::Init(const CPointEx3D& p0, const CPointEx3D& p1, const CPointEx3D& p2, const CPointEx3D& n,
	double u_interval[2], double v_interval[2]) {

	if (!Init(p0, p1, p2, n))
		return FALSE;
	m_u_interval[0] = u_interval[0]; m_u_interval[1] = u_interval[1];
	m_v_interval[0] = v_interval[0]; m_v_interval[1] = v_interval[1];
	return TRUE;
}


BOOL CPlane::Init(const CPlane& plane) {
	m_a = plane.m_a; m_b = plane.m_b; m_c = plane.m_c; m_r = plane.m_r;
	m_A = plane.m_A; m_B = plane.m_B; m_C = plane.m_C; m_D = plane.m_D;
	m_u_interval[0] = plane.m_u_interval[0]; m_u_interval[1] = plane.m_u_interval[1];
	m_v_interval[0] = plane.m_v_interval[0]; m_v_interval[1] = plane.m_v_interval[1];
	ASSERT(IsValid());
	return TRUE;
}

// Init Plane from a point P0 and the normal direction N at P0.
// Write the plane in the form: (x-x0)Nx + (y-y0)Ny + (z-z0)Nz=0, and find two points
// on the plane (e.g., P1, P2). Then use Init plane to construct the Plane.
BOOL CPlane::InitPointNormal(const CPointEx3D& P0, const CPointEx3D& N)
{
	const double D = P0.x * N.x + P0.y * N.y + P0.z * N.z;
	if (fabs(N.z) > 1e-6) // compute z=f(x,y)
	{
		CPointEx3D P1(1, 0, (D - N.x) / N.z);
		CPointEx3D d1 = P1 - P0;
		CPointEx3D P2 = P0 + 10.0 * cross(d1, N);
		return Init(P0, P1, P2, N);
	}
	else
		if (fabs(N.x) > 1e-6) // compute x=f(y,z)
		{
			CPointEx3D P1((D - N.y) / N.x, 1, 0);
			CPointEx3D d1 = P1 - P0;
			CPointEx3D P2 = P0 + 10.0 * cross(d1, N);
			return Init(P0, P1, P2, N);
		}
		else
			if (fabs(N.y) > 1e-6) // compute y=f(x,z)
			{
				CPointEx3D P1(1, (D - N.x) / N.y, 0);
				CPointEx3D d1 = P1 - P0;
				CPointEx3D P2 = P0 + 10.0 * cross(d1, N);
				return Init(P0, P1, P2, N);
			}
			else
			{
				ASSERT(FALSE); // illegal N; e.g. ||n||-->0
				return FALSE;
			}
}

// Init plane by computing the best fit to selected data points
BOOL CPlane::InitBestFit(const CPointEx3DVector& data)
{
	const int NC = data.GetSize();
	ASSERT(NC > 0);
	// get mean point
	CPointEx3D x(0, 0, 0);
	for (int i = 0; i < NC; i++)
		x += data[i];
	const double fInvNC = 1.0 / double(NC);
	x.x *= fInvNC;
	x.y *= fInvNC;
	x.z *= fInvNC;

	// Compute SIGMA matrix: NOT using EQ.(4), which is wrong, but utilizing Eq.(9) and 
	// the solution based on my notes
	double a1 = 0, a2 = 0, a3 = 0, a1a2 = 0, a1a3 = 0, a2a3 = 0;
	CPointEx3D ak(0, 0, 0);
	for (int k = 0; k < NC; k++)
	{
		ak = data[k] - x;
		a1 += SQR(ak.x);
		a2 += SQR(ak.y);
		a3 += SQR(ak.z);
		a1a2 += ak.x * ak.y;
		a1a3 += ak.x * ak.z;
		a2a3 += ak.y * ak.z;
	}
	a1 *= fInvNC; a2 *= fInvNC; a3 *= fInvNC;
	a1a2 *= fInvNC; a1a3 *= fInvNC; a2a3 *= fInvNC;

	static CDoubleMatrix SIGMA(3, 3);
	SIGMA(0, 0) = a1;	SIGMA(0, 1) = a1a2;	SIGMA(0, 2) = a1a3;
	SIGMA(1, 0) = a1a2;	SIGMA(1, 1) = a2;	SIGMA(1, 2) = a2a3;
	SIGMA(2, 0) = a1a3;	SIGMA(2, 1) = a2a3;	SIGMA(2, 2) = a3;
	double eval = 0;
	CPointEx3D plane_normal;
	BOOL bSuccess = FALSE;
	// compute smallest eigenvalue/vector
	bSuccess = LIN_SYS_MKL::ComputeSmallestEigenVector(SIGMA, plane_normal);
	//bSuccess = LIN_SYS_EIGEN::ComputeSmallestEigenVector(SIGMA, plane_normal);

	if (bSuccess)
		return InitPointNormal(x, plane_normal);
	else
		return false;
}

// Init plane by computing the best fit to selected data points
BOOL CPlane::BestFitPlaneNormal(const CPointEx3DVector& data, CPointEx3D& plane_normal)
{
    const int NC = data.GetSize();
    ASSERT(NC > 0);
    // get mean point
    CPointEx3D x(0, 0, 0);
    for (int i = 0; i < NC; i++)
        x += data[i];
    const double fInvNC = 1.0 / double(NC);
    x.x *= fInvNC;
    x.y *= fInvNC;
    x.z *= fInvNC;

    // Compute SIGMA matrix: NOT using EQ.(4), which is wrong, but utilizing Eq.(9) and
    // the solution based on my notes
    double a1 = 0, a2 = 0, a3 = 0, a1a2 = 0, a1a3 = 0, a2a3 = 0;
    CPointEx3D ak(0, 0, 0);
    for (int k = 0; k < NC; k++)
    {
        ak = data[k] - x;
        a1 += SQR(ak.x);
        a2 += SQR(ak.y);
        a3 += SQR(ak.z);
        a1a2 += ak.x * ak.y;
        a1a3 += ak.x * ak.z;
        a2a3 += ak.y * ak.z;
    }
    a1 *= fInvNC; a2 *= fInvNC; a3 *= fInvNC;
    a1a2 *= fInvNC; a1a3 *= fInvNC; a2a3 *= fInvNC;

    static CDoubleMatrix SIGMA(3, 3);
    SIGMA(0, 0) = a1;	SIGMA(0, 1) = a1a2;	SIGMA(0, 2) = a1a3;
    SIGMA(1, 0) = a1a2;	SIGMA(1, 1) = a2;	SIGMA(1, 2) = a2a3;
    SIGMA(2, 0) = a1a3;	SIGMA(2, 1) = a2a3;	SIGMA(2, 2) = a3;
    double eval = 0;
    BOOL bSuccess = FALSE;
    // compute smallest eigenvalue/vector
    bSuccess = LIN_SYS_MKL::ComputeSmallestEigenVector(SIGMA, plane_normal);
    //bSuccess = LIN_SYS_EIGEN::ComputeSmallestEigenVector(SIGMA, plane_normal);
    return bSuccess;
}

// Analytic Geometry; Hliadis Book, pp.106-107
inline BOOL CPlane::InitAnalytic() {
	double x0 = m_a.x; double y0 = m_a.y; double z0 = m_a.z;
	double a1 = m_b.x; double b1 = m_b.y; double c1 = m_b.z;
	double a2 = m_c.x; double b2 = m_c.y; double c2 = m_c.z;
	m_A = b1 * c2 - b2 * c1; m_B = c1 * a2 - c2 * a1; m_C = a1 * b2 - a2 * b1; m_D = -m_A * x0 - m_B * y0 - m_C * z0;
	return (m_A != 0 || m_B != 0 || m_C != 0); // at least one of the 3 coefficents should be non-zero
}

void CPlane::Rotate(double rotX, double rotY, double rotZ) {
	m_b = rotateX(m_b, rotX);
	m_b = rotateY(m_b, rotY);
	m_b = rotateZ(m_b, rotZ);

	m_c = rotateX(m_c, rotX);
	m_c = rotateY(m_c, rotY);
	m_c = rotateZ(m_c, rotZ);

	VERIFY(InitAnalytic());
}

void CPlane::Translate(double dx, double dy, double dz) {
	m_a += CPointEx3D(dx, dy, dz);
	VERIFY(InitAnalytic());
}

void CPlane::Translate(const CPointEx3D& t) {
	m_a += t;
	VERIFY(InitAnalytic());
}

BOOL CPlane::GetInterval(double& u0, double& u1, double& v0, double& v1) const
{
	if (!IsValid()) return FALSE;
	u0 = m_u_interval[0]; u1 = m_u_interval[1];
	v0 = m_v_interval[0]; v1 = m_v_interval[1];
	return TRUE;
}

void CPlane::SetInterval(double u0, double u1, double v0, double v1) {
	m_u_interval[0] = u0;
	m_u_interval[1] = u1;
	m_v_interval[0] = v0;
	m_v_interval[1] = v1;
}

void CPlane::SetInterval(double u_interval[2], double v_interval[2]) {
	m_u_interval[0] = u_interval[0];
	m_u_interval[1] = u_interval[1];
	m_v_interval[0] = v_interval[0];
	m_v_interval[1] = v_interval[1];
}

inline CPointEx3D CPlane::GetPoint(double u, double v) const {
	ASSERT(IsValid());
	if (u < m_u_interval[0]) u = m_u_interval[0];
	if (u > m_u_interval[1]) u = m_u_interval[1];
	if (v < m_v_interval[0]) v = m_v_interval[0];
	if (v > m_v_interval[1]) v = m_v_interval[1];
	return (m_a + u * m_b + v * m_c);
}

// Project a point p on the plane.
// Return TRUE is succesful and the projection on the plane
inline BOOL CPlane::ProjectOnPlane(const CPointEx3D& p, CPointEx3D& p_on_plane) {
	ASSERT(IsValid());
	if (!IsValid())	return FALSE;
	double w = dot(m_r, (m_a - p)); // this is eq.(6.119) for <b,c,r> = orthonormal frame

	p_on_plane = p + w * m_r;
	return TRUE;
}

// Project a point p on the plane. Return BOTH p_on_plane and the corresponding parameters
inline BOOL CPlane::ProjectOnPlane(const CPointEx3D& p, CPointEx3D& p_on_plane, double& u, double& v) {
	ASSERT(IsValid());
	if (!IsValid())	return FALSE;
	CPointEx3D ap = m_a - p;
	v = dot(m_c, ap);
	u = dot(m_b, ap);
	p_on_plane = GetPoint(u, v);

	return TRUE;
}

// Project a point p on the plane. Return ONLY the corresponding parameters
inline BOOL CPlane::ProjectOnPlane(const CPointEx3D& p, double& u, double& v) {
	ASSERT(IsValid());
	if (!IsValid())	return FALSE;
	CPointEx3D ap = m_a - p;
	v = dot(m_c, ap);
	u = dot(m_b, ap);

	return TRUE;
}

// Project orthogonally a point p on the plane. Return ONLY the corresponding parameters
inline BOOL CPlane::ProjectOnPlaneOrtho(const CPointEx3D& p, double& u, double& v) {
	ASSERT(IsValid());
	if (!IsValid())	return FALSE;
	CPointEx3D c_cross_r = cross(m_c, m_r);
	CPointEx3D b_cross_r = cross(m_b, m_r);
	u = c_cross_r * (p - m_a) / (c_cross_r * m_b);
	v = b_cross_r * (p - m_a) / (b_cross_r * m_c);

	return TRUE;
}

// Project orthogonally a point p on the plane. Return proj point and the corresponding parameters
inline BOOL CPlane::ProjectOnPlaneOrtho(const CPointEx3D& p, CPointEx3D& p_on_plane, double& u, double& v) {
	ASSERT(IsValid());
	if (!IsValid())	return FALSE;
	CPointEx3D c_cross_r = cross(m_c, m_r);
	CPointEx3D b_cross_r = cross(m_b, m_r);
	u = c_cross_r * (p - m_a) / (c_cross_r * m_b);
	v = b_cross_r * (p - m_a) / (b_cross_r * m_c);
	p_on_plane = GetPoint(u, v);

	return TRUE;
}

// Project orthogonally a point p on the plane. Return proj point only
inline BOOL CPlane::ProjectOnPlaneOrtho(const CPointEx3D& p, CPointEx3D& p_on_plane) {
	ASSERT(IsValid());
	if (!IsValid())	return FALSE;
	CPointEx3D c_cross_r = cross(m_c, m_r);
	CPointEx3D b_cross_r = cross(m_b, m_r);
	double u = c_cross_r * (p - m_a) / (c_cross_r * m_b);
	double v = b_cross_r * (p - m_a) / (b_cross_r * m_c);
	p_on_plane = GetPoint(u, v);

	return TRUE;
}

// Find the distance between a point and a plane
inline double CPlane::GetDistance(const CPointEx3D& p) const {
	ASSERT(IsValid());
	if (!IsValid())	return FALSE;
	double w = dot(m_r, (m_a - p)); // this is eq.(6.119) for <b,c,r> = orthonormal frame
	CPointEx3D p_on_plane = p + w * m_r;
	return distance(p, p_on_plane);
}

// Find the square distance between a point and a plane. Similar like previous.
inline double CPlane::GetSqrDistance(const CPointEx3D& p) const {
	ASSERT(IsValid());
	if (!IsValid())	return FALSE;
	double w = dot(m_r, (m_a - p)); // this is eq.(6.119) for <b,c,r> = orthonormal frame
	CPointEx3D p_on_plane = p + w * m_r;
	return sqr_distance(p, p_on_plane);
}

inline CPointEx3D CPlane::GetPlaneNormal() const {
	ASSERT(IsValid());
	if (!IsValid())	return no_point;
	ASSERT(fabs(m_r.sqr_norm() - 1) < 1e-6);
	return m_r;
}

// return 1 if p lies in the 'positive' semispace
// return -1 if p lies in the 'negative' semispace
// return 0 if p lies on the plane
// Analytic Geometry; Hliadis Book, pp.114
inline int CPlane::GetSemiSpace(const CPointEx3D& p) const {
	double r = ZERO_INTERVAL_0(m_A * p.x + m_B * p.y + m_C * p.z + m_D);
	return (r > 0 ? 1 : (r < 0) ? -1 : 0);
}

// Is this plane initialized correctly?
inline const BOOL CPlane::IsValid() const {
	return ((m_a != m_b) && (m_a != m_c) && (m_b != m_c) && (m_r.sqr_norm() > 0)
		&& (m_A != 0 || m_B != 0 || m_C != 0));
}


// Find the intersection of a line P0P1 with the plane.
// Return the point of intersection p and the parameter t of line
// if 0<=t<=1 the line segment P0P1 intersects the plane at p, otherwise the intersection lies 
// at the extension of P0P1.
BOOL CPlane::IntersectWithLine(const CPointEx3D& p0, const CPointEx3D& p1, CPointEx3D& p, double& t)
{
	CPointEx3D e = p1 - p0;
	CPointEx3D c = cross(m_b, m_c);
	double denom = c * e;
	if (fabs(denom) < 1e-9)
		return FALSE; // plane and line are parralel

	t = (c * m_a - c * p0) / denom;
	p = p0 + t * e;

	return (t >= 0 && t <= 1);
}

BOOL CPlane::IntersectWithTriangle(const CTriangle& tr, CPointEx3D& p0, CPointEx3D& p1)
{
	double t0, t1;
	return (
		(IntersectWithLine(tr[0], tr[1], p0, t0) && (IntersectWithLine(tr[1], tr[2], p1, t1) || IntersectWithLine(tr[2], tr[0], p1, t1)))
		||
		(IntersectWithLine(tr[1], tr[2], p0, t0) && (IntersectWithLine(tr[2], tr[0], p1, t1) || IntersectWithLine(tr[0], tr[1], p1, t1)))
		||
		(IntersectWithLine(tr[2], tr[0], p0, t0) && (IntersectWithLine(tr[0], tr[1], p1, t1) || IntersectWithLine(tr[1], tr[2], p1, t1)))
		);
}

