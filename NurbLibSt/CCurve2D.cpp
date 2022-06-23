#include "pch.h"
#include "NurbsLibSt.h"
#include <tchar.h>


/////////////////////////////////////////////////////////////////////////////
BOOL CCurve2D::Init(const CCurve2D& curve)
{
	ASSERT(curve.m_p >= 1); // invalid curve degree
	ASSERT(curve.m_n >= 1); // invalid ctrl points
	ASSERT(curve.m_U.GetData() != NULL); // surf not initialized!!

	return Init(curve.m_p, curve.m_ctrPoints, curve.m_U.GetData());
}

// initialize basic structures for curve
// return true if input values are OK!
BOOL CCurve2D::Init(int p, const CPointExVector& P, const double* U)
{
	// Curves' degrees
	m_p = p; // set degree of u-curve
	m_n = P.GetSize() - 1; // set num of ctrl points

	ASSERT(m_p >= 1); // invalid curve degree
	ASSERT(m_n >= 1); // invalid ctrl points

	if (m_n <= 1 || m_p < 1)
		return FALSE;

	// *** Control points
	m_ctrPoints.ReSize(m_n + 1);
	for (int i = 0; i <= m_n; i++)
		m_ctrPoints[i] = P[i];

	// *** Knot vector
		// compute knot dimensions
	int mu = m_n + m_p + 1;
	m_U.ReSize(mu + 1); // allocate knots

	if (U != NULL)
		memcpy(m_U.GetData(), U, sizeof(double) * (mu + 1));
	else
	{
		CDoubleVector uk; uk.ReSize(m_n + 1);
		// by chord length method
		ComputeChordLengthParams(m_n, m_ctrPoints, uk);

		// Compute Knot vector by averaging
		CompCurveParamsAverag(m_p, m_n + m_p + 1, uk, m_U);
	}

	/*	// Compute an equally spaced knot vector in the intreval [0...1]
		for(int i = 0; i <= mu; i++)
		{
			if(i >= 0 && i <= m_p)
				m_U[i] = 0;
			else
			if(i > m_p && i <= (mu - m_p - 1))
				m_U[i] = double(i - m_p) / double(mu - m_p - m_p);
			else
			if(i > (mu - m_p - 1) )
				m_U[i] =  1;
		}
		*/

	return TRUE;
}

BOOL CCurve2D::Init(int p, int nCtrlPoints, const CPointEx* P, const double* U)
{
	ASSERT(P != NULL);

	// Curves' degrees
	m_p = p; // set degree of u-curve
	m_n = nCtrlPoints - 1; // set num of ctrl points

	ASSERT(m_p >= 1); // invalid curve degree
	ASSERT(m_n >= 1); // invalid ctrl points

	if (m_n <= 1 || m_p < 1)
		return FALSE;

	// *** Control points
	m_ctrPoints.ReSize(m_n + 1);
	for (int i = 0; i <= m_n; i++)
		m_ctrPoints[i] = P[i];

	// *** Knot vector
		// compute knot dimensions
	int mu = m_n + m_p + 1;
	m_U.ReSize(mu + 1); // allocate knots

	if (U != NULL)
		memcpy(m_U.GetData(), U, sizeof(double) * (mu + 1));
	else
	{
		CDoubleVector uk; uk.ReSize(m_n + 1);
		// by chord length method
		ComputeChordLengthParams(m_n, m_ctrPoints, uk);

		// Compute Knot vector by averaging
		CompCurveParamsAverag(m_p, m_n + m_p + 1, uk, m_U);
	}

	/*
		// Compute an equally spaced knot vector in the intreval [0...1]
		for(int i = 0; i <= mu; i++)
		{
			if(i >= 0 && i <= m_p)
				m_U[i] = 0;
			else
			if(i > m_p && i <= (mu - m_p - 1))
				m_U[i] = double(i - m_p) / double(mu - m_p - m_p);
			else
			if(i > (mu - m_p - 1) )
				m_U[i] =  1;
		}
		*/

	return TRUE;
}

BOOL CCurve2D::InitLSQApprox(int p, int& n, const CPointExVector& Q, double TOL)
{
	CPointExVector P;
	CDoubleVector knot, Qparams;
	double maxSqrResidual = -1;

	BOOL bStop = FALSE;
	if (n >= Q.GetSize() - 4)
	{
		n = Q.GetSize() - 4;
		bStop = TRUE;
	}

	if (!LeastSquaresCurvApprox(p, n, Q, P, knot, Qparams, maxSqrResidual) || !Init(p, P, knot.GetData()))
		return FALSE;
	double SQRTOL = TOL * TOL;

	if (TOL <= 0 || maxSqrResidual <= SQRTOL)
		return TRUE;
	else
		if (bStop)
			return maxSqrResidual <= SQRTOL;

	BOOL bCont;
	do {
		int nIncr = max(2, int(1.5 * (maxSqrResidual + .5 - SQRTOL))); // increase number of control points
		ASSERT(nIncr > 0);
		n += nIncr;
		if (n >= Q.GetSize() - 4)
		{
			n = Q.GetSize() - 4;
			bStop = TRUE;
		}
		_tprintf((LPCSTR)L"NURBSLIB:InitLSQApprox number of control points is increased by %d. n = %d\n", nIncr, n);

		// repeat LQ fit
		if (!LeastSquaresCurvApprox(p, n + 1, Q, P, knot, Qparams, maxSqrResidual) || !Init(p, P, knot.GetData()))
			return FALSE;
		bCont = maxSqrResidual > SQRTOL;

	} while (bCont && !bStop);

	return !bCont;
}


// get a curve point at u; ASSERTS if curve not initialized!
CPointEx CCurve2D::GetPoint(double u) const
{
	ASSERT(m_p >= 1); // invalid curve degree
	ASSERT(m_n >= 1); // invalid ctrl points
	ASSERT(m_U.GetData() != NULL); // surf not initialized!!
	ASSERT(0 <= u && u <= 1);

	//	CPointEx CurvePoint1(int n, int p, const double* U, const CPointExVector& P, double u, double* N);
	return CurvePoint(m_n, m_p, m_U.GetData(), m_ctrPoints, u);
}

void CCurve2D::GetCurveDerivs(double u, int k, CPointExVector& derivs) const
{
	ASSERT(m_p >= 1); // invalid curve degree
	ASSERT(m_n >= 1); // invalid ctrl points
	ASSERT(m_U.GetData() != NULL); // surf not initialized!!
	ASSERT(0 <= u && u <= 1);

	CurveDerivs(m_n, m_p, m_U.GetData(), m_ctrPoints, u, k, derivs);
}


CPointEx CCurve2D::GetNormal(double u) const
{
	ASSERT(m_p >= 1); // invalid curve degree
	ASSERT(m_n >= 1); // invalid ctrl points
	ASSERT(m_U.GetData() != NULL); // surf not initialized!!
	ASSERT(0 <= u && u <= 1);

	return CurveNormal(m_n, m_p, m_U.GetData(), m_ctrPoints.GetData(), u);
}

CPointEx CCurve2D::GetTangent(double u) const
{
	CPointEx T = CurveFirstDeriv(m_n, m_p, m_U.GetData(), m_ctrPoints.GetData(), u);
	T.normalize();
	return T;
}

CPointEx CCurve2D::GetCurveFirstDeriv(double u) const
{
	ASSERT(m_p >= 1); // invalid curve degree
	ASSERT(m_n >= 1); // invalid ctrl points
	ASSERT(m_U.GetData() != NULL); // curve not initialized!!
	ASSERT(0 <= u && u <= 1);
	return CurveFirstDeriv(m_n, m_p, m_U.GetData(), m_ctrPoints.GetData(), u);
}

