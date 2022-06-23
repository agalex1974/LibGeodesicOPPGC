#include "pch.h"
#include <tchar.h>

#ifndef MIN
#define MIN(a,b) (a<b)?a:b
#endif
#ifndef MAX
#define MAX(a,b) (a<b)?b:a
#endif

#ifdef IMSL_V6
#define __stdcall
#endif


void SetErrorHandleIMSL()
{
	// set error handling parameters
	long lzero = 0, lone = 1, lthree1 = -3, fd = 3;
	UMACH(&lthree1, &fd);
	ERSET(&lzero, &lone, &lzero);
}

void GetStrLastErrorIMSL(std::string& strErr)
{
	// get last error
	long ITYPE, ICODE;
	ICODE = IERCD();
	long lone = 1;
	ITYPE = N1RTY(&lone);

	std::string strCode = string_format("  --  ICODE = %d", ICODE);

	switch (ITYPE) {

	case 1:
		strErr = "Error - Level = 1" + strCode; break;
	case 2:
		strErr = "Error - Level = 2" + strCode; break;
	case 4:
		strErr = "Error - Level = 4" + strCode; break;
	case 5:
		strErr = "Error - Level = 5" + strCode; break;

	case 3:
		strErr = "Warning - Level = 3" + strCode; break;
	case 6:
		strErr = "Warning - Level = 6" + strCode; break;

	default:
		strErr = "No errors or warnings" + strCode; break;
	}
}


/////////////////////////////////////////////////////////////////////////////////
CCurve2D* g_pCurve2D = NULL; // static variable used as a helper with calls to IMSL
static CPointEx3D* g_p;
static CPointEx3D* g_dir;
CallBackFuncOptim g_fncCallBackOptim = NULL;
/////////////////////////////////////////////////////////////////////////////////


CallBackFunc g_fncCallBack = NULL;
BOOL g_bPLInterpolation = TRUE;
//const double _ERRABS = MACHINE_EPSILON;
double _ERRREL = 1e-3;
const int _MAX_RES = 400; // the discretization of curves into linear pieces

/////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////
// USED by FindMinDistOfCurves
static CURVE_PARAMS3D* g_curve1 = NULL;
static CURVE_PARAMS3D* g_curve2 = NULL;
static CURVE_PARAMS2D* g_curve1_2D = NULL;
static CURVE_PARAMS2D* g_curve2_2D = NULL;
static double g_usol1 = 0, g_usol2 = 0;
static CPointEx g_dif2D;
/////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////
// USED by CSurf::Project
static double g_u = 0, g_v = 0;
static CPointEx3D g_p_on_surf;
static CPointEx3D g_dif;
static CPointEx3DMatrix g_derivs;
/////////////////////////////////////////////////////////////////////////////////

extern void __stdcall ENERGY_FOR_CURVE_DIST(long* N, double* X, double* F)
{
	ASSERT(g_curve1 != NULL);
	ASSERT(g_curve2 != NULL);

	g_usol1 = X[0], g_usol2 = X[1];
	VERIFY(g_usol1 >= 0 && g_usol1 <= 1);
	VERIFY(g_usol2 >= 0 && g_usol2 <= 1);

	// the solution for the curve1
	CPointEx3D p1 = CurvePoint(g_curve1->P.GetSize() - 1, g_curve1->p, g_curve1->U, g_curve1->P, g_usol1);
	CPointEx3D p2 = CurvePoint(g_curve2->P.GetSize() - 1, g_curve2->p, g_curve2->U, g_curve2->P, g_usol2);
	g_dif = p1 - p2;

	*F = g_dif.sqr_norm();
	if (g_fncCallBackOptim != NULL)
		g_fncCallBackOptim(*F);
}

extern void __stdcall ENERGY_GRD_FOR_CURVE_DIST(long* N, double* X, double* G)
{
	ASSERT(g_curve1 != NULL);
	ASSERT(g_curve2 != NULL);

	double usol1 = X[0], usol2 = X[1];
	VERIFY(usol1 >= 0 && usol1 <= 1);
	VERIFY(usol2 >= 0 && usol2 <= 1);

	// the solution for the surf
	if (usol1 != g_usol1 || usol2 != g_usol2)
	{
		// the solution for the 2 curves
		CPointEx3D p1 = CurvePoint(g_curve1->n, g_curve1->p, g_curve1->U, g_curve1->P, usol1);
		CPointEx3D p2 = CurvePoint(g_curve2->n, g_curve2->p, g_curve2->U, g_curve2->P, usol2);
		g_dif = p1 - p2;
	}

	CPointEx3D dusol1 = CurveFirstDeriv(g_curve1->n, g_curve1->p, g_curve1->U, g_curve1->P.GetData(), usol1);
	CPointEx3D dusol2 = CurveFirstDeriv(g_curve2->n, g_curve2->p, g_curve2->U, g_curve2->P.GetData(), usol2);

	G[0] = 2. * g_dif * dusol1;
	G[1] = -2. * g_dif * dusol2;
}

extern void __stdcall ENERGY_HESSIAN_FOR_CURVE_DIST(long* N, double* X, double* H, long* LDH)
{
	ASSERT(g_curve1 != NULL);
	ASSERT(g_curve2 != NULL);

	double usol1 = X[0], usol2 = X[1];
	VERIFY(usol1 >= 0 && usol1 <= 1);
	VERIFY(usol2 >= 0 && usol2 <= 1);

	// the solution for the 2 curves
	// the solution for the surf
	if (usol1 != g_usol1 || usol2 != g_usol2)
	{
		// the solution for the 2 curves
		CPointEx3D p1 = CurvePoint(g_curve1->n, g_curve1->p, g_curve1->U, g_curve1->P, usol1);
		CPointEx3D p2 = CurvePoint(g_curve2->n, g_curve2->p, g_curve2->U, g_curve2->P, usol2);
		g_dif = p1 - p2;
	}

	static CPointEx3DVector CK1;
	static CPointEx3DVector CK2;
	CK1.ReSize(3);	CK2.ReSize(3);
	CurveDerivs(g_curve1->n, g_curve1->p, g_curve1->U, g_curve1->P, usol1, 2, CK1);
	CurveDerivs(g_curve2->n, g_curve2->p, g_curve2->U, g_curve2->P, usol2, 2, CK2);

	CPointEx3D Cu1 = CK1[1];
	CPointEx3D Cuu1 = CK1[2];
	CPointEx3D Cu2 = CK2[1];
	CPointEx3D Cuu2 = CK2[2];

	static CDoubleMatrix HESS;
	HESS.SetData(H, 2, 2);

	HESS(0, 0) = 2.0 * (Cu1*Cu1 + g_dif * Cuu1);
	HESS(1, 1) = 2.0 * (Cu2*Cu2 - g_dif * Cuu2);
	HESS(1, 0) = -2.0 * (Cu1*Cu2);
	HESS(0, 1) = HESS(1, 0);

	HESS.FreeData();
}

// Compute the parameter values usol1 and usol2 which correspond to 
// the points of two curves with the closest distance
// Input : p1,p2 = curve degree
//		 : U1, U2 = knot vector
//		 : P1,P2 = control points
//		 : usol1, usol2 = solution parameter value
// Output: TRUE = succedded, FALSE = failed
int FindMinDistOfCurves(int p1, const double* U1, const CPointEx3DVector& P1,
	int p2, const double* U2, const CPointEx3DVector& P2,
	double& usol1, double& usol2)
{
	ASSERT(U1 != NULL);
	ASSERT(U2 != NULL);
	ASSERT(P1.GetSize() > 1);
	ASSERT(P2.GetSize() > 1);
	ASSERT(p1 >= 1); // invalid curve degree
	ASSERT(p2 >= 1); // invalid curve degree

	// store curve params
	CURVE_PARAMS3D curve1, curve2;
	curve1.n = P1.GetSize() - 1; curve1.p = p1; curve1.U = (double*)U1; curve1.P = P1;
	curve2.n = P2.GetSize() - 1; curve2.p = p2; curve2.U = (double*)U2; curve2.P = P2;
	g_curve1 = &curve1;
	g_curve2 = &curve2;


	long numVars = 2; // the usol1 and usol2
	// set error handling parameters first
	long lzero = 0, lone = 1, lthree1 = -3, fd = 3;
	UMACH(&lthree1, &fd);
	ERSET(&lzero, &lone, &lzero);

	double XSCALE[2]; XSCALE[0] = XSCALE[1] = 1.0;
	double XSOL[2]; XSOL[0] = XSOL[1] = 0.0;

	// assing vars for optim method
	double GRADTL = 1e-6;
	long MAXFCN = 1000;
	double DFPRED = 1.0;
	double FVALUE = 0.0, FSCALE = 1.0;
	long IBTYPE = 0; // I will supply all the bounds
	double XLB[2], XUB[2];
	XLB[0] = XLB[1] = 0.0000001;
	XUB[0] = XUB[1] = 0.9999999;

	long IPARAM[7];
	double RPARAM[7];
	IPARAM[0] = 0;

	double XGUESS[2];
	XGUESS[0] = usol1; XGUESS[1] = usol2;
	if (XGUESS[0] == 0) { XGUESS[0] = 0.1; }
	if (XGUESS[0] == 1) { XGUESS[0] = 0.9; }
	if (XGUESS[1] == 0) { XGUESS[1] = 0.1; }
	if (XGUESS[1] == 1) { XGUESS[1] = 0.9; }

	////////////////////////////////////////////////////////////////////////////////////
	// For testing the User Gradient, Hessian
	////////////////////////////////////////////////////////////////////////////////////
	//	TestUserGrd(XGUESS, numVars, ENERGY_FOR_CURVE_DIST, ENERGY_GRD_FOR_CURVE_DIST);
	//	TestUserHessian(XGUESS, numVars, numVars, ENERGY_GRD_FOR_CURVE_DIST, ENERGY_HESSIAN_FOR_CURVE_DIST);
	//	return 0;
	////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////

		// user supplied Hessian with simply bounds
	DBCOAH(ENERGY_FOR_CURVE_DIST, ENERGY_GRD_FOR_CURVE_DIST, ENERGY_HESSIAN_FOR_CURVE_DIST,
		(long*)&numVars, XGUESS, &IBTYPE, XLB, XUB, XSCALE,
		&FSCALE, IPARAM, RPARAM, XSOL, &FVALUE);
	// Finite difference Hessian
//	DBCONG(ENERGY_FOR_CURVE_DIST, ENERGY_GRD_FOR_CURVE_DIST, (long*) &numVars, 
//		XGUESS, &IBTYPE, XLB, XUB, XSCALE,
//		&FSCALE, IPARAM, RPARAM, XSOL, &FVALUE);

	usol1 = XSOL[0]; usol2 = XSOL[1];

	g_curve1 = NULL;
	g_curve2 = NULL;
	g_fncCallBackOptim = NULL;

	// get last error
	long ITYPE, ICODE;
	ICODE = IERCD();
	ITYPE = N1RTY(&lone);

	switch (ITYPE) {
	case 1:
	case 2:
	case 4:
	case 5:
	case 3:
	case 6:
		return ITYPE;

	default:
		return 0;
	}

	return 0;
}

extern void __stdcall ENERGY_FOR_CURVE_DIST2D(long* N, double* X, double* F)
{
	ASSERT(g_curve1_2D != NULL);
	ASSERT(g_curve2_2D != NULL);

	g_usol1 = X[0], g_usol2 = X[1];
	VERIFY(g_usol1 >= 0 && g_usol1 <= 1);
	VERIFY(g_usol2 >= 0 && g_usol2 <= 1);

	// the solution for the curve1
	CPointEx p1 = CurvePoint(g_curve1_2D->P.GetSize() - 1, g_curve1_2D->p, g_curve1_2D->U, g_curve1_2D->P, g_usol1);
	CPointEx p2 = CurvePoint(g_curve2_2D->P.GetSize() - 1, g_curve2_2D->p, g_curve2_2D->U, g_curve2_2D->P, g_usol2);
	g_dif2D = p1 - p2;

	*F = g_dif2D.sqr_norm();

	if (g_fncCallBackOptim != NULL)
		g_fncCallBackOptim(*F);
}

extern void __stdcall ENERGY_GRD_FOR_CURVE_DIST2D(long* N, double* X, double* G)
{
	ASSERT(g_curve1_2D != NULL);
	ASSERT(g_curve2_2D != NULL);

	double usol1 = X[0], usol2 = X[1];
	VERIFY(usol1 >= 0 && usol1 <= 1);
	VERIFY(usol2 >= 0 && usol2 <= 1);

	// the solution for the 2 curves
	if (usol1 != g_usol1 || usol2 != g_usol2)
	{
		// the solution for the 2 curves
		CPointEx p1 = CurvePoint(g_curve1_2D->n, g_curve1_2D->p, g_curve1_2D->U, g_curve1_2D->P, usol1);
		CPointEx p2 = CurvePoint(g_curve2_2D->n, g_curve2_2D->p, g_curve2_2D->U, g_curve2_2D->P, usol2);
		g_dif2D = p1 - p2;
	}

	CPointEx dusol1 = CurveFirstDeriv(g_curve1_2D->n, g_curve1_2D->p, g_curve1_2D->U, g_curve1_2D->P.GetData(), usol1);
	CPointEx dusol2 = CurveFirstDeriv(g_curve2_2D->n, g_curve2_2D->p, g_curve2_2D->U, g_curve2_2D->P.GetData(), usol2);

	G[0] = 2. * g_dif2D * dusol1;
	G[1] = -2. * g_dif2D * dusol2;
}




// Numerical Integration

const int N = 8;

double val1(double a, double b)
{												 // Numerical, titas math series //
	double h = (b - a) / 2.0;                        // page no 312                          //
	int n;                                      // Romberg integration             //
	double arry[N + 1], x0;                       // using Trapezoidal rule         //
	for (n = 0; n <= 2; n++) {
		x0 = a;
		arry[n] = FCURVE2D(&x0);//(cos(x0)*log(sin(x0)))/(1+sin(x0));
		a = x0 + h;
	}
	return (h / 2.0)*((arry[0] + arry[2]) + 2 * arry[1]);
}

double val2(double a, double b)
{
	double h = (b - a) / 4.0;
	int n;
	double arry[N + 1], x0;
	for (n = 0; n <= 4; n++) {
		x0 = a;
		arry[n] = FCURVE2D(&x0);//( cos(x0)*log(sin(x0) ) ) / (1+sin(x0));
		a = x0 + h;
	}

	return (h / 2.0)*((arry[0] + arry[4]) + 2 * (arry[1] + arry[2] + arry[3]));
}

double val3(double a, double b)
{
	double h = (b - a) / 8.0;
	int n;
	double arry[N + 1], x0;
	for (n = 0; n <= 8; n++) {
		x0 = a;

		arry[n] = FCURVE2D(&x0);//( cos(x0)*log(sin(x0))) / (1+sin(x0));
		a = x0 + h;
	}
	return (h / 2.0)*((arry[0] + arry[8]) + 2 * (arry[1] + arry[2] + arry[3] + arry[4] + arry[5] + arry[6] + arry[7]));
}

// this is fast but not accurate. I need to find the source book to improve the code; TODO!!
double CCurve2D::Integrate2() const
{
	ASSERT(m_p >= 1); // invalid curve degree
	ASSERT(m_n >= 1); // invalid ctrl points
	ASSERT(m_U.GetData() != NULL); // curve not initialized!!

	g_pCurve2D = (CCurve2D*)this;
	double A = 0., B = 1.0;// limit

	double pal, val;
	double I1, I2, I3;

	I1 = val1(A, B);
	I2 = val3(A, B);
	I3 = val3(A, B);

	val = I2 + (1.0 / 3.0)*(I2 - I1);
	pal = I3 + (1.0 / 3.0)*(I3 - I2);

	double RESULT = 0;

	if (fabs(val - pal) <= 0.001)
	{
		RESULT = pal;
	}
	else
		RESULT = pal - (1.0 / 3.0)*(pal - val);

	g_pCurve2D = NULL;
	return RESULT;
}

double CCurve2D::Integrate3(int numOfSamples) const
{
	double d = 0;
	int nPoints = numOfSamples > 0 ? numOfSamples : 2000;
	double dt = 1. / double(nPoints - 1);
	double t = 0;
	CPointEx pi = GetPoint(t);
	for (int i = 1; i <= nPoints; i++)
	{
		t += dt;
		t = min(t, 1.0);
		CPointEx pii = GetPoint(t);
		d += DISTANCE_EX(pi, pii);
		pi = pii;
	}
	ASSERT(t == 1);
	return d;
}
