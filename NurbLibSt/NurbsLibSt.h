#pragma once

#ifndef DECLSPEC
# if defined (BUILD_DLL_EXPORTS)
#     define DECLSPEC __declspec( dllexport ) 
# else
#     define DECLSPEC __declspec( dllimport ) 
# endif
#endif

#include "MCLSEXST.H"


/////////////////////////////////////////////////////////////////////////////////
// storage structures
typedef CMatrix<CPointEx3D> CPointEx3DMatrix;
typedef CMatrix<CPointEx> CPointExMatrix;
typedef CVector<CPointEx3D> CPointEx3DVector;
typedef CVector<CPointEx> CPointExVector;
typedef CVector<CPointEx3D> CPointEx3DVector;
typedef double ddouble[4][4];
/////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////
// USED for storing the curve parameters
typedef struct _CURVE_PARAMS2D {
	int n; //num of ctrl points-1
	int p; // curve degree
	double* U; // knot vector
	CPointExVector P; // control points
} CURVE_PARAMS2D;
/////////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////////////


// for use with CPointNDVector
#define SET_DIM(x,n) {for(int i=0;i<(x).GetSize();i++) (x)[i].SetDim(n);}

// for use with CPointNDMatrix
#define SET_DIM2(x,n) {for(int i=0;i<(x).GetRows();i++)	\
	for(int j=0;j<(x).GetCols();j++) (x)(i,j).SetDim(n);}

// for use with CPointND[]
#define SET_DIM3(x,n,nDim) {for(int i=0;i<n;i++) (x)[i].SetDim(nDim);}

///////////////////////////////////////////////////////////////////////////////////
// HERMITE Polynomials
inline double HERM_F1(double u)
{
	double uu = u * u;
	return 2.*uu*u - 3.*uu + 1.;
}

inline double HERM_F2(double u)
{
	double uu = u * u;
	return -2.*uu*u + 3.*uu;
}

inline double HERM_F3(double u)
{
	double uu = u * u;
	return uu * u - 2.*uu + u;
}

inline double HERM_F4(double u)
{
	double uu = u * u;
	return uu * u - uu;
}



/////////////////////////////////////////////////////////////////////////////
// DLLs export functions are here

/////////////////////////////////////////////////////////////////////////////
//*** Declarations of export functions.
//*** The NURBS API
/////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////
//*** NURBS SPECIFIC FUNCTIONS ///////////////////////////////////////////////////////
//*** Used to Evaluate a B-spline basis and their derivatives at a given control point

extern "C" {
	// NURBS computes the basis of a B-spline at the i control point
	// This function computes Ni,p recursivelly. So, for p > 4 it is
	// not efficient and NURBS_K should be used instead.
	// Input : i = the control points
	//		 : p = curve degree
	//		 : U = knot vector
	//		 : u = parameter value
	// Output: the value of the B-spline basis
	DECLSPEC double NURBS(int i, int p, double u, double* U);

	// NURBS_K computes the basis of a B-spline at the k control point
	// Input : k = the control points
	//		 : n = number of control points
	//		 : p = curve degree
	//		 : U = knot vector
	//		 : u = parameter value
	//		 : N = should have enought storage for the basis functions N[0],...,N[p]
	// Output: the value of the B-spline basis
	DECLSPEC double NURBS_K(int k, int n, int p, const double* U, double u, double* N);


	// NURBS_K computes the derivatives of the basis of a B-spline 
	// at the K control point.
	// Input : K = the control points
	//		 : n = number of control points
	//		 : p = curve degree
	//		 : U = knot vector
	//		 : u = parameter value
	//		 : d = required derivative's degree
	//		 : P = control points
	// Output: DERIV = the value of the dth derivative of the basis function at the K control point
	//		 : CK[] = the value of the derivative of B-spline at the K control point
	DECLSPEC void NURBS_K_DERIV(int K, int n, int p, const double* U, double u, int d,
		const CPointExVector& P,
		double& DERIV, CPointExVector& CK);
	// Output: DERIV = the value of the derivatives (0,1,2,..d) of the basis function at the K control point
	DECLSPEC void NURBS_K_DERIVS(int K, int n, int p, const double* U, double u, int d,
		const CPointExVector& P,
		CDoubleVector& DERIV, CPointExVector& CK);
}

/////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////
extern "C++" {

	// FindSpan returns the span index where u belongs using binary search.
	// Input : n = m - p - 1; m+1 = number of knots
	//		 : p = current degree of basis functions
	//		 : u = the parameter value
	//		 : U = the knot vector
	//		 : gspan = a guess span (used when calling requersively FindSpan)
	// Return: the knot span index
	DECLSPEC int FindSpan(int n, int p, double u, const double* U, int gspan = -1);

	// BasisFuns computes all the nonvanishing basis functions and stores
	// them in the array N[0], ... , N[p].
	// Input : i = Span index
	//		 : u = parameter
	//		 : p = basis function degree
	//		 : U = knot vector
	// Output: N = the basis functions N[0],...,N[p]
	DECLSPEC void BasisFuns(int i, double u, int p, const double* U, double* N);

	// OneBasisFun computes only the Ni,p(u) basis function
	// Input : i = Span index
	//		 : u = parameter
	//		 : p = basis function degree
	//		 : U = knot vector
	//		 : m = the high index of U (m+1 knots)
	// Output: Nip = the basis function Ni,p(u)
	DECLSPEC void OneBasisFun(int p, int m, double* U, int i, double u, double& Nip);

	// CurvePoint computes a point of the actual B-spline which corresponds
	// to the parameter value u.
	// Input : n = number of control points
	//		 : p = curve degree
	//		 : U = knot vector
	//		 : P = control points
	//		 : u = parameter value
	// Output: The corresponding curve point
	DECLSPEC CPointEx CurvePoint(int n, int p, const double* U, const CPointEx* P, double u, double* N);
	DECLSPEC CPointEx CurvePoint(int n, int p, const double* U, const CPointExVector& P, double u);
	DECLSPEC CPointEx3D CurvePoint(int n, int p, const double* U, const CPointEx3D* P, double u, double* N);
	DECLSPEC CPointEx3D CurvePoint(int n, int p, const double* U, const CPointEx3DVector& P, double u);

	// DersBasisFuns computes nonzero basis functions and their
	// derivates. The derivates are stored in a two dimensional array
	// of double values.
	// Input : i = Span index
	//		 : u = parameter
	//		 : p = basis function degree
	//		 : n = derivate degree; n <= p
	//		 : U = knot vector
	// Output: ders = two dimensional array of the computed derivates
	DECLSPEC void DersBasisFuns(int i, double u, int p, int n, const double* U, ddouble ders);
	DECLSPEC void DersBasisFunsGen(int i, double u, int p, int n, const double* U, CDoubleMatrix& ders); // general implementation for p > 3

	// CurveDerivs computes the dth derivative of a B-spline with respect
	// to the parameter value u.
	// Input : n = number of control points
	//		 : p = curve degree
	//		 : U = knot vector
	//		 : P = control points
	//		 : u = parameter value
	//		 : d = required derivative's degree
	// Output: C[k] is the kth derivative; 0 <= k <= d
	DECLSPEC CPointEx* CurveDerivs(int n, int p, const double* U, const CPointEx* P, double u, int d);
	DECLSPEC void CurveDerivs(int n, int p, const double* U, const CPointExVector& P, double u, int d, CPointExVector& CK);
	DECLSPEC void CurveDerivs(int n, int p, const double* U, const CPointEx* P, double u, int d, CPointEx* CK);
	DECLSPEC CPointEx3D* CurveDerivs(int n, int p, const double* U, const CPointEx3D* P, double u, int d);
	DECLSPEC void CurveDerivs(int n, int p, const double* U, const CPointEx3DVector& P, double u, int d, CPointEx3DVector& CK);
	DECLSPEC void CurveDerivs(int n, int p, const double* U, const CPointEx3D* P, double u, int d, CPointEx3D* CK);

	// CurveFirstDeriv computes the 1st derivative of a B-spline with respect
	// to the parameter value u.
	// Input : n = number of control points
	//		 : p = curve degree
	//		 : U = knot vector
	//		 : P = control points
	//		 : u = parameter value
	// Output: C is the 1st derivative
	DECLSPEC CPointEx CurveFirstDeriv(int n, int p, const double* U, const CPointEx* P, double u);
	DECLSPEC CPointEx3D CurveFirstDeriv(int n, int p, const double* U, const CPointEx3D* P, double u);
	DECLSPEC CPointEx CurveFirstDerivGen(int n, int p, const double* U, const CPointEx* P, double u); // for p > 3

	// RefineKnotVectCurve inserts the elements of X into the knot vector
	// and returns the new knot vector and the new control points of the curve.
	// Input : n = number of control points;
	//		 : p = curve degree
	//		 : U = knot vector
	//		 : P = control polygon
	//		 : X = the vector of knots to insert
	//		 : r = the dimension of X (zero based)
	// Output: Ubar = the new knot vector
	//		 : Q = the new control polygon
	DECLSPEC void RefineKnotVectCurve(int n, int p, const double* U, const CPointEx* P, const double* X, int r, double* Ubar, CPointEx* Q);
	DECLSPEC void RefineKnotVectCurve(int n, int p, const double* U, const CPointEx3D* P, const double* X, int r, double* Ubar, CPointEx3D* Q);

	// RemoveCurveKnot removes any un-necessary knots and control points of the curve.
	// Input : n = number of control points;
	//		 : p = curve degree
	//		 : U = knot vector
	//		 : P = control polygon
	//		 : u = the knot param
	//		 : r = knot index
	//		 : s = knot multiplicity
	//		 : num = num of times to remove u
	// Output: U = the new knot vector (if algorithm succeeds)
	//		 : P = the new control polygon (if algorithm succeeds)
	//		 : t = num of the actually times u was removed
	DECLSPEC void RemoveCurveKnot(int n, int p, double* U, CPointEx* P,
		double u, int r, int s, int num, int& t, double error = TOL);
	DECLSPEC void RemoveCurveKnot(int n, int p, double* U, CPointEx3D* P,
		double u, int r, int s, int num, int& t, double error = TOL);

	DECLSPEC CPointEx CurveNormal(int n, int p, const double* U, const CPointEx* P, double u);

	// **** Local interpolation using cubic B-splines ****
	DECLSPEC void FindCtrlPointsK_K1(const CPointEx& Q0, const CPointEx& Q3,
		const CPointEx& T0, const CPointEx& T3,
		CPointEx& P1, CPointEx& P2);
	DECLSPEC void FindCtrlPointsK_K1(const CPointEx3D& Q0, const CPointEx3D& Q3,
		const CPointEx3D& T0, const CPointEx3D& T3,
		CPointEx3D& P1, CPointEx3D& P2);

	DECLSPEC void FindTans(int nPoints, const CPointEx* Q, LPPOINTEX* T);
	DECLSPEC void FindTans(int nPoints, const CPointEx3D* Q, LPPOINTEX3D* T);

	DECLSPEC void FindCtrlPoints(int nPoints, const CPointEx* Q, const CPointEx* T,
		int nCtrlPoints, int m,
		LPPOINTEX* P, CDoubleVector& U);
	DECLSPEC void FindCtrlPoints(int nPoints, const CPointEx3D* Q, const CPointEx3D* T,
		int nCtrlPoints, int m,
		LPPOINTEX3D* P, CDoubleVector& U);

	// Compute the parameter values usol1 and usol2 which correspond to 
	// the points of two curves with the closest distance
	// Input : p1,p2 = curve degree
	//		 : U1, U2 = knot vector
	//		 : P1,P2 = control points
	//		 : usol1, usol2 = solution parameter value
	// Output: the result of the optimization (ITYPE)
	DECLSPEC int FindMinDistOfCurves(int p1, const double* U1, const CPointEx3DVector& P1,
		int p2, const double* U2, const CPointEx3DVector& P2,
		double& usol1, double& usol2);
	DECLSPEC int FindMinDistOfCurves(int p1, const double* U1, const CPointExVector& P1,
		int p2, const double* U2, const CPointExVector& P2,
		double& usol1, double& usol2);


	// ***********************************************************************
	// Parameterization ******************************************************
	// ***********************************************************************

	// Compute Curve parameter for a given point
	DECLSPEC BOOL CurveParam(int n, int p, const double* U, const CPointEx* P,
		double* N, const CPointEx& point, int nDataPoints, double& u, double ug = -1);
	DECLSPEC BOOL CurveParam(int n, int p, const double* U, const CPointEx3D* P,
		double* N, const CPointEx3D& point, int nDataPoints, double& u, double ug = -1);
	DECLSPEC BOOL CurveParam(int n, int p, const double* U, const CPointEx* P,
		const CPointEx& point, int nDataPoints, double& u, double ug = -1);
	DECLSPEC BOOL CurveParamGen(int n, int p, const double* U, const CPointEx* P,
		const CPointEx& point, int nDataPoints, double& u, double ug = -1);

	// ComputeKnotVectorChordLength computes a knot vector based on chord lengths
	// Input : n = m - p - 1; num of ctrl points-zero based; (m=n+p+2 = number of knots-unit based)
	//		 : p = current degree of basis functions
	//		 : P = control points
	//		 : U = the computed knot vector (Input: must have enough space for m = n+p+2 values)
	// Output: U = the computed knot vector
	DECLSPEC BOOL ComputeKnotVectorChordLength(int n, int p, const CPointEx3DVector& P, double* U);

	// Compute curve parametric values throught averaging
	// Input : m = n + p + 1; num of knots zero based;
	//		 : p = current degree of basis functions
	//		 : uj = initial parameters
	// Output: u = the computed parameters by averaging
	DECLSPEC void CompCurveParamsAverag(int p, int m, const CDoubleVector& uj, CDoubleVector& u);

	// Compute the parametric values U for global curve fitting by chord length method
	// Input : P = fit data
	// Output: U = parametric values of fit data
	// Return: TRUE = success, FALSE = length = 0!!
	DECLSPEC BOOL ComputeChordLengthParams(int n, const CPointEx3DVector& P,
		CDoubleVector& U);
	DECLSPEC BOOL ComputeChordLengthParams(int n, const CPointExVector& P,
		CDoubleVector& U);

	// ***********************************************************************
	// ***********************************************************************
	// Linear Least Squares  *************************************************
	// ***********************************************************************

	// Linear least Squares
	// Input : A = the coefficients matrix (NRAxNCA with NRA >= NCA)
	//		 : b = right-hand coefficients of system
	// Output: x = solution of the least squares
	//		 : RES = the residual vector
	DECLSPEC void SolveLeastSquares(CDoubleMatrix& A, CDoubleVector& b, CDoubleVector& x, double TOL = 1e-6);
	DECLSPEC void SolveLeastSquares(CDoubleMatrix& A, CDoubleVector& b, CDoubleVector& x, CDoubleVector& RES, double TOL = 1e-6);
	DECLSPEC void SolveLeastSquaresIter(CDoubleMatrix& A, CDoubleVector& b, CDoubleVector& x, CDoubleVector& RES, double TOL = 1e-6);
	namespace LIN_SYS_MKL {
		// Linear least Squares
		// Input : A = the coefficients matrix (NRAxNCA)
		//		 : b = right-hand matrix with the right-hand system vectors stored columnwise
		// Output: A = the QR factorization of A if NRA >= NCA, otherwise the LQ
		//		 : b = solution of the system in columnwise order
		DECLSPEC bool SolveLeastSquares(CDoubleMatrix& A, CDoubleMatrix& b);
	}

	namespace LIN_SYS_EIGEN {
		// This will compute the solution to the system (AtA)x = (At)b
		DECLSPEC bool SolveLeastSquares(const CDoubleMatrix& A, const CDoubleMatrix& b, CDoubleMatrix& x);
	}


	// ***********************************************************************
	// General Matrix Inverse ************************************************
	// ***********************************************************************

	// Compute the Inverse of a General Real Matrix
	// Input : A = NxN matrix
	// Output: A = the inverse of A 
	//		 : AINV = the inverse of A			
	DECLSPEC int MINV(CDoubleMatrix& A);
	DECLSPEC int MINV(const CDoubleMatrix& A, CDoubleMatrix& AINV);
	// ***********************************************************************


	// ***********************************************************************
	// General Matrix Mult ***************************************************
	// ***********************************************************************

	// Compute the Product C = A * B
	// Input : A = MxK matrix
	// Input : B = KxN matrix
	// Output: C = A*B
	//		 : C = MxN
	DECLSPEC void MULT(const CDoubleMatrix& A, const CDoubleMatrix& B, CDoubleMatrix& C);
	DECLSPEC void MULT(const CFloatMatrix& A, const CFloatMatrix& B, CFloatMatrix& C);
	// ***********************************************************************

	// ***********************************************************************
	// General Matrix-Vector Mult ********************************************
	// ***********************************************************************

	// Compute the Product y = A * x
	// Input : A = MxN matrix
	// Input : x = N vector
	// Output: y = A*x
	//		 : y = M
	DECLSPEC void MULT(const CDoubleMatrix& A, const CDoubleVector& x, CDoubleVector& y);
	DECLSPEC void MULT(const CDoubleMatrix& A, const CPointExVector& x, CPointExVector& y);
	DECLSPEC void MULT(const CDoubleMatrix& A, const CPointEx3DVector& x, CPointEx3DVector& y);


	// ***********************************************************************
	// General Linear System Solution ****************************************
	// ***********************************************************************

	// Linear System Solution of a General Real Matrix
	// Input : A = NxN matrix
	//		 : X = right-hand coefficients of system
	// Output: X = solution of the system
	DECLSPEC int SolveLinSys(const CDoubleMatrix& A, CDoubleVector& X);
	DECLSPEC int SolveLinSys(const CDoubleMatrix& A, CPointEx3DVector& X, BOOL bOMP = false);
	DECLSPEC int SolveLinSys(const CDoubleMatrix& A, CPointExVector& X, BOOL bOMP = false);
	// ***********************************************************************


	// ***********************************************************************
	// EIGEN VECTORS Computation *********************************************
	// ***********************************************************************
	// Input : SIGMA = real symmetric matrix 3x3
	// Output: evec = the eigenvector that corresponds to the smallest eigenvalue system
	namespace LIN_SYS_MKL {
		DECLSPEC bool ComputeSmallestEigenVector(const CDoubleMatrix& SIGMA, CPointEx3D& evec);
	}
	namespace LIN_SYS_EIGEN {
		DECLSPEC bool ComputeSmallestEigenVector(const CDoubleMatrix& SIGMA, CPointEx3D& evec);
	}


	namespace LIN_SYS_IMSL {
		// ***********************************************************************
		// General Sparse Linear System Solution *********************************
		// ***********************************************************************

		// Linear System Solution of a General Real Matrix
		// Input : a = vector of size NZ holding all the non-zero values of A
		//		 : b = right-hand coefficients of system of size N
		//		 : IROW
		//		 : JROW:  such as Airow(i), jcol(i) = a(i)
		// Output: X = solution of the system
		DECLSPEC int SolveLinSysSparse(const CDoubleVector& a, const CDoubleVector& b, const CLongVector& IROW, const CLongVector& JCOL, CDoubleVector& X);

		// ***********************************************************************
		// Linear System Solution after LU Factorization *************************
		// ***********************************************************************

		// Linear System Solution after LU Factorization 
		// Input : FAC = LU matrix coefficients
		//		 : IPVT = LU matrix coefficients for pivoting
		//		 : N = dimension of system
		//		 : sbw = semibandwidth
		//		 : NUCA = non-zero upper co-diagonals
		//		 : NLCA = non-zero lower co-diagonals
		//		 : B = right-hand coefficients of system
		// Output: B = solution of the system
		DECLSPEC int SolveLULinSys(const CDoubleMatrix& FAC, const CLongVector& IPVT,
			long N, int sbw, CDoubleVector& B);
		DECLSPEC int SolveLULinSys(const CDoubleMatrix& FAC, const CLongVector& IPVT,
			long N, long NUCA, long NLCA, CDoubleVector& B);
		DECLSPEC int SolveLULinSys(const CDoubleMatrix& FAC, const CLongVector& IPVT,
			long N, long NUCA, long NLCA, CPointEx3DVector& B);
		DECLSPEC int SolveLULinSys(const CDoubleMatrix& FAC, const CLongVector& IPVT,
			long N, int sbw, CPointEx3DVector& B);
		DECLSPEC int SolveLULinSys(const CDoubleMatrix& FAC, const CLongVector& IPVT,
			long N, int sbw, CPointExVector& B);

		// Only for general squared NxN matrices
		DECLSPEC int SolveLULinSysNxN(const CDoubleMatrix& FAC, const CLongVector& IPVT,
			long N, CPointEx3DVector& B);
		DECLSPEC int SolveLULinSysNxN(const CDoubleMatrix& FAC, const CLongVector& IPVT,
			long N, CPointExVector& B);
		// ***********************************************************************

		// ***********************************************************************
		// Linear Algebra*********************************************************
		// ***********************************************************************

		// Multiply M(mxn) * V(n*1) ->vec(mx1)
		// bOMP: TRUE --> uses parallel OMP implementation (no ckecked yet!)
		DECLSPEC void Mult(const CDoubleMatrix& matr, const CPointEx3DVector& vec, CPointEx3DVector& res);

		// Uniform Variational subdivision minimizing 2nd Differences
		// Input : f = the fixed points at i-step
		// Output: q = the computed points at i-step
		//		 : curv = the total curve
		DECLSPEC BOOL VariationalSubdivision(const CPointEx3DVector& f, CPointEx3DVector& q, CPointEx3DVector& curv);
		DECLSPEC BOOL VariationalSubdivision(const CPointEx3DVector& f, CPointEx3DVector& curv, int nIters);
	}

	// ***********************************************************************
	// Global Curve Interpolation ********************************************
	// ***********************************************************************

	// Global Curve Interpolation [using LU Decomposition (should be faster)]
	// Input : n = m - p - 1; num of ctrl points-zero based;
	//		 : p = degree of basis functions
	//		 : uk = paameters of fixed data
	//		 : U = the knot vector
	//		 : P = fixed-data for interpolation
	// Output: P = the control points of the solution
	DECLSPEC void GlobalCurveInterLU(int n, int p, const CDoubleVector& uk,
		const CDoubleVector& U, CPointEx3DVector& P);
	DECLSPEC void GlobalCurveInterLU(int n, int p, const CDoubleVector& uk,
		const CDoubleVector& U, CPointExVector& P);
	DECLSPEC void GlobalCurveInter(int n, int p, const CDoubleVector& uk,
		const CDoubleVector& U, CPointEx3DVector& P);
	DECLSPEC void GlobalCurveInter(int n, int p, const CDoubleVector& uk,
		const CDoubleVector& U, CPointExVector& P);


	// Compute a Cubic Spline Interpolation to Q fixed data
	// Input : U = knot vector
	//		 : Q = the fixed interpolation data
	//		 : P[0], P[1], P[n+1], P[n+2] = the control points at the two ends
	// Output: P = the rest of the control points
	DECLSPEC void SolveTridiagonal(const CPointEx3DVector& Q, double* U, CPointEx3DVector& P);
	// ***********************************************************************


	// ***********************************************************************
	// Global Least Squares Curve Approximation ******************************
	// ***********************************************************************

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
	DECLSPEC BOOL LeastSquaresCurvApprox(int p, int n, const CPointExVector& Q, CPointExVector& P, CDoubleVector& knot,
		CDoubleVector& Qparams, double& maxSqrResidual);
	// ***********************************************************************



	// ***********************************************************************
	// Utility functions
	//
	// This function finds a point in the curve which is the
	// closest to the given point. It is used in the GetParam
	// function to give a good intitial value for performing 
	// the necessary optimization.
	DECLSPEC double GetClosestPoint(const CPointEx& point, int nDataPoints,
		int n, int p, const double* U, const CPointEx* P,
		double* N);

	
	// ***********************************************************************
	// Discrete Geometry******************************************************
	// ***********************************************************************

	// Compute discrete curvature for a piecewise linear curve
	// Input : pnts = points on curve
	// Output: curv = the computed curvature excluding the first and last points
	// Output: curvVec = the computed curvature vector excludin the first and last points
	// SqrDiscreteCurvature is for computing the square of curvature
	DECLSPEC void DiscreteCurvature(const CPointExVector& p, CDoubleVector& curv);
	DECLSPEC void DiscreteCurvature(const CPointEx3DVector& p, CDoubleVector& curv);
	DECLSPEC void DiscreteCurvature(const CPointExVector& p, CDoubleVector& curv, CPointExVector& curvVec);
	DECLSPEC void DiscreteCurvature(const CPointEx3DVector& p, CDoubleVector& curv, CPointEx3DVector& curvVec);
	DECLSPEC void SqrDiscreteCurvature(const CPointExVector& p, CDoubleVector& curv);
	DECLSPEC void SqrDiscreteCurvature(const CPointEx3DVector& p, CDoubleVector& curv);
	DECLSPEC void SqrDiscreteCurvature(const CPointExVector& p, CDoubleVector& curv, CPointEx3DVector& curvVec);
	DECLSPEC void SqrDiscreteCurvature(const CPointEx3DVector& p, CDoubleVector& curv, CPointEx3DVector& curvVec);
}
/////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////
// The classes below encapsulate some of the above functions to
// provide a C++ interface for surfaces and curves
// planar curve; WARNING the InitInter has NOT BEEN FULLY TESTED
// CCurve2D approximates m_ctrPoints in the plane
class CCurve2D {

protected:
	int m_p; // curves degree
	int m_n; // num of ctrl points
	CDoubleVector m_U; // the knot vector

public:
	CPointExVector m_ctrPoints;
	CPointExVector m_Fixed; // used only when InitInter is used

public:
	CCurve2D() { m_p = 0; m_n = 0; m_U.ReSize(0); }
	CCurve2D(const CCurve2D& curve) { Init(curve); }
	~CCurve2D() { Destroy(); }

	void Destroy() {
		m_ctrPoints.ReSize(0);
		m_p = 0; m_n = 0;
		m_U.ReSize(0);
	}

	BOOL IsValid() const {
		return (m_ctrPoints.GetData() != NULL && m_p != 0 && m_n != 0 && m_U.GetData() != NULL);
	}

	// clone curve
	virtual BOOL Init(const CCurve2D& curve);

	// initialize basic structures for curve
	// return true if input values are OK!
	virtual BOOL Init(int p, const CPointExVector& P, const double* U = NULL);

	// initialize basic structures for curve
	// return true if input values are OK!
	virtual BOOL Init(int p, int nCtrlPoints, const CPointEx* P, const double* U = NULL);

	// initialize curve using Least Squares approximation
	// p = degree (input)
	// Q  = m+1 fit data (input)
	// TOL = fit accuracy; the routine will increase control points to achieve TOL when TOL > 0 (input)
	// n = number of desired control points is n+1 (input and output) 
	//	   if TOL > 0 on output n will hold the actual number of control points with n<m
	virtual BOOL InitLSQApprox(int p, int& n, const CPointExVector& Q, double TOL = -1);

	// get a curve point at u; ASSERTS if curve not initialized!
	CPointEx GetPoint(double u) const;

	// get curve derivs up to and including k at u
	void GetCurveDerivs(double u, int k, CPointExVector& derivs) const;
	// special case; return first deriv
	CPointEx GetCurveFirstDeriv(double u) const;

	CPointEx GetNormal(double u) const;
	CPointEx GetTangent(double u) const;

	int GetDegree() const { return m_p; } // one minus order
	int GetKnotDim() const { return m_n + m_p + 1 + 1; } // unit based
	int GetCtrlDim() const { return m_n + 1; } // unit based

	double* GetKnotU() { return m_U.GetData(); }
	const double* GetKnotU() const { return m_U.GetData(); } // the knot vector of u-curves

	// calculate the length of this curve through numerical integration
	// decrease ERRREL for higher accuracy / higher execution time
	double Integrate(double ERRREL = 1e-3) const;
	double Integrate2() const; // numerical integration without the use of IMSL; fast but not accurate!!
	double Integrate3(int numOfSamples = -1) const; // piecewise linear approximation; fast and relatively accurate
	double Integrate4() const; // uses alglib; very accurate

	// Move
	void Move(const CPointEx& d)
	{
		for (int i = 0; i < m_ctrPoints.GetSize(); i++)
			m_ctrPoints[i] += d;
	}
};



// *****************************************************************
// ** CTriangle
// ** ~~~~~~~~~~~
// ** Operations with triangles
// **
// *****************************************************************

typedef CPointEx3D CTriangle[3];
typedef CPointEx CTriangle2D[3];

// This function computes the Metacenter of a triagle defined by
// the three points entered in the pointer 'points'
inline
CPointEx3D GetMetacenter(const CTriangle& points)
{
	return (points[0] + points[1] + points[2]) / 3.0;
}

// return the normal of the triangle given in tr
inline
CPointEx3D GetPolyNormal(const CTriangle& tr)
{
	CPointEx3D q = tr[1] - tr[0];
	CPointEx3D r = tr[2] - tr[1];
	CPointEx3D n;
	n.x = (q.y * r.z - q.z * r.y);
	n.y = (q.z * r.x - q.x * r.z);
	n.z = (q.x * r.y - q.y * r.x);

	double d = n.norm();
	if (d > 1e-9)
		n *= (1. / d);
	else
		n.x = n.y = n.z = 0.0;

	return n;
}

// returns the normalized normal of a planar triangle
inline
CPointEx3D GetNormal(const CTriangle& points)
{
	return GetPolyNormal(points);
}

inline
double SignedTriangleArea(const CTriangle2D& tr, BOOL bZero = FALSE)
{
	if (!bZero)
		return(tr[0].x * (tr[1].y - tr[2].y) +
			tr[1].x * (tr[2].y - tr[0].y) +
			tr[2].x * (tr[0].y - tr[1].y));
	else
		return(tr[1].x * tr[2].y - tr[2].x * tr[1].y);
}

inline
double SignedTriangleArea(const CPointEx& tr0, const CPointEx& tr1, const CPointEx& tr2, BOOL bZero = FALSE)
{
	if (!bZero)
		return(tr0.x * (tr1.y - tr2.y) +
			tr1.x * (tr2.y - tr0.y) +
			tr2.x * (tr0.y - tr1.y));
	else
		return(tr1.x * tr2.y - tr2.x * tr1.y);
}

// compute the area of a triangle in 3D or a triangle in the plane
inline
double TriangleArea(const CTriangle2D& tr)
{
	return fabs(SignedTriangleArea(tr, (tr[0] == CPointEx(0, 0))));
}

// this returns 1 if the triangle has positive orientation
// -1 if it has negative orientation and 0 is it has no area
// note: positive orientation = anti-clockwise
inline
int turn(const CTriangle2D& tr, BOOL bZero = FALSE)
{
	double s_a = SignedTriangleArea(tr, bZero);
	return s_a > 0 ? 1 : (s_a < 0 ? -1 : 0);
}

// return the signed area of a triangle given if its first apex is zero or not
// care has to be taken if the first apex is zero
inline double SignedTriangleArea(const CTriangle& tr, BOOL bZero = FALSE)
{
	if (!bZero)
		return(tr[0].x * (tr[1].y - tr[2].y) +
			tr[1].x * (tr[2].y - tr[0].y) +
			tr[2].x * (tr[0].y - tr[1].y));
	else
		return(tr[1].x * tr[2].y - tr[2].x * tr[1].y);
}

// this returns 1 if the triangle has positive orientation
// -1 if it has negative orientation and 0 is it has no area
// note: positive orientation = anti-clockwise
inline int turn(const CTriangle& tr, BOOL bZero = FALSE)
{
	double s_a = SignedTriangleArea(tr, bZero);
	return s_a > 0 ? 1 : (s_a < 0 ? -1 : 0);
}

// compute the area of a triangle in 3D or a triangle in the plane
inline double TriangleArea(const CTriangle& tr, const BOOL b3D = TRUE)
{
	if (b3D)
	{
		CPointEx3D Mij = tr[1] - tr[0];
		CPointEx3D Mik = tr[2] - tr[0];
		return (cross(Mij, Mik)).norm();
	}
	else
	{
		/*		CPointEx mij,  mik;
				mij.x = tr[1].x - tr[0].x;
				mij.y = tr[1].y - tr[0].y;
				mik.x = tr[2].x - tr[0].x;
				mik.y = tr[2].y - tr[0].y;
				return fabs( mij.x * mik.y - mik.x * mij.y );*/
		return fabs(SignedTriangleArea(tr, (tr[0] == CPointEx3D(0, 0, 0))));
	}
}

// Test if a 3D point p is inside a spatial triangle
// Two tests are performed which ensure the TRUE or FALSE
// of the question within a given tolerence.
// A. Ensure that v1,v2,v3 and p are coplanar
// B. Ensure that (v1,v2,p), (v2,v3,p), (v3,v1,p) have the same orientation
//	  with respect to the trangles normal vector
inline
BOOL IsInside(const CTriangle& points, const CPointEx3D& p)
{
	const double e = 1e-6;

	// compute plane's equation: Ax + By + Cz = D = 0;
	CPointEx3D u = points[1] - points[0];
	CPointEx3D v = points[2] - points[1];
	u.normalize();
	v.normalize();

	double a1 = u[0],
		b1 = u[1],
		c1 = u[2],
		a2 = v[0],
		b2 = v[1],
		c2 = v[2],
		x0 = points[0][0],
		y0 = points[0][1],
		z0 = points[0][2];

	double A = (b1 * c2) - (b2 * c1),
		B = (c1 * a2) - (c2 * a1),
		C = (a1 * b2) - (a2 * b1),
		D = -A * x0 - B * y0 - C * z0;

	if (fabs(A * p.x + B * p.y + C * p.z + D) > e)
		return FALSE; // not coplanar

	CPointEx3D v1 = points[0];
	CPointEx3D v2 = points[1];
	CPointEx3D v3 = points[2];

	// compute normals

	// (v1, v2, p)
	u = v2 - v1;
	v = p - v2;
	CPointEx3D n1 = cross(u, v);
	n1.normalize();

	// (v2, v3, p)
	u = v3 - v2;
	v = p - v3;
	CPointEx3D n2 = cross(u, v);
	n2.normalize();

	// (v3, v1, p)
	u = v1 - v3;
	v = p - v1;
	CPointEx3D n3 = cross(u, v);
	n3.normalize();

	// check that all normals point to the same direction

	double n12 = n1 * n2 - 1;
	double n13 = n1 * n3 - 1;
	double n23 = n2 * n3 - 1;

	if (fabs(n12) < e &&
		fabs(n13) < e &&
		fabs(n23) < e)
		return TRUE;
	else
	{
		return FALSE;
	}

	/*	if( (n1 * n2 > 0) && (n1 * n3 > 0) && (n2 * n3 > 0) )
			return TRUE;
		else
			return FALSE;*/
}

// returns the weighed normalized normal of a planar triangle at point p
inline
CPointEx3D GetNormal(const CTriangle& points, const CPointEx3D& p,
	const CPointEx3D& N0, const CPointEx3D& N1, const CPointEx3D& N2)
{
	double area = TriangleArea(points, TRUE);

	// (v1, v2, p)
	CTriangle tr3 = { points[0], points[1], p };
	double area3 = TriangleArea(tr3, TRUE);

	// (v2, v3, p)
	CTriangle tr1 = { points[1], points[2], p };
	double area1 = TriangleArea(tr1, TRUE);

	// (v3, v1, p)
	CTriangle tr2 = { points[2], points[0], p };
	double area2 = TriangleArea(tr2, TRUE);

	CPointEx3D N = (area1 / area) * N0 + (area2 / area) * N1 + (area3 / area) * N2;
	N.normalize();

	return N;

	/*	double l0 = distance(points[0], p);
		if(l0 <= 1e-6)
			return N0;

		double l1 = distance(points[1], p);
		if(l1 <= 1e-6)
			return N1;

		double l2 = distance(points[2], p);
		if(l2 <= 1e-6)
			return N2;

		CPointEx3D N = ( 1. / l0 ) * N0 + ( 1. / l1 ) * N1 + ( 1. / l2 ) * N2;
		N.normalize();

		return N;*/
}

// returns the weighed normalized normal of a planar triangle at point p
inline
void GetWeights(const CTriangle& points, const CPointEx3D& p,
	double& w0, double& w1, double& w2)
{
	double area = TriangleArea(points, TRUE);

	// (v1, v2, p)
	CTriangle tr3 = { points[0], points[1], p };
	double area3 = TriangleArea(tr3, TRUE);

	// (v2, v3, p)
	CTriangle tr1 = { points[1], points[2], p };
	double area1 = TriangleArea(tr1, TRUE);

	// (v3, v1, p)
	CTriangle tr2 = { points[2], points[0], p };
	double area2 = TriangleArea(tr2, TRUE);

	w0 = area1 / area;
	w1 = area2 / area;
	w2 = area3 / area;
}

// returns a weighted point from a triangle: w0 + w1 + w2 = 1.0;
inline
CPointEx3D GetPoint(const CTriangle& points, double w0, double w1, double w2)
{
	return	CPointEx3D(w0 * points[0] + w1 * points[1] + w2 * points[2]);
}

inline static void
Move_Triangle_To_Zero(CTriangle& T)
{
	T[1] = T[1] - T[0];
	T[2] = T[2] - T[0];
	T[0] = T[0] - T[0];
}

inline static void
Move_Triangle_To_Zero(CTriangle& T, int dim)
{
	CPointEx3D p(T[0]);
	for (int i = 0; i < dim; i++)
		T[i] = T[i] - p;
}




// *****************************************************************
// ** CPlane
// ** ~~~~~~~~~~~
// ** Implementation of plane geometry
// **
// *****************************************************************
class CPlane {

protected:
	// 1st representation: parametric P(u,v) = m_a + u*m_b + v*m_c
	CPointEx3D m_a, m_b, m_c, m_r; // the directions of the plane and its normal vector; SUCH AS ||b|| = ||c|| = ||r|| = 1
	double m_u_interval[2], m_v_interval[2]; // by default -1<=u,v<=1

	// 2nd representation: analytic Ax+By+Cz+D=0
	double m_A, m_B, m_C, m_D;

public:
	// default constructor. need to use Init(..) afterwards to initialize the plane
	// and use the functions of the interface
	CPlane();
	CPlane(const CPlane& plane);

	// construct plane from three district points
	CPlane(const CPointEx3D& p0, const CPointEx3D& p1, const CPointEx3D& p2);

	void Destroy();

	// init plane after serialization
	BOOL Init();

	// init plane from three district points
	BOOL Init(const CPointEx3D& p0, const CPointEx3D& p1, const CPointEx3D& p2);

	BOOL Init(const CPointEx3D& p0, const CPointEx3D& p1, const CPointEx3D& p2,
		double u_interval[2], double v_interval[2]);


	// init plane from three district points and its normal
	BOOL Init(const CPointEx3D& p0, const CPointEx3D& p1, const CPointEx3D& p2, const CPointEx3D& n);

	BOOL Init(const CPointEx3D& p0, const CPointEx3D& p1, const CPointEx3D& p2, const CPointEx3D& n,
		double u_interval[2], double v_interval[2]);

	BOOL Init(const CPlane& plane);

	// Init Plane from a point P0 and the normal direction N at P0.
	// Write the plane in the form: (x-x0)Nx + (y-y0)Ny + (z-z0)Nz=0, and find two points
	// on the plane (e.g., P1, P2). Then use Init plane to construct the Plane.
	BOOL InitPointNormal(const CPointEx3D& P0, const CPointEx3D& N);

	// Find the plane that fits given point data
	BOOL InitBestFit(const CPointEx3DVector& data);

	// Analytic Geometry; Hliadis Book, pp.106-107
	BOOL InitAnalytic();

	void Rotate(double rotX, double rotY, double rotZ);

	void Translate(double dx, double dy, double dz);

	void Translate(const CPointEx3D& t);

	BOOL GetInterval(double& u0, double& u1, double& v0, double& v1) const;

	void SetInterval(double u0, double u1, double v0, double v1);

	void SetInterval(double u_interval[2], double v_interval[2]);

	// return a point on plane
	inline CPointEx3D GetPoint(double u, double v) const;

	// return plane normal vector (normalized)
	DECLSPEC inline CPointEx3D GetPlaneNormal() const;

	// return the constant factor of the plane
	inline double GetOrigin() const { return -m_D; }
	CPointEx3D GetTopLeft() const { return GetPoint(m_u_interval[0], m_v_interval[1]); }
	CPointEx3D GetTopRight() const { return GetPoint(m_u_interval[1], m_v_interval[1]); }
	CPointEx3D GetBottomLeft() const { return GetPoint(m_u_interval[0], m_v_interval[0]); }
	CPointEx3D GetBottomRight() const { return GetPoint(m_u_interval[1], m_v_interval[0]); }
    static BOOL BestFitPlaneNormal(const CPointEx3DVector& data, CPointEx3D& plane_normal);
	// Project a point p on the plane.
	// Return TRUE is succesful and the projection on the plane
	inline BOOL ProjectOnPlane(const CPointEx3D& p, CPointEx3D& p_on_plane);

	// Project a point p on the plane. Return BOTH p_on_plane and the corresponding parameters
	inline BOOL ProjectOnPlane(const CPointEx3D& p, CPointEx3D& p_on_plane, double& u, double& v);

	// Project a point p on the plane. Return ONLY the corresponding parameters
	inline BOOL ProjectOnPlane(const CPointEx3D& p, double& u, double& v);

	// Project orthogonally a point p on the plane. Return ONLY the corresponding parameters
	inline BOOL ProjectOnPlaneOrtho(const CPointEx3D& p, double& u, double& v);

	// Project orthogonally a point p on the plane. Return proj point and the corresponding parameters
	inline BOOL ProjectOnPlaneOrtho(const CPointEx3D& p, CPointEx3D& p_on_plane, double& u, double& v);

	// Project orthogonally a point p on the plane. Return proj point only
	inline BOOL ProjectOnPlaneOrtho(const CPointEx3D& p, CPointEx3D& p_on_plane);

	// Find the distance between a point and a plane
	inline double GetDistance(const CPointEx3D& p) const;

	// Find the square distance between a point and a plane. Similar like previous.
	inline double GetSqrDistance(const CPointEx3D& p) const;


	// return 1 if p lies in the 'positive' semispace
	// return -1 if p lies in the 'negative' semispace
	// return 0 if p lies on the plane
	// Analytic Geometry; Hliadis Book, pp.114
	inline int GetSemiSpace(const CPointEx3D& p) const;

	// Is this plane initialized correctly?
	inline const BOOL IsValid() const;


	// Find the intersection of a line P0P1 with the plane.
	// Return the point of intersection p and the parameter t of line
	// if 0<=t<=1 the line segment P0P1 intersects the plane at p, otherwise the intersection lies 
	// at the extension of P0P1.
	BOOL IntersectWithLine(const CPointEx3D& p0, const CPointEx3D& p1, CPointEx3D& p, double& t);

	BOOL IntersectWithTriangle(const CTriangle& tr, CPointEx3D& p0, CPointEx3D& p1);

};


// *****************************************************************
// ** CBox3D
// ** ~~~~~~~~~~~
// *****************************************************************
// ** A 3D box used for bounding testing
// ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class CBox3D {

protected:
	CPointEx3D m_near, m_rear; // the most closer and most distant point of box
	CPointEx3D m_center; // the center of the box

public:
	// default constructor. need to can Init(..) afterwards to initialize the box
	// and use the functions of the interface
	CBox3D() {
		m_near = CPointEx3D(0, 0, 0); m_rear = CPointEx3D(0, 0, 0);
		m_center = CPointEx3D(0, 0, 0);
	}

	CBox3D(const CBox3D& box) {
		m_near = box.m_near; m_rear = box.m_rear;
		m_center = box.m_center; 
	}

	CBox3D(const CPointEx3D& pcenter, double dx, double dy, double dz) {
		VERIFY(Init(pcenter, dx, dy, dz));
	}

	CBox3D(const CPointEx3D& pnear, const CPointEx3D& prear) {
		VERIFY(Init(pnear, prear));
	}

	BOOL Init(const CPointEx3D& pcenter, double dx, double dy, double dz) {
		m_near.x = pcenter.x - dx / 2.; m_near.y = pcenter.y - dy / 2.;	m_near.z = pcenter.z - dz / 2.;
		m_rear.x = pcenter.x + dx / 2.; m_rear.y = pcenter.y + dy / 2.;	m_rear.z = pcenter.z + dz / 2.;
		m_center = pcenter; return TRUE;
	}

	BOOL Init(CBox3D& box) {
		m_near = box.m_near; m_rear = box.m_rear;
		m_center = box.m_center; return TRUE;
	}

	BOOL Init(const CPointEx3D& pnear, const CPointEx3D& prear) {
		m_near = pnear; m_rear = prear;
		m_center = (pnear + prear) / 2.0; return TRUE;
	}

	inline const BOOL IsInBox(const CPointEx3D& p) const {
		return ((p.x >= m_near.x) && (p.y >= m_near.y) && (p.z >= m_near.z) &&
			(p.x < m_rear.x) && (p.y < m_rear.y) && (p.z < m_rear.z));
	}
};

