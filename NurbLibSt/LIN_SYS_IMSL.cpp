// LIN_SY.cpp : Defines the basic Linear Systems Utilities
//
#include "NurbsLibSt.h"
#include <assert.h>
#include "IMSL_LINK.h"

#ifdef IMSL_V6
#define __stdcall
#endif

#ifdef IMSL_V6
extern "C"
void	SGEMM(char*, char*, long*, long*,
	long*, float*, float*, long*, float*, long*, float*,
	float*, long*);
#else
extern "C"
void	__stdcall SGEMM(char*, unsigned int, char*, unsigned int, long*, long*,
	long*, float*, float*, long*, float*, long*, float*,
	float*, long*);
#endif


namespace LIN_SYS_IMSL {

	//#define DGEMM _MKL_BLAS_psc_dgemm
	//#define DGEMM cblas_dgemm

	void SetErrorHandleIMSL()
	{
		// set error handling parameters
		long lzero = 0, lone = 1, lthree1 = -3, fd = 3;
		UMACH(&lthree1, &fd);
		ERSET(&lzero, &lone, &lzero);
	}


	// Multiply M(mxn) * V(n*1) ->vec(mx1)
	void Mult(const CDoubleMatrix& matr, const CPointEx3DVector& vec, CPointEx3DVector& res)
	{
		int m = matr.GetRows();
		int n = matr.GetCols();
		ASSERT(vec.GetSize() == n);
		res.SetSize(m);

		for (int i = 0; i < m; i++)
		{
			CPointEx3D data = CPointEx3D(0, 0, 0);
			for (int j = 0; j < n; j++)
				data += (matr(i, j) * vec[j]);
			res[i] = data;
		}
	}

	void MULT(const CDoubleMatrix& A, const CDoubleMatrix& B, CDoubleMatrix& C)
	{
		long M = A.GetRows(),
			K = B.GetRows(),
			N = C.GetCols();
		ASSERT(A.GetCols() == B.GetRows());
		ASSERT(C.GetCols() == N && C.GetRows() == M);

		// set error handling parameters first
		SetErrorHandleIMSL();

		double ALPHA = 1.0, BETA = 0.;
		long LDA = M,
			LDB = K,
			LDC = M;
		char TRANSA = 'N';
		char TRANSB = 'N';

#ifdef IMSL_V6
		DGEMM(&TRANSA, &TRANSB, &M, &N, &K, &ALPHA, A.GetData(), &LDA, B.GetData(), &LDB, &BETA, C.GetData(), &LDC);
#else
		DGEMM(&TRANSA, sizeof(TRANSA), &TRANSB, sizeof(TRANSA), &M, &N, &K, &ALPHA, A.GetData(), &LDA, B.GetData(), &LDB, &BETA, C.GetData(), &LDC);
#endif
	}

	void MULT(const CFloatMatrix& A, const CFloatMatrix& B, CFloatMatrix& C)
	{
		long M = A.GetRows(),
			K = B.GetRows(),
			N = C.GetCols();
		ASSERT(A.GetCols() == B.GetRows());
		ASSERT(C.GetCols() == N && C.GetRows() == M);

		// set error handling parameters first
		SetErrorHandleIMSL();

		float ALPHA = 1.0f, BETA = 0.f;
		long LDA = M,
			LDB = K,
			LDC = M;
		char TRANSA = 'N';
		char TRANSB = 'N';

#ifdef IMSL_V6
		SGEMM(&TRANSA, &TRANSB, &M, &N, &K, &ALPHA, A.GetData(), &LDA, B.GetData(), &LDB, &BETA, C.GetData(), &LDC);
#else
		SGEMM(&TRANSA, sizeof(TRANSA), &TRANSB, sizeof(TRANSB), &M, &N, &K, &ALPHA, A.GetData(), &LDA, B.GetData(), &LDB, &BETA, C.GetData(), &LDC);
#endif
	}

	void MULT(const CDoubleMatrix& A, const CDoubleVector& x, CDoubleVector& y)
	{
		long M = A.GetRows(),
			N = A.GetCols();
		ASSERT(A.GetCols() == x.GetSize());
		ASSERT(A.GetRows() == y.GetSize());

		// set error handling parameters first
		SetErrorHandleIMSL();

		double ALPHA = 1.0, BETA = 0.;
		long LDA = M,
			INCX = 1,
			INCY = 1;
		char TRANSA = 'N';

		//CALL DGEMV (TRANS, M, N, DALPHA, DA, LDA, DX, INCX, DBETA,DY, INCY)

#ifdef IMSL_V6
		DGEMV(&TRANSA, &M, &N, &ALPHA, A.GetData(), &LDA, (double*)x.GetData(), &INCX, &BETA, y.GetData(), &INCY);
#else
		DGEMV(&TRANSA, sizeof(TRANSA), &M, &N, &ALPHA, A.GetData(), &LDA, (double*)x.GetData(), &INCX, &BETA, y.GetData(), &INCY);
#endif
	}

	DECLSPEC void MULT(const CDoubleMatrix& A, const CPointExVector& x, CPointExVector& y)
	{
		int N = A.GetCols();
		int M = A.GetRows();
		ASSERT(x.GetSize() == N);
		y.ReSize(M);
		static CDoubleVector X; X.ReSize(N);
		static CDoubleVector Y; Y.ReSize(M);

		for (int i = 0; i < 2; i++)
		{
			for (int j = 0; j < N; j++)
				X[j] = i == 0 ? x[j].x : x[j].y;

			LIN_SYS_IMSL::MULT(A, X, Y);

			for (int j = 0; j < M; j++)
			{
				if (i == 0)
					y[j].x = Y[j];
				else
					y[j].y = Y[j];
			}
		}
	}

	DECLSPEC void MULT(const CDoubleMatrix& A, const CPointEx3DVector& x, CPointEx3DVector& y)
	{
		int N = x.GetSize();
		y.ReSize(N);
		static CDoubleVector X; X.ReSize(N);
		static CDoubleVector Y; Y.ReSize(N);

		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < N; j++)
			{
				if (i == 0)
					X[j] = x[j].x;
				else
					if (i == 1)
						X[j] = x[j].y;
					else
						X[j] = x[j].z;
			}

			LIN_SYS_IMSL::MULT(A, X, Y);

			for (int j = 0; j < N; j++)
			{
				if (i == 0)
					y[j].x = Y[j];
				else
					if (i == 1)
						y[j].y = Y[j];
					else
						y[j].z = Y[j];
			}
		}
	}

	int LUDecompositionB_NxN(const CDoubleMatrix& A, long N, long NUCA, long NLCA,
		CDoubleMatrix& FAC, CLongVector& IPVT)
	{
		long LDA = N;
		long LDB = NUCA + 1 + NLCA;
		static CDoubleMatrix A_BANDED; A_BANDED.ReSize(LDB, N);

		// set error handling parameters first
		SetErrorHandleIMSL();

		// make A(NxN) a band matrix A_BANDED
		DCRGRB(&N, A.GetData(), &LDA, &NLCA, &NUCA, A_BANDED.GetData(), &LDB);

		// now perform LU decomposition to A_BANDED
		LDA = NUCA + NLCA + 1;
		FAC.ReSize(2 * NLCA + NUCA + 1, N);
		long LDFAC = 2 * NLCA + NUCA + 1;
		IPVT.ReSize(N);

		DLFTRB(&N, A_BANDED.GetData(), &LDA, &NLCA, &NUCA, FAC.GetData(), &LDFAC, IPVT.GetData());

		return TRUE;
	}

	int LUDecompositionNxN(const CDoubleMatrix& A,
		CDoubleMatrix& FAC, CLongVector& IPVT)
	{
		long N = A.GetRows();
		long LDA = A.GetRows();
		long LDFAC = N;
		ASSERT(A.GetCols() == LDA);

		// set error handling parameters first
		SetErrorHandleIMSL();

		// now perform LU decomposition to A_BANDED
		FAC.ReSize(N, N);
		IPVT.ReSize(N);

		DLFTRG(&N, A.GetData(), &LDA, FAC.GetData(), &LDFAC, IPVT.GetData());

		return TRUE;
	}

	int SolveLinSys(const CDoubleMatrix& A, CPointEx3DVector& X)
	{
		ASSERT(A.GetRows() == A.GetCols());
		long N = A.GetRows();
		ASSERT(N == X.GetSize());
		//	long LDA = N;
		//	long IPATH = 1;

		//	static CPointEx3DVector XX;XX.ReSize(N);
		//	static CDoubleMatrix AINV; AINV.ReSize(N,N);

			// set error handling parameters first
		SetErrorHandleIMSL();

		/*	for(int i=0; i<3;i++)
			{
				for(int j=0; j<X.GetSize(); j++)
					XX[j] = X[j][i];
				DLSLRG(&N, A.GetData(), &LDA, XX.GetData(), &IPATH, XX.GetData());
				for(j=0; j<X.GetSize(); j++)
					X[j][i] = XX[j];
			}*/

			// first compute the inverse of A ---> AINV
		/*	DLINRG(&N, A.GetData(), &LDA, AINV.GetData(), &LDA);

			// multiply to get solution
			CPointEx3D* p_bb = XX.GetData();
			for(int i = 0; i < N; i++)
			{
				(*p_bb) = CPointEx3D(0,0,0);
				for(int j = 0; j < N; j++)
					(*p_bb) += (AINV(i,j) * X[j]);
				p_bb++;
			}
			X = XX;*/

			// LU decomposition of A
		static CDoubleMatrix FAC;	static CLongVector IPVT;

		LIN_SYS_IMSL::LUDecompositionNxN(A, FAC, IPVT);
		LIN_SYS_IMSL::SolveLULinSysNxN(FAC, IPVT, N, X);

		return TRUE;
	}

	int SolveLinSys(const CDoubleMatrix& A, CPointExVector& X)
	{
		ASSERT(A.GetRows() == A.GetCols());
		long N = A.GetRows();
		ASSERT(N == X.GetSize());


		// set error handling parameters first
		SetErrorHandleIMSL();

		// LU decomposition of A
		static CDoubleMatrix FAC;	static CLongVector IPVT;

		LIN_SYS_IMSL::LUDecompositionNxN(A, FAC, IPVT);
		LIN_SYS_IMSL::SolveLULinSysNxN(FAC, IPVT, N, X);

		return TRUE;
	}

	int SolveLULinSys(const CDoubleMatrix& FAC, const CLongVector& IPVT,
		long N, long NUCA, long NLCA, CPointEx3DVector& B)
	{
		ASSERT(B.GetSize() == N);
		long LDA = NUCA + NLCA + 1;
		long LDFAC = 2 * NLCA + NUCA + 1;
		long IPATH = 1;

		// set error handling parameters first
		SetErrorHandleIMSL();

		static CDoubleVector X; X.ReSize(N);
		int i, j;

		for (i = 0; i < 3; i++)
		{
			for (j = 0; j < N; j++)
				X[j] = B[j][i];
			DLFSRB(&N, FAC.GetData(), &LDFAC, &NLCA, &NUCA, (long*)IPVT.GetData(),
				X.GetData(), &IPATH, X.GetData());
			for (j = 0; j < N; j++)
				B[j][i] = X[j];
		}

		return TRUE;
	}

	int SolveLULinSysNxN(const CDoubleMatrix& FAC, const CLongVector& IPVT,
		long N, CPointEx3DVector& B)
	{
		ASSERT(N == B.GetSize());
		ASSERT(FAC.GetRows() == N);
		long LDFAC = N;
		long IPATH = 1;
		static CDoubleVector X; X.ReSize(N);

		int i, j;

		// set error handling parameters first
		SetErrorHandleIMSL();

		for (i = 0; i < 3; i++)
		{
			for (j = 0; j < N; j++)
				X[j] = B[j][i];
			DLFSRG(&N, FAC.GetData(), &LDFAC, (long*)IPVT.GetData(), X.GetData(), &IPATH, X.GetData());
			for (j = 0; j < N; j++)
				B[j][i] = X[j];
		}

		return TRUE;
	}

	int SolveLULinSysNxN(const CDoubleMatrix& FAC, const CLongVector& IPVT,
		long N, CPointExVector& B)
	{
		ASSERT(N == B.GetSize());
		ASSERT(FAC.GetRows() == N);
		long LDFAC = N;
		long IPATH = 1;
		static CDoubleVector X; X.ReSize(N);

		int i, j;

		// set error handling parameters first
		SetErrorHandleIMSL();

		for (i = 0; i < 2; i++)
		{
			for (j = 0; j < N; j++)
				X[j] = B[j][i];
			DLFSRG(&N, FAC.GetData(), &LDFAC, (long*)IPVT.GetData(), X.GetData(), &IPATH, X.GetData());
			for (j = 0; j < N; j++)
				B[j][i] = X[j];
		}

		return TRUE;
	}


	DECLSPEC BOOL VariationalSubdivision(const CPointEx3DVector& f, CPointEx3DVector& q, CPointEx3DVector& curv)
	{
		int N = f.GetSize(), i = 0;
		ASSERT(N > 2);
		// 1. Formulate system Aq = Bf, A=(N-1)x(N-1), B = (N-1)xN, q=(N-1), f=N

		// 1. Matrix A, B
		static CDoubleMatrix A;	A.ReSize(N - 1, N - 1); A.SetVal(0);
		static CDoubleMatrix B;	B.ReSize(N - 1, N); B.SetVal(0);
		for (i = 0; i < N - 1; i++)
		{
			if (i == 0)
			{
				A(0, 0) = 5;	A(0, 1) = 1;
				B(0, 0) = 2;	B(0, 1) = 4;
			}
			else
				if (i == N - 2)
				{
					A(N - 2, N - 3) = 1; A(N - 2, N - 2) = 5;
					B(N - 2, N - 2) = 4; B(N - 2, N - 1) = 2;
				}
				else
				{
					A(i, i - 1) = 1; A(i, i) = 6; A(i, i + 1) = 1;
					B(i, i) = 4; B(i, i + 1) = 4;
				}
		}

		// 2. Compute q=Bf --> vector (N-1)
		LIN_SYS_IMSL::Mult(B, f, q);

		// 3. compute the solution
		long NUC = 1;
		long NLC = 1;

		// LU decomposition of A
		static CDoubleMatrix FAC;
		static CLongVector IPVT;

		LIN_SYS_IMSL::LUDecompositionB_NxN(A, N - 1, NUC, NLC, FAC, IPVT);
		LIN_SYS_IMSL::SolveLULinSys(FAC, IPVT, N - 1, NUC, NLC, q);

		curv.SetSize(2 * N - 1);
		for (int j = 1; j <= N; j++)
		{
			int index = 2 * j - 1;
			curv[/*2*j-1 -1*/index - 1] = f[j - 1];
			if (j != N)
				curv[/*2*j -1*/index] = q[j - 1];
		}

		return true;
	}

	DECLSPEC BOOL VariationalSubdivision(const CPointEx3DVector& f, CPointEx3DVector& curv, int nIters)
	{
		int N = f.GetSize(), i = 0;
		ASSERT(N > 2);
		if (N <= 2 || nIters <= 0)
			return false;
		CPointEx3DVector iF = f;
		int nIter = 0;
		do {
			CPointEx3DVector q;
			if (!LIN_SYS_IMSL::VariationalSubdivision(iF, q, curv))
				return false;
			iF = curv;
		} while (++nIter < nIters);
		return true;
	}
} // namespace

#undef __stdcall 
