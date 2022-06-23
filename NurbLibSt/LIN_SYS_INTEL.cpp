// LIN_SY.cpp : Defines the basic Linear Systems Utilities
//
#include "pch.h"
#include "NurbsLibSt.h"


// $(ICPP_COMPILER19)\mkl\include
#include <mkl.h> // INTEL MKL
#include <assert.h>

namespace LIN_SYS_MKL {

	//#define DGEMM _MKL_BLAS_psc_dgemm
	//#define DGEMM cblas_dgemm

	void MULT(const CDoubleMatrix& A, const CDoubleMatrix& B, CDoubleMatrix& C)
	{
		//	dfsIntegrate1D(

		int M = A.GetRows(),
			K = B.GetRows(),
			N = C.GetCols();
		ASSERT(A.GetCols() == B.GetRows());
		ASSERT(C.GetCols() == N && C.GetRows() == M);

		double ALPHA = 1.0, BETA = 0.;
		int LDA = M,
			LDB = K,
			LDC = M;
		char TRANSA = 'N';
		char TRANSB = 'N';

		DGEMM(&TRANSA, &TRANSB, &M, &N, &K, &ALPHA, A.GetData(), &LDA, B.GetData(), &LDB, &BETA, C.GetData(), &LDC);
	}

	void MULT(const CFloatMatrix& A, const CFloatMatrix& B, CFloatMatrix& C)
	{
		int M = A.GetRows(),
			K = B.GetRows(),
			N = C.GetCols();
		ASSERT(A.GetCols() == B.GetRows());
		ASSERT(C.GetCols() == N && C.GetRows() == M);

		float ALPHA = 1.0f, BETA = 0.f;
		int LDA = M,
			LDB = K,
			LDC = M;
		char TRANSA = 'N';
		char TRANSB = 'N';

		SGEMM(&TRANSA, &TRANSB, &M, &N, &K, &ALPHA, A.GetData(), &LDA, B.GetData(), &LDB, &BETA, C.GetData(), &LDC);
	}

	void MULT(const CDoubleMatrix& A, const CDoubleVector& x, CDoubleVector& y)
	{
		int M = A.GetRows(),
			N = A.GetCols();
		ASSERT(A.GetCols() == x.GetSize());
		ASSERT(A.GetRows() == y.GetSize());

		double ALPHA = 1.0, BETA = 0.;
		int LDA = M,
			INCX = 1,
			INCY = 1;
		char TRANSA = 'N';

		DGEMV(&TRANSA, &M, &N, &ALPHA, A.GetData(), &LDA, (double*)x.GetData(), &INCX, &BETA, y.GetData(), &INCY);
	}



	DECLSPEC bool SolveLeastSquares(CDoubleMatrix& A, CDoubleMatrix& b)
	{
		int M = A.GetRows(),
			N = A.GetCols(),
			nrhs = b.GetCols();
		char TRANSA = 'N';
		int LDA = M, LDB = b.GetRows();

		int info = 0;
		// use C interface
		info = LAPACKE_dgels(LAPACK_COL_MAJOR, TRANSA, M, N, nrhs, A.GetData(), LDA, b.GetData(), LDB);
		return(info == 0);
	}

	DECLSPEC bool ComputeSmallestEigenVector(const CDoubleMatrix& SIGMA, CPointEx3D& evec)
	{
		// basic parameters; SIGMAX = 3x3; col-major storage
		const MKL_INT N = 3, NSELECT = 1, LDA = N, LDZ = NSELECT;

		/* Locals */
		/* Set il, iu to compute NSELECT smallest eigenvalues */
		const MKL_INT n = N, il = 1, iu = NSELECT, lda = LDA, ldz = LDZ;
		double abstol, vl=0., vu=0.;
		/* Local arrays */
		MKL_INT ifail[N], info, m;
		double w[N], z[LDZ * N];
		double* a = SIGMA.GetData();

		/* Executable statements */
		/* Negative abstol means using the default value */
		abstol = -1.0;
		/* Solve eigenproblem */
		info = LAPACKE_dsyevx(LAPACK_ROW_MAJOR, 'V', 'I', 'U', n, a, lda, vl, vu, il, iu, abstol, &m, w, z, ldz, ifail);

		/* Check for convergence */
		if (info > 0)
			return false;

		// copy solution and return succesful
		evec.x = z[0];
		evec.y = z[1];
		evec.z = z[2];
		return true;
	}

} // namespace

