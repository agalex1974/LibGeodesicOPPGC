// LIN_SY.cpp : Defines the basic Linear Systems Utilities
//
#include "pch.h"
#include "NurbsLibSt.h"

#include <assert.h>

#undef min
#undef max
#include <Eigen/Core>
#include <Spectra/SymEigsSolver.h>
#include <Eigen/Dense>

using namespace Eigen; 
using namespace Spectra;

namespace LIN_SYS_EIGEN {


	bool SolveLeastSquares(const CDoubleMatrix& A, const CDoubleMatrix& b, CDoubleMatrix& x)
	{
		// This will compute the solution to the system (AtA)x = (At)b

		// prepare matrices
		MatrixXd AA = Map<MatrixXd>(A.GetData(), A.GetRows(), A.GetCols());
		MatrixXd bb = Map<MatrixXd>(b.GetData(), b.GetRows(), b.GetCols());

		// solve the system (AtA)x = (At)b
		MatrixXd xx = (AA.transpose() * AA).ldlt().solve(AA.transpose() * bb);

		// copy solution
		x.CopyData(xx.data(), (int)xx.rows(), (int)xx.cols());

		// check for solution accuracy (extra cpu time, but good for robust development)
		const bool bSolved = (AA.transpose() * bb).isApprox(AA.transpose() * AA * xx);
		return bSolved;
	}

	bool ComputeSmallestEigenVector(const CDoubleMatrix& SIGMA, CPointEx3D& evec)
	{
		// We are going to calculate the eigenvalues of M
		const int	NR = SIGMA.GetRows(),
					NC = SIGMA.GetCols();
		Eigen::MatrixXd M(NR, NC);
		for (int i = 0; i < NR; i++)
			for (int j = 0; j < NC; j++)
				M(i, j) = SIGMA(i, j);

		// Construct matrix operation object using the wrapper class DenseSymMatProd
		DenseSymMatProd<double> op(M);

		// Construct eigen solver object, requesting the smallest eigenvalue
		SymEigsSolver< double, SMALLEST_ALGE, DenseSymMatProd<double> > eigs(&op, 1, NR);

		// Initialize and compute
		eigs.init();
		const int nconv = (int)eigs.compute();

		// Retrieve results
		Eigen::VectorXd evalues;
		if (eigs.info() == SUCCESSFUL)
			evalues = eigs.eigenvalues();
		else
			return false;

		// retrieve eigenvector
		Eigen::MatrixXd evectors = eigs.eigenvectors();

		evec = CPointEx3D(ZERO_INTERVAL(evectors(0, 0)),
			ZERO_INTERVAL(evectors(1, 0)),
			ZERO_INTERVAL(evectors(2, 0)));
		evec.normalize();
		return true;
	}
}