#ifndef CMATHUTILITIES_H
#define CMATHUTILITIES_H

#include "Arrays.h"

class CMathUtilities
{
public:
    CMathUtilities();
    static void PreconditionedConjugateGradient(CMatrix<double> A, CMatrix<double> P,
                                                const CVector<double>& xo,
                                                const CVector<double>& b,
                                                CVector<double>& x);
    static void ConjugateGradientMKL(double* A, double* x, double* b, int N);
};

#endif // CMATHUTILITIES_H
