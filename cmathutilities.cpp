#include "pch.h"
#include <mkl.h>
#include "cmathutilities.h"
#include "Arrays.h"
#include <math.h>


static CDoubleVector Add(const CDoubleVector& v1, const CDoubleVector& v2){
    CDoubleVector v(v1.GetSize());
    for (int i = 0; i < v.GetSize(); i++){
        v[i] = v1[i] + v2[i];
    }
    return v;
}

static CDoubleVector Subtract(const CDoubleVector& v1, const CDoubleVector& v2){
    CDoubleVector v(v1.GetSize());
    for (int i = 0; i < v.GetSize(); i++){
        v[i] = v1[i] - v2[i];
    }
    return v;
}

static double Dot(const CDoubleVector& v1, const CDoubleVector& v2){
    double sum = 0.0;
    for (int i = 0; i < v1.GetSize(); i++){
        sum += v1[i] * v2[i];
    }
    return sum;
}

static CDoubleVector Mult(double alpha, const CDoubleVector& vin){
    CDoubleVector v(vin.GetSize());
    for (int i = 0; i < v.GetSize(); i++){
        v[i] = alpha * vin[i];
    }
    return v;
}

static double Norm(const CDoubleVector& v){
    double sum = 0.0;
    for (int i = 0; i < v.GetSize(); i++){
        sum += v[i] * v[i];
    }
    return sqrt(sum);
}

CMathUtilities::CMathUtilities()
{

}

void CMathUtilities::PreconditionedConjugateGradient(CMatrix<double> A, CMatrix<double> P,
                                            const CVector<double>& xo,
                                            const CVector<double>& b,
                                            CVector<double>& x){
    x = xo;
    auto r = Subtract(b, Mult(A, xo));
    auto z = Mult(P, r);
    auto p = z;
    while (true){
        auto Apk = Mult(A, p);
        double ak = Dot(r, z) / Dot(p, Apk);
        x = Add(x, Mult(ak, p));
        auto rk = r;
        auto zk = z;
        r = Subtract(r, Mult(ak, Apk));
        //qDebug() << "norm:" << Norm(r);
        if (Norm(r) < 1e-9) break;
        z = Mult(P, r);
        double bk = Dot(r, z) / Dot(rk, zk);
        p = Add(z, Mult(bk, p));
    }
}

void CMathUtilities::ConjugateGradientMKL(double* A, double* x, double* b, int N){
   double* r = (double*)mkl_malloc(N * sizeof(double), 64);
   cblas_dcopy(N, b, 1, r, 1);
   cblas_dgemv(CblasRowMajor, CblasNoTrans, N, N, -1.0, A, N, x, 1, 1.0, r, 1);
   double* p = (double*)mkl_malloc(N * sizeof(double), 64);
   cblas_dcopy(N, r, 1, p, 1);
   double* Ap = (double*)mkl_malloc(N * sizeof(double), 64);
   const int max_iterations = 1000;
   int k = 0;
   while (k < max_iterations){
       cblas_dgemv(CblasRowMajor, CblasNoTrans, N, N, 1.0, A, N, p, 1, 0.0, Ap, 1);
       double norm2rp = cblas_ddot(N, r, 1, r, 1);
       double ak = norm2rp / cblas_ddot(N, p, 1, Ap, 1);
       cblas_daxpy(N, ak, p, 1, x, 1);
       cblas_daxpy(N, -ak, Ap, 1, r, 1);
       double norm2r = cblas_ddot(N, r, 1, r, 1);
       if (sqrt(norm2r) < 1e-8) break;
       double bk = norm2r / norm2rp;
       cblas_daxpy(N, 1.0/bk, r, 1, p, 1);
       cblas_dscal(N, bk, p, 1);
       k++;
   }
   //if (k == max_iterations) qDebug() << "CG did not converge!";
   mkl_free(r);
   mkl_free(p);
   mkl_free(Ap);
}