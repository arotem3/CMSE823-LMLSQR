#include "lmlsqr.hpp"

extern "C" void dormqr_(char * side, char * trans, int * m, int * n, int * k, double * a, int * lda, double * tau, double * c, int * ldc, double * work, int * lwork, int * info);
extern "C" void dtrtrs_(char * uplo, char * trans, char * diag, int * n, int * nrhs, double * a, int * lda, double * b, int * ldb, int * info);

Matrix tausolve_qr(tausolve_helper& helper, double tau, Matrix& b)
{
    if (helper.work1 != 1)
    {
        char side = 'L';
        char trans = 'T';
        int m = helper.X1.n_rows;
        int n = 1;
        int k = helper.X1.n_cols;
        double * a = helper.X1.data();
        int lda = m;
        double * tau_qr = helper.X2.data();
        double * c = b.data();
        int ldc = m;
        int lwork = -1;
        double * work = new double[n];
        int info = 0;

        // query workspace
        dormqr_(&side, &trans, &m, &n, &k, a, &lda, tau_qr, c, &ldc, work, &lwork, &info);

        lwork = int(work[0]);
        if (lwork > n)
        {
            delete[] work;
            work = new double[lwork];
        }

        dormqr_(&side, &trans, &m, &n, &k, a, &lda, tau_qr, c, &ldc, work, &lwork, &info);

        delete[] work;
        helper.work1 = 1;
    }

    int n = helper.X1.n_cols;
    Matrix b1 = b;
    Matrix b2 = zeros(n,1);
    Matrix G = tau * eye(n);
    Matrix R = zeros(n, n);

    for (int i=0; i < n; ++i)
    {
        for (int j=i; j < n; ++j)
        {
            R(i, j) = helper.X1(i, j);
        }
    }

    // givens rotation
    for (int i=n-1; i >= 0; --i)
    {
        for (int j=i; j < n; ++j)
        {
            double h = std::hypot(G(i,j), R(j,j));
            double c = R(j, j) / h;
            double s = -G(i, j) / h;

            double w1, w2;
            for (int k=0; k < n; ++k)
            {
                w1 = c * R(j, k) - s * G(i, k);
                w2 = s * R(j, k) + c * G(i, k);
                R(j, k) = w1;
                G(i, k) = w2;
            }

            w1 = c * b1[j] - s * b2[i];
            w2 = s * b1[j] + c * b2[i];
            b1[j] = w1;
            b2[i] = w2;
        }
    }
    
    // solve x = R \ b1
    char uplo = 'U';
    char trans = 'N';
    char diag = 'N';
    int nrhs = 1;
    double * a = R.data();
    int lda = n;
    double * rhs = b1.data();
    int ldb = n;
    int info = 0;

    dtrtrs_(&uplo, &trans, &diag, &n, &nrhs, a, &lda, rhs, &ldb, &info);

    if (info != 0)
        throw std::runtime_error("Failed to solve R*x = b, R may be singular?");
    
    return b1;
}