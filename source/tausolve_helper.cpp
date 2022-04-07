#include "lmlsqr.hpp"

extern "C" void dgesvd_(char* jobu, char * jobvt, int * m, int * n, double * a, int * lda, double * s, double * u, int * ldu, double * vt, int * ldvt, double * work, int * lwork, int * info);

tausolve_helper::tausolve_helper(Matrix& J, int solver)
{
    if (solver == 0)
    {
        X1 = J; // this is for cholesky
    }
    else if (solver == 1) // compute svd of J, X1 = V^T, X2 = S
    {
        char jobu = 'N';
        char jobvt = 'A';
        int m = J.n_rows;
        int n = J.n_cols;
        double * a = J.data();
        int lda = m;

        X2 = zeros(n, 1);
        double * s = X2.data();
        double * u = nullptr;
        int ldu = 1;

        X1 = zeros(n, n);
        double * vt = X1.data();
        int ldvt = n;

        int lwork = -1;
        double * work = new double[5*m];

        int info = 0;

        // query for work variables
        dgesvd_(&jobu, &jobvt, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, &info);

        lwork = work[0];
        if (lwork > 5*m)
        {
            delete[] work;
            work = new double[lwork];
        }

        // actually compute svd
        dgesvd_(&jobu, &jobvt, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, &info);

        if (info != 0)
            throw std::runtime_error("failed to compute singular value decomposition.");

        delete[] work;
    }
    else if (solver  == 2) // compute QR of J
    {
        // fill X1 with Q, fill X2 with R
    }
    work1 = 0;
}