#include "lmlsqr.hpp"

extern "C" void dgesvd_(char* jobu, char * jobvt, int * m, int * n, double * a, int * lda, double * s, double * u, int * ldu, double * vt, int * ldvt, double * work, int * lwork, int * info);
extern "C" void dgeqrf_(int * m, int * n, double * a, int * lda, double * tau, double * work, int * lwork, int * info);

tausolve_helper::tausolve_helper(Matrix& J, int solver) : _J(J)
{
    if (solver == 0)
    {
        X1 = J; // this is for cholesky
    }
    else if (solver == 1) // compute svd of J, X1 = V^T, X2 = S
    {
        char jobu = 'N';
        char jobvt = 'O';
        int m = J.n_rows;
        int n = J.n_cols;
        X1 = J;
        double * a = X1.data();
        int lda = m;

        X2 = zeros(n, 1);
        double * s = X2.data();
        double * u = nullptr;
        double * vt = nullptr;
        int ldu = 1;
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

        Matrix Vt = zeros(n,n);
        for (int i=0; i < n; ++i)
            for (int j=0; j < n; ++j)
                Vt(i,j) = X1(i,j); // extract the top n,n submatrix of X2

        X1 = std::move(Vt);

        delete[] work;
    }
    else if (solver  == 2) // compute QR of J, X1 is Q and R compacted, X2 is tau from lapack dgeqrf
    {
        X1 = J;
        int m = J.n_rows;
        int n = J.n_cols;
        double * a = X1.data();
        int lda = m;
        X2 = zeros(n, 1);
        double * tau_qr = X2.data();
        int lwork = -1;
        double * work = new double[n];
        int info = 0;

        // workspace query
        dgeqrf_(&m, &n, a, &lda, tau_qr, work, &lwork, &info);

        lwork = work[0];
        if (lwork > n)
        {
            delete[] work;
            work = new double[lwork];
        }

        // actual decomp
        dgeqrf_(&m, &n, a, &lda, tau_qr, work, &lwork, &info);

        if (info != 0)
            throw std::runtime_error("failed to compute QR decomposition.");

        delete[] work;
    }
    work1 = 0;
}