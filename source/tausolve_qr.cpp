#include "lmlsqr.hpp"

extern "C" void dgeqrf_(int * m, int * n, double * a, int * lda, double * tau, double * work, int * lwork, int * info);
extern "C" void dormqr_(char * side, char * trans, int * m, int * n, int * k, double * a, int * lda, double * tau, double * c, int * ldc, double * work, int * lwork, int * info);
extern "C" void dtrtrs_(char * uplo, char * trans, char * diag, int * n, int * nrhs, double * a, int * lda, double * b, int * ldb, int * info);

TauSolverQR::TauSolverQR(const Matrix& J, const Matrix& b)
{
    QR = J;
    int m = J.n_rows;
    int n = J.n_cols;
    int lda = m;
    qr_tau = zeros(n, 1);
    int lwork = -1;
    double * work = new double[n];
    int info = 0;

    // workspace query
    dgeqrf_(&m, &n, QR.data(), &lda, qr_tau.data(), work, &lwork, &info);

    lwork = work[0];
    if (lwork > n)
    {
        delete[] work;
        work = new double[lwork];
    }

    // actual decomp
    dgeqrf_(&m, &n, QR.data(), &lda, qr_tau.data(), work, &lwork, &info);

    if (info != 0)
        throw std::runtime_error("failed to compute QR decomposition.");

    // apply Q.t() * b
    char side = 'L';
    char trans = 'T';
    n = 1;
    int k = QR.n_cols;
    int ldc = m;
    int nwork = lwork;
    lwork = -1;
    info = 0;

    Matrix bb = b;

    // query workspace
    dormqr_(&side, &trans, &m, &n, &k, QR.data(), &lda, qr_tau.data(), bb.data(), &ldc, work, &lwork, &info);

    lwork = int(work[0]);
    if (lwork > nwork)
    {
        delete[] work;
        work = new double[lwork];
    }

    dormqr_(&side, &trans, &m, &n, &k, QR.data(), &lda, qr_tau.data(), bb.data(), &ldc, work, &lwork, &info);

    delete[] work;

    // resize b to correct size.
    _b = zeros(QR.n_cols,1);
    std::copy_n(bb.data(), QR.n_cols, _b.data());

    // fill initialize R, G
    n = QR.n_cols;
    R = zeros(n, n);
    G = zeros(n, n);
    b2 = zeros(n, 1);
}

Matrix TauSolverQR::operator()(double tau)
{
    int n = QR.n_cols;
    b1 = _b;
    b2.fill(0.0);
    G.fill(0.0);
    G.diag().fill(std::sqrt(std::abs(tau)));

    for (int i=0; i < n; ++i)
        for (int j=i; j < n; ++j)
            R(i, j) = QR(i, j);

    // givens rotation
    for (int i=n-1; i >= 0; --i)
    {
        for (int j=i; j < n; ++j)
        {
            double h = std::hypot(G(i,j), R(j,j));
            double c = R(j, j) / h;
            double s = -G(i, j) / h;

            double w1, w2;
            for (int k=i; k < n; ++k)
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
    int lda = n;
    int ldb = n;
    int info = 0;

    dtrtrs_(&uplo, &trans, &diag, &n, &nrhs, R.data(), &lda, b1.data(), &ldb, &info);

    if (info != 0)
        throw std::runtime_error("Failed to solve R*x = b, R may be singular?");
    
    return b1;
}
