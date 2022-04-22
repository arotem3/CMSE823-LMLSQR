#include "lmlsqr.hpp"

extern "C" void dgeqrf_(int * m, int * n, double * a, int * lda, double * tau, double * work, int * lwork, int * info);
extern "C" void dormqr_(char * side, char * trans, int * m, int * n, int * k, double * a, int * lda, double * tau, double * c, int * ldc, double * work, int * lwork, int * info);
extern "C" void dtrtrs_(char * uplo, char * trans, char * diag, int * n, int * nrhs, double * a, int * lda, double * b, int * ldb, int * info);
extern "C" void drot_(int * n, double * x, int * incx, double * y, int * incy, double * c, double * s);
extern "C" void drotg_(double * a, double * b, double * c, double * s);
extern "C" void dgels_(char * trans, int * m, int * n, int * nrhs, double * a, int * lda, double * b, int * ldb, double * work, int * lwork, int * info);

TauSolverQR::TauSolverQR(const Matrix& J, const Matrix& b)
{
    Matrix QR = J;
    int m = QR.n_rows;
    int n = QR.n_cols;
    _n = QR.n_cols;
    int lda = m;
    Matrix qr_tau(n, 1);
    int lwork = -1;
    double wt;
    int info = 0;

    // workspace query
    dgeqrf_(&m, &n, QR.data(), &lda, qr_tau.data(), &wt, &lwork, &info);

    lwork = wt;
    Matrix work(lwork, 1);

    // actual decomp
    dgeqrf_(&m, &n, QR.data(), &lda, qr_tau.data(), work.data(), &lwork, &info);

    // apply Q.t() * b
    char side = 'L';
    char trans = 'T';
    n = 1;
    int k = QR.n_cols;
    int ldc = m;
    lwork = -1;
    info = 0;

    Matrix bb = b;

    // // query workspace
    dormqr_(&side, &trans, &m, &n, &k, QR.data(), &lda, qr_tau.data(), bb.data(), &ldc, work.data(), &lwork, &info);

    lwork = int(work[0]);
    if (lwork > work.size())
        work = zeros(lwork, 1);

    dormqr_(&side, &trans, &m, &n, &k, QR.data(), &lda, qr_tau.data(), bb.data(), &ldc, work.data(), &lwork, &info);

    n = J.n_cols;
    R = zeros(n+1, n); // copy R from QR, as transposed for more efficient givens rotations.
    for (int j=0; j < n; ++j)
        for (int i=0; i <= j; ++i)
            R(j, i) = QR(i, j);
    
    for (int i=0; i < n; ++i)
        R(n, i) = -bb[i];

    G = zeros(n+1, n);
}

Matrix TauSolverQR::operator()(double tau)
{
    double c, s, r, z;
    int len = 1, one = 1;
    R1 = R;
    G.fill(0);
    G.diag().fill( std::sqrt(tau) );

    // givens rotations applied to R stored in row-major for efficient caching
    for (int i=_n-1; i >= 0; --i)
    {
        for (int j=i; j < _n; ++j)
        {
            r = R1(j, j);
            z = G(j, i);
            drotg_(&r, &z, &c, &s);

            len = _n + 1 - j;
            drot_(&len, &R1(j, j), &one, &G(j, i), &one, &c, &s);
        }
    }

    // solve x = R \ b. R1 is stored in row major, so it looks like a lower
    // triangular matrix in memory, thus we tell the solver to solve a
    // transposed (trans = 'T') lower triangular (uplo = 'L') system.
    char uplo = 'L';
    char trans = 'T';
    char diag = 'N';
    int nrhs = 1;
    int lda = _n+1;
    int ldb = _n;
    int info = 0;

    Matrix x(_n, 1); // get Q'*b from last row of R1
    for (double * xi = x.data(), * ri = &R1(_n, 0); xi != x.data() + _n; ++xi, ri += _n+1)
        (*xi) = (*ri);

    dtrtrs_(&uplo, &trans, &diag, &_n, &nrhs, R1.data(), &lda, x.data(), &ldb, &info);

    if (info != 0)
        throw std::runtime_error("Failed to solve R*x = b, R may be singular?");

    return x;
}

TauSolverQRfull::TauSolverQRfull(const Matrix& J, const Matrix& b) : _J(J), _b(b)
{
    QR = zeros(J.n_rows + J.n_cols, J.n_cols);
    b1 = zeros(J.n_rows + J.n_cols, 1);
}

Matrix TauSolverQRfull::operator()(double tau)
{
    int m = _J.n_rows, n = _J.n_cols;

    // construct QR = [J; sqrt(tau)*I]
    for (int j=0; j < n; ++j)
        for (int i=0; i < m; ++i)
            QR(i, j) = _J(i, j);

    for (int j=0; j < n; ++j)
        for (int i=0; i < n; ++i)
            QR(m + i, j) = 0.0;

    for (int j=0; j < n; ++j)
        QR(m + j, j) = std::sqrt(tau);

    std::copy_n(_b.data(), m, b1.data());
    std::fill_n(&b1[m], n, 0.0);

    char trans = 'N';
    int nrhs = 1;
    int M = m + n;
    int lda = M;
    int ldb = M;
    int info = 0;

    if (work.size() == 0)
    {
        lwork = -1;
        double wt;

        dgels_(&trans, &M, &n, &nrhs, QR.data(), &lda, b1.data(), &ldb, &wt, &lwork, &info);

        lwork = wt;
        work = zeros(lwork, 1);
    }

    dgels_(&trans, &M, &n, &nrhs, QR.data(), &lda, b1.data(), &ldb, work.data(), &lwork, &info);

    if (info != 0)
        throw std::runtime_error("Could not solve least squares problem using QR.");

    Matrix x(n, 1);
    for (double * xi = x.data(), * bi = b1.data(); xi != x.data() + n; ++xi, ++bi)
        (*xi) = -(*bi);
    return x;
}