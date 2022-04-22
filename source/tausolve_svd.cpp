#include "lmlsqr.hpp"
#include "Matrix.hpp"

extern "C" void dgesvd_(char* jobu, char * jobvt, int * m, int * n, double * a, int * lda, double * s, double * u, int * ldu, double * vt, int * ldvt, double * work, int * lwork, int * info);

TauSolverSVD::TauSolverSVD(const Matrix& J, const Matrix& b)
{
    char jobu = 'O';
    char jobvt = 'A';
    int m = J.n_rows;
    int n = J.n_cols;
    Matrix U = J;
    int lda = m;

    s = zeros(n, 1);
    Vt = zeros(n, n);
    int ldu = m;
    int ldvt = n;

    int lwork = -1;
    double wt;

    int info = 0;

    // query for work variables
    dgesvd_(&jobu, &jobvt, &m, &n, U.data(), &lda, s.data(), nullptr, &ldu, Vt.data(), &ldvt, &wt, &lwork, &info);

    lwork = wt;
    Matrix work(lwork, 1);

    // actually compute svd
    dgesvd_(&jobu, &jobvt, &m, &n, U.data(), &lda, s.data(), nullptr, &ldu, Vt.data(), &ldvt, work.data(), &lwork, &info);

    if (info != 0)
        throw std::runtime_error("failed to compute singular value decomposition.");

    // compute -s * U.t() * b
    _b = U.t() * b;
    for (double * bi = _b.data(), * si =  s.data(); bi != _b.data() + n; ++bi, ++si)
        (*bi) *= -(*si);

    // square s
    s.for_each([](double& si) -> void {si = square(si);});
}

// solves diag(s + tau)*x = b where s is a 1d vector, b is a 1d vector.
Matrix diagonal_solve(const Matrix& s, double tau, const Matrix& b)
{
    Matrix x(b.size(), 1);

    const double * si = s.data();
    const double * bi = b.data();
    double * xi = x.data();

    for (int i=0; i < x.size(); ++i, ++xi, ++si, ++bi)
        (*xi) = (*bi) / ( (*si) + tau );

    return x;
}

Matrix TauSolverSVD::operator()(double tau)
{
    Matrix x = diagonal_solve(s, tau, _b);

    x = Vt.t() * x;

    return x;
}
