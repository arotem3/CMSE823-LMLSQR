#include "lmlsqr.hpp"
#include "Matrix.hpp"

// solves diag(s.^2 + tau^2)*x = b where s is a 1d vector, b is a 1d vector.
Matrix diagonal_solve(const Matrix& s, double tau, const Matrix& b)
{
    Matrix x = b;

    const double * si = s.data();
    double * xi = x.data();

    for (int i=0; i < x.size(); ++i, ++xi, ++si)
        (*xi) /= square(*si) + square(tau);

    return x;
}

Matrix tausolve_svd(tausolve_helper& helper, double tau, Matrix& b)
{
    Matrix & Vt = helper.X1;
    Matrix & s = helper.X2; // for clarity

    if (helper.work1 == 0)
        b = Vt * (helper._J.t() * b); // set b <- V.t() * b

    Matrix x = diagonal_solve(s, tau, b);

    x = Vt.t() * x;

    return x;
}