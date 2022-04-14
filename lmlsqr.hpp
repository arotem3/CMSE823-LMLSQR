#ifndef LMLSQR_HPP
#define LMLSQR_HPP

#include "Matrix.hpp"

// computes the solution to (J'*J + tau^2 * I) * x = b using Cholesky
// decomposition
class TauSolverChol
{
    public:
    explicit TauSolverChol(const Matrix& J, const Matrix& b);

    Matrix operator()(double tau);

    private:
    const Matrix & _J;
    Matrix _b;
};

// compute the solution to (J'*J + tau^2 *I) * x = b using SVD
class TauSolverSVD
{
    public:
    explicit TauSolverSVD(const Matrix& J, const Matrix& b);

    Matrix operator()(double tau);

    private:
    Matrix Vt;
    Matrix s;
    Matrix _b;
};

// compute solution to (J'*J + tau^2 * I) * x = b using modified QR
class TauSolverQR
{
    public:
    explicit TauSolverQR(const Matrix& J, const Matrix& b);

    Matrix operator()(double tau);

    private:
    Matrix QR;
    Matrix qr_tau;
    Matrix R;
    Matrix G;
    Matrix b1;
    Matrix b2;
    Matrix _b;
};

#endif