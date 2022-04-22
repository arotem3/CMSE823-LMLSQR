#ifndef LMLSQR_HPP
#define LMLSQR_HPP

#include <iostream>
#include <limits>
#include <cmath>
#include <chrono>
#include "Matrix.hpp"

// computes the solution to (J'*J + |tau| * I) * x = -J'*b using Cholesky
// decomposition
class TauSolverChol
{
    public:
    explicit TauSolverChol(const Matrix& J, const Matrix& b);

    Matrix operator()(double tau);

    private:
    Matrix _A;
    Matrix _b;
};

// compute the solution to (J'*J + |tau| *I) * x = -J'*b using SVD
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

// compute solution to (J'*J + |tau| * I) * x = -J'*b using modified QR
class TauSolverQR
{
    public:
    explicit TauSolverQR(const Matrix& J, const Matrix& b);

    Matrix operator()(double tau);

    private:
    Matrix R;
    Matrix R1;
    Matrix G;
    int _n;
};

// computes solution to (J'*J + |tau| * I) * x = -J'*b by naively recomputing
// the QR decomposition for all tau.
class TauSolverQRfull
{
    public:
    TauSolverQRfull(const Matrix& J, const Matrix& b);

    Matrix operator()(double tau);

    private:
    const Matrix& _J;
    const Matrix& _b;
    Matrix QR;
    Matrix b1;
    Matrix work;
    int lwork;
};

// computes solution to trust-region problem minimize{ 0.5*||J*s + r||^2 } subject
// to { ||s|| == delta } using the secant method on the Lagrangian dual variable.
template <class tau_solver>
Matrix trust_solve(Matrix &J, Matrix &r, double delta, double tolerance = 1e-8, int max_iter = 10)
{
    tau_solver my_tau_solve(J, r);

    double tau1 = 0, tau2 = 1e-5;
    double g1, g2, dg;
    double ns;

    Matrix s1 = my_tau_solve( std::abs(tau1) );
    ns = norm(s1);
    
    if (ns <= delta)
        return s1;

    g1 = (1.0/ns) - (1.0/delta);

    for (int i=0; i < max_iter; ++i)
    {
        Matrix s2 = my_tau_solve( std::abs(tau2) );
        ns = norm(s2);

        if (std::abs(ns - delta) < tolerance)
            return s2;

        g2 = (1.0/ns) - (1.0/delta);
        dg = (g2 - g1) / (tau2 - tau1);

        // If tau2 is better than tau1, keep tau2 and throw away tau1.
        // Otherwise, keep tau1 and try a new tau2.
        if (std::abs(g2) < std::abs(g1))
        {
            tau1 = tau2;
            s1 = std::move(s2);
            g1 = g2;
        }
        tau2 -= (g2 / dg);
    }

    return s1;
}

#endif
