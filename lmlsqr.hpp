#ifndef LMLSQR_HPP
#define LMLSQR_HPP

#include <iostream>
#include <limits>
#include <cmath>
#include <chrono>
#include "Matrix.hpp"

// computes the solution to (J'*J + tau^2 * I) * x = b using Cholesky
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

//
template <class tau_solver>
Matrix trust_solve(Matrix &J, Matrix &r, double delta, double tolerance)
{
    tau_solver my_tau_solve(J, r);

    double tau = 0;
    double tau_old;
    //double eps = std::sqrt(std::numeric_limits<double>::epsilon());
    double eps = 1e-5;

    double g;
    double g_forward;
    double g_derivative;

    Matrix s;
    Matrix s_forward;

    int max_iter = 10;
    int i = 0;
    while (i < max_iter) { 
        s = my_tau_solve(tau);
        // if (norm(s) <= delta)
            // break;
        s_forward = my_tau_solve(tau + eps);

        g = (1/norm(s)) - (1/delta);
        g_forward = (1/norm(s_forward)) - (1/delta);
        g_derivative = (g_forward-g)/eps;

        //std::cout << "Iter " << i << "; tau = " << tau << "; g' = " << g_derivative << "\n";

        tau_old = tau;
        tau = tau_old - (g/g_derivative);

        // if (std::abs(/tau-tau_old) < tolerance)
            // break;

        ++i;
    }
    // std::cout << i << '\n';
    return s;

}

class Timer
{
    public:
    inline Timer()
    {
        t0 = std::chrono::high_resolution_clock::now();
    }

    inline double elapsed_time()
    {
        auto t1 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> diff = t1 - t0;
        return diff.count();
    }

    private:
    std::chrono::_V2::high_resolution_clock::time_point t0;
};

#endif
