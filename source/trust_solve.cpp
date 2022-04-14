#include "lmlsqr.hpp"
#include "Matrix.hpp"
#include <limits>
#include <cmath>

template <class tau_solver>
Matrix trust_solve(Matrix &J, Matrix &r, double delta, double tolerance)
{
    tau_solver my_tau_solve(J, r);

    double tau = 0;
    double tau_old;
    double eps = std::sqrt(std::numeric_limits<double>::epsilon());

    double g;
    double g_forward;
    double g_derivative;

    Matrix s;
    Matrix s_forward;

    int max_iter = 1000;
    int i = 0;
    while (i < max_iter) { 
        s = my_tau_solve(tau);
        s_forward = my_tau_solve(tau + eps);

        g = (1/norm(s)) - (1/delta);r
        g_forward = (1/norm(s_forward)) - (1/delta);r
        g_derivative = (g_forward-g)/eps;

        tau_old = tau;
        tau = tau_old - (g/g_derivative);

        if (std::abs(tau-tau_old) < tolerance)
            break;

        ++i;
    }
    return s;

}
