#include <iostream>
#include "lmlsqr.hpp"
#include "Matrix.hpp"

#include "all_tests.hpp"

int TEST_tausolve_chol()
{
    Matrix J = {{0.5, -std::sqrt(2.0), 0.5},
                {0.5, std::sqrt(2.0), 0.5}};

    double tau = 0.1;

    Matrix b = {{1., 1., 1.}};

    Matrix xe = {{0.637427298586004, 1.342770721215976}};

    // Matrix xe = {{1.0, 3.14}};

    // Matrix b = (J.t() * J + square(tau)*eye(2)) * xe;

    tausolve_helper helper(J, 0);

    Matrix x = tausolve_chol(helper, tau, b);

    bool success = norm(x - xe) < 1e-10;

    if (not success)
        std::cout << "tausolve_chol() failed to compute accurate solution to test problem.\n";

    return success;
}