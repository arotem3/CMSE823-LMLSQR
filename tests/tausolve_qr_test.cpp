#include <iostream>
#include "lmlsqr.hpp"
#include "Matrix.hpp"

#include "all_tests.hpp"

int TEST_tausolve_qr()
{
    Matrix J = {{0.5, -std::sqrt(2.0), 0.5},
                {0.5, std::sqrt(2.0), 0.5}};

    double tau = 0.01;

    Matrix b = {{1., 1., 1.}};

    Matrix xe = {{-0.637427298586004, -1.342770721215976}};

    TauSolverQR solver(J, b);
    Matrix x = solver(tau);

    bool success = norm(x - xe) < 1e-10;

    if (not success)
        std::cout << "tausolve_qr() failed to compute accurate solution to test problem.\n";

    return success;
}
