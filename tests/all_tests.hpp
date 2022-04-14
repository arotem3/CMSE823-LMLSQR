#ifndef ALL_TESTS_HPP
#define ALL_TESTS_HPP

#include <iostream>
#include "lmlsqr.hpp"
#include "Matrix.hpp"

int TEST_tausolve_svd();
int TEST_tausolve_chol();
int TEST_tausolve_qr();

template <class tau_solver>
int TEST_trust_solve()
{
    Matrix J = {{0.5, -std::sqrt(2.0), 0.5},
                {0.5, std::sqrt(2.0), 0.5}};
    Matrix b = {{1., 1., 1.}};

    double delta = 0.1;
    double tolerance = 1e-4;


    Matrix x = trust_solve<tau_solver>(J, b, delta, tolerance);
    
    Matrix xe = {{0.637427298586004, 1.342770721215976}};
    bool success = norm(x - xe) < 1e-10;

    if (not success)
        std::cout << "trust_solve() failed to compute accurate solution to test problem.\n";

    return success;

}


#endif
