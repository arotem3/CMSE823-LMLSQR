#ifndef ALL_TESTS_HPP
#define ALL_TESTS_HPP

#include <algorithm>
#include <cmath>
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
    double tolerance = 1e-11;


    Matrix s = trust_solve<tau_solver>(J, b, delta, tolerance);
    
    Matrix se = {{-0.01099024468912938,0.09939423786957201}};

    double tau = 4.755029833554435;

    bool success = norm(s - se) < 1e-10;
    //bool success = std::abs(norm(s) - delta) < 1e-8;
    //std::cout << "||s|| = " << norm(s) << "\n";

    if (not success) {
        std::cout << "s: " << s.print();
        std::cout << "se: " << se.print();
        std::cout << "trust_solve() failed to compute accurate solution to test problem.\n";
    }

    return success;

}


#endif
