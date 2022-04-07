#include <iostream>
#include "lmlsqr.hpp"
#include "Matrix.hpp"

bool TEST_tausolve_svd()
{
    Matrix J = {{0.5, -std::sqrt(2.0), 0.5},
                {0.5, std::sqrt(2.0), 0.5}};

    double tau = 0.1;

    Matrix xe = {{1.0, 3.14}};

    Matrix b = (J.t() * J + square(tau)*eye(2)) * xe;

    tausolve_helper helper(J, 1);

    Matrix x = tausolve_svd(helper, tau, b);

    bool success = norm(x - xe) < 1e-10;

    return success;
}

// g++ ./tests/tausolve_svd.cpp ./source/tausolve_helper.cpp ./source/tausolve_svd.cpp ./source/Matrix.cpp -I ./ -llapack

int main()
{

    if (TEST_tausolve_svd())
        std::cout << "tausolve_svd() passed test\n";
    else
        std::cout << "tausolve_svd() failed test\n";

    return 0;
}
