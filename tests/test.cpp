#include "all_tests.hpp"
#include "lmlsqr.hpp"

int main()
{
    int n_tests = 0, n_passed = 0;

    std::cout << "Testing Tausolve SVD\n";
    n_passed += TEST_tausolve_svd(); ++n_tests;
    std::cout << "Testing Tausolve Cholesky\n";
    n_passed += TEST_tausolve_chol(); ++n_tests;
    std::cout << "Testing Tausolve QR\n";
    n_passed += TEST_tausolve_qr(); ++n_tests;
    std::cout << "Testing Trust Solve SVD\n";
    n_passed += TEST_trust_solve<TauSolverSVD>(); ++n_tests;
    std::cout << "Testing Trust Solve Cholesky\n";
    n_passed += TEST_trust_solve<TauSolverChol>(); ++n_tests;
    std::cout << "Testing Trust Solve QR\n";
    n_passed += TEST_trust_solve<TauSolverQR>(); ++n_tests;

    std::cout << n_passed << " / " << n_tests << " tests completed successfully.\n";

    return 0;
}
