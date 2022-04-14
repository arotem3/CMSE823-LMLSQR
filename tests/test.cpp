#include "all_tests.hpp"
#include "lmlsqr.hpp"

int main()
{
    int n_tests = 0, n_passed = 0;

    n_passed += TEST_tausolve_svd(); ++n_tests;
    n_passed += TEST_tausolve_chol(); ++n_tests;
    n_passed += TEST_tausolve_qr(); ++n_tests;
    n_passed += TEST_trust_solve<TauSolverSVD>(); ++n_tests;
    n_passed += TEST_trust_solve<TauSolverChol>(); ++n_tests;
    n_passed += TEST_trust_solve<TauSolverQR>(); ++n_tests;

    std::cout << n_passed << " / " << n_tests << " tests completed successfully.\n";

    return 0;
}
