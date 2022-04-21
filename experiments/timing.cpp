#include "timing.hpp"
#include <unordered_map>
#include <fstream>

int main()
{
    int ms[] = {100, 500, 1000};//, 2500};
    double ratio_n_to_m[] = {0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.};
    const int num_tests = 5;

    std::ofstream fm("ms.txt");
    for (int m : ms)
        fm << m << ' ';
    fm << '\n';

    std::ofstream fratios("ratios.txt");
    for (double ratio : ratio_n_to_m)
        fratios << ratio << ' ';
    fratios << '\n';

    std::ofstream fchol("timing_chol.txt");
    std::ofstream fqr("timing_qr.txt");
    std::ofstream fsvd("timing_svd.txt");

    for (int m : ms)
    {
        for (double ratio : ratio_n_to_m)
        {
            int n = m * ratio;
            
            double tchol = 0, tsvd = 0, tqr = 0;

            for (int i=0; i < num_tests; ++i)
            {
                tqr += TIMING_trust_solve<TauSolverQR>(m, n) / num_tests;
                tchol += TIMING_trust_solve<TauSolverChol>(m, n) / num_tests;
                tsvd += TIMING_trust_solve<TauSolverSVD>(m, n) / num_tests;
            }

            fchol << tchol << ' ';
            fqr << tqr << ' ';
            fsvd << tsvd << ' ';

            std::cout << "m: " << m << " | n: " << n << '\n';
        }
        fchol << std::endl;
        fqr << std::endl;
        fsvd << std::endl;
    }

    return 0;
}