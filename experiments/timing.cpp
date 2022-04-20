#include "timing.hpp"
#include <unordered_map>
#include <fstream>

int main()
{
    int ms[] = {100, 500, 1000, 2500, 5000};
    double ratio_n_to_m[] = {0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.};

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
            fchol << TIMING_trust_solve<TauSolverChol>(m, n) << ' ';
            fqr << TIMING_trust_solve<TauSolverQR>(m, n) << ' ';
            fsvd << TIMING_trust_solve<TauSolverSVD>(m, n) << ' ';

            std::cout << "m: " << m << " | n: " << n << '\n';
        }
        fchol << std::endl;
        fqr << std::endl;
        fsvd << std::endl;
    }

    return 0;
}