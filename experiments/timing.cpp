#include "timing.hpp"
#include <unordered_map>
#include <fstream>

void TIMING_fixed_max_iter(int k)
{
    constexpr int ms[] = {100, 500, 1000};
    constexpr double ratio_n_to_m[] = {0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.};
    constexpr int num_tests[] = {10, 10, 5};

    std::ofstream fm("ms.txt");
    for (int m : ms)
        fm << m << ' ';
    fm << std::endl;

    std::ofstream fratios("ratios.txt");
    for (double ratio : ratio_n_to_m)
        fratios << ratio << ' ';
    fratios << std::endl;

    std::string kstr = std::to_string(k);
    std::ofstream fchol("timing_chol_" + kstr + ".txt");
    std::ofstream fqr("timing_qr_" + kstr + ".txt");
    std::ofstream fqrf("timing_qrf_" + kstr + ".txt");
    std::ofstream fsvd("timing_svd_" + kstr + ".txt");

    int j = 0;
    for (int m : ms)
    {
        for (double ratio : ratio_n_to_m)
        {
            int n = m * ratio;
            
            double tchol = 0, tsvd = 0, tqr = 0, tqrf = 0;

            for (int i=0; i < num_tests[j]; ++i)
            {
                tqr += TIMING_trust_solve<TauSolverQR>(m, n, k) / num_tests[j];
                tqrf += TIMING_trust_solve<TauSolverQRfull>(m, n, k) / num_tests[j];
                tchol += TIMING_trust_solve<TauSolverChol>(m, n, k) / num_tests[j];
                tsvd += TIMING_trust_solve<TauSolverSVD>(m, n, k) / num_tests[j];
            }

            fchol << tchol << ' ' << std::flush;
            fqr << tqr << ' ' << std::flush;
            fsvd << tsvd << ' ' << std::flush;
            fqrf << tqrf << ' ' << std::flush;

            std::cout << "m: " << std::setw(5) << m << " | n: " << std::setw(5) << n << '\r';
        }
        std::cout << '\n';
        fchol << std::endl;
        fqr << std::endl;
        fsvd << std::endl;
        fqrf << std::endl;
    }
    ++j;
}

void TIMING_large_m_small_n()
{
    constexpr int m = 5000;
    constexpr int ns[] = {10, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500};
    constexpr int num_tests = 5;
    constexpr int max_iter = 5;

    std::ofstream fn("lmsn_n.txt");
    for (int n : ns)
        fn << n << ' ';
    fn << std::endl;

    std::ofstream fchol("lmsn_chol.txt");
    std::ofstream fsvd("lmsn_svd.txt");
    std::ofstream fqr("lmsn_qr.txt");
    std::ofstream fqrf("lmsn_qrf.txt");
    
    for (int n : ns)
    {
        double tchol=0, tsvd=0, tqr=0, tqrf=0;
        for (int i=0; i < num_tests; ++i)
        {
            tchol = TIMING_trust_solve<TauSolverChol>(m, n, max_iter) / num_tests;
            tsvd = TIMING_trust_solve<TauSolverSVD>(m, n, max_iter) / num_tests;
            tqr = TIMING_trust_solve<TauSolverQR>(m, n, max_iter) / num_tests;
            tqrf = TIMING_trust_solve<TauSolverQRfull>(m, n, max_iter) / num_tests;
        }

        fchol << tchol << ' ' << std::flush;
        fsvd << tsvd << ' ' << std::flush;
        fqr << tqr << ' ' << std::flush;
        fqrf << tqrf << ' ' << std::flush;

        std::cout << "n: " << std::setw(5) << n << '\r';
    }
}

int main()
{
    TIMING_fixed_max_iter(2);
    TIMING_fixed_max_iter(5);
    TIMING_fixed_max_iter(10);

    TIMING_large_m_small_n();

    return 0;
}