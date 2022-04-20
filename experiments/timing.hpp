#include <lmlsqr.hpp>

template <typename TauSolver>
long TIMING_trust_solve(int m, int n)
{
    Matrix J = randn(m, n);
    Matrix r = randn(m, 1);

    double delta = 0.1;
    double tolerance = 1e-10;

    Timer clock;
    trust_solve<TauSolver>(J, r, delta, tolerance);
    long ns = clock.elapsed_time();

    return ns;
}