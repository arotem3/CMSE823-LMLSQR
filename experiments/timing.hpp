#include <lmlsqr.hpp>

template <typename TauSolver>
double TIMING_trust_solve(int m, int n)
{
    Matrix J = randn(m, n);
    Matrix r = J * randn(n, 1) + 0.01*randn(m,1); // r is approximately in column space of J

    double delta = 0.5;
    double tolerance = 1e-10;

    Timer clock;
    trust_solve<TauSolver>(J, r, delta, tolerance);
    double seconds = clock.elapsed_time();

    return seconds;
}