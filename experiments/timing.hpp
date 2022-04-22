#include <lmlsqr.hpp>

class Timer
{
    public:
    inline Timer()
    {
        t0 = std::chrono::high_resolution_clock::now();
    }

    inline double elapsed_time()
    {
        auto t1 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> diff = t1 - t0;
        return diff.count();
    }

    private:
    std::chrono::_V2::high_resolution_clock::time_point t0;
};

template <typename TauSolver>
double TIMING_trust_solve(int m, int n, int max_iter)
{
    Matrix J = randn(m, n);
    Matrix r = J * randn(n, 1) + 0.01*randn(m,1); // r is approximately in column space of J

    double delta = 0.5;
    double tolerance = 0;

    Timer clock;
    trust_solve<TauSolver>(J, r, delta, tolerance, max_iter);
    double seconds = clock.elapsed_time();

    return seconds;
}