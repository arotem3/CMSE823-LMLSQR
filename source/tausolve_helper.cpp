#include "lmlsqr.hpp"

tausolve_helper::tausolve_helper(Matrix& J, int solver)
{
    if (solver == 0)
    {
        X1 = J; // this is for cholesky
    }
    else if (solver == 1) // compute svd of J
    {
        // fill X1 with V, fill X2 with singular values
    }
    else if (solver  == 2) // compute QR of J
    {
        // fill X1 with Q, fill X2 with R
    }
}