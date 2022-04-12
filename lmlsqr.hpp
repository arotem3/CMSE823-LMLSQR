#ifndef LMLSQR_HPP
#define LMLSQR_HPP

#include "Matrix.hpp"

// initialize this in lm function by passing jacobian matrix J, and it computes
// factorization if necessary
struct tausolve_helper
{
    Matrix X1, X2;
    Matrix & _J;
    
    // for svd, work1 == 1 if b <- V.t() * b. for qr, work1 == 1 if b <- Q.t() * b
    int work1;

    tausolve_helper(Matrix& J, int solver);
    // do something
};

// inline void check_problem(Matrix& J, Matrix& b)
// {
//     if (J.n_rows != b.n_rows)
//         throw std::invalid_argument("matrix dimensions ")
// }

// computes the solution to (J'*J + tau^2 * I) * x = b using Cholesky
// decomposition
Matrix tausolve_chol(tausolve_helper& helper, double tau, Matrix& b);

// compute the solution to (J'*J + tau^2 *I) * x = b using SVD
Matrix tausolve_svd(tausolve_helper& helper, double tau, Matrix& b);

// compute solution to (J'*J + tau^2 * I) * x = b using modified QR
Matrix tausolve_qr(tausolve_helper& helper, double tau, Matrix& b);

#endif