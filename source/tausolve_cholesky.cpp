#include "lmlsqr.hpp"
#include "Matrix.hpp"

extern "C" void dposv_(char *uplo, int *n, int *nrhs, double *A, int *lda, double *B, int *ldb, int *info);
// solves Ax = B for B, using cholesky, decomposition, where A is
// symmetric/hermitian positive definite

// Solving (J'J + tau^2I)x = b
Matrix tausolve_chol(tausolve_helper& helper, double tau, Matrix& b) {
    // helper.X1 = J for cholesky
    Matrix J = helper.X1;
    Matrix RHS = J.t()*b;
    int n = RHS.n_rows;
    Matrix LHS = (J.t()*J) + (tau*tau*eye(n));

    int nrhs = RHS.n_cols;
    char uplo = 'U'; // Store upper triangular part (although we won't use it)
    // Leading dimensions should be number of rows, unless I am operating on a submatrix or something
    int lda = n;
    int ldb = n;
    int info;
    
    dposv_(&uplo, &n, &nrhs, LHS.data(), &lda, RHS.data(), &ldb, &info);
    return RHS;
}
