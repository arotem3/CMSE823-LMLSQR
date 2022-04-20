#include "lmlsqr.hpp"
#include "Matrix.hpp"

extern "C" void dposv_(char *uplo, int *n, int *nrhs, double *A, int *lda, double *B, int *ldb, int *info);
// solves Ax = B for B, using cholesky, decomposition, where A is
// symmetric/hermitian positive definite


TauSolverChol::TauSolverChol(const Matrix& J, const Matrix& b) : _J(J)
{
    _b = J.t() * b;
}

Matrix TauSolverChol::operator()(double tau)
{
    Matrix RHS = _b;
    int n = RHS.n_rows;
    Matrix LHS = (_J.t()*_J) + (std::abs(tau)*eye(n));

    int nrhs = RHS.n_cols;
    char uplo = 'U'; // Store upper triangular part (although we won't use it)
    // Leading dimensions should be number of rows, unless I am operating on a submatrix or something
    int lda = n;
    int ldb = n;
    int info;
    
    dposv_(&uplo, &n, &nrhs, LHS.data(), &lda, RHS.data(), &ldb, &info);
    return RHS;
}
