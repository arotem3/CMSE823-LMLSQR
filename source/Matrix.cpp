#include "Matrix.hpp"

#define SAFE_MATRIX_MULT

extern "C" void dgemm_(char * transa, char * transb, int * m, int * n, int * k, double * alpha, double * a, int * lda, double * b, int * ldb, double * beta, double * c, int * ldc);
extern "C" void dgemv_(char * trans, int * m, int * n, double * alpha, double * a, int * lda, double * x, int * incx, double * beta, double * y, int * incy);
extern "C" double dnrm2_(int * n, double * x, int * inc);

double dot(unsigned long n, unsigned long stride_x, unsigned long stride_y, const double * x, const double * y)
{
    double p = 0.0;
    for (unsigned long i=0; i < n; ++i, x+=stride_x, y+=stride_y)
        p += (*x) * (*y);
    return p;
}

void check_mult(int c_rows, int c_cols, int a_rows, int a_cols, int b_rows, int b_cols)
{
    if (c_rows != a_rows || c_cols != b_cols || a_cols != b_rows)
        throw std::invalid_argument("cannot multiply matrices of size ("
                + std::to_string(a_rows) + ", " + std::to_string(a_cols)
                + ") and (" + std::to_string(b_rows) + ", " + std::to_string(b_cols)
                + ") and produce a matrix of size (" + std::to_string(c_rows) + ", " + std::to_string(c_cols)
                + ").");
}

// c = a * b, with specialization when b is a vector
void mult(Matrix& c, const Matrix& a, const Matrix& b)
{
    #ifdef SAFE_MATRIX_MULT
    check_mult(c.n_rows, c.n_cols, a.n_rows, a.n_cols, b.n_rows, b.n_cols);
    #endif
    
    char transa = 'N';
    char transb = 'N';
    int m = a.n_rows;
    int n = b.n_cols;
    int k = a.n_cols;
    double alpha = 1.0;
    int lda = m;
    int ldb = k;
    double beta = 0.0;
    int ldc = c.n_rows;
    
    if (n == 1)
        dgemv_(&transa, &m, &k, &alpha, const_cast<double*>(a.data()), &m, const_cast<double*>(b.data()), &n, &beta, c.data(), &n);
    else
        dgemm_(&transa, &transb, &m, &n, &k, &alpha, const_cast<double*>(a.data()), &lda, const_cast<double*>(b.data()), &ldb, &beta, c.data(), &ldc);
}

// c = a.t() * b, with specialization when b is a vector
void mult(Matrix& c, const Matrix::MatrixTranspose& a, const Matrix& b)
{
    #ifdef SAFE_MATRIX_MULT
    check_mult(c.n_rows, c.n_cols, a.n_rows, a.n_cols, b.n_rows, b.n_cols);
    #endif

    char transa = 'T';
    char transb = 'N';
    int m = a.n_rows;
    int n = b.n_cols;
    int k = a.n_cols;
    double alpha = 1.0;
    int lda = k;
    int ldb = k;
    double beta = 0.0;
    int ldc = c.n_rows;

    if (n == 1)
        dgemv_(&transa, &k, &m, &alpha, const_cast<double*>(a.data()), &k, const_cast<double*>(b.data()), &n, &beta, c.data(), &n);
    else
        dgemm_(&transa, &transb, &m, &n, &k, &alpha, const_cast<double*>(a.data()), &lda, const_cast<double*>(b.data()), &ldb, &beta, c.data(), &ldc);
}

// c = a * b.t()
void mult(Matrix& c, const Matrix& a, const Matrix::MatrixTranspose& b)
{
    #ifdef SAFE_MATRIX_MULT
    check_mult(c.n_rows, c.n_cols, a.n_rows, a.n_cols, b.n_rows, b.n_cols);
    #endif
    
    char transa = 'N';
    char transb = 'T';
    int m = a.n_rows;
    int n = b.n_cols;
    int k = a.n_cols;
    double alpha = 1.0;
    int lda = m;
    int ldb = n;
    double beta = 0.0;
    int ldc = c.n_rows;
    dgemm_(&transa, &transb, &m, &n, &k, &alpha, const_cast<double*>(a.data()), &lda, const_cast<double*>(b.data()), &ldb, &beta, c.data(), &ldc);
}

// c = a.t() * b.t()
void mult(Matrix& c, const Matrix::MatrixTranspose& a, const Matrix::MatrixTranspose& b)
{
    #ifdef SAFE_MATRIX_MULT
    check_mult(c.n_rows, c.n_cols, a.n_rows, a.n_cols, b.n_rows, b.n_cols);
    #endif
    
    char transa = 'T';
    char transb = 'T';
    int m = a.n_rows;
    int n = b.n_cols;
    int k = a.n_cols;
    double alpha = 1.0;
    int lda = k;
    int ldb = n;
    double beta = 0.0;
    int ldc = c.n_rows;
    dgemm_(&transa, &transb, &m, &n, &k, &alpha, const_cast<double*>(a.data()), &lda, const_cast<double*>(b.data()), &ldb, &beta, c.data(), &ldc);
}

void add(Matrix& c, const Matrix& a, const Matrix& b)
{
    if (c.size() != a.size() || a.size() != b.size())
        throw std::invalid_argument("cannot add matrices of incompatible lengths.");

    double * ci = c.data();
    const double * ai = a.data(), * bi = b.data();

    for (int i=0; i < c.size(); ++i, ++ci, ++bi, ++ai)
        (*ci) = (*ai) + (*bi);
}

void subtract(Matrix& c, const Matrix& a, const Matrix& b)
{
    if (c.size() != a.size() || a.size() != b.size())
        throw std::invalid_argument("cannot add matrices of incompatible lengths.");

    double * ci = c.data();
    const double * ai = a.data(), * bi = b.data();

    for (int i=0; i < c.size(); ++i, ++ci, ++bi, ++ai)
        (*ci) = (*ai) - (*bi);
}

Matrix::MatrixDiagonalView& Matrix::MatrixDiagonalView::operator+=(const Matrix& B)
{
    int n = std::min(_nr, _nc);

    double * a = arr;
    const double * b = B.data();
    for (int i=0; i < n; ++i, a += _nr, ++b)
        *(a + i) += (*b);

    return *this;
}

Matrix::MatrixDiagonalView& Matrix::MatrixDiagonalView::operator-=(const Matrix& B)
{
    int n = std::min(_nr, _nc);

    double * a = arr;
    const double * b = B.data();
    for (int i=0; i < n; ++i, a += _nr, ++b)
        *(a + i) -= (*b);

    return *this;
}

Matrix vcat(const Matrix& a, const Matrix& b)
{
    if (a.n_cols != b.n_cols)
        throw std::invalid_argument(
            "cannot vertically concatenate arrays of size ("
                + std::to_string(a.n_rows) + ", " + std::to_string(a.n_cols)
                + ") and (" + std::to_string(b.n_rows) + ", " + std::to_string(b.n_cols) + ")."
            );

    Matrix c(a.n_rows + b.n_rows, a.n_cols);

    double * cij = c.data();
    const double * aij = a.data();
    const double * bij = b.data();

    for (int i=0; i < c.n_cols; ++i)
    {
        for (int j=0; j < a.n_rows; ++j, ++cij, ++aij)
            (*cij) = (*aij);

        for (int j=0; j < b.n_rows; ++j, ++cij, ++bij)
            (*cij) = (*bij);
    }

    return c;
}

Matrix randn(int m, int n)
{
    Matrix A(m,n);

    std::default_random_engine generator;
    std::normal_distribution<double> distribution(0., 1.);

    A.for_each([&](double& a) -> void {a = distribution(generator);});

    return A;
}

double norm(const Matrix& x)
{
    int n = x.size();
    int one = 1;
    return dnrm2_(&n, const_cast<double*>(x.data()), &one);
}

double dot(const Matrix& a, const Matrix& b)
{
    return dot(a.size(), 1, 1, a.data(), b.data());
}