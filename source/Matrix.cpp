#include "Matrix.hpp"

extern "C" void dgemm_(char * transa, char * transb, int * m, int * n, int * k, double * alpha, double * a, int * lda, double * b, int * ldb, double * beta, double * c, int * ldc);

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

void mult(Matrix& c, const Matrix& a, const Matrix& b)
{
    check_mult(c.n_rows, c.n_cols, a.n_rows, a.n_cols, b.n_rows, b.n_cols);
    
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
    dgemm_(&transa, &transb, &m, &n, &k, &alpha, const_cast<double*>(a.data()), &lda, const_cast<double*>(b.data()), &ldb, &beta, c.data(), &ldc);
    // matmult(c.data(), a.n_rows, a.n_cols, b.n_cols, a.data(), b.data());
}

void mult(Matrix& c, const Matrix::MatrixTranspose& a, const Matrix& b)
{
    check_mult(c.n_rows, c.n_cols, a.n_rows, a.n_cols, b.n_rows, b.n_cols);
    
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
    dgemm_(&transa, &transb, &m, &n, &k, &alpha, const_cast<double*>(a.data()), &lda, const_cast<double*>(b.data()), &ldb, &beta, c.data(), &ldc);
}

void mult(Matrix& c, const Matrix& a, const Matrix::MatrixTranspose& b)
{
    check_mult(c.n_rows, c.n_cols, a.n_rows, a.n_cols, b.n_rows, b.n_cols);
    
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

void mult(Matrix& c, const Matrix::MatrixTranspose& a, const Matrix::MatrixTranspose& b)
{
    check_mult(c.n_rows, c.n_cols, a.n_rows, a.n_cols, b.n_rows, b.n_cols);
    
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

Matrix operator*(const Matrix& a, const Matrix& b)
{
    Matrix c(a.n_rows, b.n_cols);
    mult(c, a, b);

    return c;
}

Matrix operator*(const Matrix::MatrixTranspose& a, const Matrix& b)
{
    Matrix c(a.n_rows, b.n_cols);
    mult(c, a, b);
    return c;
}

Matrix operator*(const Matrix& a, const Matrix::MatrixTranspose& b)
{
    Matrix c(a.n_rows, b.n_cols);
    mult(c, a, b);
    return c;
}

Matrix operator*(const Matrix::MatrixTranspose& a, const Matrix::MatrixTranspose& b)
{
    Matrix c(a.n_rows, b.n_cols);
    mult(c, a, b);
    return c;
}

Matrix operator*(double c, const Matrix& a)
{
    Matrix b = a;
    b *= c;
    // double * bij = b.data();
    // for (int i=0; i < b.size(); ++i, ++bij)
    //     (*bij) *= c;
    
    return b;
}

Matrix operator-(const Matrix& a, const Matrix& b)
{
    Matrix c(a.n_rows, a.n_cols);
    subtract(c, a, b);
    // Matrix c = a;
    // double * ci = c.data();
    // const double * bi = b.data();
    // for (int i=0; i < c.size(); ++i, ++ci, ++bi)
    //     (*ci) -= (*bi);
    
    return c;
}

Matrix operator+(const Matrix& a, const Matrix& b)
{
    Matrix c(a.n_rows, a.n_cols);
    add(c, a, b);
    // Matrix c = a;
    // double * ci = c.data();
    // const double * bi = b.data();
    // for (int i=0; i < c.size(); ++i, ++ci, ++bi)
    //     (*ci) += (*bi);
    
    return c;
}

Matrix operator+(const Matrix& a, double c)
{
    Matrix b = a;
    b += c;

    // for (double * bi = b.data(); bi != b.data() + b.size(); ++bi)
    //     (*bi) += c;
    
    return b;
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

Matrix eye(int n)
{
    Matrix I(n,n);
    I.eye();
    // I.diag().fill(1.0);

    // double * Iij = I.data();

    // for (int i=0; i < n; ++i, Iij += n)
    //     *(Iij + i) = 1.0;
    
    return I;
}

Matrix ones(int m, int n)
{
    Matrix x(m,n);
    x.fill(1.0);
    // std::fill_n(x.data(), x.size(), 1.0);

    return x;   
}

Matrix zeros(int m, int n)
{
    return Matrix(m,n);
}

Matrix randn(int m, int n)
{
    Matrix A(m,n);

    std::default_random_engine generator;
    std::normal_distribution<double> distribution(0., 1.);

    for (int i=0; i < A.size(); ++i)
        A[i] = distribution(generator);

    return A;
}

double norm(const Matrix& x)
{
    double s = 0;
    for (const double * xi = x.data(); xi != x.data() + x.size(); ++xi)
        s += square(*xi);
    
    return std::sqrt(s);
}

double dot(const Matrix& a, const Matrix& b)
{
    return dot(a.size(), 1, 1, a.data(), b.data());
}