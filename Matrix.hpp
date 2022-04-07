#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <algorithm>
#include <sstream>
#include <iomanip>
#include <cmath>

inline double square(double x)
{
    return x*x;
}

// Kahan stable sum: https://en.wikipedia.org/wiki/Kahan_summation_algorithm
// double dot(unsigned long n, unsigned long stride_x, unsigned long stride_y, const double * x, const double * y)
// {
//     double p = 0.0;
//     double c = 0.0;

//     for (unsigned long i=0; i < n; ++i, x += stride_x, y += stride_y)
//     {
//         double q = (*x)*(*y) - c;
//         double t = p + q;
//         c = (t - p); c -= q;
//         p = t;
//     }

//     return p;
// }

double dot(unsigned long n, unsigned long stride_x, unsigned long stride_y, const double * x, const double * y)
{
    double p = 0.0;
    for (unsigned long i=0; i < n; ++i, x+=stride_x, y+=stride_y)
        p += (*x) * (*y);
    return p;
}

void matmult(double * z, unsigned long m, unsigned long n, unsigned long k, const double * x, const double * y)
{
    double * zz = z;
    for (unsigned long c=0; c < k; ++c)
        for (unsigned long r=0; r < m; ++r, ++zz)
            (*zz) = dot(n, m, 1ul, x + r, y + n*c);
}

class Matrix
{
    public:
    const unsigned long& n_rows = _nr;
    const unsigned long& n_cols = _nc;

    Matrix()
    {
        _nr = 0;
        _nc = 0;
    }

    explicit Matrix(unsigned long m, unsigned long n)
    {
        _nr = m;
        _nc = n;
        arr = new double[m*n]{0.0};
    }

    ~Matrix()
    {
        delete[] arr;
    }

    Matrix(const Matrix& mat)
    {
        _nr = mat.n_rows;
        _nc = mat.n_cols;
        arr = new double[mat.n_cols*mat.n_rows];
        std::copy(mat.arr, mat.arr+mat.n_cols*mat.n_rows, arr);
    }

    Matrix(Matrix&& mat)
    {
        arr = mat.arr;
        _nc = mat.n_cols;
        _nr = mat.n_rows;

        mat.arr = nullptr;
    }

    Matrix(const std::initializer_list<std::initializer_list<double>>& x)
    {
        _nc = x.size();
        _nr = x.begin()->size();
        arr = new double[_nr*_nc];

        double * y = arr;
        for (auto xi = x.begin(); xi != x.end(); ++xi)
            for (auto xij = xi->begin(); xij != xi->end(); ++xij, ++y)
                (*y) = *xij;
    }

    Matrix& operator=(const Matrix& mat)
    {
        delete[] arr;
        _nr = mat.n_rows;
        _nc = mat.n_cols;
        arr = new double[mat.n_cols*mat.n_rows];
        std::copy(mat.arr, mat.arr+mat.n_cols*mat.n_rows, arr);

        return *this;
    }

    Matrix& operator=(Matrix&& mat)
    {
        delete[] arr;
        arr = mat.arr;
        mat.arr = nullptr;
        _nc = mat.n_cols;
        _nr = mat.n_rows;

        return *this;
    }

    int size() const
    {
        return _nr * _nc;
    }

    double* data()
    {
        return arr;
    }

    const double* data() const
    {
        return arr;
    }

    double& operator[](unsigned long i)
    {
        return arr[i];
    }

    double operator[](unsigned long i) const
    {
        return arr[i];
    }

    double& operator()(unsigned long i, unsigned long j)
    {
        return arr[i + _nr*j];
    }

    double operator()(unsigned long i, unsigned long j) const
    {
        return arr[i + _nr*j];
    }

    Matrix t() const
    {
        Matrix At(_nc, _nr);
        for (unsigned long i=0; i < _nr; ++i)
            for (unsigned long j=0; j < _nc; ++j)
                At(j, i) = (*this)(i,j);

        return At;
    }

    Matrix rows(int first, int last) const
    {
        Matrix a(last-first, _nc);
        int k = 0;
        for (int i=first; i < last; ++i, ++k)
            for (int j=0; j < _nc; ++j)
                a(k, j) = (*this)(i,j);

        return a;
    }

    std::string print()
    {
        std::stringstream out;
        out << std::fixed << std::setprecision(6);
        for (unsigned long i=0; i < n_rows; ++i)
            {
                for (unsigned long j=0; j < n_cols; ++j)
                    out << '\t' << (*this)(i,j);
                out << '\n';
            }

        return out.str();
    }

    private:
    double * arr;
    unsigned long _nr, _nc;
};

Matrix operator*(const Matrix& a, const Matrix& b)
{
    Matrix c(a.n_rows, b.n_cols);
    matmult(c.data(), a.n_rows, a.n_cols, b.n_cols, a.data(), b.data());

    return c;
}

Matrix operator*(double c, const Matrix& a)
{
    Matrix b = a;
    double * bij = b.data();
    for (int i=0; i < b.size(); ++i, ++bij)
        (*bij) *= c;
    
    return b;
}

Matrix operator-(const Matrix& a, const Matrix& b)
{
    Matrix c = a;
    double * ci = c.data();
    const double * bi = b.data();
    for (int i=0; i < c.size(); ++i, ++ci, ++bi)
        (*ci) -= (*bi);
    
    return c;
}

Matrix operator+(const Matrix& a, const Matrix& b)
{
    Matrix c = a;
    double * ci = c.data();
    const double * bi = b.data();
    for (int i=0; i < c.size(); ++i, ++ci, ++bi)
        (*ci) += (*bi);
    
    return c;
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

    double * Iij = I.data();

    for (int i=0; i < n; ++i, Iij += n)
        *(Iij + i) = 1.0;
    
    return I;
}

Matrix ones(int m, int n)
{
    Matrix x(m,n);
    std::fill_n(x.data(), x.size(), 1.0);

    return x;   
}

Matrix zeros(int m, int n)
{
    return Matrix(m,n);
}

double fro(const Matrix& x)
{
    double norm = 0;
    for (const double * xi = x.data(); xi != x.data() + x.size(); ++xi)
        norm += square(*xi);
    
    return std::sqrt(norm);
}

#endif