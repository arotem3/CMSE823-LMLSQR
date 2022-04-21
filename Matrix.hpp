#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <algorithm>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <iostream>
#include <vector>
#include <random>

inline double square(double x)
{
    return x*x;
}

double dot(unsigned long n, unsigned long stride_x, unsigned long stride_y, const double * x, const double * y);
void matmult(double * z, unsigned long m, unsigned long n, unsigned long k, const double * x, const double * y);

class Matrix
{
    public:
    const unsigned long& n_rows = _nr;
    const unsigned long& n_cols = _nc;

    Matrix() : _nr(0), _nc(0) {}

    explicit Matrix(unsigned long m, unsigned long n) : _nr(m), _nc(n)
    {
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
        arr = new double[mat.size()];
        std::copy_n(mat.arr, mat.size(), arr);
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
        arr = new double[_nr*_nc];

        double * y = arr;
        for (auto xi = x.begin(); xi != x.end(); ++xi)
            for (auto xij = xi->begin(); xij != xi->end(); ++xij, ++y)
                (*y) = *xij;
    }

    Matrix& operator=(const Matrix& mat)
    {
        // delete[] arr;
        if (size() < mat.size())
        {
            delete[] arr;
            arr = new double[mat.size()];
        }
        _nr = mat.n_rows;
        _nc = mat.n_cols;
        // arr = new double[mat.size()];
        std::copy_n(mat.arr, mat.size(), arr);

        return *this;
    }

    Matrix& operator=(Matrix&& mat)
    {
        delete[] arr;
        arr = mat.arr;
        _nc = mat.n_cols;
        _nr = mat.n_rows;
        mat.arr = nullptr;

        return *this;
    }

    inline int size() const
    {
        return _nr * _nc;
    }

    inline double* data()
    {
        return arr;
    }

    inline const double* data() const
    {
        return arr;
    }

    inline double& operator[](unsigned long i)
    {
        return arr[i];
    }

    inline double operator[](unsigned long i) const
    {
        return arr[i];
    }

    inline double& operator()(unsigned long i, unsigned long j)
    {
        return arr[i + _nr*j];
    }

    inline double operator()(unsigned long i, unsigned long j) const
    {
        return arr[i + _nr*j];
    }

    class MatrixTranspose
    {
        public:
        const unsigned long & n_rows = _nr;
        const unsigned long & n_cols = _nc;

        explicit MatrixTranspose(const double* a, unsigned long nr, unsigned long nc) : arr(a), _nr(nr), _nc(nc) {}

        inline int size() const
        {
            return _nr * _nc;
        }

        inline const double * data() const
        {
            return arr;
        }

        inline double operator()(unsigned long i, unsigned long j) const
        {
            return arr[j + _nc*i];
        }

        private:
        const double * arr;
        unsigned long _nr, _nc;
    };

    Matrix(const MatrixTranspose& At)
    {
        _nr = At.n_rows;
        _nc = At.n_cols;
        arr = new double[size()];
        for (unsigned long i=0; i < _nr; ++i)
            for (unsigned long j=0; j < _nc; ++j)
                arr[i + _nr*j] = At(i, j);
    }

    Matrix& operator=(const MatrixTranspose& At)
    {
        if (size() < At.size())
        {
            delete[] arr;
            arr = new double[size()];
        }

        _nr = At.n_rows;
        _nc = At.n_cols;
        for (unsigned long i=0; i < _nr; ++i)
            for (unsigned long j=0; j < _nc; ++j)
                arr[i + _nr*j] = At(i, j);

        return *this;
    }

    MatrixTranspose t() const
    {
        return MatrixTranspose(arr, n_cols, n_rows);
        // Matrix At(_nc, _nr);
        // for (unsigned long i=0; i < _nr; ++i)
        //     for (unsigned long j=0; j < _nc; ++j)
        //         At(j, i) = (*this)(i,j);

        // return At;
    }

    std::string print()
    {
        std::stringstream out;
        out << std::setprecision(4);
        for (unsigned long i=0; i < n_rows; ++i)
            {
                for (unsigned long j=0; j < n_cols; ++j)
                    out << std::setw(10) << (*this)(i,j);
                out << '\n';
            }

        return out.str();
    }

    inline void fill(double x = 0.0)
    {
        std::fill_n(arr, size(), x);
    }

    inline void eye()
    {
        fill(0.0);
        diag().fill(1.0);
    }

    template <typename Func>
    inline void for_each(Func f)
    {
        std::for_each_n(arr, size(), f);
    }

    inline Matrix& operator*=(double c)
    {
        std::for_each_n(arr, size(), [&](double& a) -> void {a *= c;});
        // for (double * a = arr; a != arr + size(); ++a)
            // (*a) *= c;

        return *this;
    }

    inline Matrix& operator+=(double c)
    {
        std::for_each_n(arr, size(), [&](double& a) -> void {a += c;});
        return *this;
    }

    inline Matrix& operator-=(double c)
    {
        std::for_each_n(arr, size(), [&](double& a) -> void {a -= c;});
        return *this;
    }

    inline Matrix& operator+=(const Matrix& b)
    {
        const double * bi = b.data();
        for (double * ai = arr; ai != arr + size(); ++ai, ++bi)
            (*ai) += (*bi);
        return *this;
    }

    inline Matrix& operator-=(const Matrix& b)
    {
        const double * bi = b.data();
        for (double * ai = arr; ai != arr + size(); ++ai, ++bi)
            (*ai) -= (*bi);
        return *this;
    }

    class MatrixDiagonalView
    {
        public:
        const unsigned long& n_rows = _nr;
        const unsigned long& n_cols = _nc;

        explicit MatrixDiagonalView(double * a, unsigned long nr, unsigned long nc) : arr(a), _nr(nr), _nc(nc) {}

        inline MatrixDiagonalView& operator*=(double c)
        {
            for_each([&](double& a) -> void {a *= c;});
            return *this;
        }
        inline MatrixDiagonalView& operator+=(double c)
        {
            for_each([&](double& a) -> void {a += c;});
            return *this;
        }
        inline MatrixDiagonalView& operator-=(double c)
        {
            for_each([&](double& a) -> void {a -= c;});
            return *this;
        }
        MatrixDiagonalView& operator+=(const Matrix&);
        MatrixDiagonalView& operator-=(const Matrix&);
        
        inline void fill(double c = 0.0)
        {
            for_each([&](double& a) -> void {a = c;});
        }

        template <typename Func>
        void for_each(Func f)
        {
            int n = std::min(_nr, _nc);

            double * a = arr;
            for (int i=0; i < n; ++i, a += _nr)
                f( *(a + i) );
        }
        
        private:
        double * arr;
        unsigned long _nr, _nc;
    };

    inline MatrixDiagonalView diag()
    {
        return MatrixDiagonalView(arr, _nr, _nc);
    }

    private:
    double * arr = nullptr;
    unsigned long _nr, _nc;
};

void mult(Matrix&, const Matrix&, const Matrix&);
void mult(Matrix&, const Matrix::MatrixTranspose&, const Matrix&);
void mult(Matrix&, const Matrix&, const Matrix::MatrixTranspose&);
void mult(Matrix&, const Matrix::MatrixTranspose&, const Matrix::MatrixTranspose&);
void add(Matrix&, const Matrix&, const Matrix&);
void subtract(Matrix&, const Matrix&, const Matrix&);

inline Matrix operator*(const Matrix& a, const Matrix& b)
{
    Matrix c(a.n_rows, b.n_cols);
    mult(c, a, b);
    return c;
}

inline Matrix operator*(const Matrix::MatrixTranspose& a, const Matrix& b)
{
    Matrix c(a.n_rows, b.n_cols);
    mult(c, a, b);
    return c;
}

inline Matrix operator*(const Matrix& a, const Matrix::MatrixTranspose& b)
{
    Matrix c(a.n_rows, b.n_cols);
    mult(c, a, b);
    return c;
}

inline Matrix operator*(const Matrix::MatrixTranspose& a, const Matrix::MatrixTranspose& b)
{
    Matrix c(a.n_rows, b.n_cols);
    mult(c, a, b);
    return c;
}

inline Matrix operator*(double c, const Matrix& a)
{
    Matrix b = a;
    b *= c;
    return b;
}

inline Matrix operator-(const Matrix& a, const Matrix& b)
{
    Matrix c(a.n_rows, a.n_cols);
    subtract(c, a, b);
    return c;
}

inline Matrix operator+(const Matrix& a, const Matrix& b)
{
    Matrix c(a.n_rows, a.n_cols);
    add(c, a, b);
    return c;
}

inline Matrix operator+(const Matrix& a, double c)
{
    Matrix b = a;
    b += c;
    return b;
}

double norm(const Matrix& x);
double dot(const Matrix& a, const Matrix& b);

Matrix vcat(const Matrix& a, const Matrix& b);

inline Matrix eye(int n)
{
    Matrix I(n,n);
    I.eye();
    return I;
}

inline Matrix ones(int m, int n)
{
    Matrix A(m,n);
    A.fill(1.0);
    return A;
}

inline Matrix zeros(int m, int n)
{
    return Matrix(m, n);
}

Matrix randn(int m, int n);

#endif