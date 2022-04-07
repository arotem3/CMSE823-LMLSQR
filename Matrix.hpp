#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <algorithm>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <iostream>
#include <vector>

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
        delete[] arr;
        _nr = mat.n_rows;
        _nc = mat.n_cols;
        arr = new double[mat.size()];
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

    private:
    double * arr = nullptr;
    unsigned long _nr, _nc;
};

Matrix operator*(const Matrix& a, const Matrix& b);
Matrix operator*(double c, const Matrix& a);
Matrix operator-(const Matrix& a, const Matrix& b);
Matrix operator+(const Matrix& a, const Matrix& b);
Matrix operator+(const Matrix& a, double c);

double norm(const Matrix& x);
double dot(const Matrix& a, const Matrix& b);

Matrix vcat(const Matrix& a, const Matrix& b);

Matrix eye(int n);
Matrix ones(int m, int n);
Matrix zeros(int m, int n);

#endif