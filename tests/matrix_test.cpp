#include "all_tests.hpp"

int TEST_matrix()
{
    Matrix a = {{1.0, 3.0, 5.0},{2.0, 4.0, 6.0}};
    Matrix b = {{0.0, -1.0}, {0.5, 0.5}};

    Matrix ab = a * b;
    Matrix ab_exact = {{-2., -4., -6.},{1.5, 3.5, 5.5}};

    bool test1 = norm(ab - ab_exact) < 1e-10;

    if (not test1)
        std::cout << "failed to compute a*b\n";
    
    Matrix abt = a * b.t();
    Matrix abt_exact = {{1., 2., 3.}, {0., -1., -2.}};

    bool test2 = norm(abt - abt_exact) < 1e-10;

    if (not test2)
        std::cout << "failed to compute a*b.t()\n";

    Matrix c = {{8., 3., 4.},{1., 5., 9.}, {6., 7., 2.}};

    Matrix atc = a.t() * c;
    Matrix atc_exact = {{37., 52.},{61., 76.}, {37., 52.}};

    bool test3 = norm(atc - atc_exact) < 1e-10;

    if (not test3)
        std::cout << "failed to compute a.t()*c\n";

    Matrix atct = a.t() * c.t();
    Matrix atct_exact = {{41., 56.}, {53., 68.}, {41., 56.}};

    bool test4 = norm(atct - atct_exact) < 1e-10;

    if (not test4)
        std::cout << "failed to compute a.t()*c.t()\n";

    return test1 && test2 && test3 && test4;
}