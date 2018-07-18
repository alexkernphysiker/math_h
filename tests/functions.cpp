// this file is distributed under
// LGPLv3 license
#include <gtest/gtest.h>
#include <math_h/functions.h>
#include <math_h/integrate.h>
using namespace std;
using namespace MathTemplates;
//This file contains unit tests for functions.h
template<class numt>
void Test_Peak_Shape(const function<numt(numt)> f, const numt &from, const numt &peak, const numt &to)
{
    for (numt x = from; x <= to; x += 0.06) {
        EXPECT_TRUE(f(x) >= 0);
        if (x <= peak)
            EXPECT_TRUE(f(x) >= f(x - 0.02));
        else
            EXPECT_TRUE(f(x) >= f(x + 0.02));
    }
}
template<class numt>
void Test_Norm(const function<numt(numt)> f)
{
    const auto S = Sympson(f, -300., +300., 0.001);
    EXPECT_TRUE(pow(S - 1., 2) < 0.0001);
}

TEST(FunctionsTest,Gaussian)
{
    for (double sigma = 0.5; sigma < 5; sigma += 0.5)
        for (double X = -5; X <= 5; X += 1) {
            Test_Peak_Shape<double>([sigma, X](double x) {
                return Gaussian(x, X, sigma);
            }, X - sigma * 5, X, X + sigma * 5);
            Test_Norm<double>([sigma, X](double x) {
                return Gaussian(x, X, sigma);
            });
        }
}
TEST(FunctionsTest,Lorentzian)
{
    for (double gamma = 0.5; gamma < 5; gamma += 0.5)
        for (double X = -5; X <= 5; X += 1) {
            Test_Peak_Shape<double>([gamma, X](double x) {
                return Lorentzian(x, X, gamma);
            }, X - gamma * 5, X, X + gamma * 5);
            Test_Norm<double>([gamma, X](double x) {
                return Lorentzian(x, X, gamma);
            });
        }
}
TEST(FunctionsTest,BreitWigner)
{
    for (double gamma = 0.5; gamma < 5; gamma += 0.5)
        for (double X = -5; X <= 5; X += 1) {
            Test_Peak_Shape<double>([gamma, X](double x) {
                return BreitWigner(x, X, gamma);
            }, X - gamma * 2, X, X + gamma * 2);
            Test_Norm<double>([gamma, X](double x) {
                return BreitWigner(x, X, gamma);
            });
        }
}
TEST(FunctionsTest,Novosibirsk)
{
    for (double sigma = 0.5; sigma < 10; sigma += 0.5)
        for (double X = -5; X <= 5; X += 1)for (double asym = -1; asym <= 1; asym += 0.1)
                Test_Peak_Shape<double>([sigma, X, asym](double x) {
                return Novosibirsk(x, X, sigma, asym);
            }, X - sigma * 2, X, X + sigma * 2);
}

template<class numt>
void Test_stair_Shape(function<numt(numt)> f, numt from, numt middle, numt to)
{
    numt a = f(from);
    numt b = f(to);
    numt m = f(middle);
    EXPECT_TRUE(a != b);
    EXPECT_TRUE(pow((m - a) * (b - a) - 0.5, 2) < 0.4);
    for (numt x = from; x <= to; x += 0.05) {
        EXPECT_TRUE((f(x) - f(x - 0.01)) * (a - b) <= 0);
    }
}
TEST(FunctionsTest,FermiFunc)
{
    for (double diffuse = 0.5; diffuse < 10; diffuse += 0.5)
        for (double X = -5; X <= 5; X += 1)
            Test_stair_Shape<double>([diffuse, X](double x) {
            return FermiFunc(x, X, diffuse);
        }, X - diffuse * 5, X, X + diffuse * 5);
}

template<size_t P=10>
void test_polynom(function<double(int)> f)
{
    double C[P + 1];
    for (int i = 0; i <= P; i++)C[i] = f(i);
    for (double x = -2; x <= 2; x += 0.01) {
        double V = 0;
        for (int p = 0; p <= P; p++) {
            V += pow(x, p) * C[p];
        }
        EXPECT_EQ(true, pow(V - Polynom<P>(x, C), 2) < 0.0001);
    }
#ifdef ____full_version_of_math_h_____
    if constexpr(P>0) test_polynom<P-1>(f);
#endif
}

TEST(FunctionsTest,Polynom_zeros)
{
    test_polynom([](int) {
        return 0;
    });
}
TEST(FunctionsTest,Polynom_ones)
{
    test_polynom([](int) {
        return 1;
    });
}
TEST(FunctionsTest,Polynom_minus_ones)
{
    test_polynom([](int) {
        return -1;
    });
}
TEST(FunctionsTest,Polynom_variable_coefs)
{
    test_polynom([](int i) {
        return i;
    });
    test_polynom([](int i) {
        return pow(-1, i);
    });
    test_polynom([](int i) {
        return pow(-1, i + 1);
    });
    test_polynom([](int i) {
        return i * pow(-1, i);
    });
    test_polynom([](int i) {
        return i * pow(-1, i + 1);
    });
}
