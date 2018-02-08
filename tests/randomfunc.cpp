// this file is distributed under
// LGPLv3 license
#include <gtest/gtest.h>
#include <math_h/randomfunc.h>
#include <math_h/sigma.h>
#include <math_h/hists.h>
using namespace std;
using namespace MathTemplates;
#define ALMOST_EQ(a,b) EXPECT_TRUE(abs(a-b)<0.0001)
#define ALMOST_EQ2(a,b) EXPECT_TRUE(abs(a-b)<0.01)
#define ALMOST_EQ3(a,b) EXPECT_TRUE(abs(a-b)<1.0)
TEST(RandomUniform, BaseTest)
{
    RandomUniform<> R(0, 1);
    for (int i = 0; i < 1000; i++) {
        double r = R();
        EXPECT_TRUE((r >= 0) && (r <= 1));
    }
    auto R2 = R;
    for (int i = 0; i < 1000; i++) {
        double r = R2();
        EXPECT_TRUE((r >= 0) && (r <= 1));
    }
}
TEST(RandomValueTableDistr, initFromList_n_copy)
{
    RandomValueTableDistr<> R = Points<>{{0.0, 0.25}, {0.5, 0.75}, {1.0, 0.25}};
    for (int i = 0; i < 1000; i++) {
        double r = R();
        EXPECT_TRUE((r >= 0) && (r <= 1));
    }
    auto R2 = R;
    for (int i = 0; i < 1000; i++) {
        double r = R2();
        EXPECT_TRUE((r >= 0) && (r <= 1));
    }
}
TEST(RandomValueTableDistr, initfromfunc)
{
    RandomValueTableDistr<> R([](double)->double {
        return 1;
    }, ChainWithCount(10, 0.0, 1.0));
    for (int i = 0; i < 1000; i++) {
        double r = R();
        EXPECT_TRUE((r >= 0) && (r <= 1));
    }
}
TEST(RandomValueTableDistr, Throwing)
{
    auto f = [](double) {
        return 1;
    };
    auto z = [](double) {
        return 0;
    };
    EXPECT_THROW(RandomValueTableDistr<>(f, ChainWithCount(0, 0.0, 1.0)), Exception<Chain<>>);
    EXPECT_NO_THROW(RandomValueTableDistr<>(f, ChainWithCount(1, 0.0, 1.0)));
    EXPECT_THROW(RandomValueTableDistr<>(f, ChainWithCount(2, 0.0, -1.0)), Exception<Chain<>>);
    EXPECT_THROW(RandomValueTableDistr<>(f, ChainWithCount(2, 0.0, 0.0)), Exception<Chain<>>);
    EXPECT_NO_THROW(RandomValueTableDistr<>(f, ChainWithCount(2, 0.0, 1.0)));
    EXPECT_NO_THROW(RandomValueTableDistr<>(z, ChainWithCount(2, 0.0, 1.0)));
    EXPECT_NO_THROW(RandomValueTableDistr<>(f, Chain<>{0.0, 1.0}));
    EXPECT_NO_THROW(RandomValueTableDistr<>({{0.0, 1.0}, {1.0, 1.0}}));
    EXPECT_THROW(RandomValueTableDistr<>(z, ChainWithCount(0, 0.0, 0.0)), Exception<Chain<>>);
    auto n = [](double x) {
        return sin(10 * x);
    };
    EXPECT_THROW(RandomValueTableDistr<>(n, ChainWithCount(100, 0.0, 1.0)), Exception<SortedPoints<>>);
}
void TestRandomDistribution(function<double(double)> F, double from, double to, int bins, int accu = 10)
{
    RandomValueTableDistr<> R(F, ChainWithCount(bins * accu, from, to));
    const double step = (to - from) / double(bins);
    Distribution1D<> D(BinsByStep(from, step, to));
    const double norm = Sympson(F, from, to, 0.00001);
    const int N = 10000;
    //Generate random values and fill distribution
    for (int i = 0; i < N; i++)
	D.Fill(R());
    //calculate chi-square
    double S = 0;
    for (const auto &p : D)
        S += (p.Y()/step/N).NumCompare(F(p.X().val()) / norm);
    S /= D.size();
    EXPECT_TRUE(S < 2.5);
}
TEST(RandomValueTableDistr, Uniform)
{
    TestRandomDistribution([](double) {
        return 1.0;
    }, 0, 10, 20);
}
TEST(RandomValueTableDistr, Linear)
{
    TestRandomDistribution([](double x) {
        return x;
    }, 0, 10, 20);
}
TEST(RandomValueTableDistr, Parabolic)
{
    TestRandomDistribution([](double x) {
        return x * x;
    }, 0, 10, 20);
}
TEST(RandomValueTableDistr, Sin)
{
    TestRandomDistribution([](double x) {
        return sin(x);
    }, 0, 3, 30);
}
TEST(RandomValueTableDistr, Gauss)
{
    TestRandomDistribution([](double x) {
        return Gaussian(x, 5.0, 1.0);
    }, 0, 10, 50);
}
TEST(RandomValueTableDistr, Gauss2)
{
    TestRandomDistribution([](double x) {
        return Gaussian(x, 5.0, 1.5);
    }, 0, 10, 50);
}
TEST(RandomValueTableDistr, Gauss3)
{
    TestRandomDistribution([](double x) {
        return Gaussian(x, 5.0, 2.0);
    }, 0, 10, 50);
}
TEST(RandomValueTableDistr, Gauss4)
{
    TestRandomDistribution([](double x) {
        return Gaussian(x, 3.0, 1.0);
    }, -2, 10, 50);
}
TEST(RandomValueTableDistr, Gauss5)
{
    TestRandomDistribution([](double x) {
        return Gaussian(x, 3.0, 1.5);
    }, -2, 10, 50);
}
TEST(RandomValueTableDistr, Gauss6)
{
    TestRandomDistribution([](double x) {
        return Gaussian(x, 3.0, 2.0);
    }, -2, 10, 50);
}

TEST(RandomGauss, BaseTest)
{
    RandomGauss<> R(0, 1);
    StandardDeviation<> D,D2;
    for (int i = 0; i < 10000; i++)D<<R();
    ALMOST_EQ2(D.val(),0);
    ALMOST_EQ2(D.uncertainty(),1);
    auto R2 = R;
    for (int i = 0; i < 10000; i++)D2<<R2();
    ALMOST_EQ2(D2.val(),0);
    ALMOST_EQ2(D2.uncertainty(),1);
}
TEST(Poisson, BaseTest)
{
    RandomUniform<> L(1,20);
    for(size_t i=0;i<10;i++){
	const auto l=L();
	const auto R=Poisson(l);
	StandardDeviation<> D;
	for (int i = 0; i < 10000; i++)D<<R();
	cout<< l<< " ; " << D;
	ALMOST_EQ3(D.val(),l);
    }
}
