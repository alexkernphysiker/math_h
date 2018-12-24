// this file is distributed under
// LGPLv3 license
#include <gtest/gtest.h>
#include <math_h/sigma.h>
#include <math_h/randomfunc.h>
//this file contains unit tests for sigma.h
using namespace std;
using namespace MathTemplates;
#define ALMOST_EQ(a,b) EXPECT_TRUE(abs(a-b)<0.0001)
#define ALMOST_EQ2(a,b) EXPECT_TRUE(abs(a-b)<0.01)
#define ALMOST_EQ3(a,b) EXPECT_TRUE(abs(a-b)<0.1)
TEST(value, properties)
{
    RandomGauss<> G(0,6);
    for (size_t i = 0; i < 10; i++) {
        double x = G();
        value<> V(x, 0.1);
        EXPECT_EQ(x, V.val());
        EXPECT_EQ(0.1, V.uncertainty());
        EXPECT_EQ(0.1 / x, V.epsilon());
        EXPECT_EQ(x - 0.1, V.min());
        EXPECT_EQ(x + 0.1, V.max());
    }
}
TEST(value,initlist)
{
    RandomGauss<> G(0,6);
    for (size_t i = 0; i < 10; i++) {
        double x = G();
        value<> V{x, 0.1};
        EXPECT_EQ(x, V.val());
        EXPECT_EQ(0.1, V.uncertainty());
        EXPECT_EQ(0.1 / x, V.epsilon());
        EXPECT_EQ(x - 0.1, V.min());
        EXPECT_EQ(x + 0.1, V.max());
    }
}
TEST(value,compare)
{
    value<> V(0, 0.1);
    EXPECT_EQ(false, V.Contains(-0.2));
    EXPECT_EQ(true, V.Contains(-0.05));
    EXPECT_EQ(true, V.Contains(0));
    EXPECT_EQ(true, V.Contains(0.05));
    EXPECT_EQ(false, V.Contains(0.2));

    EXPECT_EQ(true, V.NotEqual(-0.2));
    EXPECT_EQ(false, V.NotEqual(-0.05));
    EXPECT_EQ(false, V.NotEqual(0));
    EXPECT_EQ(false, V.NotEqual(0.05));
    EXPECT_EQ(true, V.NotEqual(0.2));

    bool v = value<>(2, 1).Above(value<>(0, 0.9));
    EXPECT_TRUE(v);
    v = value<>(2, 1).Above(value<>(0, 1.1));
    EXPECT_FALSE(v);
    v = value<>(0, 1).Below(value<>(2, 0.9));
    EXPECT_TRUE(v);
    v = value<>(0, 1).Below(value<>(2, 1.1));
    EXPECT_FALSE(v);

    v = value<>(2, 1) .Above(0.9);
    EXPECT_TRUE(v);
    v = value<>(2, 1) .Above(1.1);
    EXPECT_FALSE(v);
    v = value<>(0, 1).Below(1.1);
    EXPECT_TRUE(v);
    v = value<>(0, 1).Below(0.9);
    EXPECT_FALSE(v);
}
TEST(value, arithmetic_actions)
{
    {
    RandomUniform<> val(1, 50), unc(0.1, 10);
    for (size_t cnt = 0; cnt < 1000; cnt++) {
        value<> A(val(), unc()), B(val(), unc());
        auto sum = A + B;
        EXPECT_EQ(A.val() + B.val(), sum.val());
        ALMOST_EQ(pow(A.uncertainty(), 2) + pow(B.uncertainty(), 2), pow(sum.uncertainty(), 2));
        auto sub = A - B;
        EXPECT_EQ(A.val() - B.val(), sub.val());
        ALMOST_EQ(pow(A.uncertainty(), 2) + pow(B.uncertainty(), 2), pow(sub.uncertainty(), 2));
        auto prod = A * B;
        EXPECT_EQ(A.val()*B.val(), prod.val());
        ALMOST_EQ(pow(A.uncertainty()*B.val(), 2) + pow(A.val()*B.uncertainty(), 2), pow(prod.uncertainty(), 2));
        auto ratio = A / B;
        EXPECT_EQ(A.val() / B.val(), ratio.val());
        ALMOST_EQ(pow(A.uncertainty() / B.val(), 2) + pow((A.val()*B.uncertainty()) / pow(B.val(), 2), 2), pow(ratio.uncertainty(), 2));
    }
    }
    {
    RandomUniform<> val(1, 50), unc(0.1, 10);
    for (size_t cnt = 0; cnt < 1000; cnt++) {
        value<> A(val(), unc());
        double B = val();
        auto sum = A + B;
        EXPECT_EQ(A.val() + B, sum.val());
        EXPECT_EQ(A.uncertainty(), sum.uncertainty());
        auto sub = A - B;
        EXPECT_EQ(A.val() - B, sub.val());
        EXPECT_EQ(A.uncertainty(), sub.uncertainty());
        auto prod = A * B;
        EXPECT_EQ(A.val()*B, prod.val());
        EXPECT_EQ(A.uncertainty()*B, prod.uncertainty());
        auto ratio = A / B;
        EXPECT_EQ(A.val() / B, ratio.val());
        EXPECT_EQ(A.uncertainty() / B, ratio.uncertainty());
    }
    }
}
TEST(value, NumCompare)
{
    {
        const double v = value<>(0, 1).NumCompare(0);
        EXPECT_EQ(v, 0);
    }
    {
        const double v = value<>(0, 1).NumCompare(0.5);
        EXPECT_EQ(v, 0.25);
    }
    {
        const double v = value<>(0, 1).NumCompare(1);
        EXPECT_EQ(v, 1);
    }
    {
        const double v = value<>(0, 1).NumCompare(2);
        EXPECT_EQ(v, 4);
    }
    {
        const double v = value<>(2, 1).NumCompare(2);
        EXPECT_EQ(v, 0);
    }
    {
        const double v = value<>(2, 1).NumCompare(1.5);
        EXPECT_EQ(v, 0.25);
    }
    {
        const double v = value<>(2, 1).NumCompare(1);
        EXPECT_EQ(v, 1);
    }
    {
        const double v = value<>(2, 1).NumCompare(0);
        EXPECT_EQ(v, 4);
    }
    {
        const double v = value<>(0, 1).NumCompare(value<>(0, 1));
        EXPECT_EQ(v, 0);
    }
    {
        const double v = value<>(0, 1).NumCompare(value<>(1, 1));
        EXPECT_EQ(v, 0.25);
    }
    {
        const double v = value<>(0, 1).NumCompare(value<>(2, 1));
        EXPECT_EQ(v, 1);
    }
    {
        const double v = value<>(0, 1).NumCompare(value<>(4, 1));
        EXPECT_EQ(v, 4);
    }
    {
        const double v = value<>(2, 1).NumCompare(value<>(2, 1));
        EXPECT_EQ(v, 0);
    }
    {
        const double v = value<>(2, 1).NumCompare(value<>(1, 1));
        EXPECT_EQ(v, 0.25);
    }
    {
        const double v = value<>(2, 1).NumCompare(value<>(0, 1));
        EXPECT_EQ(v, 1);
    }
    {
        const double v = value<>(2, 1).NumCompare(value<>(-2, 1));
        EXPECT_EQ(v, 4);
    }
}
TEST(value, func_val1)
{
     value<> A(1, 0.1),
     V=func_with_uncertainty([](double x){return 2.0 * x;},A);
     EXPECT_EQ(2 * A.val(), V.val());
     ALMOST_EQ(2 * A.uncertainty(), V.uncertainty());
}
TEST(value, func_val2)
{
     value<> A1(1, 0.1),A2(2, 0.5),
     V=func_with_uncertainty([](double x,double y){return x+y;},A1,A2);
     EXPECT_EQ((A1+A2).val(), V.val());
     ALMOST_EQ((A1+A2).uncertainty(), V.uncertainty());
}
TEST(value, func_val3)
{
     value<> A1(1, 0.1),A2(2, 0.5),A3(1.5, 0.2),
     V=func_with_uncertainty([](double x,double y,double z){return (x+y)*z;},A1,A2,A3);
     EXPECT_EQ(((A1+A2)*A3).val(), V.val());
     ALMOST_EQ(((A1+A2)*A3).uncertainty(), V.uncertainty());
}
TEST(value, wider)
{
    value<> V1(1, 0.1), V2(2, 0.1), V3(1, 0.5);
    EXPECT_EQ(V1, V1.make_wider(1));
    EXPECT_EQ(V2, V2.make_wider(1));
    EXPECT_EQ(V3, V3.make_wider(1));
    auto V11 = V1.make_wider(0.5);
    auto V12 = V2.make_wider(0.5);
    auto V13 = V3.make_wider(0.5);
    EXPECT_EQ(V11.val(), V1.val());
    EXPECT_EQ(V12.val(), V2.val());
    EXPECT_EQ(V13.val(), V3.val());
    EXPECT_EQ(V11.uncertainty(), V1.uncertainty() * 0.5);
    EXPECT_EQ(V12.uncertainty(), V2.uncertainty() * 0.5);
    EXPECT_EQ(V13.uncertainty(), V3.uncertainty() * 0.5);
    auto V21 = V1.make_wider(2.5);
    auto V22 = V2.make_wider(2.5);
    auto V23 = V3.make_wider(2.5);
    EXPECT_EQ(V21.val(), V1.val());
    EXPECT_EQ(V22.val(), V2.val());
    EXPECT_EQ(V23.val(), V3.val());
    EXPECT_EQ(V21.uncertainty(), V1.uncertainty() * 2.5);
    EXPECT_EQ(V22.uncertainty(), V2.uncertainty() * 2.5);
    EXPECT_EQ(V23.uncertainty(), V3.uncertainty() * 2.5);
}
TEST(value, error_handling)
{
    EXPECT_NO_THROW(value<>(1, -0.1));
    EXPECT_NO_THROW(value<>(1, 0));
    EXPECT_NO_THROW(value<>(1, 0.1));
    EXPECT_EQ(value<>(1, 0).uncertainty(), 0);
    EXPECT_EQ(value<>(1, 0.1).uncertainty(), 0.1);
    EXPECT_FALSE(isfinite(value<>(1, -0.1).uncertainty()));
    EXPECT_NO_THROW(value<>(1, -0.1) + value<>(1, 0.1));
    EXPECT_EQ((value<>(1, -0.1) + value<>(1, 0.1)).val(), 2);
    EXPECT_FALSE(isfinite((value<>(1, -0.1) + value<>(1, 0.1)).uncertainty()));
}

TEST(StandardDeviation, Throwing)
{
    StandardDeviation<> S;
    EXPECT_EQ(0, S.Sample().count());
    EXPECT_ANY_THROW(S.val());
    EXPECT_ANY_THROW(S.uncertainty());
    EXPECT_EQ(&S, &(S << 0.0));
    EXPECT_EQ(1, S.Sample().count());
    EXPECT_EQ(S.val(),0.0);
    EXPECT_ANY_THROW(S.uncertainty());
    EXPECT_EQ(&S, &(S << 0.0));
    EXPECT_EQ(2, S.Sample().count());
    EXPECT_EQ(0, S.val());
    EXPECT_EQ(0, S.uncertainty());
}
TEST(StandardDeviation, Base)
{
    StandardDeviation<> S;
    EXPECT_ANY_THROW(S.val());
    EXPECT_ANY_THROW(S.uncertainty());
    S << 0.0;
    EXPECT_EQ(1, S.Sample().count());
    EXPECT_EQ(S.val(),0.0);
    EXPECT_ANY_THROW(S.uncertainty());
    S << 1.0;
    EXPECT_EQ(2, S.Sample().count());
    EXPECT_EQ(0.5, S.val());
    EXPECT_EQ(sqrt(0.5), S.uncertainty());
}
TEST(StandardDeviation, Base2)
{
    StandardDeviation<> S;
    S << 0.0 << 1.0;
    size_t cnt=0;
    for(const auto&x:S.Sample())cnt++;
    EXPECT_EQ(2,cnt);
    EXPECT_EQ(2,S.Sample().size());
    EXPECT_EQ(2,S.Sample().count());
    EXPECT_EQ(0.0,S.Sample()[0].x());
    EXPECT_EQ(1.0,S.Sample()[1].x());
    EXPECT_ANY_THROW(S.Sample()[2].x());
}
TEST(StandardDeviation, WithRandomValues)
{
    StandardDeviation<> S;
    RandomGauss<> generator(1.0, 3.0);
    for (int i = 0; i < 200000; i++)
        S << generator();
    ALMOST_EQ3(1.0, S.val());
    ALMOST_EQ3(3.0, S.uncertainty());
}
TEST(StandardDeviation, Constructor2)
{
    StandardDeviation<> S(2);
    S << 0;
    ALMOST_EQ(1.0, S.val());
}

TEST(WeightedAverage, Zeros)
{
    WeightedAverage<> W;
    EXPECT_THROW(W.val(), Exception<WeightedAverage<>>);
    EXPECT_THROW(W.uncertainty(), Exception<WeightedAverage<>>);
    EXPECT_THROW(W << value<>(0, 0), Exception<WeightedAverage<>>);
    EXPECT_EQ(&W, &(W << value<>(0, 1)));
    ALMOST_EQ(0, W.val());
    ALMOST_EQ(1, W.uncertainty());
    EXPECT_EQ(&W, &(W << value<>(0, 1)));
    ALMOST_EQ(0, W.val());
    ALMOST_EQ(1.0 / sqrt(2.0), W.uncertainty());
    EXPECT_EQ(&W, &(W << value<>(0, 1)));
    ALMOST_EQ(0, W.val());
    ALMOST_EQ(1.0 / sqrt(3.0), W.uncertainty());
}
TEST(WeightedAverage, Ones)
{
    WeightedAverage<> W;
    EXPECT_THROW(W.val(), Exception<WeightedAverage<>>);
    EXPECT_THROW(W.uncertainty(), Exception<WeightedAverage<>>);
    EXPECT_THROW(W << value<>(1, 0), Exception<WeightedAverage<>>);
    EXPECT_EQ(&W, &(W << value<>(1, 1)));
    ALMOST_EQ(1, W.val());
    ALMOST_EQ(1, W.uncertainty());
    EXPECT_EQ(&W, &(W << value<>(1, 1)));
    ALMOST_EQ(1, W.val());
    ALMOST_EQ(1.0 / sqrt(2.0), W.uncertainty());
    EXPECT_EQ(&W, &(W << value<>(1, 1)));
    ALMOST_EQ(1, W.val());
    ALMOST_EQ(1.0 / sqrt(3.0), W.uncertainty());
}
TEST(WeightedAverage, Zeros_plus_Ones)
{
    WeightedAverage<> W;
    EXPECT_THROW(W.val(), Exception<WeightedAverage<>>);
    EXPECT_THROW(W.uncertainty(), Exception<WeightedAverage<>>);
    EXPECT_EQ(&W, &(W << value<>(1, 1)));
    ALMOST_EQ(1, W.val());
    ALMOST_EQ(1, W.uncertainty());
    EXPECT_EQ(&W, &(W << value<>(0, 1)));
    ALMOST_EQ(0.5, W.val());
    ALMOST_EQ(1.0 / sqrt(2.0), W.uncertainty());
    EXPECT_EQ(&W, &(W << value<>(1, 1)));
    ALMOST_EQ(2.0 / 3.0, W.val());
    ALMOST_EQ(1.0 / sqrt(3.0), W.uncertainty());
    EXPECT_EQ(&W, &(W << value<>(0, 1)));
    ALMOST_EQ(0.5, W.val());
    ALMOST_EQ(1.0 / sqrt(4.0), W.uncertainty());
}
