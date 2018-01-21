// this file is distributed under
// LGPLv3 license
#include <gtest/gtest.h>
#include <math_h/tabledata.h>
#include <math_h/sigma.h>
#include <math_h/randomfunc.h>
using namespace std;
using namespace MathTemplates;
double TestArray[] = { -2, -1, 0, 1, 2, 3, 4, 5, 6, 7};
TEST(WhereToInsert, BorderConditions)
{
    for (int beg = 0; beg < 10; beg++)
        for (double V = TestArray[beg] - 2; V <= TestArray[beg] + 2; V += 0.5) {
            int index = 0;
            if (beg > 0) {
                index = details::WhereToInsert(beg, beg - 1, TestArray, V);
                EXPECT_EQ(beg, index);
            }
            index = details::WhereToInsert(beg, beg, TestArray, V);
            if (V < TestArray[beg]) {
                EXPECT_EQ(beg, index);
            }
            if (V > TestArray[beg]) {
                EXPECT_EQ(beg + 1, index);
            }
            if (V == TestArray[beg]) {
                EXPECT_TRUE((index == beg) || (index == (beg + 1)));
            }
        }
}
TEST(WhereToInsert, NormalConditions)
{
    for (double x = TestArray[0] - 0.5; x <= TestArray[9] + 0.5; x += 0.5)
        for (int beg = 0; beg < 10; beg++)for (int end = beg + 1; end < 10; end++) {
                int index = details::WhereToInsert(beg, end, TestArray, x);
                if (x < TestArray[beg])EXPECT_EQ(beg, index);
                else if (x > TestArray[end])EXPECT_EQ(end + 1, index);
                else EXPECT_TRUE((x >= TestArray[index - 1]) && (x <= TestArray[index]));
            }
}
TEST(InsertSorted, BasicTest)
{
    Chain<int> X;
    for (int i = 0; i < 50; i++) {
        details::InsertSorted(rand() % 10, X,std_size(X),std_insert(X,int));
        for (int j = 0; j < i; j++)
            EXPECT_TRUE(X[j] <= X[j + 1]);
    }
}
TEST(SortedChain, basetest)
{
    RANDOM engine;
    SortedChain<> chain;
    RandomUniform<> get(-10, 10);
    for (size_t i = 0; i < 100; i++)
        chain << get(engine);
    for (size_t i = 1; i < 100; i++)
        EXPECT_TRUE(chain[i] > chain[i - 1]);
}
TEST(point, basetest)
{
    double x(25), y(78);
    point<double> p(x, y);
    EXPECT_EQ(x, p.X());
    EXPECT_EQ(y, p.Y());
}
TEST(point, basetest1)
{
    double x(25), y(78);
    point<> p = {x, y};
    EXPECT_EQ(x, p.X());
    EXPECT_EQ(y, p.Y());
}
TEST(point, basetest2)
{
    int x(25);
    double y(78);
    point<int, double> p = {25, 78};
    EXPECT_EQ(x, p.X());
    EXPECT_EQ(y, p.Y());
    point<int, double> p2 = {25., 78.};
    EXPECT_EQ(x, p2.X());
    EXPECT_EQ(y, p2.Y());
}
TEST(point3d, basetest)
{
    double x(25), y(78), z(55);
    point3d<double> p(x, y, z);
    EXPECT_EQ(x, p.X());
    EXPECT_EQ(y, p.Y());
    EXPECT_EQ(z, p.Z());
}
TEST(point, basetest_value)
{
    value<double> x(25), y(78);
    point<value<double>> p(x, y);
    EXPECT_EQ(x.val(), p.X().val());
    EXPECT_EQ(x.uncertainty(), p.X().uncertainty());
    EXPECT_EQ(y.val(), p.Y().val());
    EXPECT_EQ(y.uncertainty(), p.Y().uncertainty());
}
TEST(point3d, basetest_value)
{
    value<double> x(25), y(78), z(55);
    point3d<value<double>> p(x, y, z);
    EXPECT_EQ(x.val(), p.X().val());
    EXPECT_EQ(x.uncertainty(), p.X().uncertainty());
    EXPECT_EQ(y.val(), p.Y().val());
    EXPECT_EQ(y.uncertainty(), p.Y().uncertainty());
    EXPECT_EQ(z.val(), p.Z().val());
    EXPECT_EQ(z.uncertainty(), p.Z().uncertainty());
}
TEST(SortedPoints, size1)
{
    SortedPoints<double> chain;
    EXPECT_EQ(0, chain.size());
    chain << point<double>(0, 1);
    EXPECT_EQ(1, chain.size());
    chain << point<double>(2, 1);
    EXPECT_EQ(2, chain.size());
    chain << point<double>(1, 1);
    EXPECT_EQ(3, chain.size());
    chain << point<double>(3, 1);
    EXPECT_EQ(4, chain.size());
    chain << point<double>(5, 1);
    EXPECT_EQ(5, chain.size());
    chain << point<double>(4, 1);
    EXPECT_EQ(6, chain.size());
    EXPECT_EQ(0, chain.left().X());
    EXPECT_EQ(5, chain.right().X());
    for (size_t i = 0; i < chain.size(); i++) {
        EXPECT_EQ(1, chain[i].Y());
        EXPECT_EQ(i, chain[i].X());
    }
}
TEST(SortedPoints, range)
{
    SortedPoints<> chain=Points<>{{0, 1},{2, 1.1},{3, 1.2},{4, 1.3},{5, 1.4}};
    auto xcut = chain.XRange(0.5, 4.5);
    EXPECT_EQ(3, xcut.size());
    EXPECT_EQ(2, xcut.left().X());
    EXPECT_EQ(4, xcut.right().X());
    auto ycut = chain.YRange(1.15, 1.35).TransponateAndSort();
    EXPECT_EQ(2, ycut.size());
    EXPECT_EQ(1.2, ycut.left().X());
    EXPECT_EQ(1.3, ycut.right().X());
}
TEST(SortedPoints, size2)
{
    SortedPoints<value<>> chain;
    EXPECT_EQ(0, chain.size());
    chain << make_point(0, 1);
    EXPECT_EQ(1, chain.size());
    chain << make_point(2, 1);
    EXPECT_EQ(2, chain.size());
    chain << make_point(1, 1);
    EXPECT_EQ(3, chain.size());
    chain << make_point(3, 1);
    EXPECT_EQ(4, chain.size());
    chain << make_point(5, 1);
    EXPECT_EQ(5, chain.size());
    chain << make_point(4, 1);
    EXPECT_EQ(6, chain.size());
    EXPECT_EQ(0, chain.left().X().val());
    EXPECT_EQ(5, chain.right().X().val());
    for (size_t i = 0; i < chain.size(); i++) {
        EXPECT_EQ(1, chain[i].Y().val());
        EXPECT_EQ(i, chain[i].X().val());
    }
}
TEST(SortedPoints, WeightedAvr)
{
    SortedPoints<double,WeightedAverage<>> chain;
    chain
    << make_point(0, WeightedAverage<>())
    << make_point(1, WeightedAverage<>())
    << make_point(2, WeightedAverage<>());
    for (size_t i = 0; i < chain.size(); i++) {
        EXPECT_EQ(i, chain[i].X());
        EXPECT_THROW(chain[i].Y()(),Exception<WeightedAverage<>>);
    }
    SortedPoints<double,value<>> opr1=Points<double,value<>>{{0,{0,1}},{1,{1,1}},{2,{2,1}}};
    chain.leftArrow(opr1);
    for (size_t i = 0; i < chain.size(); i++) {
        EXPECT_EQ(i, chain[i].X());
        EXPECT_EQ(i, chain[i].Y()().val());
    }
    chain.leftArrow(Points<double,value<>>{{0,{0,1}},{1,{0,1}},{2,{0,1}}});
    for (size_t i = 0; i < chain.size(); i++) {
        EXPECT_EQ(i, chain[i].X());
        EXPECT_EQ(double(i)/2, chain[i].Y()().val());
    }
#ifdef ____full_version_of_tabledata_h_____
    const auto chain2=chain+2;
    for (size_t i = 0; i < chain2.size(); i++) {
        EXPECT_EQ(i, chain2[i].X());
        EXPECT_EQ(double(i)/2 + 2, chain2[i].Y().val());
    }
    const auto chain3=chain-2;
    for (size_t i = 0; i < chain3.size(); i++) {
        EXPECT_EQ(i, chain3[i].X());
        EXPECT_EQ(double(i)/2 - 2, chain3[i].Y().val());
    }
    const auto chain4=chain*2;
    for (size_t i = 0; i < chain4.size(); i++) {
        EXPECT_EQ(i, chain4[i].X());
        EXPECT_EQ(double(i), chain4[i].Y().val());
    }
    const auto chain5=chain/2;
    for (size_t i = 0; i < chain5.size(); i++) {
        EXPECT_EQ(i, chain5[i].X());
        EXPECT_EQ(double(i)/4, chain5[i].Y().val());
    }
#endif
}
TEST(SortedPoints, WeightedAvr2)
{
    SortedPoints<double,WeightedAverage<>> chain=Points<double,value<>>{{0,{0,1}},{1,{1,1}},{2,{2,1}}};
    for (size_t i = 0; i < chain.size(); i++) {
        EXPECT_EQ(i, chain[i].X());
        EXPECT_EQ(i, chain[i].Y()().val());
    }
    chain.leftArrow(Points<double,value<>>{{0,{0,1}},{1,{0,1}},{2,{0,1}}});
    for (size_t i = 0; i < chain.size(); i++) {
        EXPECT_EQ(i, chain[i].X());
        EXPECT_EQ(double(i)/2, chain[i].Y()().val());
    }
}
