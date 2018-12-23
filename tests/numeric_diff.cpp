// this file is distributed under
// LGPLv3 license
#include <gtest/gtest.h>
#include <math_h/numeric_diff.h>
//this file contains unit tests for sigma.h
using namespace std;
using namespace MathTemplates;
#define ALMOST_EQ(a,b) EXPECT_TRUE(abs(a-b)<0.0001)
#define ALMOST_EQ2(a,b) EXPECT_TRUE(abs(a-b)<0.01)
#define ALMOST_EQ3(a,b) EXPECT_TRUE(abs(a-b)<0.1)

TEST(der1,test){
    EXPECT_EQ(0,der1<>(0.1)([](const double&){return 5.0;},3.0));
    EXPECT_EQ(2,der1<>(0.1)([](const double&x){return 2.0*x;},3.0));
    EXPECT_EQ(0,der1<>(0.1)([](const double&x){return 2.0*x*x;},0.0));
}
TEST(der2,test){
    EXPECT_EQ(0,der2<>(0.1)([](const double&){return 5.0;},3.0));
    EXPECT_EQ(0,der2<>(0.1)([](const double&x){return 2.0*x;},3.0));
    EXPECT_EQ(2.0,der2<>(0.1)([](const double&x){return x*x;},0.0));
}
