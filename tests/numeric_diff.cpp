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
    EXPECT_EQ(0,num_der1(3.0,0.1)*[](const double&)->double{return 5.0;});
    EXPECT_EQ(2,num_der1(3.0,0.1)*[](const double&x)->double{return 2.0*x;});
    EXPECT_EQ(0,num_der1(0.0,0.1)*[](const double&x)->double{return 2.0*x*x;});
}
TEST(der2,test){
    EXPECT_EQ(0,num_der2(3.0,0.1)*[](const double&)->double{return 5.0;});
    EXPECT_EQ(0,num_der2(3.0,0.1)*[](const double&x)->double{return 2.0*x;});
    EXPECT_EQ(2.0,num_der2(0.0,0.1)*[](const double&x)->double{return x*x;});
}
// TEST(nabla,test){
//     const auto F=[](const Vector<>&x){return x.M_sqr();};
//     const auto G=nabla<>(Zero(),0.1);
//     EXPECT_EQ(0,(G*F).M());
// }
