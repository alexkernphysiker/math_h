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
TEST(Pder1,test){
    const auto F=[](const Vector<>&r){return r.M_sqr();};
    EXPECT_EQ(0,num_Pder1(Zero(),X()*0.1)*F);
    EXPECT_EQ(0,num_Pder1(Zero(),Y()*0.1)*F);
    EXPECT_EQ(0,num_Pder1(Zero(),Z()*0.1)*F);
    const auto F2=[](const Vector<>&r){return r.x()+r.y()-r.z();};
    EXPECT_EQ( 1,num_Pder1(Zero(),X()*0.1)*F2);
    EXPECT_EQ( 1,num_Pder1(Zero(),Y()*0.1)*F2);
    EXPECT_EQ(-1,num_Pder1(Zero(),Z()*0.1)*F2);
}
TEST(nabla,test){
    const auto F=[](const Vector<>&r){return r.M_sqr();};
    const auto G=nabla(Zero(),0.1)*F;
    EXPECT_EQ(0,G.x());
    EXPECT_EQ(0,G.y());
    EXPECT_EQ(0,G.z());
    const auto F2=[](const Vector<>&r){return r.x()+r.y()-r.z();};
    const auto G2=nabla(Zero(),0.1)*F2;
    EXPECT_EQ( 1,G2.x());
    EXPECT_EQ( 1,G2.y());
    EXPECT_EQ(-1,G2.z());
    const auto F3=[](const Vector<>&){return 0.0;};
    const auto G3=nabla(Zero(),0.1)*F3;
    EXPECT_EQ(0,G3.x());
    EXPECT_EQ(0,G3.y());
    EXPECT_EQ(0,G3.z());
    std::function<double(const Vector<>&)> f1=F,f2=F2,f3=F3;
    const auto div=nabla(Zero(),0.1)*vec(f1,f2,f3);
    EXPECT_EQ(1,div);
}
