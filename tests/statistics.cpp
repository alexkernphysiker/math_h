// this file is distributed under
// LGPLv3 license
#include <gtest/gtest.h>
#include <math_h/randomfunc.h>
#include <math_h/statistics.h>
#include <math_h/vectortransformations.h>
//this file contains unit tests for statistics.h
using namespace std;
using namespace MathTemplates;
#define ALMOST_EQ(a,b) EXPECT_TRUE(abs(a-b)<0.1)
TEST(Statistics,AverageObtainer)
{
    AverageObtainer<double> S;
    EXPECT_EQ(S.count(),0);
    S.Fill(2);
    S.Fill(3);
    EXPECT_EQ(S.Average(),2.5);
    EXPECT_EQ(S.count(),2);
    S.Fill(2);
    S.Fill(3);
    EXPECT_EQ(S.Average(),2.5);
    EXPECT_EQ(S.count(),4);
}
class TestStatistics:public AverageObtainer<double>{
public:
    TestStatistics():AverageObtainer<double>(){}
    virtual ~TestStatistics(){}
    void Test()const{
	size_t c=0;
	double res=0;
	ForEach([&c,&res](double x){
	    c++;
	    res+=x;
	});
	EXPECT_EQ(c,count());
	EXPECT_EQ(res,Average()*c);
    }
};
TEST(Statistics,abstract){
    TestStatistics S;
    S.Fill(1);S.Test();
    S.Fill(5);S.Test();
    S.Fill(7);S.Test();
    S.Fill(4);S.Test();
    S.Fill(6);S.Test();
    S.Fill(2);S.Test();
    S.Fill(7);S.Test();
    S.Fill(3);S.Test();
}
TEST(Statistics,Sampling){
    Sampling<2> S;
    RandomGauss<> X(0,1),Y(0,2);
    for(size_t i=0;i<100000;i++){
	S.Fill(X(),Y());
    }
    const auto&s11=S.Cov().element<1,1>();
    const auto&s12=S.Cov().element<1,2>();
    const auto&s21=S.Cov().element<2,1>();
    const auto&s22=S.Cov().element<2,2>();
    ALMOST_EQ(s11,1);
    ALMOST_EQ(s22,4);
    EXPECT_EQ(s12,s21);
    ALMOST_EQ(s12,0);
}
TEST(Statistics,SamplingXY){
    {
    SamplingXY<2,3> S;
    RandomGauss<> M(0,2);
    for(size_t i=0;i<100000;i++){
	S.Fill(randomIsotropic<2>()*M(),randomIsotropic<3>()*M());
    }
    EXPECT_EQ(S.CovX(),S.CovX().transponate());
    EXPECT_EQ(S.CovY(),S.CovY().transponate());
    const auto&s11=S.CovXY().element<1,1>();
    const auto&s12=S.CovXY().element<1,2>();
    const auto&s13=S.CovXY().element<1,3>();
    const auto&s21=S.CovXY().element<2,1>();
    const auto&s22=S.CovXY().element<2,2>();
    const auto&s23=S.CovXY().element<2,3>();
    ALMOST_EQ(s11,0);
    ALMOST_EQ(s12,0);
    ALMOST_EQ(s13,0);
    ALMOST_EQ(s21,0);
    ALMOST_EQ(s22,0);
    ALMOST_EQ(s23,0);
    }
    {
    SamplingXY<2,2> S;
    RandomGauss<> M(0,2);
    for(size_t i=0;i<100000;i++){
	const auto x=randomIsotropic<2>()*M();
	S.Fill(x,x);
    }
    const auto&s11=S.CovX().element<1,1>();
    const auto&s12=S.CovX().element<1,2>();
    const auto&s21=S.CovX().element<2,1>();
    const auto&s22=S.CovX().element<2,2>();
    ALMOST_EQ(s11,2);
    ALMOST_EQ(s22,2);
    EXPECT_EQ(s12,s21);
    ALMOST_EQ(s12,0);
    EXPECT_EQ(S.CovX(),S.CovY());
    EXPECT_EQ(S.CovX(),S.CovXY());
    }
}
