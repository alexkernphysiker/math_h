// this file is distributed under 
// MIT license
#include <gtest/gtest.h>
#include <vector>
#include <math.h>
#include <math_h/randomfunc.h>
using namespace std;
using namespace MathTemplates;
mt19937 rnd;
TEST(RandomValueTableDistr,BaseTest){
	RandomValueTableDistr<double> R([](double)->double{return 1;},ChainWithCount(10,0.0,1.0));
	for(int i=0;i<100;i++){
		double r=R(rnd);
		EXPECT_TRUE((r>=0)&&(r<=1));
	}
	auto R2=R;
	for(int i=0;i<100;i++){
		double r=R2(rnd);
		EXPECT_TRUE((r>=0)&&(r<=1));
	}
}
TEST(RandomValueTableDistr,Throwing){
	auto f=[](double){return 1;};
	auto z=[](double){return 0;};
	EXPECT_THROW(RandomValueTableDistr<double>(f,ChainWithCount(0,0.0,1.0)),Exception<SortedChain<double>>);
	EXPECT_NO_THROW(RandomValueTableDistr<double>(f,ChainWithCount(1,0.0,1.0)));
	EXPECT_THROW(RandomValueTableDistr<double>(f,ChainWithCount(2,0.0,-1.0)),Exception<SortedChain<double>>);
	EXPECT_THROW(RandomValueTableDistr<double>(f,ChainWithCount(2,0.0,0.0)),Exception<SortedChain<double>>);
	EXPECT_NO_THROW(RandomValueTableDistr<double>(f,ChainWithCount(2,0.0,1.0)));
	EXPECT_NO_THROW(RandomValueTableDistr<double>(z,ChainWithCount(2,0.0,1.0)));
	EXPECT_NO_THROW(RandomValueTableDistr<double>(f,{0.0,1.0}));
	EXPECT_NO_THROW(RandomValueTableDistr<double>({point<double>(0.0,1.0),point<double>(1.0,1.0)}));
	EXPECT_THROW(RandomValueTableDistr<double>(z,ChainWithCount(0,0.0,0.0)),Exception<SortedChain<double>>);
	auto n=[](double x){return sin(10*x);};
	EXPECT_THROW(RandomValueTableDistr<double>(n,ChainWithCount(100,0.0,1.0)),Exception<SortedPoints<double>>);
}
#include <math_h/hists.h>
void TestRandomDistribution(function<double(double)> F,double from,double to,int bins,int accu=10){
	RandomValueTableDistr<double> R(F,ChainWithCount(bins*accu,from,to));
	double step=(to-from)/double(bins);
	Distribution1D<double> D(BinsByStep(from,step,to));
	double norm=Sympson(F,from,to,0.0001);
	int N=500;
	for(int i=0;i<N;i++)D.Fill(R(rnd));
	double S=0;
	for(const auto&p:D){
		double n_exp=step*double(N);
		double d=p.Y().uncertainty();
		double y=p.Y().val()/n_exp;
		double f=F(p.X().val())/norm;
		S+=pow((y-f)/d,2);
	}
	S/=D.size();
	printf("chi^2=%f\n",S);
	EXPECT_TRUE(S<0.5);
}
TEST(RandomValueTableDistr,Uniform){TestRandomDistribution([](double){return 1.0;},0,10,20);}
TEST(RandomValueTableDistr,Linear){TestRandomDistribution([](double x){return x;},0,10,20);}
TEST(RandomValueTableDistr,Parabolic){TestRandomDistribution([](double x){return x*x;},0,10,20);}
TEST(RandomValueTableDistr,Sin){TestRandomDistribution([](double x){return sin(x);},0,3,30);}
TEST(RandomValueTableDistr,Gauss){TestRandomDistribution([](double x){return Gaussian(x,5.0,1.0);},0,10,50);}
TEST(RandomValueTableDistr,Gauss2){TestRandomDistribution([](double x){return Gaussian(x,5.0,1.5);},0,10,50);}
TEST(RandomValueTableDistr,Gauss3){TestRandomDistribution([](double x){return Gaussian(x,5.0,2.0);},0,10,50);}
TEST(RandomValueTableDistr,Gauss4){TestRandomDistribution([](double x){return Gaussian(x,3.0,1.0);},-2,10,50);}
TEST(RandomValueTableDistr,Gauss5){TestRandomDistribution([](double x){return Gaussian(x,3.0,1.5);},-2,10,50);}
TEST(RandomValueTableDistr,Gauss6){TestRandomDistribution([](double x){return Gaussian(x,3.0,2.0);},-2,10,50);}
