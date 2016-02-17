// this file is distributed under 
// MIT license
#include <gtest/gtest.h>
#include <vector>
#include <math.h>
#include <math_h/randomfunc.h>
#include <math_h/hist.h>
#include <math_h/functions.h>
using namespace std;
using namespace MathTemplates;
mt19937 rnd;
TEST(RandomValueGenerator,BaseTest){
	RandomValueGenerator<double> R([](double)->double{return 1;},ChainWithCount(10,0.0,1.0));
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
TEST(RandomValueGenerator,Throwing){
	auto f=[](double){return 1;};
	auto z=[](double){return 0;};
	EXPECT_THROW(RandomValueGenerator<double>(f,ChainWithCount(0,0.0,1.0)),Exception<vector<double>>);
	EXPECT_NO_THROW(RandomValueGenerator<double>(f,ChainWithCount(1,0.0,1.0)));
	EXPECT_THROW(RandomValueGenerator<double>(f,ChainWithCount(2,0.0,-1.0)),Exception<vector<double>>);
	EXPECT_THROW(RandomValueGenerator<double>(f,ChainWithCount(2,0.0,0.0)),Exception<vector<double>>);
	EXPECT_NO_THROW(RandomValueGenerator<double>(f,ChainWithCount(2,0.0,1.0)));
	EXPECT_NO_THROW(RandomValueGenerator<double>(z,ChainWithCount(2,0.0,1.0)));
	EXPECT_THROW(RandomValueGenerator<double>(z,ChainWithCount(0,0.0,0.0)),Exception<vector<double>>);
	auto n=[](double x){return sin(10*x);};
	EXPECT_THROW(RandomValueGenerator<double>(n,ChainWithCount(100,0.0,1.0)),Exception<LinearInterpolation<double>>);
}
void TestRandomDistribution(function<double(double)> F,double from,double to,int bins,int accu=10){
	RandomValueGenerator<double> R(F,ChainWithCount(bins*accu,from,to));
	double step=(to-from)/double(bins);
	Distribution1D<double> D(BinsByStep(from,step,to));
	double norm=Sympson(F,from,to,0.0001);
	int N=5000;
	for(int i=0;i<N;i++)D<<R(rnd);
	double S=0;
	for(const auto&p:D){
		double n_exp=step*double(N);
		double d=p.Y().delta();
		double y=p.Y().val()/n_exp;
		double f=F(p.X().val())/norm;
		S+=pow((y-f)/d,2);
	}
	S/=D.size();
	printf("chi^2=%f\n",S);
	EXPECT_TRUE(S<1.5);
}
TEST(RandomValueGenerator,Uniform){TestRandomDistribution([](double){return 1.0;},0,10,20);}
TEST(RandomValueGenerator,Linear){TestRandomDistribution([](double x){return x;},0,10,20);}
TEST(RandomValueGenerator,Parabolic){TestRandomDistribution([](double x){return x*x;},0,10,20);}
TEST(RandomValueGenerator,Sin){TestRandomDistribution([](double x){return sin(x);},0,3,30);}
TEST(RandomValueGenerator,Gauss){TestRandomDistribution([](double x){return Gaussian(x,5.0,1.0);},0,10,50);}
TEST(RandomValueGenerator,Gauss2){TestRandomDistribution([](double x){return Gaussian(x,5.0,1.5);},0,10,50);}
TEST(RandomValueGenerator,Gauss3){TestRandomDistribution([](double x){return Gaussian(x,5.0,2.0);},0,10,50);}
TEST(RandomValueGenerator,Gauss4){TestRandomDistribution([](double x){return Gaussian(x,3.0,1.0);},-2,10,50);}
TEST(RandomValueGenerator,Gauss5){TestRandomDistribution([](double x){return Gaussian(x,3.0,1.5);},-2,10,50);}
TEST(RandomValueGenerator,Gauss6){TestRandomDistribution([](double x){return Gaussian(x,3.0,2.0);},-2,10,50);}
