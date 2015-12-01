// this file is distributed under 
// MIT license
#include <gtest/gtest.h>
#include <vector>
#include <math.h>
#include <randomfunc.h>
#include <interpolate.h>
#include <sigma.h>
#include <functions.h>
using namespace std;
mt19937 rnd;
TEST(RandomValueGenerator,BaseTest){
	RandomValueGenerator<double> R([](double){return 1;},0,1,10);
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
	EXPECT_THROW(RandomValueGenerator<double>(f,0,1,0),math_h_error<LinearInterpolation_fixedsize<double>>);
	EXPECT_THROW(RandomValueGenerator<double>(f,0,1,-1),math_h_error<LinearInterpolation_fixedsize<double>>);
	EXPECT_THROW(RandomValueGenerator<double>(f,0,1,1),math_h_error<LinearInterpolation_fixedsize<double>>);
	EXPECT_THROW(RandomValueGenerator<double>(f,0,-1,2),math_h_error<LinearInterpolation_fixedsize<double>>);
	EXPECT_THROW(RandomValueGenerator<double>(f,0,0,2),math_h_error<LinearInterpolation_fixedsize<double>>);
	EXPECT_NO_THROW(RandomValueGenerator<double>(f,0,1,2));
	EXPECT_THROW(RandomValueGenerator<double>(z,0,1,2),math_h_error<RandomValueGenerator<double>>);
	EXPECT_THROW(RandomValueGenerator<double>(z,0,0,0),math_h_error<LinearInterpolation_fixedsize<double>>);
	auto n=[](double x){return sin(10*x);};
	EXPECT_THROW(RandomValueGenerator<double>(n,0,1,100),math_h_error<RandomValueGenerator<double>>);
}
void TestRandomDistribution(function<double(double)> F,double from,double to,int bins,int accu=10){
	RandomValueGenerator<double> R(F,from,to,bins*accu);
	Distribution<double> D(from,to,bins);
	double norm=Sympson(F,from,to,0.0001);
	int N=5000;
	for(int i=0;i<N;i++)D.AddValue(R(rnd));
	double S=0;
	for(int i=0,n=D.size();i<n;i++){
		double n_exp=D.BinWidth()*double(N);
		double d=sqrt(D.getY(i));if(d<1)d=1;d/=n_exp;
		double y=D.getY(i)/n_exp;
		double f=F(D.getX(i))/norm;
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
