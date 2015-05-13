#include <gtest/gtest.h>
#include <vector>
#include <math.h>
#include <randomfunc.h>
#include <interpolate.h>
#include <sigma.h>
#include <functions.h>
using namespace std;
TEST(RandomUniformlyI,BasicTest){
	for(int l=-10;l<=10;l++)
		for(int r=l;r<=10;r++)
			for(int c=0;c<100;c++){
				int V=RandomUniformlyI(l,r);
				EXPECT_TRUE((V>=l)&&(V<=r));
			}
}
TEST(RandomUniformlyI,Throwing){
	for(int l=-10;l<=10;l++)
		for(int r=-10;r<l;r++)
			EXPECT_THROW(RandomUniformlyI(l,r),exception);
}
TEST(RandomUniformlyR,BasicTest){
	for(double l=-10;l<=10;l+=0.5)
		for(double r=l;r<=10;r+=0.5)
			for(int c=0;c<100;c++){
				double V=RandomUniformlyR(l,r);
				EXPECT_TRUE((V>=l)&&(V<=r));
			}
}
TEST(RandomUniformlyR,Throwing){
	for(double l=-10;l<=10;l+=0.5)
		for(double r=-10;r<l;r+=0.5)
			EXPECT_THROW(RandomUniformlyR(l,r),exception);
}
TEST(RandomGauss,BasicTest){
	for(double X=-10;X<=10;X+=0.5){
		for(int c=0;c<100;c++){
			double V=RandomGauss(0.0,X);
			EXPECT_EQ(X,V);
		}
		for(double sigma=0.1;sigma<3;sigma+=0.1){
			int cg=0,cl=0;
			for(int i=0;i<100;i++){
				if(pow(RandomGauss(sigma,X)-X,2)<pow(2*sigma,2))cl++;
				else cg++;
			}
			EXPECT_TRUE(cg<cl);
		}
	}
}
TEST(RandomGauss,Throwing){
	for(double X=-1;X<=1;X+=0.5){
		auto f=[&X](){return RandomGauss(-1.0,X);};
		EXPECT_THROW(f(),exception);
	}
	EXPECT_THROW(RandomGauss(1.0,5.0,0),exception);
	EXPECT_NO_THROW(RandomGauss(1.0,5.0,1));
	EXPECT_NO_THROW(RandomGauss(1.0,5.0,2));
}
TEST(RandomValueGenerator,BaseTest){
	RandomValueGenerator<double> R([](double){return 1;},0,1,10);
	for(int i=0;i<100;i++){
		double r=R();
		EXPECT_TRUE((r>=0)&&(r<=1));
	}
	auto R2=R;
	for(int i=0;i<100;i++){
		double r=R2();
		EXPECT_TRUE((r>=0)&&(r<=1));
	}
}
TEST(RandomValueGenerator,Throwing){
	auto f=[](double){return 1;};
	auto z=[](double){return 0;};
	EXPECT_THROW(RandomValueGenerator<double>(f,0,1,0),exception);
	EXPECT_THROW(RandomValueGenerator<double>(f,0,1,-1),exception);
	EXPECT_THROW(RandomValueGenerator<double>(f,0,1,1),exception);
	EXPECT_THROW(RandomValueGenerator<double>(f,0,-1,2),exception);
	EXPECT_THROW(RandomValueGenerator<double>(f,0,0,2),exception);
	EXPECT_NO_THROW(RandomValueGenerator<double>(f,0,1,2));
	EXPECT_THROW(RandomValueGenerator<double>(z,0,1,2),exception);
	EXPECT_THROW(RandomValueGenerator<double>(z,0,0,0),exception);
}
TEST(RandomGauss,Compare){
	int N=1000000,bins=100;
	Distribution<double> D(0.0,10.0,bins);
	for(int i=0;i<N;i++)
		D.AddValue(RandomGauss(1.0,5.0));
	auto F=[](double x){return Gaussian(x,5.0,1.0);};
	double C=Sympson(F,0.0,10.0,0.0001);
	double S=0;
	int k=0;
	for(int i=0,n=D.size();i<n;i++){
		double d1=D.getY(i);S+=d1;
		double d2=F(D.getX(i))*double(N)*D.BinWidth()/C;
		printf("pract: %f\t theor: %f\n",d1,d2);
		double estim=d1;if(estim<1)estim=1;estim*=2;
		if(pow(d1-d2,2)>estim)k++;
	}
	EXPECT_EQ(N,S);
	EXPECT_TRUE(k<=(bins/10));
}
#define _EQ2(a,b) EXPECT_TRUE(pow(a-b,2)<0.02)
TEST(RandomValueGenerator,Compare){
	auto F=[](double x){return Gaussian(x,5.0,1.0);};
	double C=Sympson(F,0.0,10.0,0.0001);
	int N=1000000,bins=100;
	RandomValueGenerator<double> R(F,0.0,10.0,1000);
	Distribution<double> D(0.0,10.0,bins);
	Sigma<double> Sig;
	for(int i=0;i<N;i++){
		double d=R();
		D.AddValue(d);
		Sig.AddValue(d);
	}
	_EQ2(5.0,Sig.getAverage());
	_EQ2(1.0,Sig.getSigma());
	double S=0;
	int k=0;
	for(int i=0,n=D.size();i<n;i++){
		double d1=D.getY(i);S+=d1;
		double d2=F(D.getX(i))*double(N)*D.BinWidth()/C;
		printf("pract: %f\t theor: %f\n",d1,d2);
		double estim=d1;if(estim<1)estim=1;estim*=2;
		if(pow(d1-d2,2)>estim)k++;
	}
	EXPECT_EQ(N,S);
	EXPECT_TRUE(k<=(bins/10));
}
