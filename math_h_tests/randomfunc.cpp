#include <gtest/gtest.h>
#include <vector>
#include <math.h>
#include <randomfunc.h>
#include <interpolate.h>
#include <sigma.h>
#include <functions.h>
using namespace std;

void Test_Random_Func(function<double()> R,function<double(double)> F,double from,double to,int bins){
	Distribution<double> D(from,to,bins);
	double norm=Sympson(F,from,to,0.0001);
	int N=5000;
	for(int i=0;i<N;i++)D.AddValue(R());
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
	EXPECT_TRUE(S<=2.0);
}
TEST(RandomUniformlyI,BasicTest){
	for(int l=-10;l<=10;l++)
		for(int r=l;r<=10;r++)
			for(int c=0;c<20;c++){
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
			for(int c=0;c<20;c++){
				double V=RandomUniformlyR(l,r);
				EXPECT_TRUE((V>=l)&&(V<=r));
			}
}
TEST(RandomUniformlyR,Throwing){
	for(double l=-10;l<=10;l+=0.5)
		for(double r=-10;r<l;r+=0.5)
			EXPECT_THROW(RandomUniformlyR(l,r),exception);
}
TEST(RandomUniformlyR,Distr){
	Test_Random_Func([](){
		return RandomUniformlyR<double>(0.0,10.0);
	},[](double){return 1.0;},0.0,10.0,20);
}
TEST(RandomGauss,Zeros){
	for(double X=-10;X<=10;X+=0.5)
		for(int c=0;c<100;c++){
			double V=RandomGauss(0.0,X);
			EXPECT_EQ(X,V);
		}
}
TEST(RandomGauss,Throwing){
	for(double X=-1;X<=1;X+=0.5){
		auto f=[&X](){return RandomGauss(-1.0,X);};
		EXPECT_THROW(f(),exception);
	}
}
void TestGauss(double mean,double sigma, double from,double to,int bins){
	Test_Random_Func([sigma,mean](){
		return RandomGauss(sigma,mean);
	},[sigma,mean](double x){
		return Gaussian(x,mean,sigma);
	},from,to,bins);
}
TEST(RandomGauss,ShapeTest){TestGauss(5,0.5, 0,10,20);}
TEST(RandomGauss,ShapeTest2){TestGauss(5,1, 0,10,20);}
TEST(RandomGauss,ShapeTest3){TestGauss(5,1.5, 0,10,20);}
TEST(RandomGauss,ShapeTest4){TestGauss(5,2, 0,10,20);}
TEST(RandomGauss,ShapeTest5){TestGauss(4,1, -3,10,26);}
TEST(RandomGauss,ShapeTest6){TestGauss(3,1, -3,10,26);}
TEST(RandomGauss,ShapeTest7){TestGauss(2,1, -3,10,26);}

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
	auto n=[](double x){return sin(10*x);};
	EXPECT_THROW(RandomValueGenerator<double>(n,0,1,100),exception);
}
void TestRandomDistribution(function<double(double)> F,double from,double to,int bins,int accu=10){
	RandomValueGenerator<double> R(F,from,to,bins*accu);
	Test_Random_Func([&R](){return R();},F,from,to,bins);
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
