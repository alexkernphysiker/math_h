#include <gtest/gtest.h>
#include <sympson.h>
#include <functions.h>
using namespace std;
#define _EQ(a,b) EXPECT_TRUE(pow(a-b,2)<0.000001)
TEST(Sympson,BaseTest){
	auto F=[](double x){return x;};
	_EQ(0.5,Sympson(F,0.0,1.0,0.1));
	_EQ(0.5,Sympson(F,0.0,1.0,-0.1));
	_EQ(-0.5,Sympson(F,1.0,0.0,0.1));
	_EQ(-0.5,Sympson(F,1.0,0.0,-0.1));
}
TEST(Sympson,BaseTest2){
	_EQ(0.0,Sympson([](double x){return 0.0;},0.0,1.0,0.00001));
	_EQ(1.0,Sympson([](double x){return 1.0;},0.0,1.0,0.00001));
	_EQ(0.5,Sympson([](double x){return x;},0.0,1.0,0.00001));
	_EQ(0.333333,Sympson([](double x){return x*x;},0.0,1.0,0.00001));
	_EQ(0.25,Sympson([](double x){return x*x*x;},0.0,1.0,0.00001));
	_EQ(1.0,Sympson([](double x){return Gaussian(x,5.0,1.0);},0.0,10.0,0.0001));
}
TEST(SympsonTable,Simplest){
	double *X=nullptr;
	double *Y=SympsonTable<double,double*>([](double){return 0.0;},X,0);
	EXPECT_EQ(nullptr,Y);
}
TEST(SympsonTable,BasicTest){
	double *X=new double[11];
	for(int i=0;i<=10;i++)X[i]=0.1*i;
	auto F=[](double x){return x;};
	double S=Sympson(F,0.0,1.0,0.1);
	double *ST=SympsonTable<double,double*>(F,X,11);
	_EQ(0,ST[0]);
	_EQ(S,ST[10]);
	for(int i=1;i<=10;i++)EXPECT_TRUE(ST[i-1]<=ST[i]);
	EXPECT_NO_THROW(delete ST);
	auto F2=[](double x){return x*x;};
	S=Sympson(F2,0.0,1.0,0.1);
	ST=SympsonTable<double,double*>(F2,X,11);
	_EQ(0,ST[0]);
	_EQ(S,ST[10]);
	for(int i=1;i<=10;i++)EXPECT_TRUE(ST[i-1]<=ST[i]);
	EXPECT_NO_THROW(delete ST);
	EXPECT_NO_THROW(delete X);
}