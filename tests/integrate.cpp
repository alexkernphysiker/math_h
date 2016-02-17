// this file is distributed under 
// MIT license
#include <iostream>
#include <gtest/gtest.h>
#include <math_h/integrate.h>
#include <math_h/functions.h>
using namespace std;
using namespace MathTemplates;
#define _EQ(a,b) EXPECT_TRUE(pow(a-b,2)<0.00001)
TEST(Sympson,BaseTest){
	auto F=[](double x){return x;};
	_EQ(0.5,Sympson(F,0.0,1.0,0.1));
	_EQ(0.5,Sympson(F,0.0,1.0,-0.1));
	_EQ(-0.5,Sympson(F,1.0,0.0,0.1));
	_EQ(-0.5,Sympson(F,1.0,0.0,-0.1));
}
TEST(Sympson,BaseTest2){
	_EQ(0.0,Sympson([](double ){return 0.0;},0.0,1.0,0.00001));
	_EQ(1.0,Sympson([](double ){return 1.0;},0.0,1.0,0.00001));
	_EQ(0.5,Sympson([](double x){return x;},0.0,1.0,0.00001));
	_EQ(0.333333,Sympson([](double x){return x*x;},0.0,1.0,0.00001));
	_EQ(0.25,Sympson([](double x){return x*x*x;},0.0,1.0,0.00001));
	_EQ(1.0,Sympson([](double x){return Gaussian(x,5.0,1.0);},0.0,10.0,0.0001));
}
TEST(Int_Trapez_Table,BasicTest){
	LinearInterpolation<double> func;
	for(double x=0;x<=1;x+=0.0001)
		func<<make_pair(x,x*x);
	auto res=Int_Trapez_Table(func);
	EXPECT_EQ(func.left().first,res.left().first);
	EXPECT_EQ(func.right().first,res.right().first);
	EXPECT_EQ(0,res.left().second);
	_EQ(0.33333333,res.right().second);
	cout<<res.right().second<<endl;
}
TEST(Int_Trapez_Table_PositiveStrict,BasicTest){
	LinearInterpolation<double> func;
	for(double x=0;x<=1;x+=0.0001)
		func<<make_pair(x,x*x);
	auto res=Int_Trapez_Table_PositiveStrict(func);
	EXPECT_EQ(func.left().first,res.left().first);
	EXPECT_EQ(func.right().first,res.right().first);
	EXPECT_EQ(0,res.left().second);
	_EQ(0.333333,res.right().second);
	cout<<res.right().second<<endl;
}
TEST(Convolution,BasicTest){
	auto F1=[](double x){return x;};
	auto F2=[](double x){return x*x;};
	Convolution<double> C(F1,F2,0,1,0.01);
	double x=0.1;
	double test_value=Sympson([x,F1,F2](double k){return F1(k)*F2(x-k);},0.0,1.0,0.01);
	EXPECT_EQ(test_value,C(x));
}