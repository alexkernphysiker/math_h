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