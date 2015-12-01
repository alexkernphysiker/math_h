// this file is distributed under 
// MIT license
#include <gtest/gtest.h>
#include <singleparam.h>
#include <functions.h>
using namespace std;
TEST(SingleParam,BaseTest){
	SingleParam<double,0,double,double,double> G(&Gaussian,INFINITY,5,1);
	for(double x=0;x<=10;x+=0.1)
		EXPECT_EQ(Gaussian<double>(x,5,1),G(x));
	SingleParam<double,1,double,double,double> G2(&Gaussian,1,INFINITY,1);
	for(double x=0;x<=10;x+=0.1)
		EXPECT_EQ(Gaussian<double>(1,x,1),G2(x));
}
