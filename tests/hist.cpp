// this file is distributed under 
// MIT license
#include <gtest/gtest.h>
#include <vector>
#include <math.h>
#include <math_h/hist.h>
using namespace std;
using namespace MathTemplates;
TEST(point,basetest){
	value<double> x(25),y(78);
	point<double> p(x,y);
	EXPECT_EQ(x.val(),p.X().val());
	EXPECT_EQ(x.delta(),p.X().delta());
	EXPECT_EQ(y.val(),p.Y().val());
	EXPECT_EQ(y.delta(),p.Y().delta());
}