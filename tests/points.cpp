// this file is distributed under 
// MIT license
#include <gtest/gtest.h>
#include <math_h/points.h>
using namespace std;
using namespace MathTemplates;
TEST(point,basetest){
	double x(25),y(78);
	point<double> p(x,y);
	EXPECT_EQ(x,p.X());
	EXPECT_EQ(y,p.Y());
}
TEST(point3d,basetest){
	double x(25),y(78),z(55);
	point3d<double> p(x,y,z);
	EXPECT_EQ(x,p.X());
	EXPECT_EQ(y,p.Y());
	EXPECT_EQ(z,p.Z());
}
#include <math_h/sigma.h>
TEST(point,basetest_value){
	value<double> x(25),y(78);
	point<value<double>> p(x,y);
	EXPECT_EQ(x.val(),p.X().val());
	EXPECT_EQ(x.uncertainty(),p.X().uncertainty());
	EXPECT_EQ(y.val(),p.Y().val());
	EXPECT_EQ(y.uncertainty(),p.Y().uncertainty());
}
TEST(point3d,basetest_value){
	value<double> x(25),y(78),z(55);
	point3d<value<double>> p(x,y,z);
	EXPECT_EQ(x.val(),p.X().val());
	EXPECT_EQ(x.uncertainty(),p.X().uncertainty());
	EXPECT_EQ(y.val(),p.Y().val());
	EXPECT_EQ(y.uncertainty(),p.Y().uncertainty());
	EXPECT_EQ(z.val(),p.Z().val());
	EXPECT_EQ(z.uncertainty(),p.Z().uncertainty());
}
