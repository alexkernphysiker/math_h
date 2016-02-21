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
TEST(Distribution1D,basetest){
	Distribution1D<double> D{value<double>(-1.0,0.5),value<double>(-0.0,0.5),value<double>(1.0,0.5)};
	ASSERT_EQ(3,D.size());
	ASSERT_EQ(-1,D[0].X().val());
	ASSERT_EQ(0,D[1].X().val());
	ASSERT_EQ(1,D[2].X().val());
	ASSERT_EQ(0.5,D[0].X().delta());
	ASSERT_EQ(0.5,D[1].X().delta());
	ASSERT_EQ(0.5,D[2].X().delta());
	EXPECT_EQ(0,D[0].Y().val());
	EXPECT_EQ(1,D[0].Y().delta());
	EXPECT_EQ(0,D[1].Y().val());
	EXPECT_EQ(1,D[1].Y().delta());
	EXPECT_EQ(0,D[2].Y().val());
	EXPECT_EQ(1,D[2].Y().delta());
	D<<0.0;
	EXPECT_EQ(0,D[0].Y().val());
	EXPECT_EQ(1,D[0].Y().delta());
	EXPECT_EQ(1,D[1].Y().val());
	EXPECT_EQ(1,D[1].Y().delta());
	EXPECT_EQ(0,D[2].Y().val());
	EXPECT_EQ(1,D[2].Y().delta());
	D<<1.0;
	EXPECT_EQ(0,D[0].Y().val());
	EXPECT_EQ(1,D[0].Y().delta());
	EXPECT_EQ(1,D[1].Y().val());
	EXPECT_EQ(1,D[1].Y().delta());
	EXPECT_EQ(1,D[2].Y().val());
	EXPECT_EQ(1,D[2].Y().delta());
	D<<-0.7;
	EXPECT_EQ(1,D[0].Y().val());
	EXPECT_EQ(1,D[0].Y().delta());
	EXPECT_EQ(1,D[1].Y().val());
	EXPECT_EQ(1,D[1].Y().delta());
	EXPECT_EQ(1,D[2].Y().val());
	EXPECT_EQ(1,D[2].Y().delta());
	D<<-0.2;
	EXPECT_EQ(1,D[0].Y().val());
	EXPECT_EQ(1,D[0].Y().delta());
	EXPECT_EQ(2,D[1].Y().val());
	EXPECT_EQ(sqrt(2.0),D[1].Y().delta());
	EXPECT_EQ(1,D[2].Y().val());
	EXPECT_EQ(1,D[2].Y().delta());
}
TEST(Distribution2D,BaseTest){
	Distribution2D<double> D(
		{value<double>(-1.0,0.5),value<double>(-0.0,0.5),value<double>(1.0,0.5)},
		{value<double>(-0.0,0.5),value<double>(1.0,0.5)}
	);
	ASSERT_EQ(3,D.X().size());
	ASSERT_EQ(2,D.Y().size());
	ASSERT_EQ(3,D.size());
	ASSERT_EQ(2,D[0].size());
	ASSERT_EQ(2,D[1].size());
	ASSERT_EQ(2,D[2].size());
	EXPECT_EQ(0,D[0][0].val());
	EXPECT_EQ(0,D[0][1].val());
	EXPECT_EQ(0,D[1][0].val());
	EXPECT_EQ(0,D[1][1].val());
	EXPECT_EQ(0,D[2][0].val());
	EXPECT_EQ(0,D[2][1].val());
	D<<make_pair(0.0,0.0)<<make_pair(0.1,0.2)<<make_pair(0.9,1.1)<<make_pair(-0.9,1.1)<<make_pair(0.6,0.6)<<make_pair(0.6,0.4);
	EXPECT_EQ(0,D[0][0].val());
	EXPECT_EQ(1,D[0][1].val());
	EXPECT_EQ(2,D[1][0].val());
	EXPECT_EQ(0,D[1][1].val());
	EXPECT_EQ(1,D[2][0].val());
	EXPECT_EQ(2,D[2][1].val());
	vector<point3d<double>> dbg;
	D.FullCycle([&dbg](point3d<double>&&p){dbg.push_back(p);});
	ASSERT_EQ(6,dbg.size());
	EXPECT_EQ(-1,dbg[0].X().val());
	EXPECT_EQ( 0,dbg[0].Y().val());
	EXPECT_EQ(0,dbg[0].Z().val());
	EXPECT_EQ(-1,dbg[1].X().val());
	EXPECT_EQ( 1,dbg[1].Y().val());
	EXPECT_EQ(1,dbg[1].Z().val());
	EXPECT_EQ( 0,dbg[2].X().val());
	EXPECT_EQ( 0,dbg[2].Y().val());
	EXPECT_EQ(2,dbg[2].Z().val());
	EXPECT_EQ( 0,dbg[3].X().val());
	EXPECT_EQ( 1,dbg[3].Y().val());
	EXPECT_EQ(0,dbg[3].Z().val());
	EXPECT_EQ( 1,dbg[4].X().val());
	EXPECT_EQ( 0,dbg[4].Y().val());
	EXPECT_EQ(1,dbg[4].Z().val());
	EXPECT_EQ( 1,dbg[5].X().val());
	EXPECT_EQ( 1,dbg[5].Y().val());
	EXPECT_EQ(2,dbg[5].Z().val());
}