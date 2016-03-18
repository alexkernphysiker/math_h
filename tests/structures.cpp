// this file is distributed under 
// MIT license
#include <gtest/gtest.h>
#include <vector>
#include <math.h>
#include <math_h/structures.h>
using namespace std;
using namespace MathTemplates;
double TestArray[]={-2,-1,0,1,2,3,4,5,6,7};
TEST(WhereToInsert,BorderConditions){
	for(int beg=0;beg<10;beg++)
		for(double V=TestArray[beg]-2;V<=TestArray[beg]+2;V+=0.5){
			int index=details::WhereToInsert(beg,beg-1,TestArray,V);
			EXPECT_EQ(beg,index);
			index=details::WhereToInsert(beg,beg,TestArray,V);
			if(V<TestArray[beg])EXPECT_EQ(beg,index);
			if(V>TestArray[beg])EXPECT_EQ(beg+1,index);
			if(V==TestArray[beg])EXPECT_TRUE((index==beg)||(index==(beg+1)));
		}
}
TEST(WhereToInsert,NormalConditions){
	for(double x=TestArray[0]-0.5;x<=TestArray[9]+0.5;x+=0.5)
		for(int beg=0;beg<10;beg++)for(int end=beg+1;end<10;end++){
			int index=details::WhereToInsert(beg,end,TestArray,x);
			if(x<TestArray[beg])EXPECT_EQ(beg,index);
			else if(x>TestArray[end])EXPECT_EQ(end+1,index);
			else EXPECT_TRUE((x>=TestArray[index-1])&&(x<=TestArray[index]));
		}
}
TEST(WhereToInsert,Search){
	for(double x=TestArray[0]-0.5;x<=TestArray[9]+0.5;x+=0.5)
		for(int beg=0;beg<10;beg++)for(int end=beg+1;end<10;end++){
			int index=details::Search(beg,end,TestArray,x);
			if(index>=beg)EXPECT_EQ(x,TestArray[index]);
		}
}
TEST(InsertSorted,BasicTest){
	vector<int> X;
	for(int i=0;i<50;i++){
		InsertSorted(rand()%10,X,std_size(X),std_insert(X,int));
		for(int j=0;j<i;j++)
			EXPECT_TRUE(X[j]<=X[j+1]);
	}
}
TEST(point,basetest){
	value<double> x(25),y(78);
	point<value<double>> p(x,y);
	EXPECT_EQ(x.val(),p.X().val());
	EXPECT_EQ(x.delta(),p.X().delta());
	EXPECT_EQ(y.val(),p.Y().val());
	EXPECT_EQ(y.delta(),p.Y().delta());
}
TEST(SortedPoints,size1){
	SortedPoints<double> chain;
	EXPECT_EQ(0,chain.size());
	chain<<point<double>(0,1);
	EXPECT_EQ(1,chain.size());
	chain<<point<double>(2,1);
	EXPECT_EQ(2,chain.size());
	chain<<point<double>(1,1);
	EXPECT_EQ(3,chain.size());
	chain<<point<double>(3,1);
	EXPECT_EQ(4,chain.size());
	chain<<point<double>(5,1);
	EXPECT_EQ(5,chain.size());
	chain<<point<double>(4,1);
	EXPECT_EQ(6,chain.size());
	EXPECT_EQ(0,chain.left().X());
	EXPECT_EQ(5,chain.right().X());
	for(size_t i=0;i<chain.size();i++){
		EXPECT_EQ(1,chain[i].Y());
		EXPECT_EQ(i,chain[i].X());
	}
}
TEST(SortedPoints,size2){
	SortedPoints<value<double>> chain;
	EXPECT_EQ(0,chain.size());
	chain<<point<value<double>>(0,1);
	EXPECT_EQ(1,chain.size());
	chain<<point<value<double>>(2,1);
	EXPECT_EQ(2,chain.size());
	chain<<point<value<double>>(1,1);
	EXPECT_EQ(3,chain.size());
	chain<<point<value<double>>(3,1);
	EXPECT_EQ(4,chain.size());
	chain<<point<value<double>>(5,1);
	EXPECT_EQ(5,chain.size());
	chain<<point<value<double>>(4,1);
	EXPECT_EQ(6,chain.size());
	EXPECT_EQ(0,chain.left().X().val());
	EXPECT_EQ(5,chain.right().X().val());
	for(size_t i=0;i<chain.size();i++){
		EXPECT_EQ(1,chain[i].Y().val());
		EXPECT_EQ(i,chain[i].X().val());
	}
}
TEST(SortedPoints,range){
	SortedPoints<double> chain({point<double>(0,1),point<double>(2,1.1),point<double>(3,1.2),point<double>(4,1.3),point<double>(5,1.4)});
	auto xcut=chain.XRange(0.5,4.5);
	EXPECT_EQ(3,xcut.size());
	EXPECT_EQ(2,xcut.left().X());
	EXPECT_EQ(4,xcut.right().X());
	auto ycut=chain.YRange(1.15,1.35).Transponate();
	EXPECT_EQ(2,ycut.size());
	EXPECT_EQ(1.2,ycut.left().X());
	EXPECT_EQ(1.3,ycut.right().X());
}



TEST(hist,scale_norm){
	hist<double> H(BinsByCount(10,0.0,1.0));
	for(size_t i=0;i<H.size();i++)
		H.Bin(i).varY()=value<double>(10.0+4.0*sin(H[i].X().val()));
	double s1=0;
	for(const auto&p:H)s1+=p.Y().val();
	double s2=0;
	auto H2=H.Scale(2);
	for(const auto&p:H2)s2+=p.Y().val();
	EXPECT_TRUE(pow(s1-s2,2)<0.000001);
	for(size_t x=0;x<H2.size();x++){
		EXPECT_EQ(sqrt(H2[x].Y().val()),H2[x].Y().delta());
		EXPECT_EQ(H2[x].Y().val(),H[x*2].Y().val()+H[x*2+1].Y().val());
	}
}
TEST(hist2d,scale_norm){
	hist2d<double> H(BinsByCount(10,0.0,1.0),BinsByCount(10,0.0,1.0));
	H.FullCycleVar([](const value<double>&x,const value<double>&y,value<double>&z){
		z=value<double>(10.0+4.0*sin(x.val()+y.val()));
	});
	double s1=0;
	H.FullCycle([&s1](const value<double>&,const value<double>&,const value<double>&z){s1+=z.val();});
	double s2=0;
	auto H2=H.Scale(2,2);
	H2.FullCycle([&s2](const value<double>&,const value<double>&,const value<double>&z){s2+=z.val();});
	EXPECT_TRUE(pow(s1-s2,2)<0.000001);
	for(size_t x=0;x<H2.X().size();x++)
		for(size_t y=0;y<H2.Y().size();y++){
			EXPECT_EQ(sqrt(H2[x][y].val()),H2[x][y].delta());
			EXPECT_EQ(H2[x][y].val(),H[x*2][y*2].val()+H[x*2+1][y*2].val()+H[x*2][y*2+1].val()+H[x*2+1][y*2+1].val());
		}
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
	D.Fill(0.0);
	EXPECT_EQ(0,D[0].Y().val());
	EXPECT_EQ(1,D[0].Y().delta());
	EXPECT_EQ(1,D[1].Y().val());
	EXPECT_EQ(1,D[1].Y().delta());
	EXPECT_EQ(0,D[2].Y().val());
	EXPECT_EQ(1,D[2].Y().delta());
	D.Fill(1.0);
	EXPECT_EQ(0,D[0].Y().val());
	EXPECT_EQ(1,D[0].Y().delta());
	EXPECT_EQ(1,D[1].Y().val());
	EXPECT_EQ(1,D[1].Y().delta());
	EXPECT_EQ(1,D[2].Y().val());
	EXPECT_EQ(1,D[2].Y().delta());
	D.Fill(-0.7);
	EXPECT_EQ(1,D[0].Y().val());
	EXPECT_EQ(1,D[0].Y().delta());
	EXPECT_EQ(1,D[1].Y().val());
	EXPECT_EQ(1,D[1].Y().delta());
	EXPECT_EQ(1,D[2].Y().val());
	EXPECT_EQ(1,D[2].Y().delta());
	D.Fill(-0.2);
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
	D.Fill(make_pair(0.0,0.0)).Fill(make_pair(0.1,0.2)).Fill(make_pair(0.9,1.1)).Fill(make_pair(-0.9,1.1)).Fill(make_pair(0.6,0.6)).Fill(make_pair(0.6,0.4));
	EXPECT_EQ(0,D[0][0].val());
	EXPECT_EQ(1,D[0][1].val());
	EXPECT_EQ(2,D[1][0].val());
	EXPECT_EQ(0,D[1][1].val());
	EXPECT_EQ(1,D[2][0].val());
	EXPECT_EQ(2,D[2][1].val());
	vector<point3d<value<double>>> dbg;
	D.FullCycle([&dbg](const point3d<value<double>>&p){dbg.push_back(p);});
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