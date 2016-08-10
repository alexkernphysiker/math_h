// this file is distributed under 
// MIT license
#include <gtest/gtest.h>
#include <vector>
#include <math.h>
#include <math_h/matrices.h>
using namespace std;
using namespace MathTemplates;
typedef Exception<vector<vector<double>>,0> invalid_matr;
typedef Exception<vector<vector<double>>,1> size_mismatch;
TEST(Matrix,SimpleObjects){
	auto A=Unitary<double>(3);
	EXPECT_EQ(A.height(),A.width());
	EXPECT_EQ(3,A.height());
	EXPECT_TRUE(A==A);
	EXPECT_FALSE(A!=A);
	EXPECT_TRUE(A.HasSizeAs(A));
	EXPECT_TRUE(A==Unitary<double>(3));
	EXPECT_FALSE(Unitary<double>(3)!=A);
	EXPECT_TRUE(Unitary<double>(3)==A);
	EXPECT_FALSE(A!=Unitary<double>(3));
	EXPECT_TRUE(A.HasSizeAs(Unitary<double>(3)));
	for(size_t i=0;i<3;i++)for(size_t j=0;j<3;j++)
		if(i==j)EXPECT_EQ(A(i,j),1);
		else EXPECT_EQ(A(i,j),0);
	EXPECT_TRUE(A.HasSizeAs(Zeros<double>(3,3)));
	EXPECT_FALSE(A==Zeros<double>(3,3));
	EXPECT_TRUE(Zeros<double>(3,3).HasSizeAs(A));
	EXPECT_FALSE(Zeros<double>(3,3)==A);
	EXPECT_EQ(Zeros<double>(4,2).height(),4);
	EXPECT_EQ(Zeros<double>(4,2).width(),2);
	EXPECT_TRUE(Zeros<double>(4,2)==MatrixByFormula<double>(4,2,[](size_t,size_t){return 0.0;}));
	auto Z=Zeros<double>(3,4);
	EXPECT_TRUE(Z==Z);
	EXPECT_FALSE(Z!=Z);
	EXPECT_EQ(Z.height(),3);
	EXPECT_EQ(Z.width(),4);
	for(size_t i=0;i<3;i++)for(size_t j=0;j<4;j++)
		EXPECT_EQ(Z(i,j),0);
}
const function<double(size_t,size_t)> F=[](size_t i,size_t j){return i*30+j;};
TEST(Matrix,Formula){
	auto A=MatrixByFormula<double>(3,4,F);
	auto B=MatrixData<double>({
		{0,0,0,0},
		{0,0,0,0},
		{0,0,0,0}
	});
	B.Transform([](size_t i,size_t j,const double&){return F(i,j);});
	EXPECT_TRUE(A==B);
	for(size_t i=0;i<3;i++)for(size_t j=0;j<4;j++)
		EXPECT_EQ(A(i,j),F(i,j));
}
TEST(Matrix,Transponate){
	auto A=MatrixByFormula<double>(3,4,F);
	auto B=Transponate(A);
	EXPECT_EQ(A.width(),B.height());
	EXPECT_EQ(B.width(),A.height());
	for(size_t i=0;i<A.height();i++)for(size_t j=0;j<A.width();j++)
		EXPECT_EQ(A(i,j),B(j,i));
	EXPECT_TRUE(Transponate(A)==B);
}
TEST(Matrix,Multiply){
	EXPECT_EQ(Multiply(Zeros<double>(3,3),Zeros<double>(3,3)),Zeros<double>(3,3));
	EXPECT_EQ(Multiply(Zeros<double>(3,5),Zeros<double>(5,3)),Zeros<double>(3,3));
	EXPECT_EQ(Multiply(Zeros<double>(2,3),Zeros<double>(3,4)),Zeros<double>(2,4));
	EXPECT_EQ(MatrixByFormula<double>(4,4,F),Multiply(MatrixByFormula<double>(4,4,F),Unitary<double>(4)));
	EXPECT_EQ(MatrixByFormula<double>(4,4,F),Multiply(Unitary<double>(4),MatrixByFormula<double>(4,4,F)));
	EXPECT_EQ(Zeros<double>(4,4),Multiply(MatrixByFormula<double>(4,4,F),Zeros<double>(4,4)));
	EXPECT_EQ(Zeros<double>(4,4),Multiply(Zeros<double>(4,4),MatrixByFormula<double>(4,4,F)));
	for(size_t i=0;i<5;i++)for(size_t j=0;j<5;j++)
		EXPECT_EQ(Multiply(Transponate(RVec<double>(5,i)),RVec<double>(5,j)),(i==j)?1.0:0.0);
	auto A=MatrixData<double>({
		{1,2,3,4},
		{5,6,7,8},
		{9,0,1,2}
	});
	auto B=MatrixData<double>({
		{1,2},
		{3,4},
		{5,6},
		{7,8}
	});
	auto C=MatrixData<double>({
		{50,60},
		{114,140},
		{28,40}
	});
	EXPECT_EQ(Multiply(A,B),C);
}