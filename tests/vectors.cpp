// this file is distributed under 
// MIT license
#include <gtest/gtest.h>
#include <math_h/vectors.h>
using namespace std;
using namespace MathTemplates;
TEST(Vector3,scalar_prod){
	EXPECT_EQ(0,Vector3<double>::basis_x()*Vector3<double>::basis_y());
	EXPECT_EQ(0,Vector3<double>::basis_y()*Vector3<double>::basis_z());
	EXPECT_EQ(0,Vector3<double>::basis_x()*Vector3<double>::basis_z());
	EXPECT_EQ(1,Vector3<double>::basis_x()*Vector3<double>::basis_x());
	EXPECT_EQ(1,Vector3<double>::basis_y()*Vector3<double>::basis_y());
	EXPECT_EQ(1,Vector3<double>::basis_z()*Vector3<double>::basis_z());

	EXPECT_EQ(0,Vector3<double>::basis_x()*(-Vector3<double>::basis_y()));
	EXPECT_EQ(0,(-Vector3<double>::basis_x())*Vector3<double>::basis_z());
	
	EXPECT_EQ(-1,(-Vector3<double>::basis_x())*Vector3<double>::basis_x());
	EXPECT_EQ(-1,Vector3<double>::basis_y()*(-Vector3<double>::basis_y()));
	EXPECT_EQ(1,(-Vector3<double>::basis_z())*(-Vector3<double>::basis_z()));
}
TEST(Vector3,vector_prod){
	EXPECT_EQ(Vector3<double>::basis_x().VecP(Vector3<double>::basis_y()),Vector3<double>::basis_z());
	EXPECT_EQ(Vector3<double>::basis_x().VecP(Vector3<double>::basis_x()),Vector3<double>::zero());
	EXPECT_EQ(Vector3<double>::basis_y().VecP(Vector3<double>::basis_x()),-Vector3<double>::basis_z());
	EXPECT_EQ(Vector3<double>::basis_x().VecP(-Vector3<double>::basis_x()),Vector3<double>::zero());
}
