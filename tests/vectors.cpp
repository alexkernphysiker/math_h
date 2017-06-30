// this file is distributed under 
// MIT license
#include <gtest/gtest.h>
#include <math_h/vectors.h>
#include <math_h/sigma.h>
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

	EXPECT_EQ(1,Vector3<double>::basis_x().VecP(Vector3<double>::basis_y())*Vector3<double>::basis_z());
	EXPECT_EQ(1,Vector3<double>::basis_y().VecP(Vector3<double>::basis_z())*Vector3<double>::basis_x());
	EXPECT_EQ(1,Vector3<double>::basis_z().VecP(Vector3<double>::basis_x())*Vector3<double>::basis_y());

}
TEST(Vector3,Isotropic){
    mt19937 RG;
    StandardDeviation<double> X,Y,Z;
    for(size_t i=0;i<100;i++){
	const auto V=Vector3<double>::RandomIsotropicDirection(RG);
	X<<V.x();Y<<V.y();Z<<V.z();
    }
    EXPECT_TRUE(X().Contains(0));
    EXPECT_TRUE(Y().Contains(0));
    EXPECT_TRUE(Z().Contains(0));
}
