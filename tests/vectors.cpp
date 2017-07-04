// this file is distributed under 
// MIT license
#include <gtest/gtest.h>
#include <math_h/vectors.h>
#include <math_h/sigma.h>
using namespace std;
using namespace MathTemplates;
TEST(Vector3,scalar_prod){
	EXPECT_EQ(0,Vector3<>::basis_x()*Vector3<>::basis_y());
	EXPECT_EQ(0,Vector3<>::basis_y()*Vector3<>::basis_z());
	EXPECT_EQ(0,Vector3<>::basis_x()*Vector3<>::basis_z());
	EXPECT_EQ(1,Vector3<>::basis_x()*Vector3<>::basis_x());
	EXPECT_EQ(1,Vector3<>::basis_y()*Vector3<>::basis_y());
	EXPECT_EQ(1,Vector3<>::basis_z()*Vector3<>::basis_z());

	EXPECT_EQ(0,Vector3<>::basis_x()*(-Vector3<>::basis_y()));
	EXPECT_EQ(0,(-Vector3<>::basis_x())*Vector3<>::basis_z());
	
	EXPECT_EQ(-1,(-Vector3<>::basis_x())*Vector3<>::basis_x());
	EXPECT_EQ(-1,Vector3<>::basis_y()*(-Vector3<>::basis_y()));
	EXPECT_EQ(1,(-Vector3<>::basis_z())*(-Vector3<>::basis_z()));
}
TEST(Vector3,vector_prod){
	EXPECT_EQ(Vector3<>::basis_x().VecP(Vector3<>::basis_y()),Vector3<>::basis_z());
	EXPECT_EQ(Vector3<>::basis_x().VecP(Vector3<>::basis_x()),Vector3<>::zero());
	EXPECT_EQ(Vector3<>::basis_y().VecP(Vector3<>::basis_x()),-Vector3<>::basis_z());
	EXPECT_EQ(Vector3<>::basis_x().VecP(-Vector3<>::basis_x()),Vector3<>::zero());

	EXPECT_EQ(1,Vector3<>::basis_x().VecP(Vector3<>::basis_y())*Vector3<>::basis_z());
	EXPECT_EQ(1,Vector3<>::basis_y().VecP(Vector3<>::basis_z())*Vector3<>::basis_x());
	EXPECT_EQ(1,Vector3<>::basis_z().VecP(Vector3<>::basis_x())*Vector3<>::basis_y());

}
TEST(Vector3,Isotropic){
    mt19937 RG;
    StandardDeviation<> X,Y,Z;
    for(size_t i=0;i<100;i++){
	const auto V=Vector3<>::RandomIsotropicDirection(RG);
	X<<V.x();Y<<V.y();Z<<V.z();
    }
    EXPECT_TRUE(X().Contains(0));
    EXPECT_TRUE(Y().Contains(0));
    EXPECT_TRUE(Z().Contains(0));
    //ToDo: check if all coordinates have uniform distribution
}

TEST(Vector4,Lorentz){
    mt19937 RG;
    const double epsilon=0.0000000000001;
    RandomUniform<> mr(0,1.0-epsilon),metrr(-5,5);
    for(size_t i=0;i<10000;i++){
	const auto V0=Vector4<>::bySpaceC_and_Length4(Vector3<>::RandomIsotropicDirection(RG),metrr(RG));
	const auto V1=V0.Lorentz(Vector3<>::zero());
	EXPECT_TRUE(V1==V0);
	const auto beta=Vector3<>::RandomIsotropicDirection(RG)*mr(RG);
	const auto V2=V0.Lorentz(beta).Lorentz(-beta);
	EXPECT_TRUE(pow(V2.time_component()-V0.time_component(),2)<epsilon);
	EXPECT_TRUE(pow(V2.space_component().x()-V0.space_component().x(),2)<epsilon);
	EXPECT_TRUE(pow(V2.space_component().y()-V0.space_component().y(),2)<epsilon);
	EXPECT_TRUE(pow(V2.space_component().z()-V0.space_component().z(),2)<epsilon);
	EXPECT_THROW(V0.Lorentz(Vector3<>::RandomIsotropicDirection(RG)*(1.0+mr(RG))),Exception<Vector4<>>);
	const auto L0=V0.length4(),L1=V0.Lorentz(beta).length4(),L2=V0.Lorentz(beta).Lorentz(Vector3<>::RandomIsotropicDirection(RG)*0.2).length4();
	EXPECT_TRUE(pow(L0-L1,2)<epsilon);
	EXPECT_TRUE(pow(L2-L1,2)<epsilon);
	EXPECT_TRUE(pow(L0-L2,2)<epsilon);
	const auto V00=Vector4<>(V0.length4()).Lorentz(-V0.Beta());
	EXPECT_TRUE(pow(V00.time_component()-V0.time_component(),2)<epsilon);
	EXPECT_TRUE(pow(V00.space_component().x()-V0.space_component().x(),2)<epsilon);
	EXPECT_TRUE(pow(V00.space_component().y()-V0.space_component().y(),2)<epsilon);
	EXPECT_TRUE(pow(V00.space_component().z()-V0.space_component().z(),2)<epsilon);
    }
}
