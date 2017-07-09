// this file is distributed under 
// MIT license
#include <gtest/gtest.h>
#include <math_h/vectors.h>
#include <math_h/hists.h>
#include <math_h/sigma.h>
using namespace std;
using namespace MathTemplates;
TEST(Vector2,scalar_prod){
	EXPECT_EQ(0,Vector2<>::basis_x()*Vector2<>::basis_y());
	EXPECT_EQ(0,Vector2<>::basis_y()*Vector2<>::basis_x());
	EXPECT_EQ(1,Vector2<>::basis_x()*Vector2<>::basis_x());
	EXPECT_EQ(1,Vector2<>::basis_y()*Vector2<>::basis_y());

	EXPECT_EQ(0,Vector2<>::basis_x()*(-Vector2<>::basis_y()));
	EXPECT_EQ(0,Vector2<>::basis_y()*(-Vector2<>::basis_x()));
	EXPECT_EQ(0,-Vector2<>::basis_x()*(Vector2<>::basis_y()));
	EXPECT_EQ(0,-Vector2<>::basis_y()*(Vector2<>::basis_x()));
	
	EXPECT_EQ(-1,(-Vector2<>::basis_x())*Vector2<>::basis_x());
	EXPECT_EQ(-1,Vector2<>::basis_y()*(-Vector2<>::basis_y()));
	EXPECT_EQ(-1,Vector2<>::basis_x()*(-Vector2<>::basis_x()));
	EXPECT_EQ(-1,(-Vector2<>::basis_y())*Vector2<>::basis_y());
}
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
    RANDOM RG;
    const auto c=BinsByStep(-1.0,0.1,1.0);
    Distribution1D<> X(c),Y(c),Z(c);
    for(size_t i=0;i<100000;i++){
	const auto V=Vector3<>::RandomIsotropicDirection(RG);
	X.Fill(V.x());Y.Fill(V.y());Z.Fill(V.z());
    }
    const auto x=X.TotalSum().val()/X.size();
    for(const auto&p:X)EXPECT_TRUE(p.Y().make_wider(2.5).Contains(x));
    for(const auto&p:Y)EXPECT_TRUE(p.Y().make_wider(2.5).Contains(x));
    for(const auto&p:Z)EXPECT_TRUE(p.Y().make_wider(2.5).Contains(x));
}
TEST(Vector3,IsotropicXY){
    RANDOM RG;
    StandardDeviation<> X,Y;
    for(size_t i=0;i<100000;i++){
	const auto V=Vector3<>::RandomIsotropicXYDirection(RG);
	X<<(V.x());Y<<(V.y());
	EXPECT_EQ(V.z(),0);
    }
    EXPECT_TRUE(X().Contains(0));
    EXPECT_TRUE(Y().Contains(0));
}
TEST(Vector3,IsotropicYZ){
    RANDOM RG;
    StandardDeviation<> Y,Z;
    for(size_t i=0;i<100000;i++){
	const auto V=Vector3<>::RandomIsotropicYZDirection(RG);
	Y<<(V.y());Z<<(V.z());
	EXPECT_EQ(V.x(),0);
    }
    EXPECT_TRUE(Z().Contains(0));
    EXPECT_TRUE(Y().Contains(0));
}
TEST(Vector3,IsotropicXZ){
    RANDOM RG;
    StandardDeviation<> X,Z;
    for(size_t i=0;i<100000;i++){
	const auto V=Vector3<>::RandomIsotropicXZDirection(RG);
	X<<(V.x());Z<<(V.z());
	EXPECT_EQ(V.y(),0);
    }
    EXPECT_TRUE(X().Contains(0));
    EXPECT_TRUE(Z().Contains(0));
}
TEST(Vector2,Rotation){
    RANDOM RG;
    const double epsilon=0.0000000000001;
    RandomUniform<> Phi(0,PI<>()*2.0);
    for(size_t i=0;i<10000;i++){
	const auto I=Vector2<>::RandomIsotropicDirection(RG);
	const auto ang=Phi(RG);
	const auto F=I.Rotate(ang);
	const auto p=I*F;
	EXPECT_TRUE(pow(p-cos(ang),2)<epsilon);
    }
}

TEST(Vector3,RotationX){
    RANDOM RG;
    const double epsilon=0.0000000000001;
    RandomUniform<> Phi(0,PI<>()*2.0);
    for(size_t i=0;i<10000;i++){
	const auto I=Vector3<>::RandomIsotropicYZDirection(RG);
	const auto axis=Vector3<>::basis_x();
	const auto ang=Phi(RG);
	const auto F=I.Rotate(axis,ang);
	EXPECT_EQ(axis*I,0);
	EXPECT_EQ(axis*F,0);
	const auto p=I*F;
	EXPECT_TRUE(pow(p-cos(ang),2)<epsilon);
    }
}
TEST(Vector3,RotationY){
    RANDOM RG;
    const double epsilon=0.0000000000001;
    RandomUniform<> Phi(0,PI<>()*2.0);
    for(size_t i=0;i<10000;i++){
	const auto I=Vector3<>::RandomIsotropicXZDirection(RG);
	const auto axis=Vector3<>::basis_y();
	const auto ang=Phi(RG);
	const auto F=I.Rotate(axis,ang);
	EXPECT_EQ(axis*I,0);
	EXPECT_EQ(axis*F,0);
	const auto p=I*F;
	EXPECT_TRUE(pow(p-cos(ang),2)<epsilon);
    }
}
TEST(Vector3,RotationZ){
    RANDOM RG;
    const double epsilon=0.0000000000001;
    RandomUniform<> Phi(0,PI<>()*2.0);
    for(size_t i=0;i<10000;i++){
	const auto I=Vector3<>::RandomIsotropicXYDirection(RG);
	const auto axis=Vector3<>::basis_z();
	const auto ang=Phi(RG);
	const auto F=I.Rotate(axis,ang);
	EXPECT_EQ(axis*I,0);
	EXPECT_EQ(axis*F,0);
	const auto p=I*F;
	EXPECT_TRUE(pow(p-cos(ang),2)<epsilon);
    }
}

TEST(LorentzVector,LorentzTransform){
    RANDOM RG;
    const double epsilon=0.0000000000001;
    RandomUniform<> mr(0,1.0-epsilon),metrr(-5,5);
    for(size_t i=0;i<10000;i++){
	const auto V0=LorentzVector<>::bySpaceC_and_Length4(Vector3<>::RandomIsotropicDirection(RG),metrr(RG));
	const auto V1=V0.Lorentz(Vector3<>::zero());
	EXPECT_TRUE(V1==V0);
	const auto beta=Vector3<>::RandomIsotropicDirection(RG)*mr(RG);
	const auto V2=V0.Lorentz(beta).Lorentz(-beta);
	EXPECT_TRUE(pow(V2.time_component()-V0.time_component(),2)<epsilon);
	EXPECT_TRUE(pow(V2.space_component().x()-V0.space_component().x(),2)<epsilon);
	EXPECT_TRUE(pow(V2.space_component().y()-V0.space_component().y(),2)<epsilon);
	EXPECT_TRUE(pow(V2.space_component().z()-V0.space_component().z(),2)<epsilon);
	EXPECT_THROW(V0.Lorentz(
	    Vector3<>::RandomIsotropicDirection(RG)*(1.0+mr(RG))),
	    Exception<LorentzVector<>>
	);
	const auto L0=V0.length4(),L1=V0.Lorentz(beta).length4(),
	L2=V0.Lorentz(beta).Lorentz(Vector3<>::RandomIsotropicDirection(RG)*0.2).length4();
	EXPECT_TRUE(pow(L0-L1,2)<epsilon);
	EXPECT_TRUE(pow(L2-L1,2)<epsilon);
	EXPECT_TRUE(pow(L0-L2,2)<epsilon);
	const auto V00=LorentzVector<>(V0.length4()).Lorentz(-V0.Beta());
	EXPECT_TRUE(pow(V00.time_component()-V0.time_component(),2)<epsilon);
	EXPECT_TRUE(pow(V00.space_component().x()-V0.space_component().x(),2)<epsilon);
	EXPECT_TRUE(pow(V00.space_component().y()-V0.space_component().y(),2)<epsilon);
	EXPECT_TRUE(pow(V00.space_component().z()-V0.space_component().z(),2)<epsilon);
    }
}
TEST(LorentzVector,LorentzTransform2d){
    RANDOM RG;
    const double epsilon=0.0000000000001;
    RandomUniform<> mr(0,1.0-epsilon),metrr(-5,5);
    for(size_t i=0;i<10000;i++){
	const auto V0=LorentzVector<Vector2<>>::bySpaceC_and_Length4(Vector2<>::RandomIsotropicDirection(RG),metrr(RG));
	const auto V1=V0.Lorentz(Vector2<>::zero());
	EXPECT_TRUE(V1==V0);
	const auto beta=Vector2<>::RandomIsotropicDirection(RG)*mr(RG);
	const auto V2=V0.Lorentz(beta).Lorentz(-beta);
	EXPECT_TRUE(pow(V2.time_component()-V0.time_component(),2)<epsilon);
	EXPECT_TRUE(pow(V2.space_component().x()-V0.space_component().x(),2)<epsilon);
	EXPECT_TRUE(pow(V2.space_component().y()-V0.space_component().y(),2)<epsilon);
	EXPECT_THROW(V0.Lorentz(
	    Vector2<>::RandomIsotropicDirection(RG)*(1.0+mr(RG))),
	    Exception<LorentzVector<Vector2<>>>
	);
	const auto L0=V0.length4(),L1=V0.Lorentz(beta).length4(),
	L2=V0.Lorentz(beta).Lorentz(Vector2<>::RandomIsotropicDirection(RG)*0.2).length4();
	EXPECT_TRUE(pow(L0-L1,2)<epsilon);
	EXPECT_TRUE(pow(L2-L1,2)<epsilon);
	EXPECT_TRUE(pow(L0-L2,2)<epsilon);
	const auto V00=LorentzVector<Vector2<>>(V0.length4()).Lorentz(-V0.Beta());
	EXPECT_TRUE(pow(V00.time_component()-V0.time_component(),2)<epsilon);
	EXPECT_TRUE(pow(V00.space_component().x()-V0.space_component().x(),2)<epsilon);
	EXPECT_TRUE(pow(V00.space_component().y()-V0.space_component().y(),2)<epsilon);
    }
}
