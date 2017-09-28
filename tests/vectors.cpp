// this file is distributed under
// MIT license
#include <gtest/gtest.h>
#include <math_h/vectors.h>
#include <math_h/hists.h>
#include <math_h/sigma.h>
using namespace std;
using namespace MathTemplates;
const double epsilon=0.0000000001;
TEST(Vector, scalar_prod2)
{
    EXPECT_EQ(0, x<>()*y<>());
    EXPECT_EQ(0, y<>()*x<>());
    EXPECT_EQ(1, x<>()*x<>());
    EXPECT_EQ(1, y<>()*y<>());

    EXPECT_EQ(0, x<>() * (-y<>()));
    EXPECT_EQ(0, y<>() * (-x<>()));
    EXPECT_EQ(0, (-x<>())*y<>());
    EXPECT_EQ(0, (-y<>())*x<>());

    EXPECT_EQ(-1, x<>() * (-x<>()));
    EXPECT_EQ(-1, y<>() * (-y<>()));
    EXPECT_EQ(-1, (-x<>())*x<>());
    EXPECT_EQ(-1, (-y<>())*y<>());
}
TEST(Vector, scalar_prod3)
{
    EXPECT_EQ(0, X<>()*Y<>());
    EXPECT_EQ(0, Y<>()*Z<>());
    EXPECT_EQ(0, X<>()*Z<>());
    EXPECT_EQ(1, X<>()*X<>());
    EXPECT_EQ(1, Y<>()*Y<>());
    EXPECT_EQ(1, Z<>()*Z<>());

    EXPECT_EQ(0, X<>() * (-Y<>()));
    EXPECT_EQ(0, (-X<>())*Z<>());

    EXPECT_EQ(-1, (-X<>())*X<>());
    EXPECT_EQ(-1, Y<>() * (-Y<>()));
    EXPECT_EQ(1, (-Z<>()) * (-Z<>()));
}
TEST(Vector, vector_prod3)
{
    EXPECT_EQ(X<>()^Y<>(), Z<>());
    EXPECT_EQ(X<>()^X<>(), Zero<>());
    EXPECT_EQ(Y<>()^X<>(), -Z<>());
    EXPECT_EQ(X<>() ^ (-X<>()), Zero<>());

    EXPECT_EQ(1, (X<>()^Y<>())*Z<>());
    EXPECT_EQ(1, (Y<>()^Z<>())*X<>());
    EXPECT_EQ(1, (Z<>()^X<>())*Y<>());

}
TEST(Vector,baseAllD){
    EXPECT_TRUE(desCartes(1.0)==Vector<1>({1.0}));
    EXPECT_TRUE(desCartes(1.0,2.0)==Vector<2>({1.0,2.0}));
    EXPECT_TRUE(desCartes(1.0,2.0,3.0)==Vector<3>({1.0,2.0,3.0}));
    EXPECT_TRUE(desCartes(1.0,2.0,3.0,4.0)==Vector<4>({1.0,2.0,3.0,4.0}));
    EXPECT_TRUE(desCartes(1.0)==Vector<1>(desCartes(1.0).chain()));
    EXPECT_TRUE(desCartes(1.0,2.0)==Vector<2>(desCartes(1.0,2.0).chain()));
    EXPECT_TRUE(desCartes(1.0,2.0,3.0)==Vector<3>(desCartes(1.0,2.0,3.0).chain()));
    EXPECT_TRUE(desCartes(1.0,2.0,3.0,4.0)==Vector<4>(desCartes(1.0,2.0,3.0,4.0).chain()));
    EXPECT_EQ(desCartes(1.0).chain()[0],1.0);
    EXPECT_EQ(desCartes(1.0,2.0).chain()[0],1.0);
    EXPECT_EQ(desCartes(1.0,2.0,3.0).chain()[0],1.0);
    EXPECT_EQ(desCartes(1.0,2.0,3.0,4.0).chain()[0],1.0);
    EXPECT_EQ(desCartes(1.0,2.0).chain()[1],2.0);
    EXPECT_EQ(desCartes(1.0,2.0,3.0).chain()[1],2.0);
    EXPECT_EQ(desCartes(1.0,2.0,3.0,4.0).chain()[1],2.0);
}
TEST(Vector,base1d){
    RANDOM RG;
    RandomUniform<> X(-50.,50);
    for (size_t i = 0; i < 10; i++) {
        const auto x=X(RG);
	const auto V=desCartes(x);
	EXPECT_EQ(x,V.x());
	EXPECT_EQ(x,V.component<1>());
	const auto v=V.chain();
	EXPECT_EQ(V.Dimensions,v.size());
	EXPECT_EQ(x,v[0]);
	const Vector<1> V2=v;
	EXPECT_EQ(x,V2.x());
	EXPECT_TRUE(V==V2);
    }
}
TEST(Vector,base2d){
    RANDOM RG;
    RandomUniform<> X(-50.,50);
    for (size_t i = 0; i < 10; i++) {
        const auto x=X(RG),y=X(RG);
	const auto V=desCartes(x,y);
	EXPECT_EQ(x,V.x());
	EXPECT_EQ(x,V.component<1>());
	EXPECT_EQ(y,V.y());
	EXPECT_EQ(y,V.component<2>());
	const auto v=V.chain();
	EXPECT_EQ(V.Dimensions,v.size());
	EXPECT_EQ(x,v[0]);
	EXPECT_EQ(y,v[1]);
	const Vector<2> V2=v;
	EXPECT_EQ(x,V2.x());
	EXPECT_EQ(y,V2.y());
	EXPECT_TRUE(V==V2);
    }
}
TEST(Vector,base3d){
    RANDOM RG;
    RandomUniform<> X(-50.,50);
    for (size_t i = 0; i < 10; i++) {
        const auto x=X(RG),y=X(RG),z=X(RG);
	const auto V=desCartes(x,y,z);
	EXPECT_EQ(x,V.x());
	EXPECT_EQ(x,V.component<1>());
	EXPECT_EQ(y,V.y());
	EXPECT_EQ(y,V.component<2>());
	EXPECT_EQ(z,V.z());
	EXPECT_EQ(z,V.component<3>());
	const auto v=V.chain();
	EXPECT_EQ(V.Dimensions,v.size());
	EXPECT_EQ(x,v[0]);
	EXPECT_EQ(y,v[1]);
	EXPECT_EQ(z,v[2]);
	const Vector<3> V2=v;
	EXPECT_EQ(x,V2.x());
	EXPECT_EQ(y,V2.y());
	EXPECT_EQ(z,V2.z());
	EXPECT_TRUE(V==V2);
    }
}
TEST(Vector,base4d){
    RANDOM RG;
    RandomUniform<> X(-50.,50);
    for (size_t i = 0; i < 10; i++) {
        const auto x=X(RG),y=X(RG),z=X(RG),zz=X(RG);
	const auto V=desCartes(x,y,z,zz);
	EXPECT_EQ(x,V.x());
	EXPECT_EQ(x,V.component<1>());
	EXPECT_EQ(y,V.y());
	EXPECT_EQ(y,V.component<2>());
	EXPECT_EQ(z,V.z());
	EXPECT_EQ(z,V.component<3>());
	EXPECT_EQ(zz,V.component<4>());
	const auto v=V.chain();
	EXPECT_EQ(V.Dimensions,v.size());
	EXPECT_EQ(x,v[0]);
	EXPECT_EQ(y,v[1]);
	EXPECT_EQ(z,v[2]);
	EXPECT_EQ(zz,v[3]);
	const Vector<4> V2=v;
	EXPECT_EQ(x,V2.x());
	EXPECT_EQ(y,V2.y());
	EXPECT_EQ(z,V2.z());
	EXPECT_EQ(zz,V2.component<4>());
	EXPECT_TRUE(V==V2);
    }
}

TEST(Direction,base2d){
    RANDOM RG;
    RandomUniform<> Phi(-PI<>(), PI<>());
    for (size_t i = 0; i < 10; i++) {
        const auto phi=Phi(RG);
	const auto D=direction(phi);
	EXPECT_EQ(D.phi(),phi);
    }
}
TEST(Direction,base3d){
    RANDOM RG;
    RandomUniform<> Phi(-PI<>(), PI<>());
    RandomUniform<> Theta(0, PI<>());
    for (size_t i = 0; i < 10; i++) {
        const auto phi=Phi(RG);
        const auto theta=Theta(RG);
	const auto D=direction(phi,theta);
	EXPECT_EQ(D.phi(),phi);
	EXPECT_EQ(D.th<1>(),theta);
    }
}
TEST(Direction,base4d){
    RANDOM RG;
    RandomUniform<> Phi(-PI<>(), PI<>());
    RandomUniform<> Theta(0, PI<>());
    for (size_t i = 0; i < 10; i++) {
        const auto phi=Phi(RG);
        const auto theta1=Theta(RG);
        const auto theta2=Theta(RG);
	const auto D=direction(phi,theta1,theta2);
	EXPECT_EQ(D.phi(),phi);
	EXPECT_EQ(D.th<1>(),theta1);
	EXPECT_EQ(D.th<2>(),theta2);
    }
}
TEST(Direction,base5d){
    RANDOM RG;
    RandomUniform<> Phi(-PI<>(), PI<>());
    RandomUniform<> Theta(0, PI<>());
    for (size_t i = 0; i < 10; i++) {
        const auto phi=Phi(RG);
        const auto theta1=Theta(RG);
        const auto theta2=Theta(RG);
        const auto theta3=Theta(RG);
	const auto D=direction(phi,theta1,theta2,theta3);
	EXPECT_EQ(D.phi(),phi);
	EXPECT_EQ(D.th<1>(),theta1);
	EXPECT_EQ(D.th<2>(),theta2);
	EXPECT_EQ(D.th<3>(),theta3);
    }
}
TEST(Vector, Isotropic3)
{
    RANDOM RG;
    const auto c = BinsByStep(-1.0, 0.1, 1.0);
    Distribution1D<> X(c), Y(c), Z(c);
    for (size_t i = 0; i < 1000000; i++) {
        const auto V = randomIsotropic<3>(RG)*1.0;
        X.Fill(V.x());
        Y.Fill(V.y());
        Z.Fill(V.z());
    }
    const auto x = X.TotalSum().val() / X.size();
    for (const auto &p : X)EXPECT_TRUE(p.Y().make_wider(2.5).Contains(x));
    for (const auto &p : Y)EXPECT_TRUE(p.Y().make_wider(2.5).Contains(x));
    for (const auto &p : Z)EXPECT_TRUE(p.Y().make_wider(2.5).Contains(x));
}
TEST(Vector, Rotation2)
{
    RANDOM RG;
    RandomUniform<> Phi(0, PI<>() * 2.0);
    for (size_t i = 0; i < 10; i++) {
        const auto I = randomIsotropic<2>(RG)*1.0;
        const auto ang = Phi(RG);
        const auto F = Rotate(I, ang);
    }
}

TEST(Vector, Rotation3)
{
    RANDOM RG;
    RandomUniform<> M(0, 5.0),Phi(-PI<>(),PI<>());
    for (size_t i = 0; i < 10; i++) {
        const auto I = randomIsotropic<3>(RG)*M(RG);
        const auto ang = Phi(RG);
	const auto dir =randomIsotropic<3>(RG);
        const auto F = Rotate(I, dir, ang);
        EXPECT_TRUE(pow((dir*1.0) * I - (dir*1.0) * F, 2) < 0.0000000000001);
    }
}

TEST(LorentzVector, LorentzTransform)
{
    RANDOM RG;
    RandomUniform<> mr(0, 0.99), metrr(-5, 5);
    for (size_t i = 0; i < 10; i++) {
        const auto V0 = lorentz_byPM(randomIsotropic<3>(RG)*1.0, metrr(RG));
        const auto V1 = V0.Transform(Zero<>());
        EXPECT_TRUE(V1 == V0);
        const auto beta = randomIsotropic<3>(RG) * mr(RG);
        const auto V2 = V0.Transform(beta).Transform(-beta);
        EXPECT_TRUE(abs(V2.T() - V0.T()) < epsilon);
	EXPECT_TRUE(V2.S().CloseTo(V0.S(),epsilon));
        EXPECT_THROW(V0.Transform(randomIsotropic<3>(RG) * (1.0 + mr(RG))),Exception<LorentzVector<>>);
        const auto L0 = V0.M(), L1 = V0.Transform(beta).M(),
                   L2 = V0.Transform(beta).Transform(randomIsotropic<3>(RG) * 0.2).M();
        EXPECT_TRUE(abs(L0 - L1) < epsilon);
        EXPECT_TRUE(abs(L2 - L1) < epsilon);
        EXPECT_TRUE(abs(L0 - L2) < epsilon);
        const auto V00 = lorentz_byPM(Zero<>(),V0.M()).Transform(-V0.Beta());
        EXPECT_TRUE(abs(V00.T() - V0.T()) < epsilon);
	EXPECT_TRUE(V00.S().CloseTo(V0.S(),epsilon));
    }
}
TEST(LorentzVector, LorentzTransform2)
{
    RANDOM RG;
    RandomUniform<> M(0, 5),P(0, 5);
    for (size_t i = 0; i < 10; i++) {
        const auto V0 = lorentz_byPM(randomIsotropic<3>(RG)*P(RG), M(RG));
	const auto V1=V0.Transform(V0.Beta());
        EXPECT_TRUE(V1.S().CloseTo(Zero(),epsilon));
        EXPECT_TRUE(abs(V1.M()-V0.M()) < epsilon);
    }
}
TEST(LorentzVector, LorentzTransform2d)
{
    RANDOM RG;
    RandomUniform<> mr(0, 1.0 - epsilon), metrr(-5, 5);
    for (size_t i = 0; i < 10; i++) {
        const auto V0 = lorentz_byPM(randomIsotropic<2>(RG)*1.0, metrr(RG));
        const auto V1 = V0.Transform(zero<>());
        EXPECT_TRUE(V1 == V0);
        const auto beta = randomIsotropic<2>(RG) * mr(RG);
        const auto V2 = V0.Transform(beta).Transform(-beta);
        EXPECT_TRUE(abs(V2.T() - V0.T()) < epsilon);
        EXPECT_TRUE(V2.S().CloseTo(V0.S(),epsilon));
	typedef LorentzVector<double,Vector<2>> LV2;
        EXPECT_THROW(V0.Transform(randomIsotropic<2>(RG) * (1.0 + mr(RG))),Exception<LV2>);
        const auto L0 = V0.M(), L1 = V0.Transform(beta).M(),
                   L2 = V0.Transform(beta).Transform(randomIsotropic<2>(RG) * 0.2).M();
        EXPECT_TRUE(abs(L0 - L1) < epsilon);
        EXPECT_TRUE(abs(L2 - L1) < epsilon);
        EXPECT_TRUE(abs(L0 - L2) < epsilon);
        const auto V00 = lorentz_byPM(zero<>(),V0.M()).Transform(-V0.Beta());
        EXPECT_TRUE(abs(V00.T() - V0.T()) < epsilon);
        EXPECT_TRUE(V00.S().CloseTo(V0.S(),epsilon));
    }
}

TEST(LorentzVector, decays)
{
    RANDOM RG;
    RandomUniform<> IM(2,3),M1(0,1),M2(0,1),THETA(0,PI()),PHI(0, 2.0 * PI());
    for (size_t i = 0; i < 10; i++) {
	const double im=IM(RG),m1=M1(RG),m2=M2(RG);
	const auto C2=binaryDecay(im,m1,m2,randomIsotropic<2>(RG));
	EXPECT_TRUE(abs((C2.first+C2.second).M()-im)<epsilon);
	EXPECT_TRUE(abs((C2.first).M()-m1)<epsilon);
	EXPECT_TRUE(abs((C2.second).M()-m2)<epsilon);
	EXPECT_TRUE((C2.first.S()+C2.second.S()).M()<epsilon);
	const auto C3=binaryDecay(im,m1,m2,randomIsotropic<3>(RG));
	EXPECT_TRUE(abs((C3.first+C3.second).M()-im)<epsilon);
	EXPECT_TRUE(abs(C3.first.M()-m1)<epsilon);
	EXPECT_TRUE(abs(C3.second.M()-m2)<epsilon);
	EXPECT_TRUE((C3.first.S()+C3.second.S()).M()<epsilon);
    }
}
TEST(LorentzVector, decays2)
{
    RANDOM RG;
    RandomUniform<> IM(2,3),M1(0,1),M2(0,1),THETA(0,PI()),PHI(0, 2.0 * PI());
    for (size_t i = 0; i < 10; i++) {
	const double im=IM(RG),m1=M1(RG),m2=M2(RG);
	const auto C2=binaryDecay(im,m1,m2,randomIsotropic<2>(RG));
	EXPECT_TRUE(abs((C2.first+C2.second).M()-im)<epsilon);
	EXPECT_TRUE(abs((C2.first).M()-m1)<epsilon);
	EXPECT_TRUE(abs((C2.second).M()-m2)<epsilon);
	EXPECT_TRUE((C2.first.S()+C2.second.S()).M()<epsilon);
	const auto C3=binaryDecay(im,m1,m2,randomIsotropic<3>(RG));
	EXPECT_TRUE(abs((C3.first+C3.second).M()-im)<epsilon);
	EXPECT_TRUE(abs(C3.first.M()-m1)<epsilon);
	EXPECT_TRUE(abs(C3.second.M()-m2)<epsilon);
	EXPECT_TRUE((C3.first.S()+C3.second.S()).M()<epsilon);
    }
}

TEST(Vector,angle2d)
{
    RANDOM RG;
    RandomUniform<> PHI(-PI(),PI()),M(0.0,10.0);
    for (size_t i = 0; i < 10; i++) {
	const auto phi=PHI(RG);
	const auto phi2=direction(direction(phi)*M(RG)).phi();
	EXPECT_TRUE(abs(phi-phi2)<epsilon);
    }
}
TEST(Vector,angle3d)
{
    RANDOM RG;
    RandomUniform<> THETA(0.0,PI()),PHI(-PI(),PI()),M(0.0,10.0);
    for (size_t i = 0; i < 10; i++) {
	const auto angles=direction(PHI(RG),THETA(RG));
	const auto angles2=direction(angles*M(RG));
	EXPECT_TRUE(abs(angles.phi()-angles2.phi())<epsilon);
	EXPECT_TRUE(abs(angles.th<1>()-angles2.th<1>())<epsilon);
    }
}
TEST(Vector,angle4d)
{
    RANDOM RG;
    RandomUniform<> THETA(0.0,PI()),PHI(-PI(),PI()),M(0.0,10.0);
    for (size_t i = 0; i < 10; i++) {
	const auto angles=direction(PHI(RG),THETA(RG),THETA(RG));
	const auto angles2=direction(angles*M(RG));
	EXPECT_TRUE(abs(angles.phi()-angles2.phi())<epsilon);
	EXPECT_TRUE(abs(angles.th<1>()-angles2.th<1>())<epsilon);
	EXPECT_TRUE(abs(angles.th<2>()-angles2.th<2>())<epsilon);
    }
}
TEST(Vector,angle5d)
{
    RANDOM RG;
    RandomUniform<> THETA(0.0,PI()),PHI(-PI(),PI()),M(0.0,10.0);
    for (size_t i = 0; i < 10; i++) {
	const auto angles=direction(PHI(RG),THETA(RG),THETA(RG),THETA(RG));
	const auto angles2=direction(angles*M(RG));
	EXPECT_TRUE(abs(angles.phi()-angles2.phi())<epsilon);
	EXPECT_TRUE(abs(angles.th<1>()-angles2.th<1>())<epsilon);
	EXPECT_TRUE(abs(angles.th<2>()-angles2.th<2>())<epsilon);
	EXPECT_TRUE(abs(angles.th<3>()-angles2.th<3>())<epsilon);
    }
}

TEST(Plane3D, SimplePlanes)
{
    RANDOM RG;
    RandomUniform<> mr(0, 10), th(0, 2.0 * PI());
    for (size_t i = 0; i < 10000; i++) {
        const auto v1 = randomIsotropic<2>(RG) * mr(RG), v2 = randomIsotropic<2>(RG) * mr(RG);
        const auto p1 = v1 * v2;
        {
            const auto plane = Plane3D<>(X<>(),Y<>());
            const auto p2 = plane(v1) * plane(v2);
            EXPECT_TRUE(abs(p1 - p2) < epsilon);
            EXPECT_EQ(v1.x(), plane(v1).x());
            EXPECT_EQ(v1.y(), plane(v1).y());
            EXPECT_EQ(0, plane(v1).z());
        }
        {
            const auto plane = Plane3D<>(Y<>(),X<>());
            const auto p2 = plane(v1) * plane(v2);
            EXPECT_TRUE(abs(p1 - p2) < epsilon);
            EXPECT_EQ(v1.x(), plane(v1).y());
            EXPECT_EQ(v1.y(), plane(v1).x());
            EXPECT_EQ(0, plane(v1).z());
        }
        {
            const auto plane = Plane3D<>(X<>(),Z<>());
            const auto p2 = plane(v1) * plane(v2);
            EXPECT_TRUE(abs(p1 - p2) < epsilon);
            EXPECT_EQ(v1.x(), plane(v1).x());
            EXPECT_EQ(v1.y(), plane(v1).z());
            EXPECT_EQ(0, plane(v1).y());
        }
        {
            const auto plane = Plane3D<>(Z<>(),X<>());
            const auto p2 = plane(v1) * plane(v2);
            EXPECT_TRUE(abs(p1 - p2) < epsilon);
            EXPECT_EQ(v1.x(), plane(v1).z());
            EXPECT_EQ(v1.y(), plane(v1).x());
            EXPECT_EQ(0, plane(v1).y());
        }
        {
            const auto plane = Plane3D<>(Y<>(),Z<>());
            const auto p2 = plane(v1) * plane(v2);
            EXPECT_TRUE(abs(p1 - p2) < epsilon);
            EXPECT_EQ(v1.x(), plane(v1).y());
            EXPECT_EQ(v1.y(), plane(v1).z());
            EXPECT_EQ(0, plane(v1).x());
        }
        {
            const auto plane = Plane3D<>(Z<>(),Y<>());
            const auto p2 = plane(v1) * plane(v2);
            EXPECT_TRUE(abs(p1 - p2) < epsilon);
            EXPECT_EQ(v1.x(), plane(v1).z());
            EXPECT_EQ(v1.y(), plane(v1).y());
            EXPECT_EQ(0, plane(v1).x());
        }
    }
}
TEST(Plane3D, norm)
{
    RANDOM RG;
    for (size_t i = 0; i < 100; i++) {
	const  auto n=randomIsotropic<3>(RG);
        const auto plane = Plane3D<>::ByNormalVector(n);
	EXPECT_TRUE(abs(plane.x().M()-1)<epsilon);
	EXPECT_TRUE(abs(plane.y().M()-1)<epsilon);
        EXPECT_TRUE(abs(plane.x()*plane.y()) < epsilon);
        EXPECT_TRUE(abs(plane.x()*(n*1.0)) < epsilon);
        EXPECT_TRUE(abs(plane.y()*(n*1.0)) < epsilon);
        EXPECT_TRUE(((n*1.0)^plane.x()).CloseTo(plane.y(),epsilon));
	if ((n*1.0)==Z<>()) {
	    EXPECT_TRUE(plane.x().CloseTo(X<>(),epsilon));
        } else {
	    EXPECT_TRUE(plane.x().CloseTo(direction((n*1.0)^Z<>())*1.0,epsilon));
        }
    }
}
TEST(Plane3D, norm2)
{
    RANDOM RG;
    for (size_t i = 0; i < 100; i++) {
	const  auto n=randomIsotropic<3>(RG);
	const auto rot=randomIsotropic<2>(RG);
        const auto plane = Plane3D<>::ByNormalVector(n,rot);
	EXPECT_TRUE(abs(plane.x().M()-1)<epsilon);
	EXPECT_TRUE(abs(plane.y().M()-1)<epsilon);
        EXPECT_TRUE(abs(plane.x()*plane.y()) < epsilon);
        EXPECT_TRUE(abs(plane.x()*(n*1.0)) < epsilon);
        EXPECT_TRUE(abs(plane.y()*(n*1.0)) < epsilon);
        EXPECT_TRUE(((n*1.0)^plane.x()).CloseTo(plane.y(),epsilon));
	if (((n*1.0) ^ Z<>()).M() == 0) {
	    EXPECT_TRUE(plane.x().CloseTo(Rotate(X<>(),n,rot.phi()),epsilon));
        } else {
	    EXPECT_TRUE(plane.x().CloseTo(Rotate(direction((n*1.0)^Z<>())*1.0,n,rot.phi()),epsilon));
        }
    }
}
TEST(Plane3D, scalarProductInvariance)
{
    RANDOM RG;
    RandomUniform<> mr(0, 10);
    for (size_t i = 0; i < 10000; i++) {
        const auto v1 = randomIsotropic<2>(RG) * mr(RG), v2 = randomIsotropic<2>(RG) * mr(RG);
        const auto p1 = v1 * v2;
        const auto plane = Plane3D<>::ByNormalVector(randomIsotropic<3>(RG),randomIsotropic<2>(RG));
        const auto p2 = plane(v1) * plane(v2);
        EXPECT_TRUE(abs(p1 - p2) < epsilon);
    }
}
TEST(Plane3D, vectorProductInvariance)
{
    RANDOM RG;
    RandomUniform<> mr(0, 10);
    for (size_t i = 0; i < 10000; i++) {
        const auto v1 = randomIsotropic<2>(RG) * mr(RG), v2 = randomIsotropic<2>(RG) * mr(RG);
        const auto p1 = abs(v1 ^ v2);
        const auto plane = Plane3D<>::ByNormalVector(randomIsotropic<3>(RG), randomIsotropic<2>(RG));
        const auto p2 = (plane(v1)^plane(v2)).M();
        EXPECT_TRUE(abs(p1 - p2) < epsilon);
    }
}
TEST(Plane3D, LorentzInvariance)
{
    RANDOM RG;
    RandomUniform<> mr(0, 10);
    for (size_t i = 0; i < 10000; i++) {
        const auto v = lorentz_byPM(randomIsotropic<2>(RG) * mr(RG), mr(RG));
        const auto plane = Plane3D<>::ByNormalVector(randomIsotropic<3>(RG), randomIsotropic<2>(RG));
        const auto v2 = lorentzVector(v.T(), plane(v.S()));
        EXPECT_TRUE(abs(v.M() - v2.M()) < epsilon);
    }
}
