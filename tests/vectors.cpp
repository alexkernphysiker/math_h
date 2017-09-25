// this file is distributed under
// MIT license
#include <gtest/gtest.h>
#include <math_h/vectors.h>
#include <math_h/hists.h>
#include <math_h/sigma.h>
using namespace std;
using namespace MathTemplates;
TEST(Vector2, scalar_prod)
{
    EXPECT_EQ(0, I<>()*J<>());
    EXPECT_EQ(0, J<>()*I<>());
    EXPECT_EQ(1, I<>()*I<>());
    EXPECT_EQ(1, J<>()*J<>());

    EXPECT_EQ(0, I<>() * (-J<>()));
    EXPECT_EQ(0, J<>() * (-I<>()));
    EXPECT_EQ(0, (-I<>())*J<>());
    EXPECT_EQ(0, (-J<>())*I<>());

    EXPECT_EQ(-1, I<>() * (-I<>()));
    EXPECT_EQ(-1, J<>() * (-J<>()));
    EXPECT_EQ(-1, (-I<>())*I<>());
    EXPECT_EQ(-1, (-J<>())*J<>());
}
TEST(Vector3, scalar_prod)
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
TEST(Vector3, vector_prod)
{
    EXPECT_EQ(X<>()^Y<>(), Z<>());
    EXPECT_EQ(X<>()^X<>(), Zero<>());
    EXPECT_EQ(Y<>()^X<>(), -Z<>());
    EXPECT_EQ(X<>() ^ (-X<>()), Zero<>());

    EXPECT_EQ(1, (X<>()^Y<>())*Z<>());
    EXPECT_EQ(1, (Y<>()^Z<>())*X<>());
    EXPECT_EQ(1, (Z<>()^X<>())*Y<>());

}
TEST(Vector3, Isotropic)
{
    RANDOM RG;
    const auto c = BinsByStep(-1.0, 0.1, 1.0);
    Distribution1D<> X(c), Y(c), Z(c);
    for (size_t i = 0; i < 100000; i++) {
        const auto V = RandomIsotropicDirection3<>(RG);
        X.Fill(V.x());
        Y.Fill(V.y());
        Z.Fill(V.z());
    }
    const auto x = X.TotalSum().val() / X.size();
    for (const auto &p : X)EXPECT_TRUE(p.Y().make_wider(2.5).Contains(x));
    for (const auto &p : Y)EXPECT_TRUE(p.Y().make_wider(2.5).Contains(x));
    for (const auto &p : Z)EXPECT_TRUE(p.Y().make_wider(2.5).Contains(x));
}
TEST(Vector2, Rotation)
{
    RANDOM RG;
    RandomUniform<> Phi(0, PI<>() * 2.0);
    for (size_t i = 0; i < 10000; i++) {
        const auto I = RandomIsotropicDirection2<>(RG);
        const auto ang = Phi(RG);
        const auto F = Rotate(I, ang);
    }
}

TEST(Vector3, RotationX)
{
    RANDOM RG;
    RandomUniform<> Phi(0, PI<>() * 2.0);
    for (size_t i = 0; i < 10000; i++) {
        const auto I = RandomIsotropicDirection3<>(RG);
        const auto axis = X<>();
        const auto ang = Phi(RG);
        const auto F = Rotate(I, axis, ang);
        EXPECT_TRUE(pow(axis * I - axis * F, 2) < 0.0000000000001);
    }
}
TEST(Vector3, RotationY)
{
    RANDOM RG;
    RandomUniform<> Phi(0, PI<>() * 2.0);
    for (size_t i = 0; i < 10000; i++) {
        const auto I = RandomIsotropicDirection3<>(RG);
        const auto axis = Y<>();
        const auto ang = Phi(RG);
        const auto F = Rotate(I, axis, ang);
        EXPECT_TRUE(pow(axis * I - axis * F, 2) < 0.0000000000001);
    }
}
TEST(Vector3, RotationZ)
{
    RANDOM RG;
    RandomUniform<> Phi(0, PI<>() * 2.0);
    for (size_t i = 0; i < 10000; i++) {
        const auto I = RandomIsotropicDirection3<>(RG);
        const auto axis = Z<>();
        const auto ang = Phi(RG);
        const auto F = Rotate(I, axis, ang);
        EXPECT_TRUE(pow(axis * I - axis * F, 2) < 0.0000000000001);
    }
}

TEST(LorentzVector, LorentzTransform)
{
    RANDOM RG;
    const double epsilon = 0.0000000000001;
    RandomUniform<> mr(0, 1.0 - epsilon), metrr(-5, 5);
    for (size_t i = 0; i < 10000; i++) {
        const auto V0 = lorentz_byPM(RandomIsotropicDirection3<>(RG), metrr(RG));
        const auto V1 = V0.Lorentz(Zero<>());
        EXPECT_TRUE(V1 == V0);
        const auto beta = RandomIsotropicDirection3<>(RG) * mr(RG);
        const auto V2 = V0.Lorentz(beta).Lorentz(-beta);
        EXPECT_TRUE(pow(V2.time_component() - V0.time_component(), 2) < epsilon);
        EXPECT_TRUE(pow(V2.space_component().x() - V0.space_component().x(), 2) < epsilon);
        EXPECT_TRUE(pow(V2.space_component().y() - V0.space_component().y(), 2) < epsilon);
        EXPECT_TRUE(pow(V2.space_component().z() - V0.space_component().z(), 2) < epsilon);
        EXPECT_THROW(V0.Lorentz(
                         RandomIsotropicDirection3<>(RG) * (1.0 + mr(RG))),
                     Exception<LorentzVector<>>
                    );
        const auto L0 = V0.length4(), L1 = V0.Lorentz(beta).length4(),
                   L2 = V0.Lorentz(beta).Lorentz(RandomIsotropicDirection3<>(RG) * 0.2).length4();
        EXPECT_TRUE(pow(L0 - L1, 2) < epsilon);
        EXPECT_TRUE(pow(L2 - L1, 2) < epsilon);
        EXPECT_TRUE(pow(L0 - L2, 2) < epsilon);
        const auto V00 = LorentzVector<>(V0.length4()).Lorentz(-V0.Beta());
        EXPECT_TRUE(pow(V00.time_component() - V0.time_component(), 2) < epsilon);
        EXPECT_TRUE(pow(V00.space_component().x() - V0.space_component().x(), 2) < epsilon);
        EXPECT_TRUE(pow(V00.space_component().y() - V0.space_component().y(), 2) < epsilon);
        EXPECT_TRUE(pow(V00.space_component().z() - V0.space_component().z(), 2) < epsilon);
    }
}
TEST(LorentzVector, LorentzTransform2)
{
    RANDOM RG;
    const double epsilon = 0.0000000000001;
    RandomUniform<> M(0, 5),P(0, 5);
    for (size_t i = 0; i < 10000; i++) {
        const auto V0 = lorentz_byPM(RandomIsotropicDirection3<>(RG)*P(RG), M(RG));
	const auto V1=V0.Lorentz(V0.Beta());
        EXPECT_TRUE(pow(V1.space_component().x(), 2) < epsilon);
        EXPECT_TRUE(pow(V1.space_component().y(), 2) < epsilon);
        EXPECT_TRUE(pow(V1.space_component().z(), 2) < epsilon);
        EXPECT_TRUE(pow(V1.length4()-V0.length4(), 2) < epsilon);
    }
}
TEST(LorentzVector, LorentzTransform2d)
{
    RANDOM RG;
    const double epsilon = 0.0000000000001;
    RandomUniform<> mr(0, 1.0 - epsilon), metrr(-5, 5);
    for (size_t i = 0; i < 10000; i++) {
        const auto V0 = lorentz_byPM(RandomIsotropicDirection2<>(RG), metrr(RG));
        const auto V1 = V0.Lorentz(zero<>());
        EXPECT_TRUE(V1 == V0);
        const auto beta = RandomIsotropicDirection2<>(RG) * mr(RG);
        const auto V2 = V0.Lorentz(beta).Lorentz(-beta);
        EXPECT_TRUE(pow(V2.time_component() - V0.time_component(), 2) < epsilon);
        EXPECT_TRUE(pow(V2.space_component().x() - V0.space_component().x(), 2) < epsilon);
        EXPECT_TRUE(pow(V2.space_component().y() - V0.space_component().y(), 2) < epsilon);
        typedef LorentzVector<double, Vector<2>> LVtype;
        EXPECT_THROW(V0.Lorentz(RandomIsotropicDirection2<>(RG) * (1.0 + mr(RG))),
                     Exception<LVtype>);
        const auto L0 = V0.length4(), L1 = V0.Lorentz(beta).length4(),
                   L2 = V0.Lorentz(beta).Lorentz(RandomIsotropicDirection2<>(RG) * 0.2).length4();
        EXPECT_TRUE(pow(L0 - L1, 2) < epsilon);
        EXPECT_TRUE(pow(L2 - L1, 2) < epsilon);
        EXPECT_TRUE(pow(L0 - L2, 2) < epsilon);
        const auto V00 = LorentzVector<double, Vector<2>>(V0.length4()).Lorentz(-V0.Beta());
        EXPECT_TRUE(pow(V00.time_component() - V0.time_component(), 2) < epsilon);
        EXPECT_TRUE(pow(V00.space_component().x() - V0.space_component().x(), 2) < epsilon);
        EXPECT_TRUE(pow(V00.space_component().y() - V0.space_component().y(), 2) < epsilon);
    }
}
TEST(Plane3D, SimplePlanes)
{
    RANDOM RG;
    const double epsilon = 0.0000000000001;
    RandomUniform<> mr(0, 10), th(0, 2.0 * PI());
    for (size_t i = 0; i < 10000; i++) {
        const auto v1 = RandomIsotropicDirection2<>(RG) * mr(RG), v2 = RandomIsotropicDirection2<>(RG) * mr(RG);
        const auto p1 = v1 * v2;
        {
            const auto plane = Plane3D<>::basis<1, 2>();
            const auto p2 = plane(v1) * plane(v2);
            EXPECT_TRUE(pow(p1 - p2, 2) < epsilon);
            EXPECT_EQ(v1.x(), plane(v1).x());
            EXPECT_EQ(v1.y(), plane(v1).y());
            EXPECT_EQ(0, plane(v1).z());
        }
        {
            const auto plane = Plane3D<>::basis<2, 1>();
            const auto p2 = plane(v1) * plane(v2);
            EXPECT_TRUE(pow(p1 - p2, 2) < epsilon);
            EXPECT_EQ(v1.x(), plane(v1).y());
            EXPECT_EQ(v1.y(), plane(v1).x());
            EXPECT_EQ(0, plane(v1).z());
        }
        {
            const auto plane = Plane3D<>::basis<1, 3>();
            const auto p2 = plane(v1) * plane(v2);
            EXPECT_TRUE(pow(p1 - p2, 2) < epsilon);
            EXPECT_EQ(v1.x(), plane(v1).x());
            EXPECT_EQ(v1.y(), plane(v1).z());
            EXPECT_EQ(0, plane(v1).y());
        }
        {
            const auto plane = Plane3D<>::basis<3, 1>();
            const auto p2 = plane(v1) * plane(v2);
            EXPECT_TRUE(pow(p1 - p2, 2) < epsilon);
            EXPECT_EQ(v1.x(), plane(v1).z());
            EXPECT_EQ(v1.y(), plane(v1).x());
            EXPECT_EQ(0, plane(v1).y());
        }
        {
            const auto plane = Plane3D<>::basis<2, 3>();
            const auto p2 = plane(v1) * plane(v2);
            EXPECT_TRUE(pow(p1 - p2, 2) < epsilon);
            EXPECT_EQ(v1.x(), plane(v1).y());
            EXPECT_EQ(v1.y(), plane(v1).z());
            EXPECT_EQ(0, plane(v1).x());
        }
        {
            const auto plane = Plane3D<>::basis<3, 2>();
            const auto p2 = plane(v1) * plane(v2);
            EXPECT_TRUE(pow(p1 - p2, 2) < epsilon);
            EXPECT_EQ(v1.x(), plane(v1).z());
            EXPECT_EQ(v1.y(), plane(v1).y());
            EXPECT_EQ(0, plane(v1).x());
        }
    }
}
TEST(Plane3D, scalarProductInvariance)
{
    RANDOM RG;
    const double epsilon = 0.0000000000001;
    RandomUniform<> mr(0, 10), th(0, 2.0 * PI());
    for (size_t i = 0; i < 10000; i++) {
        const auto v1 = RandomIsotropicDirection2<>(RG) * mr(RG), v2 = RandomIsotropicDirection2<>(RG) * mr(RG);
        const auto p1 = v1 * v2;
        const auto plane = Plane3D<>::ByNormalVectorAndTheta(RandomIsotropicDirection3<>(RG), th(RG));
        const auto p2 = plane(v1) * plane(v2);
        EXPECT_TRUE(pow(p1 - p2, 2) < epsilon);
    }
}
TEST(Plane3D, vectorProductInvariance)
{
    RANDOM RG;
    const double epsilon = 0.0000000000001;
    RandomUniform<> mr(0, 10), th(0, 2.0 * PI());
    for (size_t i = 0; i < 10000; i++) {
        const auto v1 = RandomIsotropicDirection2<>(RG) * mr(RG), v2 = RandomIsotropicDirection2<>(RG) * mr(RG);
        const auto p1 = pow(v1 ^ v2, 2);
        const auto plane = Plane3D<>::ByNormalVectorAndTheta(RandomIsotropicDirection3<>(RG), th(RG));
        const auto p2 = (plane(v1)^plane(v2)).mag_sqr();
        EXPECT_TRUE(pow(p1 - p2, 2) < epsilon);
    }
}
TEST(Plane3D, LorentzInvariance)
{
    RANDOM RG;
    const double epsilon = 0.0000000000001;
    RandomUniform<> mr(0, 10), th(0, 2.0 * PI());
    for (size_t i = 0; i < 10000; i++) {
        const auto v = lorentz_byPM(RandomIsotropicDirection2<>(RG) * mr(RG), mr(RG));
        const auto plane = Plane3D<>::ByNormalVectorAndTheta(RandomIsotropicDirection3<>(RG), th(RG));
        const auto v2 = lorentzVector(v.time_component(), plane(v.space_component()));
        EXPECT_TRUE(pow(v.length4() - v2.length4(), 2) < epsilon);
    }
}
TEST(LorentzVector, decays)
{
    RANDOM RG;
    const double epsilon = 0.0000000000001;
    RandomUniform<> IM(2,3),M1(0,1),M2(0,1),THETA(0,PI()),PHI(0, 2.0 * PI());
    for (size_t i = 0; i < 10000; i++) {
	const double im=IM(RG),m1=M1(RG),m2=M2(RG);
	const auto C1=binaryDecay(im,m1,m2);
	EXPECT_TRUE(pow((C1.first+C1.second).length4()-im,2)<epsilon);
	EXPECT_TRUE(pow((C1.first).length4()-m1,2)<epsilon);
	EXPECT_TRUE(pow((C1.second).length4()-m2,2)<epsilon);
	EXPECT_TRUE((C1.first.space_component()+C1.second.space_component()).mag_sqr()<epsilon);
	const auto C2=binaryDecay(im,m1,m2,PHI(RG));
	EXPECT_TRUE(pow((C2.first+C2.second).length4()-im,2)<epsilon);
	EXPECT_TRUE(pow((C2.first).length4()-m1,2)<epsilon);
	EXPECT_TRUE(pow((C2.second).length4()-m2,2)<epsilon);
	EXPECT_TRUE((C2.first.space_component()+C2.second.space_component()).mag_sqr()<epsilon);
	const auto C3=binaryDecay(im,m1,m2,THETA(RG),PHI(RG));
	EXPECT_TRUE(pow((C3.first+C3.second).length4()-im,2)<epsilon);
	EXPECT_TRUE(pow(C3.first.length4()-m1,2)<epsilon);
	EXPECT_TRUE(pow(C3.second.length4()-m2,2)<epsilon);
	EXPECT_TRUE((C3.first.space_component()+C3.second.space_component()).mag_sqr()<epsilon);
    }
}
TEST(LorentzVector, decays2)
{
    RANDOM RG;
    const double epsilon = 0.0000000000001;
    RandomUniform<> IM(2,3),M1(0,1),M2(0,1),THETA(0,PI()),PHI(0, 2.0 * PI());
    for (size_t i = 0; i < 10000; i++) {
	const double im=IM(RG),m1=M1(RG),m2=M2(RG);
	const auto C2=binaryDecay(im,m1,m2,RandomIsotropicDirection2<>(RG));
	EXPECT_TRUE(pow((C2.first+C2.second).length4()-im,2)<epsilon);
	EXPECT_TRUE(pow((C2.first).length4()-m1,2)<epsilon);
	EXPECT_TRUE(pow((C2.second).length4()-m2,2)<epsilon);
	EXPECT_TRUE((C2.first.space_component()+C2.second.space_component()).mag_sqr()<epsilon);
	const auto C3=binaryDecay(im,m1,m2,RandomIsotropicDirection3<>(RG));
	EXPECT_TRUE(pow((C3.first+C3.second).length4()-im,2)<epsilon);
	EXPECT_TRUE(pow(C3.first.length4()-m1,2)<epsilon);
	EXPECT_TRUE(pow(C3.second.length4()-m2,2)<epsilon);
	EXPECT_TRUE((C3.first.space_component()+C3.second.space_component()).mag_sqr()<epsilon);
    }
}

