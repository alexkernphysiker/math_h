
// this file is distributed under
// LGPLv3 license
#include <gtest/gtest.h>
#include <math_h/vectortransformations.h>
#include <math_h/hists.h>
#include <math_h/sigma.h>
using namespace std;
using namespace MathTemplates;
const double epsilon = 0.0000000001;
TEST(VectorTransformation, pseudoscalar_prod)
{
    EXPECT_EQ(x()^y(), vec(1.));
    EXPECT_EQ(y()^x(), vec(-1.));
    EXPECT_EQ(x()^x(), vec(0.));
    EXPECT_EQ(y()^y(), vec(0.));
    EXPECT_EQ((-x())^x(), vec(0.));
    EXPECT_EQ((-y())^y(), vec(0.));
}
TEST(VectorTransformation, vector_prod3_basis)
{
    EXPECT_EQ(X()^Y(), Z());
    EXPECT_EQ(X()^X(), Zero());
    EXPECT_EQ(Y()^X(), -Z());
    EXPECT_EQ(X() ^ (-X()), Zero());

    EXPECT_EQ(1, (X()^Y())*Z());
    EXPECT_EQ(1, (Y()^Z())*X());
    EXPECT_EQ(1, (Z()^X())*Y());

}

TEST(VectorTransformation, pseudoscalar_prod_algebra)
{
    
    RandomUniform<> M(-50, 50);
    for (size_t i = 0; i < 50; i++) {
        const auto
        A = randomIsotropic<2>() * M(),
        B = randomIsotropic<2>() * M(),
        C = randomIsotropic<2>() * M();
        const auto a = M();
        EXPECT_TRUE((A ^ A).M() < epsilon);
        EXPECT_TRUE((B ^ B).M() < epsilon);
        EXPECT_TRUE(((A ^ B) + (B ^ A)).M() < epsilon);
        EXPECT_TRUE(((A ^ B)*a - ((A * a)^B)).M() < epsilon);
        EXPECT_TRUE(((A ^ B)*a - (A ^ (B * a))).M() < epsilon);
        EXPECT_TRUE((((A + B)^C) - ((A ^ C) + (B ^ C))).M() < epsilon);
    }
}

TEST(VectorTransformation, vector_prod3_algebra)
{
    
    RandomUniform<> M(-50, 50);
    for (size_t i = 0; i < 50; i++) {
        const auto
        A = randomIsotropic<3>() * M(),
        B = randomIsotropic<3>() * M(),
        C = randomIsotropic<3>() * M();
        const auto a = M();
        EXPECT_TRUE((A ^ A).M() < epsilon);
        EXPECT_TRUE((B ^ B).M() < epsilon);
        EXPECT_TRUE((C ^ C).M() < epsilon);
        EXPECT_TRUE(((A ^ B) + (B ^ A)).M() < epsilon);
        EXPECT_TRUE(((A ^ C) + (C ^ A)).M() < epsilon);
        EXPECT_TRUE(((C ^ B) + (B ^ C)).M() < epsilon);
        EXPECT_TRUE(((A ^ B)*a - ((A * a)^B)).M() < epsilon);
        EXPECT_TRUE(((A ^ B)*a - (A ^ (B * a))).M() < epsilon);
        EXPECT_TRUE((((A + B)^C) - ((A ^ C) + (B ^ C))).M() < epsilon);
        EXPECT_TRUE((((A ^ B)^C) + ((C ^ A)^B) + ((B ^ C)^A)).M() < epsilon);
        EXPECT_TRUE(((A ^ (B ^ C)) - (B * (A * C)) + (C * (A * B))).M() < epsilon);
        EXPECT_TRUE(abs(((A ^ B)*C) - (A * (B ^ C))) < epsilon);
    }
}
TEST(VectorTransformation, vector_prod7_algebra)
{
    
    RandomUniform<> M(-50, 50);
    for (size_t i = 0; i < 50; i++) {
        const auto
        A = randomIsotropic<7>() * M(),
        B = randomIsotropic<7>() * M(),
        C = randomIsotropic<7>() * M();
        const auto a = M();
        EXPECT_TRUE((A ^ A).M() < epsilon);
        EXPECT_TRUE((B ^ B).M() < epsilon);
        EXPECT_TRUE((C ^ C).M() < epsilon);
        EXPECT_TRUE(((A ^ B) + (B ^ A)).M() < epsilon);
        EXPECT_TRUE(((A ^ C) + (C ^ A)).M() < epsilon);
        EXPECT_TRUE(((C ^ B) + (B ^ C)).M() < epsilon);
        EXPECT_TRUE(((A ^ B)*a - ((A * a)^B)).M() < epsilon);
        EXPECT_TRUE(((A ^ B)*a - (A ^ (B * a))).M() < epsilon);
        EXPECT_TRUE((((A + B)^C) - ((A ^ C) + (B ^ C))).M() < epsilon);
        EXPECT_TRUE(abs(((A ^ B)*C) - (A * (B ^ C))) < epsilon);
    }
}
TEST(Direction, base2d)
{
    
    RandomUniform<> Phi(-PI(), PI());
    RandomUniform<> A(0,10);
    for (size_t i = 0; i < 50; i++) {
        const auto phi = Phi();
        const auto D = direction(phi);
        EXPECT_EQ(D.phi(), phi);
	const auto a=A();
	EXPECT_EQ(D*a,vec(
	    a*cos(phi),
	    a*sin(phi)
	));
    }
}
TEST(Direction, base3d)
{
    
    RandomUniform<> Phi(-PI(), PI());
    RandomUniform<> Theta(0, PI());
    RandomUniform<> A(0,10);
    for (size_t i = 0; i < 50; i++) {
        const auto phi = Phi();
        const auto theta = Theta();
        const auto D = direction(phi, theta);
        EXPECT_EQ(D.phi(), phi);
        EXPECT_EQ(D.th(), theta);
        EXPECT_EQ(D.th<1>(), theta);
	const auto a=A();
	EXPECT_EQ(D*a,vec(
	    a*sin(theta)*cos(phi),
	    a*sin(theta)*sin(phi),
	    a*cos(theta)
	));
    }
}
TEST(Direction, base4d)
{
    
    RandomUniform<> Phi(-PI(), PI());
    RandomUniform<> Theta(0, PI<>());
    RandomUniform<> A(0,10);
    for (size_t i = 0; i < 50; i++) {
        const auto phi = Phi();
        const auto theta1 = Theta();
        const auto theta2 = Theta();
        const auto D = direction(phi, theta1, theta2);
        EXPECT_EQ(D.phi(), phi);
        EXPECT_EQ(D.th<1>(), theta1);
        EXPECT_EQ(D.th<2>(), theta2);
	const auto a=A();
	EXPECT_EQ(D*a,vec(
	    a*sin(theta2)*sin(theta1)*cos(phi),
	    a*sin(theta2)*sin(theta1)*sin(phi),
	    a*sin(theta2)*cos(theta1),
	    a*cos(theta2)
	));
    }
}
TEST(Direction, base5d)
{
    
    RandomUniform<> Phi(-PI(), PI());
    RandomUniform<> Theta(0, PI());
    RandomUniform<> A(0,10);
    for (size_t i = 0; i < 50; i++) {
        const auto phi = Phi();
        const auto theta1 = Theta();
        const auto theta2 = Theta();
        const auto theta3 = Theta();
        const auto D = direction(phi, theta1, theta2, theta3);
        EXPECT_EQ(D.phi(), phi);
        EXPECT_EQ(D.th<1>(), theta1);
        EXPECT_EQ(D.th<2>(), theta2);
        EXPECT_EQ(D.th<3>(), theta3);
	const auto a=A();
	EXPECT_EQ(D*a,vec(
	    a*sin(theta3)*sin(theta2)*sin(theta1)*cos(phi),
	    a*sin(theta3)*sin(theta2)*sin(theta1)*sin(phi),
	    a*sin(theta3)*sin(theta2)*cos(theta1),
	    a*sin(theta3)*cos(theta2),
	    a*cos(theta3)
	));
    }
}

TEST(VectorTransformation, Isotropic1)
{
    
    const auto c = BinsByCount(2, -2.0, 2.0);
    Distribution1D<> D(c);
    for (size_t i = 0; i < 1000000; i++) {
        D.Fill(randomIsotropic<1>().dir());
    }
    const auto x = D.TotalSum().val() / D.size();
    for (const auto &p : D)EXPECT_TRUE(p.Y().make_wider(2.5).Contains(x));
}
TEST(VectorTransformation, Isotropic2)
{
    
    const auto c = BinsByCount(10, -PI(), PI());
    Distribution1D<> Phi(c);
    for (size_t i = 0; i < 1000000; i++) {
        Phi.Fill(randomIsotropic<2>().phi());
    }
    const auto x = Phi.TotalSum().val() / Phi.size();
    for (const auto &p : Phi)EXPECT_TRUE(p.Y().make_wider(2.5).Contains(x));
}
TEST(VectorTransformation, Isotropic3)
{
    
    const auto c = BinsByStep(-1.0, 0.1, 1.0);
    Distribution1D<> X(c), Y(c), Z(c);
    for (size_t i = 0; i < 1000000; i++) {
        const auto V = randomIsotropic<3>() * 1.0;
        X.Fill(V.x());
        Y.Fill(V.y());
        Z.Fill(V.z());
    }
    const auto x = X.TotalSum().val() / X.size();
    for (const auto &p : X)EXPECT_TRUE(p.Y().make_wider(2.5).Contains(x));
    for (const auto &p : Y)EXPECT_TRUE(p.Y().make_wider(2.5).Contains(x));
    for (const auto &p : Z)EXPECT_TRUE(p.Y().make_wider(2.5).Contains(x));
}
TEST(VectorTransformation, Rotation2)
{
    
    RandomUniform<> Phi(0, PI() * 2.0);
    for (size_t i = 0; i < 50; i++) {
        const auto I = randomIsotropic<2>() * 1.0;
        const auto ang = Phi();
        const auto F = Rotation(ang) * I;
        EXPECT_TRUE(abs(I.M() - F.M()) < epsilon);
    }
}

TEST(VectorTransformation, Rotation3)
{
    
    RandomUniform<> M(0, 5.0), Phi(-PI(), PI());
    for (size_t i = 0; i < 50; i++) {
        const auto I = randomIsotropic<3>() * M();
        const auto ang = Phi();
        const auto dir = randomIsotropic<3>();
        const auto F = Rotation(dir, ang) * I;
        EXPECT_TRUE(abs((dir * 1.0) * I - (dir * 1.0) * F) < epsilon);
    }
}

TEST(VectorTransformation, decompose2)
{
    
    RandomUniform<> M(0, 5.0);
    for (size_t i = 0; i < 50; i++) {
        const auto chosen_direction = randomIsotropic<2>();
        const auto vector_to_decompose = randomIsotropic<2>() * M();
        const auto decomposition = decompose_by_direction(vector_to_decompose, chosen_direction);
        EXPECT_TRUE(abs(decomposition.n * (chosen_direction * 1.0)) < epsilon);
        EXPECT_TRUE(vector_to_decompose.CloseTo(decomposition.tau + decomposition.n, epsilon));
    }
}
TEST(VectorTransformation, decompose3)
{
    
    RandomUniform<> M(0, 5.0);
    for (size_t i = 0; i < 50; i++) {
        const auto chosen_direction = randomIsotropic<3>();
        const auto vector_to_decompose = randomIsotropic<3>() * M();
        const auto decomposition = decompose_by_direction(vector_to_decompose, chosen_direction);
        EXPECT_TRUE(abs(decomposition.n * (chosen_direction * 1.0)) < epsilon);
        EXPECT_TRUE(vector_to_decompose.CloseTo(decomposition.tau + decomposition.n, epsilon));
    }
}
TEST(VectorTransformation, decompose4)
{
    
    RandomUniform<> M(0, 5.0);
    for (size_t i = 0; i < 50; i++) {
        const auto chosen_direction = randomIsotropic<4>();
        const auto vector_to_decompose = randomIsotropic<4>() * M();
        const auto decomposition = decompose_by_direction(vector_to_decompose, chosen_direction);
        EXPECT_TRUE(abs(decomposition.n * (chosen_direction * 1.0)) < epsilon);
        EXPECT_TRUE(vector_to_decompose.CloseTo(decomposition.tau + decomposition.n, epsilon));
    }
}

TEST(VectorTransformation, direction1d)
{
    
    RandomUniform<> V(-10.0, 10.0);
    for (size_t i = 0; i < 50; i++) {
        const auto d = direction();
        const auto v = vec(V());
        EXPECT_TRUE((v - d * v.x()).M() < epsilon);
    }
}
TEST(VectorTransformation, direction2d)
{
    
    RandomUniform<> PHI(-PI(), PI()), M(0.0, 10.0);
    for (size_t i = 0; i < 50; i++) {
        const auto phi = PHI();
        const auto phi2 = direction(direction(phi) * M()).phi();
        EXPECT_TRUE(abs(phi - phi2) < epsilon);
    }
}
TEST(VectorTransformation, direction3d)
{
    
    RandomUniform<> THETA(0.0, PI()), PHI(-PI(), PI()), M(0.0, 10.0);
    for (size_t i = 0; i < 50; i++) {
        const auto angles = direction(PHI(), THETA());
        const auto angles2 = direction(angles * M());
        EXPECT_TRUE(abs(angles.phi() - angles2.phi()) < epsilon);
        EXPECT_TRUE(abs(angles.th() - angles2.th()) < epsilon);
    }
}
TEST(VectorTransformation, direction4d)
{
    
    RandomUniform<> THETA(0.0, PI()), PHI(-PI(), PI()), M(0.0, 10.0);
    for (size_t i = 0; i < 50; i++) {
        const auto angles = direction(PHI(), THETA(), THETA());
        const auto angles2 = direction(angles * M());
        EXPECT_TRUE(abs(angles.phi() - angles2.phi()) < epsilon);
        EXPECT_TRUE(abs(angles.th<1>() - angles2.th<1>()) < epsilon);
        EXPECT_TRUE(abs(angles.th<2>() - angles2.th<2>()) < epsilon);
    }
}
TEST(VectorTransformation, direction5d)
{
    
    RandomUniform<> THETA(0.0, PI()), PHI(-PI(), PI()), M(0.0, 10.0);
    for (size_t i = 0; i < 50; i++) {
        const auto angles = direction(PHI(), THETA(), THETA(), THETA());
        const auto angles2 = direction(angles * M());
        EXPECT_TRUE(abs(angles.phi() - angles2.phi()) < epsilon);
        EXPECT_TRUE(abs(angles.th<1>() - angles2.th<1>()) < epsilon);
        EXPECT_TRUE(abs(angles.th<2>() - angles2.th<2>()) < epsilon);
        EXPECT_TRUE(abs(angles.th<3>() - angles2.th<3>()) < epsilon);
    }
}

TEST(VectorTransformation, rotations1d)
{
    
    RandomUniform<> M(0.0, 10.0);
    for (size_t i = 0; i < 50; i++) {
        const auto dir = randomIsotropic<1>();
        const auto m = M();
        const auto V1 = dir * m;
        const auto V2 = dir.Rotations() * (axis<1, 1>() * m);
        EXPECT_TRUE(V1.CloseTo(V2, epsilon));
        const auto V3 = dir.AntiRotations() * V2;
        EXPECT_TRUE(V3.CloseTo(axis<1, 1>()*m, epsilon));
    }
}
TEST(VectorTransformation, rotations2d)
{
    
    RandomUniform<> M(0.0, 10.0);
    for (size_t i = 0; i < 50; i++) {
        const auto dir = randomIsotropic<2>();
        const auto m = M();
        const auto V1 = dir * m;
        const auto V2 = dir.Rotations() * (x() * m);
        EXPECT_TRUE(V1.CloseTo(V2, epsilon));
        const auto V3 = dir.AntiRotations() * V2;
        EXPECT_TRUE(V3.CloseTo(x()*m, epsilon));
    }
}
TEST(VectorTransformation, rotations3d)
{
    
    RandomUniform<> M(0.0, 10.0);
    for (size_t i = 0; i < 50; i++) {
        const auto dir = randomIsotropic<3>();
        const auto m = M();
        const auto V1 = dir * m;
        const auto V2 = dir.Rotations() * (Z() * m);
        EXPECT_TRUE(V1.CloseTo(V2, epsilon));
        const auto V3 = dir.AntiRotations() * V2;
        EXPECT_TRUE(V3.CloseTo(Z()*m, epsilon));
    }
}
TEST(VectorTransformation, rotations4d)
{
    
    RandomUniform<> M(0.0, 10.0);
    for (size_t i = 0; i < 50; i++) {
        const auto dir = randomIsotropic<4>();
        const auto m = M();
        const auto V1 = dir * m;
        const auto V2 = dir.Rotations() * (axis<4, 4>() * m);
        EXPECT_TRUE(V1.CloseTo(V2, epsilon));
        const auto V3 = dir.AntiRotations() * V2;
        EXPECT_TRUE(V3.CloseTo(axis<4, 4>()*m, epsilon));
    }
}
TEST(VectorTransformation, rotations5d)
{
    
    RandomUniform<> M(0.0, 10.0);
    for (size_t i = 0; i < 50; i++) {
        const auto dir = randomIsotropic<5>();
        const auto m = M();
        const auto V1 = dir * m;
        const auto V2 = dir.Rotations() * (axis<5, 5>() * m);
        EXPECT_TRUE(V1.CloseTo(V2, epsilon));
        const auto V3 = dir.AntiRotations() * V2;
        EXPECT_TRUE(V3.CloseTo(axis<5, 5>()*m, epsilon));
    }
}

TEST(VectorTransformation, rotations3d_plane)
{
    
    RandomUniform<> M(-10.0, 10.0), TH(-PI<>(), PI<>());
    for (size_t i = 0; i < 50; i++) {
        const auto R = randomIsotropic<3>().Rotations();
        const auto v1 = R * vec(M(), M(), 0.);
        const auto v2 = R * vec(M(), M(), 0.);
        const auto v3 = R * vec(M(), M(), 0.);
        EXPECT_TRUE(abs((v1 ^ v2)*v3) < epsilon);
    }
    for (size_t i = 0; i < 50; i++) {
        const auto R = randomIsotropic<3>().Rotations();
        const auto v1 = R * vec(0., M(), M());
        const auto v2 = R * vec(0., M(), M());
        const auto v3 = R * vec(0., M(), M());
        EXPECT_TRUE(abs((v1 ^ v2)*v3) < epsilon);
    }
    for (size_t i = 0; i < 10; i++) {
        const auto R = randomIsotropic<3>().Rotations();
        const auto v1 = R * vec(M(), 0., M());
        const auto v2 = R * vec(M(), 0., M());
        const auto v3 = R * vec(M(), 0., M());
        EXPECT_TRUE(abs((v1 ^ v2)*v3) < epsilon);
    }
    for (size_t i = 0; i < 50; i++) {
        const auto R1 = randomIsotropic<3>().Rotations();
        const auto v1 = vec(M(), 0., M());
        const auto R2 = Rotation(direction(v1), TH());
        const auto R = R1 * R2;
        const auto v2 = vec(M(), 0., M());
        const auto v3 = vec(M(), 0., M());
        EXPECT_TRUE(abs(((R * v1) ^ (R * v2)) * (R * v3)) < epsilon);
    }
}
