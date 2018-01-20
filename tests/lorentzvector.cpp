
// this file is distributed under
// LGPLv3 license
#include <gtest/gtest.h>
#include <math_h/lorentzvector.h>
#include <math_h/hists.h>
#include <math_h/sigma.h>
using namespace std;
using namespace MathTemplates;
const double epsilon = 0.0000000001;

TEST(LorentzVector, LorentzTransform1d)
{
    RANDOM RG;
    RandomUniform<> mr(0, 0.99), metrr(-5, 5);
    for (size_t i = 0; i < 50; i++) {
        const auto V0 = lorentz_byPM(randomIsotropic<1>(RG) * 1.0, metrr(RG));
	const auto V0_copy=V0;
	EXPECT_EQ(V0,V0_copy);
        const auto V1 = V0.Transform(desCartes(0.));
        EXPECT_TRUE(V1 == V0);
        const auto beta = randomIsotropic<1>(RG) * mr(RG);
        const auto V2 = V0.Transform(beta).Transform(-beta);
        EXPECT_TRUE(abs(V2.E() - V0.E()) < epsilon);
        EXPECT_TRUE(V2.P().CloseTo(V0.P(), epsilon));
        EXPECT_ANY_THROW(V0.Transform(randomIsotropic<1>(RG) * (1.0 + mr(RG))));
        const auto L0 = V0.M(), L1 = V0.Transform(beta).M(),
                   L2 = V0.Transform(beta).Transform(randomIsotropic<1>(RG) * 0.2).M();
        EXPECT_TRUE(abs(L0 - L1) < epsilon);
        EXPECT_TRUE(abs(L2 - L1) < epsilon);
        EXPECT_TRUE(abs(L0 - L2) < epsilon);
        const auto V00 = lorentz_byPM(desCartes(0.), V0.M()).Transform(-V0.Beta());
        EXPECT_TRUE(abs(V00.E() - V0.E()) < epsilon);
        EXPECT_TRUE(V00.P().CloseTo(V0.P(), epsilon));
    }
}
TEST(LorentzVector, LorentzTransform2d)
{
    RANDOM RG;
    RandomUniform<> mr(0, 0.99), metrr(-5, 5);
    for (size_t i = 0; i < 50; i++) {
        const auto V0 = lorentz_byPM(randomIsotropic<2>(RG) * 1.0, metrr(RG));
	const auto V0_copy=V0;
	EXPECT_EQ(V0,V0_copy);
        const auto V1 = V0.Transform(zero());
        EXPECT_TRUE(V1 == V0);
        const auto beta = randomIsotropic<2>(RG) * mr(RG);
        const auto V2 = V0.Transform(beta).Transform(-beta);
        EXPECT_TRUE(abs(V2.E() - V0.E()) < epsilon);
        EXPECT_TRUE(V2.P().CloseTo(V0.P(), epsilon));
        EXPECT_ANY_THROW(V0.Transform(randomIsotropic<2>(RG) * (1.0 + mr(RG))));
        const auto L0 = V0.M(), L1 = V0.Transform(beta).M(),
                   L2 = V0.Transform(beta).Transform(randomIsotropic<2>(RG) * 0.2).M();
        EXPECT_TRUE(abs(L0 - L1) < epsilon);
        EXPECT_TRUE(abs(L2 - L1) < epsilon);
        EXPECT_TRUE(abs(L0 - L2) < epsilon);
        const auto V00 = lorentz_byPM(zero(), V0.M()).Transform(-V0.Beta());
        EXPECT_TRUE(abs(V00.E() - V0.E()) < epsilon);
        EXPECT_TRUE(V00.P().CloseTo(V0.P(), epsilon));
    }
}
TEST(LorentzVector, LorentzTransform3d)
{
    RANDOM RG;
    RandomUniform<> mr(0, 0.99), metrr(-5, 5);
    for (size_t i = 0; i < 50; i++) {
        const auto V0 = lorentz_byPM(randomIsotropic<3>(RG) * 1.0, metrr(RG));
	const auto V0_copy=V0;
	EXPECT_EQ(V0,V0_copy);
        const auto V1 = V0.Transform(Zero());
        EXPECT_TRUE(V1 == V0);
        const auto beta = randomIsotropic<3>(RG) * mr(RG);
        const auto V2 = V0.Transform(beta).Transform(-beta);
        EXPECT_TRUE(abs(V2.E() - V0.E()) < epsilon);
        EXPECT_TRUE(V2.P().CloseTo(V0.P(), epsilon));
        EXPECT_THROW(V0.Transform(randomIsotropic<3>(RG) * (1.0 + mr(RG))), Exception<LorentzVector<>>);
        const auto L0 = V0.M(), L1 = V0.Transform(beta).M(),
                   L2 = V0.Transform(beta).Transform(randomIsotropic<3>(RG) * 0.2).M();
        EXPECT_TRUE(abs(L0 - L1) < epsilon);
        EXPECT_TRUE(abs(L2 - L1) < epsilon);
        EXPECT_TRUE(abs(L0 - L2) < epsilon);
        const auto V00 = lorentz_byPM(Zero(), V0.M()).Transform(-V0.Beta());
        EXPECT_TRUE(abs(V00.E() - V0.E()) < epsilon);
        EXPECT_TRUE(V00.P().CloseTo(V0.P(), epsilon));
    }
}
TEST(LorentzVector, LorentzTransform4d)
{
    RANDOM RG;
    RandomUniform<> mr(0, 0.99), metrr(-5, 5);
    for (size_t i = 0; i < 50; i++) {
        const auto V0 = lorentz_byPM(randomIsotropic<4>(RG) * 1.0, metrr(RG));
	const auto V0_copy=V0;
	EXPECT_EQ(V0,V0_copy);
        const auto V1 = V0.Transform(desCartes(0., 0., 0., 0.));
        EXPECT_TRUE(V1 == V0);
        const auto beta = randomIsotropic<4>(RG) * mr(RG);
        const auto V2 = V0.Transform(beta).Transform(-beta);
        EXPECT_TRUE(abs(V2.E() - V0.E()) < epsilon);
        EXPECT_TRUE(V2.P().CloseTo(V0.P(), epsilon));
        EXPECT_ANY_THROW(V0.Transform(randomIsotropic<4>(RG) * (1.0 + mr(RG))));
        const auto L0 = V0.M(), L1 = V0.Transform(beta).M(),
                   L2 = V0.Transform(beta).Transform(randomIsotropic<4>(RG) * 0.2).M();
        EXPECT_TRUE(abs(L0 - L1) < epsilon);
        EXPECT_TRUE(abs(L2 - L1) < epsilon);
        EXPECT_TRUE(abs(L0 - L2) < epsilon);
        const auto V00 = lorentz_byPM(desCartes(0., 0., 0., 0.), V0.M()).Transform(-V0.Beta());
        EXPECT_TRUE(abs(V00.E() - V0.E()) < epsilon);
        EXPECT_TRUE(V00.P().CloseTo(V0.P(), epsilon));
    }
}
TEST(LorentzVector, LorentzTransform1d_more)
{
    RANDOM RG;
    RandomUniform<> M(0, 5), P(0, 5);
    for (size_t i = 0; i < 50; i++) {
        const auto V0 = lorentz_byPM(randomIsotropic<1>(RG) * P(RG), M(RG));
	const auto V0_copy=V0;
	EXPECT_EQ(V0,V0_copy);
        const auto V1 = V0.Transform(V0.Beta());
        EXPECT_TRUE(V1.P().CloseTo(desCartes(0.), epsilon));
        EXPECT_TRUE(abs(V1.M() - V0.M()) < epsilon);
    }
}
TEST(LorentzVector, LorentzTransform2d_more)
{
    RANDOM RG;
    RandomUniform<> M(0, 5), P(0, 5);
    for (size_t i = 0; i < 50; i++) {
        const auto V0 = lorentz_byPM(randomIsotropic<2>(RG) * P(RG), M(RG));
	const auto V0_copy=V0;
	EXPECT_EQ(V0,V0_copy);
        const auto V1 = V0.Transform(V0.Beta());
        EXPECT_TRUE(V1.P().CloseTo(zero(), epsilon));
        EXPECT_TRUE(abs(V1.M() - V0.M()) < epsilon);
    }
}
TEST(LorentzVector, LorentzTransform3d_more)
{
    RANDOM RG;
    RandomUniform<> M(0, 5), P(0, 5);
    for (size_t i = 0; i < 50; i++) {
        const auto V0 = lorentz_byPM(randomIsotropic<3>(RG) * P(RG), M(RG));
	const auto V0_copy=V0;
	EXPECT_EQ(V0,V0_copy);
        const auto V1 = V0.Transform(V0.Beta());
        EXPECT_TRUE(V1.P().CloseTo(Zero(), epsilon));
        EXPECT_TRUE(abs(V1.M() - V0.M()) < epsilon);
    }
}
TEST(LorentzVector, LorentzTransform4d_more)
{
    RANDOM RG;
    RandomUniform<> M(0, 5), P(0, 5);
    for (size_t i = 0; i < 50; i++) {
        const auto V0 = lorentz_byPM(randomIsotropic<4>(RG) * P(RG), M(RG));
	const auto V0_copy=V0;
	EXPECT_EQ(V0,V0_copy);
        const auto V1 = V0.Transform(V0.Beta());
        EXPECT_TRUE(V1.P().CloseTo(desCartes(0., 0., 0., 0.), epsilon));
        EXPECT_TRUE(abs(V1.M() - V0.M()) < epsilon);
    }
}

TEST(LorentzVector, decays)
{
    RANDOM RG;
    RandomUniform<> IM(2, 3), M1(0, 1), M2(0, 1), THETA(0, PI()), PHI(0, 2.0 * PI());
    for (size_t i = 0; i < 50; i++) {
        const double im = IM(RG), m1 = M1(RG), m2 = M2(RG);
        const auto C2 = binaryDecay(im, m1, m2, randomIsotropic<2>(RG));
        EXPECT_TRUE(abs((C2.first + C2.second).M() - im) < epsilon);
        EXPECT_TRUE(abs((C2.first).M() - m1) < epsilon);
        EXPECT_TRUE(abs((C2.second).M() - m2) < epsilon);
        EXPECT_TRUE((C2.first.P() + C2.second.P()).M() < epsilon);
        const auto C3 = binaryDecay(im, m1, m2, randomIsotropic<3>(RG));
        EXPECT_TRUE(abs((C3.first + C3.second).M() - im) < epsilon);
        EXPECT_TRUE(abs(C3.first.M() - m1) < epsilon);
        EXPECT_TRUE(abs(C3.second.M() - m2) < epsilon);
        EXPECT_TRUE((C3.first.P() + C3.second.P()).M() < epsilon);
    }
}
TEST(LorentzVector, decays2)
{
    RANDOM RG;
    RandomUniform<> IM(2, 3), M1(0, 1), M2(0, 1), THETA(0, PI()), PHI(0, 2.0 * PI());
    for (size_t i = 0; i < 50; i++) {
        const double im = IM(RG), m1 = M1(RG), m2 = M2(RG);
        const auto C2 = binaryDecay(im, m1, m2, randomIsotropic<2>(RG));
        EXPECT_TRUE(abs((C2.first + C2.second).M() - im) < epsilon);
        EXPECT_TRUE(abs((C2.first).M() - m1) < epsilon);
        EXPECT_TRUE(abs((C2.second).M() - m2) < epsilon);
        EXPECT_TRUE((C2.first.P() + C2.second.P()).M() < epsilon);
        const auto C3 = binaryDecay(im, m1, m2, randomIsotropic<3>(RG));
        EXPECT_TRUE(abs((C3.first + C3.second).M() - im) < epsilon);
        EXPECT_TRUE(abs(C3.first.M() - m1) < epsilon);
        EXPECT_TRUE(abs(C3.second.M() - m2) < epsilon);
        EXPECT_TRUE((C3.first.P() + C3.second.P()).M() < epsilon);
    }
}