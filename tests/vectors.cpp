// this file is distributed under
// LGPLv3 license
#include <gtest/gtest.h>
#include <math_h/vectors.h>
#include <math_h/hists.h>
#include <math_h/sigma.h>
using namespace std;
using namespace MathTemplates;
TEST(Vector, scalar_prod2)
{
    EXPECT_EQ(0, x()*y());
    EXPECT_EQ(0, y()*x());
    EXPECT_EQ(1, x()*x());
    EXPECT_EQ(1, y()*y());

    EXPECT_EQ(0, x() * (-y()));
    EXPECT_EQ(0, y() * (-x()));
    EXPECT_EQ(0, (-x())*y());
    EXPECT_EQ(0, (-y())*x());

    EXPECT_EQ(-1, x() * (-x()));
    EXPECT_EQ(-1, y() * (-y()));
    EXPECT_EQ(-1, (-x())*x());
    EXPECT_EQ(-1, (-y())*y());
}
TEST(Vector, scalar_prod3)
{
    EXPECT_EQ(0, X()*Y());
    EXPECT_EQ(0, Y()*Z());
    EXPECT_EQ(0, X()*Z());
    EXPECT_EQ(1, X()*X());
    EXPECT_EQ(1, Y()*Y());
    EXPECT_EQ(1, Z()*Z());

    EXPECT_EQ(0, X() * (-Y()));
    EXPECT_EQ(0, (-X())*Z());

    EXPECT_EQ(-1, (-X())*X());
    EXPECT_EQ(-1, Y() * (-Y()));
    EXPECT_EQ(1, (-Z()) * (-Z()));
}
TEST(Vector, base1d)
{
    RANDOM RG;
    RandomUniform<> X(-50., 50);
    for (size_t i = 0; i < 50; i++) {
        const auto x = X(RG);
        const auto V = desCartes(x);
        EXPECT_EQ(x, V.x());
        EXPECT_EQ(x, V.component<1>());
    }
}
TEST(Vector, base2d)
{
    RANDOM RG;
    RandomUniform<> X(-50., 50);
    for (size_t i = 0; i < 50; i++) {
        const auto x = X(RG), y = X(RG);
        const auto V = desCartes(x, y);
        EXPECT_EQ(x, V.x());
        EXPECT_EQ(x, V.component<1>());
        EXPECT_EQ(y, V.y());
        EXPECT_EQ(y, V.component<2>());
    }
}
TEST(Vector, base3d)
{
    RANDOM RG;
    RandomUniform<> X(-50., 50);
    for (size_t i = 0; i < 50; i++) {
        const auto x = X(RG), y = X(RG), z = X(RG);
        const auto V = desCartes(x, y, z);
        EXPECT_EQ(x, V.x());
        EXPECT_EQ(x, V.component<1>());
        EXPECT_EQ(y, V.y());
        EXPECT_EQ(y, V.component<2>());
        EXPECT_EQ(z, V.z());
        EXPECT_EQ(z, V.component<3>());
    }
}
TEST(Vector, base4d)
{
    RANDOM RG;
    RandomUniform<> X(-50., 50);
    for (size_t i = 0; i < 50; i++) {
        const auto x = X(RG), y = X(RG), z = X(RG), zz = X(RG);
        const auto V = desCartes(x, y, z, zz);
        EXPECT_EQ(x, V.x());
        EXPECT_EQ(x, V.component<1>());
        EXPECT_EQ(y, V.y());
        EXPECT_EQ(y, V.component<2>());
        EXPECT_EQ(z, V.z());
        EXPECT_EQ(z, V.component<3>());
        EXPECT_EQ(zz, V.component<4>());
    }
}

