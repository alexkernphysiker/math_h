// this file is distributed under
// LGPLv3 license
#include <gtest/gtest.h>
#include <math_h/numeric_diff.h>
#include <math_h/vectortransformations.h>
//this file contains unit tests for sigma.h
using namespace std;
using namespace MathTemplates;

TEST(numeric_diff, der1) {
    EXPECT_EQ(0, (num_der1(0.1) * [](const double&)->double {return 5.0;})(3.0));
    EXPECT_EQ(2, (num_der1(0.1) * [](const double& x)->double {return 2.0 * x;})(3.0));
    EXPECT_EQ(0, (num_der1(0.1) * [](const double& x)->double {return 2.0 * x * x;})(0.0));
}
TEST(numeric_diff, der2) {
    EXPECT_EQ(0, (num_der2(0.1) * [](const double&)->double {return 5.0;})(3.0));
    EXPECT_EQ(0, (num_der2(0.1) * [](const double& x)->double {return 2.0 * x;})(3.0));
    EXPECT_EQ(2.0, (num_der2(0.1) * [](const double& x)->double {return x * x;})(0.0));
}
TEST(numeric_diff, Pder1) {
    const auto F = [](const Vector<>& r) {return r.length_sqr();};
    EXPECT_EQ(0, (num_Pder1<1, 3>(0.1) * F)(Zero()));
    EXPECT_EQ(0, (num_Pder1<2, 3>(0.1) * F)(Zero()));
    EXPECT_EQ(0, (num_Pder1<3, 3>(0.1) * F)(Zero()));
    const auto F2 = [](const Vector<>& r) {return r.x() + r.y() - r.z();};
    EXPECT_EQ(1, (num_Pder1<1, 3>(0.1) * F2)(Zero()));
    EXPECT_EQ(1, (num_Pder1<2, 3>(0.1) * F2)(Zero()));
    EXPECT_EQ(-1, (num_Pder1<3, 3>(0.1) * F2)(Zero()));
}
TEST(numeric_diff, nabla) {
    const auto F1 = [](const Vector<>& r) {return r.length_sqr();};
    const auto G = nabla<3>(0.1) * F1;
    EXPECT_EQ(0, G.x()(Zero()));
    EXPECT_EQ(0, G.y()(Zero()));
    EXPECT_EQ(0, G.z()(Zero()));
    const auto F2 = [](const Vector<>& r) {return r.x() + 2.0 * r.y() - r.z();};
    const auto G2 = nabla<3>(0.1) * F2;
    EXPECT_EQ(1, G2.x()(Zero()));
    EXPECT_EQ(2, G2.y()(Zero()));
    EXPECT_EQ(-1, G2.z()(Zero()));
    const auto F3 = [](const Vector<>&) {return 0.0;};
    const auto G3 = nabla<3>(0.1) * F3;
    EXPECT_EQ(0, G3.x()(Zero()));
    EXPECT_EQ(0, G3.y()(Zero()));
    EXPECT_EQ(0, G3.z()(Zero()));
    const FunctionWrap<double, const Vector<>&> f1 = F1, f2 = F2, f3 = F3;
    const auto div = nabla<3>(0.1) * vec(f1, f2, f3);
    EXPECT_EQ(2, div(Zero()));
}
TEST(numeric_diff, nabla2) {
    const FunctionWrap<double, const Vector<>&>
        f1 = [](const Vector<>& r) {return exp(-r.x() * r.x());},
        f2 = [](const Vector<>& r) {return exp(-r.y() * r.y());},
        f3 = [](const Vector<>& r) {return exp(-r.z() * r.z());};
    const auto x = vec(1., 1., 0.);
    const auto d = nabla<3>(0.1) * vec(f1, f2, f3);
    const auto d2 = num_Pder1<1, 3>(0.1) * f1 + num_Pder1<2, 3>(0.1) * f2 + num_Pder1<3, 3>(0.1) * f3;
    EXPECT_EQ(d(x), d2(x));
}
TEST(numeric_diff, nabla3) {
    const FunctionWrap<double, const Vector<>&>
        f1 = [](const Vector<>& r) {return exp(-r.x() * r.x());},
        f2 = [](const Vector<>& r) {return exp(-r.y() * r.y());},
        f3 = [](const Vector<>& r) {return exp(-r.z() * r.z());};
    EXPECT_EQ(0, (nabla<3>(0.1) ^ vec(f1, f2, f3)).x(vec(1., 1., 0.)));
    EXPECT_EQ(0, (nabla<3>(0.1) ^ vec(f1, f2, f3)).y(vec(1., 1., 0.)));
    EXPECT_EQ(0, (nabla<3>(0.1) ^ vec(f1, f2, f3)).z(vec(1., 1., 0.)));
    EXPECT_EQ(0, (nabla<3>(0.1) ^ vec(f1, f2, f3)).x(vec(1., 0., 1.)));
    EXPECT_EQ(0, (nabla<3>(0.1) ^ vec(f1, f2, f3)).y(vec(1., 0., 1.)));
    EXPECT_EQ(0, (nabla<3>(0.1) ^ vec(f1, f2, f3)).z(vec(1., 0., 1.)));
    EXPECT_EQ(0, (nabla<3>(0.1) ^ vec(f1, f2, f3)).x(vec(0., 1., 1.)));
    EXPECT_EQ(0, (nabla<3>(0.1) ^ vec(f1, f2, f3)).y(vec(0., 1., 1.)));
    EXPECT_EQ(0, (nabla<3>(0.1) ^ vec(f1, f2, f3)).z(vec(0., 1., 1.)));
}
TEST(numeric_diff, nabla4) {
    const FunctionWrap<double, const Vector<>&>
        f1 = [](const Vector<>& r) {return exp(-r.x() * r.x());},
        f2 = [](const Vector<>& r) {return exp(-r.y() * r.y());},
        f3 = [](const Vector<>& r) {return exp(-r.z() * r.z());};
    EXPECT_EQ(0, (nabla<3>(0.1) ^ vec(f1, f2, f3))(vec(1., 1., 0.)).length());
    EXPECT_EQ(0, (nabla<3>(0.1) ^ vec(f1, f2, f3))(vec(1., 0., 1.)).length());
    EXPECT_EQ(0, (nabla<3>(0.1) ^ vec(f1, f2, f3))(vec(0., 1., 1.)).length());
}
