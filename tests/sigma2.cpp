// this file is distributed under
// LGPLv3 license
#include <gtest/gtest.h>
#include <math_h/sigma2.h>
//this file contains unit tests for sigma2.h
using namespace std;
using namespace MathTemplates;
#define ALMOST_EQ(a,b) EXPECT_TRUE(abs(a-b)<0.0001)
#define ALMOST_EQ2(a,b) EXPECT_TRUE(abs(a-b)<0.01)
#define ALMOST_EQ3(a,b) EXPECT_TRUE(abs(a-b)<0.1)
class mock1 :public abstract_value_with_uncertainty_numeric<> {
private:
    double m_val;
    static size_t m_counter;
public:
    static inline void ResetCounter() { m_counter = 0; }
    static inline size_t Counter() { return m_counter; }
    mock1(double v) :abstract_value_with_uncertainty_numeric<>(), m_val(v) {}
    virtual ~mock1() {}
protected:
    virtual double mmeasure()const override {
        m_counter++;
        return m_val;
    }
};
size_t mock1::m_counter = 0;
TEST(ext_value, abstract_test)
{
    mock1 V(1);
    EXPECT_EQ(V.val(), 1);
    EXPECT_EQ(V.uncertainty(), 0);
    EXPECT_EQ(mock1::Counter(), 10000);
}
class mock2 :public abstract_value_with_uncertainty_numeric<> {
private:
    RandomGauss<> m_val;
    static size_t m_counter;
public:
    static inline void ResetCounter() { m_counter = 0; }
    static inline size_t Counter() { return m_counter; }
    mock2(double X, double s) :abstract_value_with_uncertainty_numeric<>(), m_val(X, s) {}
    virtual ~mock2() {}
protected:
    virtual double mmeasure()const override {
        m_counter++;
        return m_val();
    }
};
size_t mock2::m_counter = 0;
TEST(ext_value, abstract2)
{
    mock2 V(1, 0.1);
    ALMOST_EQ2(V.val(), 1);
    ALMOST_EQ2(V.uncertainty(), 0.1);
    EXPECT_EQ(mock2::Counter(), 10000);
}

TEST(ext_value, v_const)
{
    value_numeric_const<> V(1);
    EXPECT_EQ(V.val(), 1);
    EXPECT_EQ(V.uncertainty(), 0);
}
TEST(ext_value, v_gauss)
{
    value_numeric_distr<> V(1, 0.1);
    ALMOST_EQ2(V.val(), 1);
    ALMOST_EQ2(V.uncertainty(), 0.1);
    value_numeric_distr<> V2(RandomGauss<>(1, 0.1));
    ALMOST_EQ2(V2.val(), 1);
    ALMOST_EQ2(V2.uncertainty(), 0.1);
}
TEST(ext_value, v_func_chain)
{
    const value_numeric_distr<> V(12, 0.1), W(8, 0.12), X(2, 0.12);
    const Chain<value_numeric_distr<>> chain = { V,W,X };
    const function<double(const Chain<double>&)> f = [](const Chain<double>& x) {return x[0] * x[1] + x[2];};
    const auto F = FUNC(f, chain);
    const value<> v = V, w = W, x = X, r = (v * w + x);
    EXPECT_TRUE(r.Contains(F));
    ALMOST_EQ2(r.epsilon(), F.epsilon());
}
TEST(ext_value, v_func1)
{
    const value_numeric_distr<> V(12, 0.1);
    const value_f<1> F([](double x) {return x * x;}, V);
    const value<> v = V;
    EXPECT_TRUE((v * v).Contains(F));
    ALMOST_EQ2((v * v).epsilon(), F.epsilon());
}
TEST(ext_value, v_func1_sin)
{
    const value_numeric_distr<> V(12, 0.1);
    const auto F = SIN(V);
    EXPECT_TRUE(F.Contains(sin(V.val())));
}
TEST(ext_value, v_func1_cos)
{
    const value_numeric_distr<> V(12, 0.1);
    const auto F = COS(V);
    EXPECT_TRUE(F.Contains(cos(V.val())));
}
TEST(ext_value, v_func1_exp)
{
    const value_numeric_distr<> V(12, 0.1);
    const auto F = EXP(V);
    EXPECT_TRUE(F.Contains(exp(V.val())));
}
TEST(ext_value, v_func1_log)
{
    const value_numeric_distr<> V(12, 0.1);
    const auto F = LOG(V);
    EXPECT_TRUE(F.Contains(log(V.val())));
}
TEST(ext_value, v_func1_sqr)
{
    const value_numeric_distr<> V(12, 0.1);
    const auto F = SQR(V);
    EXPECT_TRUE(F.Contains(pow(V.val(), 2)));
}
TEST(ext_value, v_func1_SQRT)
{
    const value_numeric_distr<> V(12, 0.1);
    const auto F = SQRT(V);
    EXPECT_TRUE(F.Contains(sqrt(V.val())));
}
TEST(ext_value, v_func1_tan)
{
    const value_numeric_distr<> V(12, 0.1);
    const auto F = TAN(V);
    EXPECT_TRUE(F.Contains(tan(V.val())));
}
TEST(ext_value, v_func1_atan)
{
    const value_numeric_distr<> V(12, 0.1);
    const auto F = ATAN(V);
    EXPECT_TRUE(F.Contains(atan(V.val())));
}
TEST(ext_value, v_func1_exp_cos)
{
    const value_numeric_distr<> V(12, 0.1);
    const auto F = EXP(COS(V));
    EXPECT_TRUE(F.Contains(exp(cos(V.val()))));
}
TEST(ext_value, v_func2)
{
    const value_numeric_distr<> A(5, 0.1), B(3, 0.2);
    const value_f<2> F([](double x, double y) {return x * y;}, A, B);
    const value<> a = A, b = B;
    EXPECT_TRUE((a * b).Contains(F));
    ALMOST_EQ2((a * b).epsilon(), F.epsilon());
}
TEST(ext_value, v_func1_atan_2)
{
    const value_numeric_distr<> A(12, 0.1), B(-7, 0.2);
    const auto F = ATAN2(A, B);
    EXPECT_TRUE(F.Contains(atan2(A.val(), B.val())));
}
TEST(ext_value, v_func1_pow_2)
{
    const value_numeric_distr<> A(12, 0.1), B(3.0, 0.01);
    const auto F = POW(A, B);
    EXPECT_TRUE(F.Contains(pow(A.val(), B.val())));
}
TEST(ext_value, v_func2_opr1)
{
    const value_numeric_distr<> A(5, 0.1), B(3, 0.2);
    const auto F = A + B;
    const value<> a = A, b = B;
    EXPECT_TRUE((a + b).Contains(F));
    ALMOST_EQ2((a + b).epsilon(), F.epsilon());
}
TEST(ext_value, v_func2_opr1_ext)
{
    const value<> A(5, 0.1), B(3, 0.2);
    const auto F = A + B;
    const auto a = value_numeric(A);
    const auto b = value_numeric(B);
    EXPECT_TRUE((a + b).Contains(F));
    ALMOST_EQ2((a + b).epsilon(), F.epsilon());
}
TEST(ext_value, v_func2_opr2)
{
    const value_numeric_distr<> A(5, 0.1), B(3, 0.2);
    const auto F = A - B;
    const value<> a = A, b = B;
    EXPECT_TRUE((a - b).Contains(F));
    ALMOST_EQ2((a - b).epsilon(), F.epsilon());
}
TEST(ext_value, v_func2_opr3)
{
    const value_numeric_distr<> A(5, 0.1), B(3, 0.2);
    const auto F = A * B;
    const value<> a = A, b = B;
    EXPECT_TRUE((a * b).Contains(F));
    ALMOST_EQ2((a * b).epsilon(), F.epsilon());
}
TEST(ext_value, v_func2_opr4)
{
    const value_numeric_distr<> A(5, 0.1), B(3, 0.2);
    const auto F = A / B;
    const value<> a = A, b = B;
    EXPECT_TRUE((a / b).Contains(F));
    ALMOST_EQ2((a / b).epsilon(), F.epsilon());
}

TEST(ext_value, v_func1_e)
{
    const value_numeric_distr<> V(12, 0.1);
    const auto F = FUNC([](double x) {return x * x;}, V);
    const value<> v = V;
    EXPECT_TRUE((v * v).Contains(F));
    ALMOST_EQ2((v * v).epsilon(), F.epsilon());
}
TEST(ext_value, v_func1_e2)
{
    const value_numeric_distr<> V(12, 0.1);
    const auto F = LOG(V);
    const value<> v = V;
    EXPECT_TRUE(func_with_uncertainty([](double x) {return log(x);}, v).Contains(F));
    ALMOST_EQ2(func_with_uncertainty([](double x) {return log(x);}, v).epsilon(), F.epsilon());
}
TEST(ext_value, v_func2_e)
{
    const value_numeric_distr<> A(5, 0.1), B(3, 0.2);
    const auto F = FUNC([](double x, double y) {return x * y;}, A, B);
    const value<> a = A, b = B;
    EXPECT_TRUE((a * b).Contains(F));
    ALMOST_EQ2((a * b).epsilon(), F.epsilon());
}
TEST(ext_value, v_func2_e2)
{
    const value_numeric_distr<> A(5, 0.1), B(3, 0.2);
    const auto F = ATAN2(A, B);
    const value<> a = A, b = B;
    EXPECT_TRUE(func_with_uncertainty([](double x, double y) {return atan2(x, y);}, a, b).Contains(F));
    ALMOST_EQ2(func_with_uncertainty([](double x, double y) {return atan2(x, y);}, a, b).epsilon(), F.epsilon());
}
TEST(ext_value, v_func3_e)
{
    const value_numeric_distr<> A(5, 0.1), B(3, 0.2), C(0, 0.2);
    const auto F = FUNC([](double x, double y, double z) {return x * y + z;}, A, B, C);
    const value<> a = A, b = B, c = C;
    EXPECT_TRUE((a * b + c).Contains(F));
    ALMOST_EQ2((a * b + c).epsilon(), F.epsilon());
}
TEST(ext_value, v_func4_e)
{
    const value_numeric_distr<> A(5, 0.1), B(3, 0.2), C(0, 0.2), D(5, 0.02);
    const auto F = FUNC([](double x, double y, double z, double zz) {return (x * y + z) / zz;}, A, B, C, D);
    const value<> a = A, b = B, c = C, d = D;
    EXPECT_TRUE(((a * b + c) / d).Contains(F));
    ALMOST_EQ2(((a * b + c) / d).epsilon(), F.epsilon());
}
