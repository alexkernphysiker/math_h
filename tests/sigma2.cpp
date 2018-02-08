
// this file is distributed under
// LGPLv3 license
#include <gtest/gtest.h>
#include <math_h/sigma2.h>
using namespace std;
using namespace MathTemplates;
#define ALMOST_EQ(a,b) EXPECT_TRUE(abs(a-b)<0.0001)
#define ALMOST_EQ2(a,b) EXPECT_TRUE(abs(a-b)<0.01)
#define ALMOST_EQ3(a,b) EXPECT_TRUE(abs(a-b)<0.1)
class mock1:public abstract_value_with_extended_uncertainty<>{
private:
    double m_val;
    static size_t m_counter;
public:
    static inline void ResetCounter(){m_counter=0;}
    static inline size_t Counter(){return m_counter;}
    mock1(double v):abstract_value_with_extended_uncertainty<>(),m_val(v){}
    virtual ~mock1(){}
protected:
    virtual double mmeasure()const override{
	m_counter++;
	return m_val;
    }
};
size_t mock1::m_counter=0;
TEST(ext_value, base1)
{
    mock1 V(1);
    EXPECT_EQ(V.val(),1);
    EXPECT_EQ(V.uncertainty(),0);
    EXPECT_EQ(mock1::Counter(),1000);
}
class mock2:public abstract_value_with_extended_uncertainty<>{
private:
    RandomGauss<> m_val;
    static size_t m_counter;
public:
    static inline void ResetCounter(){m_counter=0;}
    static inline size_t Counter(){return m_counter;}
    mock2(double X,double s):abstract_value_with_extended_uncertainty<>(),m_val(X,s){}
    virtual ~mock2(){}
protected:
    virtual double mmeasure()const override{
	m_counter++;
	return m_val();
    }
};
size_t mock2::m_counter=0;
TEST(ext_value, base2)
{
    mock2 V(1,0.1);
    ALMOST_EQ2(V.val(),1);
    ALMOST_EQ2(V.uncertainty(),0.1);
    EXPECT_EQ(mock2::Counter(),1000);
}

TEST(ext_value, v_const)
{
    value_const<> V(1);
    EXPECT_EQ(V.val(),1);
    EXPECT_EQ(V.uncertainty(),0);
}
TEST(ext_value, v_gauss)
{
    value_ext<> V(1,0.1);
    ALMOST_EQ2(V.val(),1);
    ALMOST_EQ2(V.uncertainty(),0.1);
    value_ext<> V2(RandomGauss<>(1,0.1));
    ALMOST_EQ2(V2.val(),1);
    ALMOST_EQ2(V2.uncertainty(),0.1);
}
TEST(ext_value, v_func1)
{
    const value_ext<> V(12,0.1);
    const value_f<1> F([](double x){return x*x;},V);
    const value<> v=V;
    EXPECT_TRUE((v*v).Contains(F));
    cout<<(v*v)<<" ; "<<F<<endl;
    ALMOST_EQ2((v*v).epsilon(),F.epsilon());
}
TEST(ext_value, v_func2)
{
    const value_ext<> A(5,0.1),B(3,0.2);
    const value_f<2> F([](double x,double y){return x*y;},A,B);
    const value<> a=A,b=B;
    EXPECT_TRUE((a*b).Contains(F));
    cout<<(a*b)<<" ; "<<F<<endl;
    ALMOST_EQ2((a*b).epsilon(),F.epsilon());
}
TEST(ext_value, v_func2_opr1)
{
    const value_ext<> A(5,0.1),B(3,0.2);
    const auto F=A+B;
    const value<> a=A,b=B;
    EXPECT_TRUE((a+b).Contains(F));
    ALMOST_EQ2((a+b).epsilon(),F.epsilon());
}
TEST(ext_value, v_func2_opr1_ext)
{
    const value<> A(5,0.1),B(3,0.2);
    const auto F=A+B;
    const auto a=ext_value(A);
    const auto b=ext_value(B);
    EXPECT_TRUE((a+b).Contains(F));
    ALMOST_EQ2((a+b).epsilon(),F.epsilon());
}
TEST(ext_value, v_func2_opr2)
{
    const value_ext<> A(5,0.1),B(3,0.2);
    const auto F=A-B;
    const value<> a=A,b=B;
    EXPECT_TRUE((a-b).Contains(F));
    ALMOST_EQ2((a-b).epsilon(),F.epsilon());
}
TEST(ext_value, v_func2_opr3)
{
    const value_ext<> A(5,0.1),B(3,0.2);
    const auto F=A*B;
    const value<> a=A,b=B;
    EXPECT_TRUE((a*b).Contains(F));
    ALMOST_EQ2((a*b).epsilon(),F.epsilon());
}
TEST(ext_value, v_func2_opr4)
{
    const value_ext<> A(5,0.1),B(3,0.2);
    const auto F=A/B;
    const value<> a=A,b=B;
    EXPECT_TRUE((a/b).Contains(F));
    ALMOST_EQ2((a/b).epsilon(),F.epsilon());
}
