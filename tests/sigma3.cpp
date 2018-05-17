// this file is distributed under
// LGPLv3 license
#include <gtest/gtest.h>
#include <math_h/sigma3.h>
#include <math_h/randomfunc.h>
using namespace std;
using namespace MathTemplates;
#define ALMOST_EQ(a,b) EXPECT_TRUE(abs(a-b)<0.0001)
#define ALMOST_EQ2(a,b) EXPECT_TRUE(abs(a-b)<0.01)
#define ALMOST_EQ3(a,b) EXPECT_TRUE(abs(a-b)<0.1)

TEST(Uncertainties,test1)
{
    RandomUniform<> A(0.1,10);
    RandomUniform<> a(0.1,1);
    for (size_t i = 0; i < 50; i++) {
	const double v=A(),u1=a(),u2=a();
	const double w=A(),j1=a(),j2=a();
	const value<> V1(v,u1),V2(v,u2);
	const value<> W1(w,j1),W2(w,j2);
	const auto M=uncertainties(v,u1,u2);
	const auto N=uncertainties(w,j1,j2);

	EXPECT_EQ((M+N).val(),(V1+W1).val());
	EXPECT_EQ((M+N).val(),(V2+W2).val());
	EXPECT_EQ((M+N).template uncertainty<1>(),(V1+W1).uncertainty());
	EXPECT_EQ((M+N).template uncertainty<2>(),(V2+W2).uncertainty());

	EXPECT_EQ((M-N).val(),(V1-W1).val());
	EXPECT_EQ((M-N).val(),(V2-W2).val());
	EXPECT_EQ((M-N).template uncertainty<1>(),(V1-W1).uncertainty());
	EXPECT_EQ((M-N).template uncertainty<2>(),(V2-W2).uncertainty());

	EXPECT_EQ((M*N).val(),(V1*W1).val());
	EXPECT_EQ((M*N).val(),(V2*W2).val());
	EXPECT_EQ((M*N).template uncertainty<1>(),(V1*W1).uncertainty());
	EXPECT_EQ((M*N).template uncertainty<2>(),(V2*W2).uncertainty());

	EXPECT_EQ((M/N).val(),(V1/W1).val());
	EXPECT_EQ((M/N).val(),(V2/W2).val());
	EXPECT_EQ((M/N).template uncertainty<1>(),(V1/W1).uncertainty());
	EXPECT_EQ((M/N).template uncertainty<2>(),(V2/W2).uncertainty());
    }
}
TEST(Uncertainties,test2)
{
    RandomUniform<> A(0.1,10);
    RandomUniform<> a(0.1,1);
    for (size_t i = 0; i < 50; i++) {
	const value<> V(A(),a());
	const auto M11=extend_value<1>(V);
	EXPECT_EQ(M11.val(),V.val());
	EXPECT_EQ(M11.template uncertainty<1>(),V.uncertainty());
	const auto M12=extend_value<1,2>(V);
	EXPECT_EQ(M12.val(),V.val());
	EXPECT_EQ(M12.template uncertainty<1>(),V.uncertainty());
	EXPECT_EQ(M12.template uncertainty<2>(),0);
	const auto M13=extend_value<1,3>(V);
	EXPECT_EQ(M13.val(),V.val());
	EXPECT_EQ(M13.template uncertainty<1>(),V.uncertainty());
	EXPECT_EQ(M13.template uncertainty<2>(),0);
	EXPECT_EQ(M13.template uncertainty<3>(),0);
	const auto M22=extend_value<2>(V);
	EXPECT_EQ(M22.val(),V.val());
	EXPECT_EQ(M22.template uncertainty<1>(),0);
	EXPECT_EQ(M22.template uncertainty<2>(),V.uncertainty());
	const auto M23=extend_value<2,3>(V);
	EXPECT_EQ(M23.val(),V.val());
	EXPECT_EQ(M23.template uncertainty<1>(),0);
	EXPECT_EQ(M23.template uncertainty<2>(),V.uncertainty());
	EXPECT_EQ(M23.template uncertainty<3>(),0);
	const auto M33=extend_value<3>(V);
	EXPECT_EQ(M33.val(),V.val());
	EXPECT_EQ(M33.template uncertainty<1>(),0);
	EXPECT_EQ(M33.template uncertainty<2>(),0);
	EXPECT_EQ(M33.template uncertainty<3>(),V.uncertainty());
    }
}
TEST(Uncertainties,test3)
{
    RandomUniform<> A(0.1,10);
    RandomUniform<> a(0.1,1);
    for (size_t i = 0; i < 50; i++) {
	const value<> V1(A(),a()),V2(A(),a());
	{
	    const auto s=V1+V2,S=(extend_value<1,2>(V1)+extend_value<2,2>(V2)).wrap();
	    EXPECT_EQ(s.val(),S.val());
	    ALMOST_EQ2(s.uncertainty(),S.uncertainty());
	}
	{
	    const auto s=V1-V2,S=(extend_value<1,2>(V1)-extend_value<2,2>(V2)).wrap();
	    EXPECT_EQ(s.val(),S.val());
	    ALMOST_EQ2(s.uncertainty(),S.uncertainty());
	}
	{
	    const auto s=V1*V2,S=(extend_value<1,2>(V1)*extend_value<2,2>(V2)).wrap();
	    EXPECT_EQ(s.val(),S.val());
	    ALMOST_EQ2(s.uncertainty(),S.uncertainty());
	}
	{
	    const auto s=V1/V2,S=(extend_value<1,2>(V1)/extend_value<2,2>(V2)).wrap();
	    EXPECT_EQ(s.val(),S.val());
	    ALMOST_EQ2(s.uncertainty(),S.uncertainty());
	}
    }
}
