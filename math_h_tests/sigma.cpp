#include <gtest/gtest.h>
#include <sigma.h>
#include <randomfunc.h>
using namespace std;
TEST(Sigma,Throwing){
	Sigma<double> S;
	EXPECT_EQ(0,S.count());
	EXPECT_THROW(S.getAverage(),exception);
	EXPECT_THROW(S.getSigmaSqr(),exception);
	EXPECT_THROW(S.getSigma(),exception);
}
#define _EQ(a,b) EXPECT_TRUE(pow(a-b,2)<0.01)
TEST(Sigma,Base){
	Sigma<double> S;
	S.AddValue(0);
	EXPECT_EQ(1,S.count());
	EXPECT_EQ(0,S.getAverage());
	EXPECT_THROW(S.getSigmaSqr(),exception);
	EXPECT_THROW(S.getSigma(),exception);
	S.AddValue(1);
	EXPECT_EQ(2,S.count());
	EXPECT_EQ(0.5,S.getAverage());
	EXPECT_EQ(0.5,S.getSigmaSqr());
	EXPECT_EQ(S.getSigma(),sqrt(S.getSigmaSqr()));
}
TEST(Sigma,WithRandomValues){
	Sigma<double> S;
	for(int i=0;i<1000;i++)
		S.AddValue(RandomGauss(1.0,1.0));
	_EQ(1.0,S.getAverage());
	_EQ(1.0,S.getSigma());
	_EQ(1.0,S.getSigmaSqr());
}

