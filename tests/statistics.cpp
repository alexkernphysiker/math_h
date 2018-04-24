// this file is distributed under
// LGPLv3 license
#include <gtest/gtest.h>
#include <math_h/randomfunc.h>
#include <math_h/statistics.h>
using namespace std;
using namespace MathTemplates;
#define ALMOST_EQ(a,b) EXPECT_TRUE(abs(a-b)<0.0001)
#define ALMOST_EQ2(a,b) EXPECT_TRUE(abs(a-b)<0.01)
#define ALMOST_EQ3(a,b) EXPECT_TRUE(abs(a-b)<0.1)
TEST(AverageObtainer,base)
{
    AverageObtainer<double> S;
    EXPECT_EQ(S.count(),0);
    EXPECT_EQ(&(S<<2<<3<<2<<3),&S);
    EXPECT_EQ(S.Average(),2.5);
    EXPECT_EQ(S.count(),4);
}
class TestStatistics:public AverageObtainer<double>{
public:
    TestStatistics():AverageObtainer<double>(){}
    virtual ~TestStatistics(){}
    void Test()const{
	size_t c=0;
	double res=0;
	ForEach([&c,&res](double x){
	    c++;
	    res+=x;
	});
	EXPECT_EQ(c,count());
	EXPECT_EQ(res,Average()*c);
    }
};
TEST(TestStatistics,base){
    TestStatistics S;
    S<<0;S.Test();
    S<<1;S.Test();
    S<<2;S.Test();
    S<<2;S.Test();
    S<<5;S.Test();
    S<<6;S.Test();
    S<<3;S.Test();
    S<<1;S.Test();
}

