#include <gtest/gtest.h>
#include <vector>
#include <math.h>
#include <randomfunc.h>
using namespace std;
TEST(RandomUniformlyI,BasicTest){
	for(int l=-10;l<=10;l++)
		for(int r=l;r<=10;r++)
			for(int c=0;c<100;c++){
				int V=RandomUniformlyI(l,r);
				EXPECT_TRUE((V>=l)&&(V<=r));
			}
}
TEST(RandomUniformlyI,Throwing){
	for(int l=-10;l<=10;l++)
		for(int r=-10;r<l;r++)
			ASSERT_ANY_THROW(RandomUniformlyI(l,r));
}
TEST(RandomUniformlyR,BasicTest){
	for(double l=-10;l<=10;l+=0.5)
		for(double r=l;r<=10;r+=0.5)
			for(int c=0;c<100;c++){
				double V=RandomUniformlyR(l,r);
				EXPECT_TRUE((V>=l)&&(V<=r));
			}
}
TEST(RandomUniformlyR,Throwing){
	for(double l=-10;l<=10;l+=0.5)
		for(double r=-10;r<l;r+=0.5)
			ASSERT_ANY_THROW(RandomUniformlyR(l,r));
}
TEST(RandomGauss,BasicTest){
	for(double X=-10;X<=10;X+=0.5){
		for(int c=0;c<100;c++){
			double V=RandomGauss(0.0,X);
			EXPECT_EQ(X,V);
		}
		for(double sigma=0.1;sigma<3;sigma+=0.1){
			int cg=0,cl=0;
			for(int i=0;i<100;i++){
				double V=RandomGauss(sigma,X);
				if(pow(V-X,2)<pow(sigma,2))
					cl++;
				else
					cg++;
			}
			EXPECT_TRUE(cg<cl);
			EXPECT_TRUE(cg*6>cl);
		}
	}
}
TEST(RandomGauss,Throwing){
	for(double X=-1;X<=1;X+=0.5){
		auto f=[&X](){return RandomGauss(-1.0,X);};
		ASSERT_ANY_THROW(f());
	}
}

