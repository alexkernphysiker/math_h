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
