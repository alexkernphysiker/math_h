#include <gtest/gtest.h>
#include <vector>
#include <interpolate.h>
using namespace std;
double TestArray[]={-2,-1,0,1,2,3,4,5,6,7};
TEST(WhereToInsert,BorderConditions){
	for(int beg=0;beg<10;beg++)
		for(double V=TestArray[beg]-2;V<=TestArray[beg]+2;V+=0.5){
			int index=WhereToInsert(beg,beg-1,TestArray,V);
			EXPECT_EQ(beg,index);
			index=WhereToInsert(beg,beg,TestArray,V);
			if(V<TestArray[beg])EXPECT_EQ(beg,index);
			if(V>TestArray[beg])EXPECT_EQ(beg+1,index);
			if(V==TestArray[beg])ASSERT_TRUE((index==beg)||(index==(beg+1)));
		}
}
TEST(WhereToInsert,NormalConditions){
	for(double x=TestArray[0]-0.5;x<=TestArray[9]+0.5;x+=0.5)
		for(int beg=0;beg<10;beg++)for(int end=beg+1;end<10;end++){
			int index=WhereToInsert(beg,end,TestArray,x);
			if(x<TestArray[beg])EXPECT_EQ(beg,index);
			else if(x>TestArray[end])EXPECT_EQ(end+1,index);
			else ASSERT_TRUE((x>=TestArray[index-1])&&(x<=TestArray[index]));
		}
}
TEST(InsertSorted,BasicTest){
	vector<int> X;
	for(int i=0;i<50;i++){
		InsertSorted(rand()%10,X,std_size(X),std_insert(X,int));
		for(int j=0;j<i;j++)
			ASSERT_TRUE(X[j]<=X[j+1]);
	}
}

