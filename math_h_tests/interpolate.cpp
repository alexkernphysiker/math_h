#include <gtest/gtest.h>
#include <vector>
#include <math.h>
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
			if(V==TestArray[beg])EXPECT_TRUE((index==beg)||(index==(beg+1)));
		}
}
TEST(WhereToInsert,NormalConditions){
	for(double x=TestArray[0]-0.5;x<=TestArray[9]+0.5;x+=0.5)
		for(int beg=0;beg<10;beg++)for(int end=beg+1;end<10;end++){
			int index=WhereToInsert(beg,end,TestArray,x);
			if(x<TestArray[beg])EXPECT_EQ(beg,index);
			else if(x>TestArray[end])EXPECT_EQ(end+1,index);
			else EXPECT_TRUE((x>=TestArray[index-1])&&(x<=TestArray[index]));
		}
}
TEST(WhereToInsert,Search){
	for(double x=TestArray[0]-0.5;x<=TestArray[9]+0.5;x+=0.5)
		for(int beg=0;beg<10;beg++)for(int end=beg+1;end<10;end++){
			int index=Search(beg,end,TestArray,x);
			if(index>=beg)EXPECT_EQ(x,TestArray[index]);
		}
}
TEST(InsertSorted,BasicTest){
	vector<int> X;
	for(int i=0;i<50;i++){
		InsertSorted(rand()%10,X,std_size(X),std_insert(X,int));
		for(int j=0;j<i;j++)
			EXPECT_TRUE(X[j]<=X[j+1]);
	}
}
typedef std::pair<double,double> Pair;
TEST(LinearInterpolation,Create){
	LinearInterpolation<double> F;
	EXPECT_EQ(0,F.size());
	EXPECT_THROW(F.min(),math_h_error<LinearInterpolation<double>>);
	EXPECT_THROW(F.max(),math_h_error<LinearInterpolation<double>>);
	EXPECT_EQ(&F,&(F<<make_pair(0,0)));
	EXPECT_EQ(1,F.size());
	EXPECT_NO_THROW(F.min());
	EXPECT_NO_THROW(F.max());
	EXPECT_THROW(F(0.5),math_h_error<Pair>);
	EXPECT_EQ(&F,&(F<<make_pair(1,0)));
	EXPECT_EQ(2,F.size());
	EXPECT_NO_THROW(F.min());
	EXPECT_NO_THROW(F.max());
	EXPECT_EQ(0,F(0.5));
}
#define _EQ(a,b) EXPECT_TRUE(pow(a-b,2)<0.0001)
TEST(LinearInterpolation,SimpleLine){
	LinearInterpolation<double> F;
	F<<make_pair(0,0)<<make_pair(1,1);
	EXPECT_EQ(0,F(0));
	EXPECT_EQ(1,F(1));
	_EQ(0.1,F(0.1));
	_EQ(0.2,F(0.2));
	_EQ(0.3,F(0.3));
	_EQ(0.4,F(0.4));
	_EQ(0.5,F(0.5));
	_EQ(0.6,F(0.6));
	_EQ(0.7,F(0.7));
	_EQ(0.8,F(0.8));
	_EQ(0.9,F(0.9));
	F<<make_pair(2,0);
	EXPECT_EQ(0,F(0));
	EXPECT_EQ(1,F(1));
	_EQ(0.1,F(0.1));
	_EQ(0.2,F(0.2));
	_EQ(0.3,F(0.3));
	_EQ(0.4,F(0.4));
	_EQ(0.5,F(0.5));
	_EQ(0.6,F(0.6));
	_EQ(0.7,F(0.7));
	_EQ(0.8,F(0.8));
	_EQ(0.9,F(0.9));
	EXPECT_EQ(0,F(2));
	_EQ(0.1,F(1.9));
	_EQ(0.2,F(1.8));
	_EQ(0.3,F(1.7));
	_EQ(0.4,F(1.6));
	_EQ(0.5,F(1.5));
	_EQ(0.6,F(1.4));
	_EQ(0.7,F(1.3));
	_EQ(0.8,F(1.2));
	_EQ(0.9,F(1.1));
	F<<make_pair(0.5,0.5)<<make_pair(1.5,0.5);
	EXPECT_EQ(0,F(0));
	EXPECT_EQ(1,F(1));
	_EQ(0.1,F(0.1));
	_EQ(0.2,F(0.2));
	_EQ(0.3,F(0.3));
	_EQ(0.4,F(0.4));
	EXPECT_EQ(0.5,F(0.5));
	_EQ(0.6,F(0.6));
	_EQ(0.7,F(0.7));
	_EQ(0.8,F(0.8));
	_EQ(0.9,F(0.9));
	EXPECT_EQ(0,F(2));
	_EQ(0.1,F(1.9));
	_EQ(0.2,F(1.8));
	_EQ(0.3,F(1.7));
	_EQ(0.4,F(1.6));
	EXPECT_EQ(0.5,F(1.5));
	_EQ(0.6,F(1.4));
	_EQ(0.7,F(1.3));
	_EQ(0.8,F(1.2));
	_EQ(0.9,F(1.1));
	for(double x=0;x<2;x+=0.1)EXPECT_EQ(F(x),F.func()(x));
}
TEST(LinearInterpolation_fixedsize,Basic){
	LinearInterpolation_fixedsize<double> F(0,2,5);
	EXPECT_EQ(0,F(0));
	EXPECT_EQ(0,F(1));
	_EQ(0,F(0.1));
	_EQ(0,F(0.2));
	_EQ(0,F(0.3));
	_EQ(0,F(0.4));
	EXPECT_EQ(0,F(0.5));
	_EQ(0,F(0.6));
	_EQ(0,F(0.7));
	_EQ(0,F(0.8));
	_EQ(0,F(0.9));
	EXPECT_EQ(0,F(2));
	_EQ(0,F(1.9));
	_EQ(0,F(1.8));
	_EQ(0,F(1.7));
	_EQ(0,F(1.6));
	EXPECT_EQ(0,F(1.5));
	_EQ(0,F(1.4));
	_EQ(0,F(1.3));
	_EQ(0,F(1.2));
	_EQ(0,F(1.1));
	double values[]={0,0.5,1,0.5,0};
	for(int i=0,n=F.size();i<n;i++)
		F.setY(i,values[i]);
	EXPECT_EQ(0,F(0));
	EXPECT_EQ(1,F(1));
	_EQ(0.1,F(0.1));
	_EQ(0.2,F(0.2));
	_EQ(0.3,F(0.3));
	_EQ(0.4,F(0.4));
	EXPECT_EQ(0.5,F(0.5));
	_EQ(0.6,F(0.6));
	_EQ(0.7,F(0.7));
	_EQ(0.8,F(0.8));
	_EQ(0.9,F(0.9));
	EXPECT_EQ(0,F(2));
	_EQ(0.1,F(1.9));
	_EQ(0.2,F(1.8));
	_EQ(0.3,F(1.7));
	_EQ(0.4,F(1.6));
	EXPECT_EQ(0.5,F(1.5));
	_EQ(0.6,F(1.4));
	_EQ(0.7,F(1.3));
	_EQ(0.8,F(1.2));
	_EQ(0.9,F(1.1));
	for(double x=0;x<2;x+=0.1)EXPECT_EQ(F(x),F.func()(x));
}
TEST(LinearInterpolation_fixedsize,Throwing){
	EXPECT_THROW(LinearInterpolation_fixedsize<double>(1,0,2),math_h_error<LinearInterpolation_fixedsize<double>>);
	EXPECT_THROW(LinearInterpolation_fixedsize<double>(0,0,2),math_h_error<LinearInterpolation_fixedsize<double>>);
	EXPECT_THROW(LinearInterpolation_fixedsize<double>(0,1,0),math_h_error<LinearInterpolation_fixedsize<double>>);
	EXPECT_THROW(LinearInterpolation_fixedsize<double>(0,1,1),math_h_error<LinearInterpolation_fixedsize<double>>);
	EXPECT_NO_THROW(LinearInterpolation_fixedsize<double>(0,1,2));
	LinearInterpolation_fixedsize<double> F(0,1,2);
	ASSERT_EQ(2,F.size());
	EXPECT_THROW(F.getX(-1),math_h_error<LinearInterpolation_fixedsize<double>>);
	EXPECT_THROW(F.getY(-1),math_h_error<LinearInterpolation_fixedsize<double>>);
	EXPECT_THROW(F.getX(2),math_h_error<LinearInterpolation_fixedsize<double>>);
	EXPECT_THROW(F.getY(2),math_h_error<LinearInterpolation_fixedsize<double>>);
	EXPECT_EQ(0,F.getX(0));
	EXPECT_EQ(0,F.getY(0));
	EXPECT_EQ(1,F.getX(1));
	EXPECT_EQ(0,F.getY(1));
	EXPECT_EQ(0,F(0));
	EXPECT_EQ(0,F(1));
}
TEST(Distribution,BasicTest){
	Distribution<double> D(0,2,2);
	_EQ(1,D.BinWidth());
	ASSERT_EQ(2,D.size());
	EXPECT_EQ(0.5,D.min());
	EXPECT_EQ(1.5,D.max());
	EXPECT_EQ(0.5,D.getX(0));
	EXPECT_EQ(0,D.getY(0));
	EXPECT_EQ(1.5,D.getX(1));
	EXPECT_EQ(0,D.getY(1));
	EXPECT_EQ(&D,&(D.AddValue(0.1)));
	ASSERT_EQ(2,D.size());
	EXPECT_EQ(0.5,D.getX(0));
	EXPECT_EQ(1,D.getY(0));
	EXPECT_EQ(1.5,D.getX(1));
	EXPECT_EQ(0,D.getY(1));
	EXPECT_EQ(&D,&(D.AddValue(1.5)));
	ASSERT_EQ(2,D.size());
	EXPECT_EQ(0.5,D.getX(0));
	EXPECT_EQ(1,D.getY(0));
	EXPECT_EQ(1.5,D.getX(1));
	EXPECT_EQ(1,D.getY(1));
	EXPECT_EQ(1,D(1));
}
TEST(Distribution,Throwing){
	EXPECT_THROW(Distribution<double>(1,0,2),math_h_error<LinearInterpolation_fixedsize<double>>);
	EXPECT_THROW(Distribution<double>(0,0,2),math_h_error<LinearInterpolation_fixedsize<double>>);
	EXPECT_THROW(Distribution<double>(0,1,0),math_h_error<LinearInterpolation_fixedsize<double>>);
	EXPECT_THROW(Distribution<double>(0,1,1),math_h_error<LinearInterpolation_fixedsize<double>>);
	EXPECT_NO_THROW(Distribution<double>(0,1,2));
}
