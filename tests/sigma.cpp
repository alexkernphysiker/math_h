// this file is distributed under 
// MIT license
#include <random>
#include <gtest/gtest.h>
#include <math_h/sigma.h>
using namespace std;
using namespace MathTemplates;
TEST(value,throwing){
	EXPECT_THROW(value<double>(1,-0.1),Exception<value<double>>);
	EXPECT_NO_THROW(value<double>(1,0));
	EXPECT_NO_THROW(value<double>(1,0.1));
}
TEST(value,base){
	mt19937 gen;
	normal_distribution<double> G;
	for(size_t i=0;i<10;i++){
		double x=G(gen);
		value<double> V(x,0.1);
		EXPECT_EQ(x,V.val());
		EXPECT_EQ(0.1,V.delta());
		EXPECT_EQ(0.1/x,V.epsilon());
		EXPECT_EQ(x-0.1,V.min());
		EXPECT_EQ(x+0.1,V.max());
	}
}
TEST(value,contains){
	value<double> V(0,0.1);
	EXPECT_EQ(false,V.contains(-0.2));
	EXPECT_EQ(true,V.contains(-0.05));
	EXPECT_EQ(true,V.contains(0));
	EXPECT_EQ(true,V.contains(0.05));
	EXPECT_EQ(false,V.contains(0.2));
}
TEST(Sigma,Throwing){
	Sigma<double> S;
	EXPECT_EQ(0,S.count());
	EXPECT_THROW(S.get(),Exception<Sigma<double>>);
	EXPECT_EQ(&S,&(S<<0.0));
	EXPECT_EQ(1,S.count());
	EXPECT_THROW(S.get(),Exception<Sigma<double>>);
	EXPECT_EQ(&S,&(S<<0.0));
	EXPECT_EQ(2,S.count());
	EXPECT_EQ(0,S.get().val());
	EXPECT_EQ(0,S.get().delta());
}
#define _EQ(a,b) EXPECT_TRUE(pow(a-b,2)<0.005)
TEST(Sigma,Base){
	Sigma<double> S;
	S<<0.0;
	EXPECT_EQ(1,S.count());
	EXPECT_THROW(S.get(),Exception<Sigma<double>>);
	S<<1.0;
	EXPECT_EQ(2,S.count());
	EXPECT_EQ(0.5,S.get().val());
	EXPECT_EQ(sqrt(0.5),S.get().delta());
}
#define _EQ2(a,b) EXPECT_TRUE(pow(a-b,2)<0.01)
TEST(Sigma,WithRandomValues){
	Sigma<double> S;
	default_random_engine generator;
	normal_distribution<double> distribution(1.0,3.0);
	for(int i=0;i<2000;i++)
		S<<distribution(generator);
	_EQ2(1.0,S.get().val());
	_EQ2(3.0,S.get().delta());
}
TEST(Sigma,WithRandomValues2){
	Sigma<double> S(2);
	default_random_engine generator;
	normal_distribution<double> distribution(1.0,3.0);
	for(int i=0;i<2000;i++)
		S<<distribution(generator);
	_EQ2(1.0,S.get().val());
	_EQ2(6.0,S.get().delta());
}

TEST(WeightedAverageCalculator,Zeros){
	WeightedAverageCalculator<double> W;
	EXPECT_THROW(W.get(),Exception<WeightedAverageCalculator<double>>);
	EXPECT_THROW(W<<value<double>(0,0),Exception<WeightedAverageCalculator<double>>);
	EXPECT_EQ(&W,&(W<<value<double>(0,1)));
	_EQ(0,W.get().val());
	_EQ(1,W.get().delta());
	EXPECT_EQ(&W,&(W<<value<double>(0,1)));
	_EQ(0,W.get().val());
	_EQ(1.0/sqrt(2.0),W.get().delta());
	EXPECT_EQ(&W,&(W<<value<double>(0,1)));
	_EQ(0,W.get().val());
	_EQ(1.0/sqrt(3.0),W.get().delta());
}
TEST(WeightedAverageCalculator,Ones){
	WeightedAverageCalculator<double> W;
	EXPECT_THROW(W.get(),Exception<WeightedAverageCalculator<double>>);
	EXPECT_THROW(W<<value<double>(1,0),Exception<WeightedAverageCalculator<double>>);
	EXPECT_EQ(&W,&(W<<value<double>(1,1)));
	_EQ(1,W.get().val());
	_EQ(1,W.get().delta());
	EXPECT_EQ(&W,&(W<<value<double>(1,1)));
	_EQ(1,W.get().val());
	_EQ(1.0/sqrt(2.0),W.get().delta());
	EXPECT_EQ(&W,&(W<<value<double>(1,1)));
	_EQ(1,W.get().val());
	_EQ(1.0/sqrt(3.0),W.get().delta());
}
TEST(WeightedAverageCalculator,Zeros_plus_Ones){
	WeightedAverageCalculator<double> W;
	EXPECT_THROW(W.get(),Exception<WeightedAverageCalculator<double>>);
	EXPECT_EQ(&W,&(W<<value<double>(1,1)));
	_EQ(1,W.get().val());
	_EQ(1,W.get().delta());
	EXPECT_EQ(&W,&(W<<value<double>(0,1)));
	_EQ(0.5,W.get().val());
	_EQ(1.0/sqrt(2.0),W.get().delta());
	EXPECT_EQ(&W,&(W<<value<double>(1,1)));
	_EQ(2.0/3.0,W.get().val());
	_EQ(1.0/sqrt(3.0),W.get().delta());
	EXPECT_EQ(&W,&(W<<value<double>(0,1)));
	_EQ(0.5,W.get().val());
	_EQ(1.0/sqrt(4.0),W.get().delta());
}
