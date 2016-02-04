// this file is distributed under 
// MIT license
#include <random>
#include <gtest/gtest.h>
#include <math_h/sigma.h>
using namespace std;
using namespace MathTemplates;
TEST(Sigma,Throwing){
	Sigma<double> S;
	EXPECT_EQ(0,S.count());
	EXPECT_THROW(S.getValue(),Exception<Sigma<double>>);
	EXPECT_EQ(&S,&(S.AddValue(0)));
	EXPECT_EQ(1,S.count());
	EXPECT_THROW(S.getValue(),Exception<Sigma<double>>);
	EXPECT_EQ(&S,&(S.AddValue(0)));
	EXPECT_EQ(2,S.count());
	EXPECT_EQ(0,S.getValue().val());
	EXPECT_EQ(0,S.getValue().delta());
}
#define _EQ(a,b) EXPECT_TRUE(pow(a-b,2)<0.005)
TEST(Sigma,Base){
	Sigma<double> S;
	S.AddValue(0);
	EXPECT_EQ(1,S.count());
	EXPECT_THROW(S.getValue(),Exception<Sigma<double>>);
	S.AddValue(1);
	EXPECT_EQ(2,S.count());
	EXPECT_EQ(0.5,S.getValue().val());
	EXPECT_EQ(sqrt(0.5),S.getValue().delta());
}
#define _EQ2(a,b) EXPECT_TRUE(pow(a-b,2)<0.01)
TEST(Sigma,WithRandomValues){
	Sigma<double> S;
	default_random_engine generator;
	normal_distribution<double> distribution(1.0,3.0);
	for(int i=0;i<2000;i++)
		S.AddValue(distribution(generator));
	_EQ2(1.0,S.getValue().val());
	_EQ2(3.0,S.getValue().delta());
}

TEST(WeightedAverageCalculator,Zeros){
	WeightedAverageCalculator<double> W;
	EXPECT_THROW(W.Average(),Exception<WeightedAverageCalculator<double>>);
	EXPECT_THROW(W.Sigma(),Exception<WeightedAverageCalculator<double>>);
	EXPECT_THROW(W.AddValue(0,0),Exception<WeightedAverageCalculator<double>>);
	EXPECT_EQ(&W,&(W.AddValue(0,1)));
	_EQ(0,W.Average());
	_EQ(1,W.Sigma());
	EXPECT_EQ(&W,&(W.AddValue(0,1)));
	_EQ(0,W.Average());
	_EQ(1.0/sqrt(2.0),W.Sigma());
	EXPECT_EQ(&W,&(W.AddValue(0,1)));
	_EQ(0,W.Average());
	_EQ(1.0/sqrt(3.0),W.Sigma());
}
TEST(WeightedAverageCalculator,Ones){
	WeightedAverageCalculator<double> W;
	EXPECT_THROW(W.Average(),Exception<WeightedAverageCalculator<double>>);
	EXPECT_THROW(W.Sigma(),Exception<WeightedAverageCalculator<double>>);
	EXPECT_THROW(W.AddValue(1,0),Exception<WeightedAverageCalculator<double>>);
	EXPECT_EQ(&W,&(W.AddValue(1,1)));
	_EQ(1,W.Average());
	_EQ(1,W.Sigma());
	EXPECT_EQ(&W,&(W.AddValue(1,1)));
	_EQ(1,W.Average());
	_EQ(1.0/sqrt(2.0),W.Sigma());
	EXPECT_EQ(&W,&(W.AddValue(1,1)));
	_EQ(1,W.Average());
	_EQ(1.0/sqrt(3.0),W.Sigma());
}
TEST(WeightedAverageCalculator,Zeros_plus_Ones){
	WeightedAverageCalculator<double> W;
	EXPECT_THROW(W.Average(),Exception<WeightedAverageCalculator<double>>);
	EXPECT_THROW(W.Sigma(),Exception<WeightedAverageCalculator<double>>);
	EXPECT_THROW(W.AddValue(0,0),Exception<WeightedAverageCalculator<double>>);
	EXPECT_EQ(&W,&(W.AddValue(1,1)));
	_EQ(1,W.Average());
	_EQ(1,W.Sigma());
	EXPECT_EQ(&W,&(W.AddValue(0,1)));
	_EQ(0.5,W.Average());
	_EQ(1.0/sqrt(2.0),W.Sigma());
	EXPECT_EQ(&W,&(W.AddValue(1,1)));
	_EQ(2.0/3.0,W.Average());
	_EQ(1.0/sqrt(3.0),W.Sigma());
	EXPECT_EQ(&W,&(W.AddValue(0,1)));
	_EQ(0.5,W.Average());
	_EQ(1.0/sqrt(4.0),W.Sigma());
}
