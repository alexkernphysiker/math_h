// this file is distributed under 
// MIT license
#include <random>
#include <math.h>
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
#define _EQ(a,b) EXPECT_TRUE(pow(a-b,2)<0.005)
TEST(value,func_val){
	value<double> V(1,0.1);
	auto V2=func_value<double>([](double x){return 2.0*x;},V);
	EXPECT_EQ(2*V.val(),V2.val());
	_EQ(2*V.delta(),V2.delta());
}
TEST(value,func_val2){
	value<double> V(1,0.1);
	auto F=value_func<double>([](double x){return 2.0*x;});
	EXPECT_EQ(2*V.val(),F(V).val());
	_EQ(2*V.delta(),F(V).delta());
}
TEST(StandardDeviation,Throwing){
	StandardDeviation<double> S;
	EXPECT_EQ(0,S.count());
	EXPECT_THROW(S(),Exception<StandardDeviation<double>>);
	EXPECT_EQ(&S,&(S<<0.0));
	EXPECT_EQ(1,S.count());
	EXPECT_THROW(S(),Exception<StandardDeviation<double>>);
	EXPECT_EQ(&S,&(S<<0.0));
	EXPECT_EQ(2,S.count());
	EXPECT_EQ(0,S().val());
	EXPECT_EQ(0,S().delta());
}
TEST(StandardDeviation,Base){
	StandardDeviation<double> S;
	S<<0.0;
	EXPECT_EQ(1,S.count());
	EXPECT_THROW(S(),Exception<StandardDeviation<double>>);
	S<<1.0;
	EXPECT_EQ(2,S.count());
	EXPECT_EQ(0.5,S().val());
	EXPECT_EQ(sqrt(0.5),S().delta());
}
#define _EQ2(a,b) EXPECT_TRUE(pow(a-b,2)<0.01)
TEST(StandardDeviation,WithRandomValues){
	StandardDeviation<double> S;
	default_random_engine generator;
	normal_distribution<double> distribution(1.0,3.0);
	for(int i=0;i<2000;i++)
		S<<distribution(generator);
	_EQ2(1.0,S().val());
	_EQ2(3.0,S().delta());
}
TEST(StandardDeviation,WithRandomValues2){
	StandardDeviation<double> S(2);
	default_random_engine generator;
	normal_distribution<double> distribution(1.0,3.0);
	for(int i=0;i<2000;i++)
		S<<distribution(generator);
	_EQ2(1.0,S().val());
	_EQ2(6.0,S().delta());
}

TEST(WeightedAverageCalculator,Zeros){
	WeightedAverageCalculator<double> W;
	EXPECT_THROW(W(),Exception<WeightedAverageCalculator<double>>);
	EXPECT_THROW(W<<value<double>(0,0),Exception<WeightedAverageCalculator<double>>);
	EXPECT_EQ(&W,&(W<<value<double>(0,1)));
	_EQ(0,W().val());
	_EQ(1,W().delta());
	EXPECT_EQ(&W,&(W<<value<double>(0,1)));
	_EQ(0,W().val());
	_EQ(1.0/sqrt(2.0),W().delta());
	EXPECT_EQ(&W,&(W<<value<double>(0,1)));
	_EQ(0,W().val());
	_EQ(1.0/sqrt(3.0),W().delta());
}
TEST(WeightedAverageCalculator,Ones){
	WeightedAverageCalculator<double> W;
	EXPECT_THROW(W(),Exception<WeightedAverageCalculator<double>>);
	EXPECT_THROW(W<<value<double>(1,0),Exception<WeightedAverageCalculator<double>>);
	EXPECT_EQ(&W,&(W<<value<double>(1,1)));
	_EQ(1,W().val());
	_EQ(1,W().delta());
	EXPECT_EQ(&W,&(W<<value<double>(1,1)));
	_EQ(1,W().val());
	_EQ(1.0/sqrt(2.0),W().delta());
	EXPECT_EQ(&W,&(W<<value<double>(1,1)));
	_EQ(1,W().val());
	_EQ(1.0/sqrt(3.0),W().delta());
}
TEST(WeightedAverageCalculator,Zeros_plus_Ones){
	WeightedAverageCalculator<double> W;
	EXPECT_THROW(W(),Exception<WeightedAverageCalculator<double>>);
	EXPECT_EQ(&W,&(W<<value<double>(1,1)));
	_EQ(1,W().val());
	_EQ(1,W().delta());
	EXPECT_EQ(&W,&(W<<value<double>(0,1)));
	_EQ(0.5,W().val());
	_EQ(1.0/sqrt(2.0),W().delta());
	EXPECT_EQ(&W,&(W<<value<double>(1,1)));
	_EQ(2.0/3.0,W().val());
	_EQ(1.0/sqrt(3.0),W().delta());
	EXPECT_EQ(&W,&(W<<value<double>(0,1)));
	_EQ(0.5,W().val());
	_EQ(1.0/sqrt(4.0),W().delta());
}
