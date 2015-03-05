#include <iostream>
#include <fstream>
#include <math.h>
#define USE_RANDOM_DEVICE
#include <functions.h>
#include <randomfunc.h>
#include <sigma.h>
#include <singleparam.h>
int main(int,char**){
	SingleParam<double,0,double,double,double> G(&Gaussian,0,3,1);
	RandomValueGenerator<double,decltype(G)> rand(G);
	rand.Init(0,6,0.01);
	Sigma<double> test;
	for(int i=0;i<10000;i++)
		test.AddValue(rand());
	printf("%f %f \n",test.getAverage(),test.getSigma());
	return 0;
}