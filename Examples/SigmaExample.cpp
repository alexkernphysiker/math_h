#include <iostream>
#include <fstream>
#include <functional>
#include <functions.h>
#include <randomfunc.h>
#include <sigma.h>
using namespace std;
int main(int,char**){
	RandomValueGenerator<double> 
		rand([](double x){return Gaussian(x,3.0,1.0);},0,6,0.01);
	Sigma<double> test;
	for(int i=0;i<10000;i++)
		test.AddValue(rand());
	printf("%f %f \n",test.getAverage(),test.getSigma());
	return 0;
}