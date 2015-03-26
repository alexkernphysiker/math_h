#include <fstream>
#include <stdio.h>
#include <math.h>
#include <sympson.h>
using namespace std;
const double step=0.0001;
const double step2=0.01;
int main(int,char**){
	printf("using double function(double):\n");
	printf("Integral(0,1) x^(1/2) dx = %f\n",Sympson(sqrt,0.0,1.0,step));
	printf("Integral(1,0) x^(1/2) dx = %f\n",Sympson(sqrt,1.0,0.0,step));
	printf("using lambda:\n");
	for(int k=0; k<10; k++)
		printf("Ittegral(0,1) x^%i dx = %f \n",k,Sympson([k](double x){return pow(x,k);},0.0,1.0,step));
	printf("double integral using lambdas:\n");
	printf("Integral ~gauss(x)*gauss(y) = %f \n",Sympson([](double x){
		return Sympson([x](double y){
			return exp(-y*y);
		},-20.0,20.0,step2)*exp(-x*x);
	},-20.0,20.0,step2));
}
