#include <fstream>
#include <stdio.h>
#include <math.h>
#include <sympson.h>
#include <singleparam.h>
using namespace std;
const double step=0.0001;
const double step2=0.01;
int main(int,char**){
	printf("using double function(double):\n");
	printf("Integral(0,1) x^(1/2) dx = %f\n",Sympson(sqrt,0.0,1.0,step));
	printf("Integral(1,0) x^(1/2) dx = %f\n",Sympson(sqrt,1.0,0.0,step));
	printf("using singleparam:\n");
	for(int k=0; k<10; k++){
		SingleParam<double,0,double,double> p(&pow,INFINITY,double(k));
		printf("Ittegral(0,1) x^%i dx = %f \n",k,Sympson(p,0.0,1.0,step));
	}
	printf("using lambda:\n");
	for(int k=0; k<10; k++)
		printf("Ittegral(0,1) x^%i dx = %f \n",k,Sympson([k](double x){return pow(x,k);},0.0,1.0,step));
	printf("double integral using lambdas:\n");
	double doubleintvalue=Sympson([](double x){
		double res= Sympson([x](double y){
			return exp(-y*y);//this may depend on x and y
		},-20.0,20.0,step2);
		res*=exp(-x*x);// this may depend only on x
		return res;
	},-20.0,20.0,step2);
	printf("Integral ~gauss(x)*gauss(y) = %f \n",doubleintvalue);
}
