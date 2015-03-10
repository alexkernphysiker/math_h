#include <fstream>
#include <stdio.h>

#include <complex.h>
using namespace std;
int main(int,char**){
	Complex<double> a(1);
	Complex<double> b(0,1);
	auto c=pow(a+b,2);
	printf("c=%f+%f i \n",c.re(),c.im());
}