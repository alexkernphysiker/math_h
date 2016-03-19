// this file is distributed under 
// MIT license
#include <iostream>
#include <random>
#include <math_h/sigma.h>
#include <math_h/structures.h>
#include <gnuplot_wrap.h>
using namespace std;
using namespace MathTemplates;
using namespace GnuplotWrap;
int main(){
	mt19937 engine;
	normal_distribution<double> d1(2,1);
	normal_distribution<double> d2(5,1);
	function<double(double)> F=[](double x){return log(x);};
	Sigma<double> S1,S2,S_sum,S_mul,S_f;
	for(size_t i=0;i<1000;i++){
		double a=d1(engine),b=d2(engine);
		S1<<a;S2<<b;
		S_sum<<(a+b);
		S_mul<<(a*b);
		S_f<<F(b);
	}
	
	cout<<"magnitude 1 = "<<S1()<<endl;
	cout<<"magnitude 2 = "<<S2()<<endl;
	cout<<"theory sum = "<<S1()+S2()<<endl;
	cout<<"experiment sum = "<<S_sum()<<endl;
	cout<<"theory mul = "<<S1()*S2()<<endl;
	cout<<"experiment mul = "<<S_mul()<<endl;
	cout<<"theory func = "<<F*S2()<<endl;
	cout<<"experiment func = "<<S_f()<<endl;
}