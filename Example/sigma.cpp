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
	normal_distribution<double> d1(2,0.1);
	normal_distribution<double> d2(5,1);
	Sigma<double> S1,S2,S_sum,S_mul,S_f1,S_f2;
	auto F=[](double x){return log(x);};
	for(size_t i=0;i<1000;i++){
		double a=d1(engine),b=d2(engine);
		S1<<a;S2<<b;
		S_sum<<(a+b);
		S_mul<<(a*b);
		S_f1<<F(a);
		S_f2<<F(b);
	}
	
	cout<<"magnitude 1 = "<<S1()<<endl;
	cout<<"magnitude 2 = "<<S2()<<endl;
	
	cout<<"theory sum = "<<S1()+S2()<<endl;
	cout<<"experiment sum = "<<S_sum()<<endl;
	cout<<"OK="<<S_sum().contains(S1()+S2())<<endl;
	
	cout<<"theory mul = "<<S1()*S2()<<endl;
	cout<<"experiment mul = "<<S_mul()<<endl;
	cout<<"OK="<<S_mul().contains(S1()*S2())<<endl;

	cout<<"theory func 1 = "<<func_value<double>(F,S1())<<endl;
	cout<<"experiment func 1 = "<<S_f1()<<endl;
	cout<<"OK="<<S_f1().contains(func_value<double>(F,S1()))<<endl;
	
	cout<<"theory func 2 = "<<func_value<double>(F,S2())<<endl;
	cout<<"experiment func 2 = "<<S_f2()<<endl;
	cout<<"OK="<<S_f2().contains(func_value<double>(F,S2()))<<endl;
}