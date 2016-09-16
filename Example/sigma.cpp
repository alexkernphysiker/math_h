// this file is distributed under 
// MIT license
#include <iostream>
#include <random>
#include <math_h/sigma.h>
#include <gnuplot_wrap.h>
using namespace std;
using namespace MathTemplates;
using namespace GnuplotWrap;
int main(){
	mt19937 engine;
	normal_distribution<double> d1(2,0.1);
	normal_distribution<double> d2(5,1);
	StandardDeviation<double> S1,S2,S_sum,S_mul,S_ratio,S_f1,S_f2;
	auto F=[](double x){return log(x);};
	for(size_t i=0;i<5000;i++){
		double a=d1(engine),b=d2(engine);
		S1<<a;S2<<b;
		S_sum<<(a+b);
		S_mul<<(a*b);
		S_ratio<<(a/b);
		S_f1<<F(a);
		S_f2<<F(b);
	}
	
	cout<<"magnitude 1 = "<<S1()<<endl;
	cout<<"magnitude 2 = "<<S2()<<endl;
	cout<<endl;
	cout<<"theory sum = "<<S1()+S2()<<endl;
	cout<<"experiment sum = "<<S_sum()<<endl;
	cout<<endl;
	cout<<"theory mul = "<<S1()*S2()<<endl;
	cout<<"experiment mul = "<<S_mul()<<endl;
	cout<<endl;
	cout<<"theory ratio = "<<S1()/S2()<<endl;
	cout<<"experiment ratio = "<<S_ratio()<<endl;
	cout<<endl;
	cout<<"theory func 1 = "<<func_value<double>(F,S1())<<endl;
	cout<<"experiment func 1 = "<<S_f1()<<endl;
	cout<<endl;
	cout<<"theory func 2 = "<<func_value<double>(F,S2())<<endl;
	cout<<"experiment func 2 = "<<S_f2()<<endl;
	cout<<endl;
}