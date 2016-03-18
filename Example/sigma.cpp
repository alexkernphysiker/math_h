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
	for(size_t i=0;i<10000;i++){
		double a=d1(engine),b=d2(engine);
		S1<<a;S2<<b;
		S_sum<<(a+b);
		S_mul<<(a*b);
		S_f<<F(a+b);
	}
	
	cout<<"magnitude 1 = "<<S1().val()<<"+/-"<<S1().delta()<<endl;
	cout<<"magnitude 2 = "<<S2().val()<<"+/-"<<S2().delta()<<endl;
	auto th=S1()+S2();
	cout<<"theory sum = "<<th.val()<<"+/-"<<th.delta()<<endl;
	cout<<"experiment sum = "<<S_sum().val()<<"+/-"<<S_sum().delta()<<endl;
	th=S1()*S2();
	cout<<"theory mul = "<<th.val()<<"+/-"<<th.delta()<<endl;
	cout<<"experiment mul = "<<S_mul().val()<<"+/-"<<S_mul().delta()<<endl;
	th=F*(S1()+S2());
	cout<<"theory func = "<<th.val()<<"+/-"<<th.delta()<<endl;
	cout<<"experiment func = "<<S_f().val()<<"+/-"<<S_f().delta()<<endl;
}