// this file is distributed under 
// MIT license
#include <math_h/functions.h>
#include <math_h/randomfunc.h>
#include <math_h/hists.h>
#include <gnuplot_wrap.h>
using namespace std;
using namespace MathTemplates;
using namespace GnuplotWrap;
int main(){
	auto F=[](double x){return Lorentzian(x,2.0,0.2);};
	auto x_values_chain=ChainWithStep(0.0,0.01,4.0);
	const double binwidth=0.1;
	const size_t count =1000;
	RandomValueGenerator<double,mt19937> distr(F,x_values_chain);
	
	mt19937 engine;
	Distribution1D<double> hist(BinsByStep(0.0,binwidth,4.0));
	for(size_t i=0;i<count;i++)hist.Fill(distr(engine));
	
	Plotter::Instance().SetOutput(".","randomfunc");
	Plot<double>().Hist(hist).Line(SortedPoints<double>(F,x_values_chain)*(binwidth*count));
	return 0;
}