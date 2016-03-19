// this file is distributed under 
// MIT license
#include <math_h/functions.h>
#include <math_h/structures.h>
#include <math_h/randomfunc.h>
#include <gnuplot_wrap.h>
using namespace std;
using namespace MathTemplates;
using namespace GnuplotWrap;
int main(){
	RandomValueGenerator<double,mt19937> X([](double x){return FermiFunc(x,2.0,0.2);},ChainWithStep(0.0,0.01,4.0));
	RandomValueGenerator<double,mt19937> Y([](double y){return FermiFunc(y,2.0,0.5);},ChainWithStep(0.0,0.01,4.0));
	
	mt19937 engine;
	Distribution2D<double> hist(BinsByCount(40,0.0,4.0),BinsByCount(40,0.0,4.0));
	for(size_t i=0;i<100000;i++)hist.Fill(make_pair(X(engine),Y(engine)));
	
	Plotter::Instance().SetOutput(".","randomfunc2");
	PlotHist2d<double>(sp2).Distr(hist);
	return 0;
}