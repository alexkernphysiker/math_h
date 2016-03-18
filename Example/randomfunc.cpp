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
	RandomValueGenerator<double,mt19937> distr([](double x){return FermiFunc(x,2.0,0.2);},ChainWithStep(0.0,0.01,4.0));
	
	mt19937 engine;
	Distribution1D<double> hist(BinsByCount(20,0.0,4.0));
	for(size_t i=0;i<5000;i++)hist.Fill(distr(engine));
	
	Plotter::Instance().SetOutput(".","randomfunc");
	Plot<double>().Hist(hist);
	return 0;
}