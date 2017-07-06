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
	RandomValueTableDistr<> generator(F,x_values_chain);
	RANDOM engine;
	Plotter::Instance().SetOutput(".","randomfunc1");
	PlotDistr1D<> dist("Test","random value",BinsByStep(0.0,0.1,4.0));
	for(size_t i=0;i<10000;i++){
	    dist.Fill(generator(engine));
	}
	return 0;
}
