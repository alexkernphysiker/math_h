// this file is distributed under 
// MIT license
#include <math_h/interpolate.h>
#include <gnuplot_wrap.h>
using namespace std;
using namespace MathTemplates;
using namespace GnuplotWrap;
int main(){
	LinearInterpolation<value<double>> test{
		point<value<double>>(value<double>(0.0,0.1),value<double>(0.0,0.4)),
		point<value<double>>(value<double>(1.0,0.1),value<double>(1.0,0.4)),
		point<value<double>>(value<double>(2.0,0.1),value<double>(4.0,0.4)),
		point<value<double>>(value<double>(3.0,0.1),value<double>(9.0,0.4))
	};
	
	Plotter::Instance().SetOutput(".","interpolation-perversion");
	Plot<double>().Hist(test).Hist(SortedPoints<value<double>>(test.func(),BinsByStep(0.0,0.1,3.0)));
	return 0;
}