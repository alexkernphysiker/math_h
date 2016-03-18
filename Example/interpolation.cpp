// this file is distributed under 
// MIT license
#include <math_h/interpolate.h>
#include <gnuplot_wrap.h>
using namespace std;
using namespace MathTemplates;
using namespace GnuplotWrap;
int main(){
	Plotter::Instance().SetOutput(".","interpolation");
	Plot<double> output;
	
	LinearInterpolation<double> test;
	test<<point<double>(0,0)<<point<double>(1,1)<<point<double>(2,4)<<point<double>(3,9);
	test<<point<double>(-1,1)<<point<double>(-2,4)<<point<double>(-3,9);
	
	auto x_chain=ChainWithStep(-3.0,0.1,3.0);
	output.Line(SortedPoints<double>(test.func(),x_chain),"y");
	output.Line(SortedPoints<double>(LinearInterpolation<double>(test*2.0).func(),x_chain),"y*2");
	
	return 0;
}