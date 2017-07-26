// this file is distributed under 
// MIT license
#include <math_h/interpolate.h>
#include <gnuplot_wrap.h>
using namespace std;
using namespace MathTemplates;
using namespace GnuplotWrap;
int main(){
    LinearInterpolation<double> test{{0,0},{1,1},{2,4},{3,9}};
    test<<point<double>(-1,1)<<point<double>(-2,4)<<point<double>(-3,9);
    const auto x_chain_for_interpolating=ChainWithStep(-3.0,0.01,3.0);
    Plotter::Instance().SetOutput(".","interpolation1");
    Plot<double>().Line(SortedPoints<double>(test.func(),x_chain_for_interpolating),"interpolation")
    .Line(SortedPoints<double>(IntegratedLinearInterpolation<>(test).func(),x_chain_for_interpolating),"integrated");
    const ReverseIntegratedLinearInterpolation<> RI=test;
    Plot<double>().Line(SortedPoints<double>(RI.func(),ChainWithStep(RI.left().X(),0.01,RI.right().X())),"reverse integrated");
    return 0;
}
