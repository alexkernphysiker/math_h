// this file is distributed under 
// MIT license
#include <math_h/interpolate.h>
#include <math_h/sigma.h>
#include <gnuplot_wrap.h>
using namespace std;
using namespace MathTemplates;
using namespace GnuplotWrap;
int main(){
    LinearInterpolation<value<double>> test{
	point<value<double>>(-0.43,{ 1.3,12.0}),
	point<value<double>>( 0.04,{160.,13.6}),
	point<value<double>>( 0.75,{390.,16.8}),
	point<value<double>>( 2.58,{412.,16.5}),
	point<value<double>>( 5.00,{410.,18.1}),
	point<value<double>>( 8.01,{390.,19.0}),
	point<value<double>>(11.01,{385.,19.9}),
	point<value<double>>(20.,{320.,20.}),
	point<value<double>>(40.,{420.,20.})
    };
    const auto x_chain_for_interpolating=BinsByStep(0.,2.5,30.0);
    Plotter::Instance().SetOutput(".","interpolation2");
    Plot<double> output;
    output.Hist(hist<double>(test.func(),x_chain_for_interpolating),"interpolation");
    output.Hist(test,"points")<<"set key on";
    return 0;
}
