// this file is distributed under
// MIT license
#include <math_h/interpolate.h>
#include <math_h/sigma.h>
#include <gnuplot_wrap.h>
using namespace std;
using namespace MathTemplates;
using namespace GnuplotWrap;
int main()
{
    //create table for interpolation that contains uncertainties
    LinearInterpolation<value<>> test=Points<value<>>{
	{-1., { 1.3, 12.}},
	{0., {160., 13.}},
	{2., {412., 16.}},
	{5., {410., 18.}},
	{8., {390., 19.}},
	{11., {385., 19.}},
	{20., {320., 20.}},
	{40., {420., 20.}}
    };
    //create chain of intervals
    const auto x_chain_for_interpolating = BinsByStep(0., 2.5, 30.0);
    //interpolating
    Plot<> output("interpolation2");
    output.Hist(hist<>(test.func(), x_chain_for_interpolating), "interpolation");
    output.Hist(test, "points") << "set key on" << "set xrange [-5:45]";
    return 0;
}
