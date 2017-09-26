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
    const auto x_chain_for_interpolating = BinsByStep(0., 2.5, 30.0);
    Plot<> output("interpolation3");
    output.Hist(hist<>(test.func(), x_chain_for_interpolating), "interpolation");
    output.Hist(test, "points") << "set key on" << "set xrange [-5:45]";
    return 0;
}
