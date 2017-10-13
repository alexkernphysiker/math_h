// this file is distributed under
// MIT license
#include <math_h/interpolate.h>
#include <gnuplot_wrap.h>
using namespace std;
using namespace MathTemplates;
using namespace GnuplotWrap;
int main()
{
    //create the table of points for interpolation
    LinearInterpolation<> test=Points<>{{0, 0}, {1, 1}, {2, 4}, {3, 9}};
    test << point<>(-1, 1) << point<>(-2, 4) << point<>(-3, 9);
    //create chain of values with lesser step
    const auto x_chain_for_interpolating = ChainWithStep(-3.0, 0.01, 3.0);
    //plot interpolated tables with lesser step
    Plot<>("interpolation1").Line(SortedPoints<>(test.func(), x_chain_for_interpolating), "interpolation")
    .Line(SortedPoints<>(IntegratedLinearInterpolation<>(test).func(), x_chain_for_interpolating), "integrated");
    return 0;
}
