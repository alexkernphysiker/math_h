// this file is distributed under
// MIT license
#include <math_h/interpolate.h>
#include <gnuplot_wrap.h>
using namespace std;
using namespace MathTemplates;
using namespace GnuplotWrap;
//This is an example how to use
// linear interpolation class
int main()
{
    //create the table of points for interpolation
    LinearInterpolation<> test=Points<>{{0, 0}, {1, 1}, {2, 4}, {3, 9}};
    //you can add more points later
    test << point<>(-1, 1) << point<>(-2, 4) << point<>(-3, 9);
    //create chain of values with lesser step
    SortedPoints<> plot_data;
    for(double x=-3.;x<=3.;x+=0.1)plot_data<<make_point(x,test(x));
    Plot<>("interpolation1").Points(plot_data);
    return 0;
}
