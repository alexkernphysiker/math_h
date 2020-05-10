//This is an example how to use the math_h libary
#include <math_h/functions.h>
#include <math_h/randomfunc.h>
#include <math_h/sigma3.h>
#include <gnuplot_wrap.h>
using namespace std;
using namespace MathTemplates;
using namespace GnuplotWrap;
//This is an example how to make calculations with
//statistical systematical uncertainties
// it's modified randomfunc.cpp program
//let's assume that random value generator generates data with statistical uncertainty
// and we need to multiply them by a coefficient known with uncertainty (systematical error)
int main()
{
    const RandomValueTableDistr<> generator = Points<>{
        {0.7, 0.4}, {1.0, 0.0}, {1.8, 1.0}, {2.2, 1.0},
        {3.0, 0.5}, {3.5, 0.3}, {4.0, 0.3}
    };

    // create distribution and fill it with data
    Distribution1D<> measurements(BinsByStep(0.0, 0.1, 5.0));
    for (size_t i = 0; i < 1000; i++)
        measurements.Fill(generator());

    //Extend our histogram to the one with two uncertainties
    // 1 - statistical; 2 - systematical
    // the first one is filled from original histogram
    // the second one is not filled (is zero)
    const auto data=extend_hist<1,2>(measurements);

    //This is the theoretical coefficient known with uncertainty
    const value<> theoretical_coefficient(3.5,0.5);

    //Multiply out histogram by the coefficient filling its uncertainty to the second "slot"
    const auto result=data * extend_value<2,2>(theoretical_coefficient);

    //Plot the result with two error bars for each point
    Plot("sigma3-example").Hist_2bars<1,2>(result, "stat", "syst", "sigma3-output") << "set key on";

    //Plot result with only the first error bars
    Plot("sigma3-example-2").Hist<1>(result, "statistical uncertainty", "sigma3-output") << "set key on";

    //Plot the result with one error bar taking into account both uncertainties
    //comparing with the calculation that uses one uncertainty
    Plot("sigma3-example-3")
        .Hist(wrap_hist(result), "total uncertainty")
	    .Hist(measurements * theoretical_coefficient, "compare")
        << "set key on";
    return 0;
}

