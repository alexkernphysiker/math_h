#include <math_h/functions.h>
#include <math_h/randomfunc.h>
#include <math_h/sigma3.h>
#include <gnuplot_wrap.h>
using namespace std;
using namespace MathTemplates;
using namespace GnuplotWrap;
//This is an example how to generate random values
//distributed by a function given in a table
int main()
{
    //Let's assume this to be measurements
    const RandomValueTableDistr<> generator = Points<>{
        {0.7, 0.4}, {1.0, 0.0}, {1.8, 1.0}, {2.2, 1.0},
        {3.0, 0.5}, {3.5, 0.3}, {4.0, 0.3}
    };
    Distribution1D<> measurements(BinsByStep(0.0, 0.1, 5.0));
    for (size_t i = 0; i < 1000; i++)measurements.Fill(generator());
    const auto data=extend_hist<1,2>(measurements);//only statistical uncertainty
    //Let's assume this to be a coefficient known with some uncertainty 
    const auto theoretical_coefficient=extend_value<2,2>(value<>(3.5,0.5));//only systematical uncertainty
    //Let's save multiplied data with both uncertainties
    const auto result=data*theoretical_coefficient;
    Plot("sigma3-example").Hist_2bars<1,2>(result,"","sigma3-output");
    Plot("sigma3-example-2").Hist(wrap_hist(result));
    return 0;
}

