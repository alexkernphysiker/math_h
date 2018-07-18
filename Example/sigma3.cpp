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
    Distribution1D<> measurements(BinsByStep(0.0, 0.1, 5.0));
    const value<> coefficient(3.5,0.5);
    for (size_t i = 0; i < 1000; i++)measurements.Fill(generator());

    //Let's assume this to be statistical data
    const auto data=extend_hist<1,2>(measurements);//only statistical uncertainty
    //Let's assume this to be a coefficient known with some uncertainty 
    const value<> theoretical_coefficient(coefficient);//only systematical uncertainty

    //Let's save multiplied data with both uncertainties
    const auto result=data*extend_value<2,2>(theoretical_coefficient);
    Plot("sigma3-example").Hist_2bars<1,2>(result,"stat","syst","sigma3-output")<<"set key on";//draws two error bars
    Plot("sigma3-example-2").Hist<1>(result,"statistical uncertainty","sigma3-output")<<"set key on";//draws one error bar
    Plot("sigma3-example-3").Hist(wrap_hist(result),"total uncertainty")//calculates one uncertainty taking statistical and systematical errors into account
	.Hist(measurements*coefficient,"compare")<<"set key on";//comparing with the calculation that uses one uncertainty
    return 0;
}

