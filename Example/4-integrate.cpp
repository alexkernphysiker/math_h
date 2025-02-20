//This is an example how to use the math_h libary
#include <iostream>
#include <math_h/functions.h>
#include <math_h/integrate.h>
#include <gnuplot_wrap.h>
using namespace std;
using namespace MathTemplates;
using namespace GnuplotWrap;
// This is an example of using 
// functions for numeric integrating
int main()
{
    //How you can integrate functions
    for (size_t k = 0; k < 5; k++)
        cout
        << "int(0,1)x^" << k << " = "
        << Sympson(
            [k](double x) {return pow(x, k);},
            0.0, 1.0, 0.001
        )
        << endl;

    //How you can integrate table data
    const SortedPoints<> table(
        [](double x) {return Gaussian(x, 0.0, 1.0);},
        ChainWithStep(-5.0, 0.1, 5.0)
    );
    const auto integrated_table = Int_Trapez_Table(table);
    Plot("integrate").Points(table, "density").Points(integrated_table, "integrated");

    //How you can calculate convolution integral
    const auto conv = make_convolution(
        // first function
        [](double ksi) { return Gaussian(ksi, 1.5, 0.5); },
        // second function
        [](double ksi) { if (ksi < 0) return 0.0; else return exp(-ksi / 1.5); },
        //integrating range and delta
        -20., 20., 0.001
    );
    //calculating the integral values and fill the table to plot
    SortedPoints<> plot_conv;
    for (double x = -2.;x <= 10.;x += 0.1) {
        const auto y = conv(x);
        plot_conv << make_point(x, y);
    }
    Plot("convolution", 1).Points(plot_conv);

}
