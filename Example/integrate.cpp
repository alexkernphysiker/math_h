// this file is distributed under
// MIT license
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
	<< Sympson([k](double x) {return pow(x, k);}, 0.0, 1.0, 0.001) << endl;

    //How you can integrate table data
    SortedPoints<> density([](double x){return Gaussian(x, 0.0, 1.0);}, ChainWithStep(-5.0, 0.1, 5.0));
    Plot<>("integrate").Points(density, "density").Points(Int_Trapez_Table(density), "integrated");

    //How you can calculate convolution integral
    const auto conv=make_convolution(
	[](double ksi){return Gaussian(ksi,1.5,0.2);},
	[](double ksi){if(ksi<0)return 0.0;return exp(-ksi/1.5);},
				     -20.,20.,0.01);
    SortedPoints<> plot_conv;
    for(double x=-0.;x<=10.;x+=0.1)plot_conv<<make_point(x,conv(x));
    Plot<>("convolution").Points(plot_conv);
}
