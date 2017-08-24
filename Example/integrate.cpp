// this file is distributed under
// MIT license
#include <iostream>
#include <math_h/functions.h>
#include <math_h/integrate.h>
#include <gnuplot_wrap.h>
using namespace std;
using namespace MathTemplates;
using namespace GnuplotWrap;
int main()
{
    for (size_t k = 0; k < 5; k++)
        cout << "int(0,1)x^" << k << " = " << Sympson([k](double x) {
        return pow(x, k);
    }, 0.0, 1.0, 0.001) << endl;

    SortedPoints<double> density([](double x) {
        return Gaussian(x, 0.0, 1.0);
    }, ChainWithStep(-5.0, 0.1, 5.0));
    Plotter<>::Instance().SetOutput(".", "integrate");
    Plot<double>().Points(density, "density").Points(Int_Trapez_Table(density), "integrated");
}
