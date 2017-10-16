// this file is distributed under
// MIT license
#include <math_h/functions.h>
#include <math_h/randomfunc.h>
#include <math_h/hists.h>
#include <gnuplot_wrap.h>
using namespace std;
using namespace MathTemplates;
using namespace GnuplotWrap;
//This is an example how to generate random values
//distributed by a function given in a table
int main()
{
    //test table with distribution density
     const LinearInterpolation<> P = Points<>{
        {0.7, 0.4}, {1.0, 0.0}, {1.8, 1.0}, {2.2, 1.0},
        {3.0, 0.5}, {3.5, 0.3}, {4.0, 0.3}
    };
    //creating random value generator
    RANDOM engine;
    const RandomValueTableDistr<> generator = P;
    //Generate values and plot histogram
    PlotDistr1D<> dist("Test", "random value", BinsByStep(0.0, 0.05, 5.0),"randomfunc");
    for (size_t i = 0; i < 1000000; i++) {
        dist.Fill(generator(engine));
    }
    return 0;
}
