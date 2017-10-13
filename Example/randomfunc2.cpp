// this file is distributed under
// MIT license
#include <math_h/functions.h>
#include <math_h/randomfunc.h>
#include <math_h/hists.h>
#include <gnuplot_wrap.h>
using namespace std;
using namespace MathTemplates;
using namespace GnuplotWrap;
int main()
{
    const LinearInterpolation<> P = Points<>{
        {0.7, 0.4}, {1.0, 0.0}, {1.8, 1.0}, {2.2, 1.0},
        {3.0, 0.5}, {3.5, 0.3}, {4.0, 0.3}
    };
    const RandomValueTableDistr<> generator = P;
    RANDOM engine;
    PlotDistr1D<> dist("Test", "random value", BinsByStep(0.0, 0.05, 5.0),"randomfunc2");
    for (size_t i = 0; i < 1000000; i++) {
        dist.Fill(generator(engine));
    }
    return 0;
}
