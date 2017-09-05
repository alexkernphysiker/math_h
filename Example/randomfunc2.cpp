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
    RandomValueTableDistr<> X([](double x) {
        return FermiFunc(x, 2.0, 0.2);
    }, ChainWithStep(0.0, 0.01, 4.0));
    RandomValueTableDistr<> Y([](double y) {
        return FermiFunc(y, 2.0, 0.5);
    }, ChainWithStep(0.0, 0.01, 4.0));

    mt19937 engine;
    Distribution2D<> H(BinsByCount(40, 0.0, 4.0), BinsByCount(40, 0.0, 4.0));
    for (size_t i = 0; i < 10000; i++)H.Fill(X(engine), Y(engine));

    Plotter<>::Instance().SetOutput(".", "randomfunc2");
    PlotHist2d<>(sp2).Distr(H);
    return 0;
}
