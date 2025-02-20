//This is an example how to use the math_h libary
#include <iostream>
#include <math_h/numeric_diff.h>
#include <math_h/vectortransformations.h>
#include <gnuplot_wrap.h>
using namespace std;
using namespace MathTemplates;
using namespace GnuplotWrap;
//This is an example how to use nabla operator
int main()
{
    const double delta = 0.00001;
    const auto Potential = [](const Vector<3>& r) { return 1.0 / r.length(); };
    const auto field = -nabla<3>(delta) * Potential;
    const auto charge_density = laplace<3>(delta) * Potential;

    Plot("nabla-example-potential").Line(
        SortedPoints<>(
            [&Potential](const double& c) {return Potential(c * X() + 0.5 * Y());},
            ChainWithStep(-2.1, 0.02, 2.1)
        )
    );
    Plot("nabla-example-field").Line(
        SortedPoints<>(
            [&field](const double& c) {return field(c * X() + 0.5 * Y()).length();},
            ChainWithStep(-2.1, 0.02, 2.1)
        )
    );
    Plot("nabla-example-charge").Line(
        SortedPoints<>(
            [&charge_density](const double& c) {return charge_density(c * X() + 0.5 * Y());},
            ChainWithStep(-2.1, 0.02, 2.1)
        )
    );
}
