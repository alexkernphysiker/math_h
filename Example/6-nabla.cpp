//This is an example how to use the math_h libary
#include <iostream>
#include <math_h/numeric_diff.h>
#include <gnuplot_wrap.h>
using namespace std;
using namespace MathTemplates;
using namespace GnuplotWrap;
//This is an example how to use nabla operator
int main()
{
    const double delta = 0.001;
    const auto Potential = [](const Vector<3>&r){ return 1.0 / r.length_sqr(); };
	const auto field = -nabla<3>(delta) * Potential;
	const auto charge_density = -nabla<3>(delta) * field;

    Plot("nabla-example-potential").Line(
        SortedPoints<>(
            [&Potential](const double&c){return Potential(X()*c);},
            ChainWithStep(0.1,0.02,2.1)
        )
    );
	Plot("nabla-example-field").Line(
        SortedPoints<>(
            [&field](const double&c){return field(c*X()).length();},
            ChainWithStep(0.1,0.02,2.1)
        )
    );
    Plot("nabla-example-charge").Line(
        SortedPoints<>(
            [&charge_density](const double&c){return charge_density(c*X());},
            ChainWithStep(0.1,0.02,2.1)
        )
    );
}
