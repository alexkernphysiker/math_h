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
    const auto Potential=[](const Vector<3>&r){return 1.0/r.M_sqr();};
	auto field=nabla<3>(0.01)*Potential;
	auto charge_density = nabla<3>(0.01)*field;

    Plot("nabla-example-potential").Line(SortedPoints<>([Potential](const double&c){return Potential(X()*c);},ChainWithStep(-10.1,0.2,10.1)));
	Plot("nabla-example-field").Line(SortedPoints<>([field](const double&c){return field.x(X()*c);},ChainWithStep(-10.1,0.2,10.1)));
    Plot("nabla-example-charge").Line(SortedPoints<>([charge_density](const double&c){return charge_density(X()*c);},ChainWithStep(-10.1,0.2,10.1)));
}
