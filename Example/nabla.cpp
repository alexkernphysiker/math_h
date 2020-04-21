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
    Plot("nabla-example-potential")
	.Line(SortedPoints<>([Potential](const double&c){
	    return Potential(X()*c);
	},ChainWithStep(-10.1,0.2,10.1)));
	Plot("nabla-example-field")
	.Line(SortedPoints<>([Potential](const double&c){
		auto field=nabla<3>(X()*c,0.01)*Potential;
	    return field.x();
	},ChainWithStep(-10.1,0.2,10.1)));
    Plot("nabla-example-charge")
	.Line(SortedPoints<>([Potential](const double&c){
		auto field = nabla<3>(0.01)*Potential;
		auto charge_density = nabla<3>(0.01)*field;
	    return charge_density(X()*c);
	},ChainWithStep(-10.1,0.2,10.1)));
}
