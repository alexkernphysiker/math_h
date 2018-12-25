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
    const auto Potential=[](const Vector<3>&r){return exp(-r.M_sqr());};
    Plot("nabla-example-potential").Line(SortedPoints<>([Potential](const double&c){return Potential(X()*c);},ChainWithStep(-10.,0.1,10.)));
    Plot("nabla-example-field").Line(SortedPoints<>([Potential](const double&c){return (nabla(X()*c,0.01)*Potential).x();},ChainWithStep(-10.,0.1,10.)));
}
