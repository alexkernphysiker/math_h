//This is an example how to use the math_h libary
#include <iostream>
#include <math_h/tabledata.h>
#include <gnuplot_wrap.h>
using namespace std;
using namespace MathTemplates;
using namespace GnuplotWrap;
// This is an example of using table data
// and plotting them
int main()
{
    const Points<> table1{
        {0.7, 0.4}, {1.0, 0.0}, {1.8, 1.0}, {2.2, 1.0},
        {3.0, 0.5}, {3.5, 0.3}, {4.0, 0.3}
    };
    const SortedPoints table2(
        [](double x){return Gaussian(x, 2.5, 1.0);}, 
	    ChainWithStep(0.0, 0.1, 5.0)
    );
    Plot("test-plot").Points(table1, "table1").Points(table2, "table2") << "set key on" << "set yrange [0:2]";
}
