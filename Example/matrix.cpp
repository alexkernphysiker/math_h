// this file is distributed under
// MIT license
#include <math_h/matrices.h>
using namespace std;
using namespace MathTemplates;
//This is an example how to solve Matrix equation A*x=b
//where A is a matrix, b is a known vector 
// and x is an unknown vector
int main()
{
    const MatrixData<double> A({
        { 1., 2., 0. },
        { 3., 4., 4. },
        { 5., 6., 3. },
    }),
    b({{3}, {7}, {8}});
    const MatrixData<double> x = Solve(A, b);
    cout << A << endl << x << endl << b << endl;
    return 0;
}
