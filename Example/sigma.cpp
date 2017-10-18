// this file is distributed under
// MIT license
#include <iostream>
#include <math_h/sigma.h>
using namespace std;
using namespace MathTemplates;
//This is an example how to use
//values with uncertainties
int main()
{
    //how to declare values with uncertainties
    const value<> A(10.,1.0),B(3.,0.2);
    //how to get value's components
    cout<<"A: value = "<<A.val()<<endl;
    cout<<"A: delta = "<<A.uncertainty()<<endl;
    cout<<"A: epsilon = "<<A.epsilon()<<endl;
    cout<<"A: min = "<<A.min()<<endl;
    cout<<"A: max = "<<A.max()<<endl;
    //how to perform arithmetic operations
    cout<<"A = "<<A<<endl;
    cout<<"B = "<<B<<endl;
    cout<<"A+B = "<<A+B<<endl;
    cout<<"A-B = "<<A-B<<endl;
    cout<<"A*B = "<<A*B<<endl;
    cout<<"A/B = "<<A/B<<endl;
    cout<<"log(A) = "<<A.Func([](double x){return log(x);})<<endl;
}
