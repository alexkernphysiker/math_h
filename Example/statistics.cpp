//This is an example how to use the math_h libary
#include <math_h/randomfunc.h>
#include <math_h/statistics.h>
using namespace std;
using namespace MathTemplates;
//This is an example how to generate random values
//distributed by a function given in a table
int main()
{
    RandomGauss<> X(5,2),Y(2,2);
    Sampling<3> S;
    for (size_t i = 0; i < 1000000; i++){
	const auto x=X(),y=Y();
	S.Fill(x,y,x-y);
    }
    cout<<"Average"<<endl;
    cout<<S.Average()<<endl;
    cout<<"Cov"<<endl;
    cout<<S.Cov()<<endl;
    return 0;
}
