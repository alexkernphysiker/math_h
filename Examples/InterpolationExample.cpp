#include <iostream>
#include <fstream>
#include <interpolate.h>
using namespace std;
int main(int,char**){
	printf("Linear interpolation example\n");
	LinearInterpolation<double> myfunc;
	myfunc<<make_pair(0,0)<<make_pair(1,1)<<make_pair(2,4);
	printf("points:\n");
	for(auto p:myfunc)
		printf("(%f,%f)\n",p.first,p.second);
	vector<double> X;
	X.push_back(0.5);
	X.push_back(1.5);
	printf("Calculated:\n");
	for(auto x:X)
		printf("myfunc(%f)=%f\n",x,myfunc(x));
	return 0;
}