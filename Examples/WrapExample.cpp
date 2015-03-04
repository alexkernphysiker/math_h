#include <fstream>
#include <stdio.h>
#include <math.h>
#include <functions.h>
#define use_num_type double
#define use_indexer_type double*
#include <wrap3.h>
#undef use_num_type
#undef use_indexer_type
double params[]={1,0.5,3};
using namespace FuncWrappers_xP;
using namespace std;
double a;
int main(int , char **){
	a=1;
	{ofstream file;
		file.open("output.txt");
		if(file.is_open()){
			for(double x=-2; x<=6; x+=0.01)
				file<<x<<" " <<
					add<
						iconst<-1>,
						add<
							func3<Gaussian,arg,par<0>,par<1>>,
							func3<Gaussian,arg,par<2>,var<a>>
						>
					>(x,params)
				<<"\n";
			file.close();
		}
	}
	{ofstream file;
		file.open(".plotscript.gp");
		if(file.is_open()){
			file << "plot ";
			file <<"\"output.txt\" w l title \"func\"";
			file << "\n";
			file<<"\npause -1";
			file.close();
		}
		system("gnuplot .plotscript.gp");
	}
	return 0;
}
