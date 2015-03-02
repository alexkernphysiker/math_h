#include <fstream>
#include <stdio.h>
#include <math.h>
#include <sympson.h>
#include <functions.h>
#include <singleparam.h>
using namespace std;
int main(int , char **){
	SingleParam<double,2,double,double,double>
			E(&KExpLX,0.5,-2,INFINITY);
	SingleParam<double,0,double,double,double>
			G(&Gaussian,INFINITY,0.5,0.2);
	Convolution<double,decltype(E),decltype(G)> conv(E,G);
	conv.Init(0,10,0.001);
	{ofstream file;
		file.open("output.txt");
		if(file.is_open()){
			for(double x=0; x<=6; x+=0.01)
				file<<x<<" " << conv(x) <<"\n";
			file.close();
		}
	}
	{ofstream file;
		file.open(".plotscript.gp");
		if(file.is_open()){
			file << "plot ";
			file <<"\"output.txt\" w l title \"func\"";
			file << "\n";
			file <<"\npause -1";
			file.close();
		}
		system("gnuplot .plotscript.gp");
	}
}

