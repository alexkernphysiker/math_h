#include <fstream>
#include <stdio.h>
#include <functions.h>
#include <sympson.h>
using namespace std;
int main(int , char **){
	auto E=[](double x){return 0.5*exp(-2.0*x);};
	auto G=[](double x){return Gaussian(x,0.5,0.2);};
	Convolution<double> conv(E,G,0,10,0.001);
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

