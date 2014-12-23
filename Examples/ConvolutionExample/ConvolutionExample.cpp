#include <QFile>
#include <QTextStream>
#include <QProcess>
#include <stdio.h>
#include <math.h>
#include <sympson.h>
#include <functions.h>
#include <singleparam.h>
int main(int , char **){
	SingleParam<double,2,double,double,double>
			E(&KExpLX<double>,0.5,-2,INFINITY);
	SingleParam<double,0,double,double,double>
			G(&Gaussian<double>,INFINITY,0.5,0.2);
	Convolution<double,decltype(E),decltype(G)> conv(E,G);
	conv.Init(0,10,0.001);
	{QFile file("output.txt");
		file.open(QFile::WriteOnly);
		if(file.isOpen()){
			QTextStream str(&file);
			for(double x=0; x<=6; x+=0.01)
				str<<x<<" " << conv(x) <<"\n";
			file.close();
		}
	}
	{QString script=".plotscript.gp";
		QFile file(script);
		file.open(QFile::WriteOnly);
		if(file.isOpen()){
			QTextStream str(&file);
			str << "plot ";
			str <<"\"output.txt\" w l title \"func\"";
			str << "\n";
			str<<"\npause -1";
			file.close();
		}
		QProcess *gnuplot=new QProcess();
		gnuplot->startDetached("gnuplot",QStringList()<<script);
	}
}

