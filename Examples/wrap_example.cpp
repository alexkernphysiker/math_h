#include <QFile>
#include <QTextStream>
#include <QProcess>
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
double a;
int main(int , char **){
	a=1;
	{QFile file("output.txt");
		file.open(QFile::WriteOnly);
		if(file.isOpen()){
			QTextStream str(&file);
			for(double x=-2; x<=6; x+=0.01)
				str<<x<<" " <<
					 add<
						func3<Gaussian,arg,par<0>,par<1>>,
						func3<Gaussian,arg,par<2>,var<a>>
					>(x,params)
				<<"\n";
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
	return 0;
}
