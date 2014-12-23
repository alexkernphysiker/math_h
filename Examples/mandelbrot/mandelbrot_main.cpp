#include "mandelbrotwindow.h"
#include <QApplication>

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);
	MandelbrotWindow w;
	w.show();
	
	return a.exec();
}
