#include "mandelbrotwindow.h"
#include "ui_mandelbrotwindow.h"
#include <QPainter>
#include <math.h>
#include <complex.h>
typedef Complex<double> complex;
//data
int data_array[600][800];
//scale
int x_re_0=400;
int y_im_0=300;
int pix_per_unit=250;
//numeric calculation parameters
uint n_iter=100;
double r_max=100000;
double K=30;
//calculation
complex from_coords(int x,int y){
	return complex(double(x-x_re_0)/double(pix_per_unit),double(y-y_im_0)/double(pix_per_unit));
}
double check(const complex c){
	complex z=0;
	for(uint i=0; i<n_iter; i++)
		 z=pow(z,2)+c;
	return abs(z);
}
void init_data_array(){
	for(uint y = 0; y<600; y++){
		for(uint x=0; x<800; x++){
			data_array[y][x]=int(log(check(from_coords(x,y))/r_max)*K);
		}
	}
}
//UI
MandelbrotWindow::MandelbrotWindow(QWidget *parent) :
	QMainWindow(parent),
	ui(new Ui::MandelbrotWindow)
{
	init_data_array();
	ui->setupUi(this);
}
MandelbrotWindow::~MandelbrotWindow(){	delete ui;}
void MandelbrotWindow::paintEvent(QPaintEvent *){
	QPainter painter(this);
	for(uint y = 0; y<600; y++)
		for(uint x=0; x<800; x++){
			QColor color;
			int v=data_array[y][x];
			if(v>0)
				color.setRgba(qRgba(v/2,128-v/2,0,255));
			else
				color.setRgba(qRgba(0,-v,-v,255));
			QRect rect(x,y,1,1);
			painter.fillRect(rect,color);
		}
}
