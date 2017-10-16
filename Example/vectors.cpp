// this file is distributed under
// MIT license
#include <iostream>
#include <math_h/tabledata.h>
#include <math_h/vectors.h>
#include <gnuplot_wrap.h>
using namespace std;
using namespace MathTemplates;
using namespace GnuplotWrap;
//This is an example how to use vectors.h
//that provides a simple Monte Carlo simulation
//of particle scattering
int main()
{
    RANDOM RG;
    //let's assume two bodies scattering
    const double M1=0.5,M2=1.0,P=0.3;
    //Plots
    PlotDistr2D<> 
    P_vs_P("P-P",BinsByCount(50,0.0,0.5),BinsByCount(50,0.0,0.5),"plot_PP"),
    Th_vs_Th("Th-Th",BinsByCount(45,0.0,180.),BinsByCount(45,0.0,180.),"plot_ThTh");
    //Simulation cycle
    for(size_t i=0;i<1000;i++){
	const auto Pr0=lorentz_byPM(Z()*P,M1);
	const auto Pt0=lorentz_byPM(Zero(),M2);
	const auto Total=Pr0+Pt0;
	//let them scatter isotropically in CM
	const auto finalCM=binaryDecay(Total.M(),M1,M2,randomIsotropic<3>(RG));
	//Converting to lab system. We put the sign minus,
	//because Total describes how CM is moving in lab
	//but we are converting from CM to lab
	const auto final1=finalCM.first.Transform(-Total.Beta());
	const auto final2=finalCM.second.Transform(-Total.Beta());
	//gather kinematic information
	Th_vs_Th.Fill(direction(final1.S()).th()*180./PI(),direction(final2.S()).th()*180./PI());
	P_vs_P.Fill(final1.S().M(),final2.S().M());
    }
}
