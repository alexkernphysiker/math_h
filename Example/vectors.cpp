// this file is distributed under
// MIT license
#include <iostream>
#include <math_h/tabledata.h>
#include <math_h/vectors.h>
#include <gnuplot_wrap.h>
using namespace std;
using namespace MathTemplates;
using namespace GnuplotWrap;
int main()
{
    RANDOM RG;
    //let's assume two bodies scattering
    const double M1=0.5,M2=1.0,E=0.3;
    //Plots
    PlotDistr2D<> 
    E_vs_E("E-E",BinsByCount(60,0.0,0.3),BinsByCount(60,0.0,0.3),"plot_EE"),
    P_vs_P("P-P",BinsByCount(50,0.0,1.0),BinsByCount(50,0.0,1.0),"plot_PP"),
    Th_vs_Th("Th-Th",BinsByCount(45,0.0,180.),BinsByCount(45,0.0,180.),"plot_ThTh");
    for(size_t i=0;i<1000;i++){
	const auto Pr0=lorentz_byEM(E+M1,M1,direction(Z<>()));
	const auto Pt0=lorentz_byPM(Zero<>(),M2);
	const auto Total=Pr0+Pt0;
	const double invariant_mass=Total.M();
	//let them scatter isotropically in CM
	const auto finalCM=binaryDecay(invariant_mass,M1,M2,randomIsotropic<3>(RG));
	//minus, because Total describes how CM is moving in lab
	//but we are converting from CM to lab
	const auto final1=finalCM.first.Transform(-Total.Beta());
	const auto final2=finalCM.second.Transform(-Total.Beta());
	//gather kinematic information
	Th_vs_Th.Fill(direction(final1.S()).th<>()*180./PI(),direction(final2.S()).th<>()*180./PI());
	E_vs_E.Fill(final1.T()-final1.M(),final2.T()-final2.M());
	P_vs_P.Fill(final1.S().M(),final2.S().M());
    }
}

