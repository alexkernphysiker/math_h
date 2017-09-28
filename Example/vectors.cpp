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
    const double M1=0.5,M2=1.0,E=0.6;
    //Plots
    PlotDistr2D<> 
    Th_vs_E1("First body",BinsByCount(45,0.0,180.),BinsByCount(40,0.0,0.2)),
    Th_vs_E2("First body",BinsByCount(45,0.0,180.),BinsByCount(40,0.0,0.2)),
    E_vs_E("E-E",BinsByCount(40,0.0,0.2),BinsByCount(40,0.0,0.2)),
    Th_vs_Th("Th-Th",BinsByCount(45,0.0,180.),BinsByCount(45,0.0,180.));
    for(size_t i=0;i<1000;i++){
	const auto Pr0=lorentz_byEM(E,M1,direction(Z<>()));
	const auto Pt0=lorentz_byPM(Zero<>(),M2);
	const auto Total=Pr0+Pt0;
	const double invariant_mass=Total.M();
	//let them scatter isotropically in CM
	const auto finalCM=binaryDecay(invariant_mass,M1,M2,randomIsotropic<3>(RG));
	const auto final1=finalCM.first.Transform(-Total.Beta());
	const auto final2=finalCM.second.Transform(-Total.Beta());
	//gather kinematic information
	const auto
	E1=final1.T()-final1.M(),
	E2=final2.T()-final2.M(),
	Th1=direction(final1.S()).th<>()*180./PI(),
	Th2=direction(final2.S()).th<>()*180./PI();
	Th_vs_E1.Fill(Th1,E1);
	Th_vs_E2.Fill(Th2,E2);
	Th_vs_Th.Fill(Th1,Th2);
	E_vs_E.Fill(E1,E2);
    }
}

