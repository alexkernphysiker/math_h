//This is an example how to use the math_h libary
#include <iostream>
#include <math_h/lorentzvector.h>
#include <gnuplot_wrap.h>
using namespace std;
using namespace MathTemplates;
using namespace GnuplotWrap;
//This is an example how to use lorentzvector.h that provides
//a simple Monte Carlo simulation of particle scattering
//The vectors must be given in unit system that assumes c=1
//and equal units for energy, momentum and mass
int main()
{
    //Particles masses
    const double M1 = 0.75, M2 = 1.0;
    //momentum distribution
    const RandomGauss<> P(0.3,0.02);

    PlotDistr2D<> //distributions that will be plotted
    P_vs_P("P-P", BinsByCount(100, 0.0, 0.5), BinsByCount(100, 0.0, 0.5), "plot_PP",3),
    E_vs_E("Ek-Ek", BinsByCount(100, 0.0, 0.1), BinsByCount(100, 0.0, 0.1), "plot_EE",3),
    Th_vs_Th("Th-Th", BinsByCount(90, 0.0, 180.), BinsByCount(90, 0.0, 90.), "plot_ThTh",3);


    for (size_t i = 0; i < 100000; i++) {
	//Create Lorentz vectors for input particles
        const auto Pr0 = lorentz_byPM(P()*Z(), M1);
        const auto Pt0 = lorentz_Rest(M2);
	// calculate total 4-momentum that is conserved
        const auto Total = Pr0 + Pt0;
        //Let the particles scatter isotropically in CM frame
        const auto finalCM = binaryDecay(Total.M(), M1, M2, randomIsotropic<3>());
        //Converting to lab system.
        const auto final1 = finalCM.first.Transform(-Total.Beta());
        const auto final2 = finalCM.second.Transform(-Total.Beta());
        //Gather kinematic statistics
        Th_vs_Th.Fill(direction(final1.P()).th() * 180. / PI(), direction(final2.P()).th() * 180. / PI());
        P_vs_P.Fill(final1.P().length(), final2.P().length());
        E_vs_E.Fill(final1.Ekin(), final2.Ekin());
    }
}

