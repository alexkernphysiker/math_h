#include <math.h>
#include <list>
#include "randcmp.h"
#include <randomfunc.h>
double getPseudo(double X, double sigma){
	return RandomGauss<double>(sigma,X);
}
