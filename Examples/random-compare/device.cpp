#include <math.h>
#include "randcmp.h"
#define USE_RANDOM_DEVICE
#include <randomfunc.h>
double getDevice(double X, double sigma){
	return RandomGauss<double>(sigma,X);
}
