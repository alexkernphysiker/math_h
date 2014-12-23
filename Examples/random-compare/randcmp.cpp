#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "randcmp.h"
#include <sigma.h>
int main(int , char **){
	Sigma<double> *pseudo=new Sigma<double>();
	Sigma<double> *device=new Sigma<double>();
	for(int i=0; i<10000;i++){
		pseudo->AddValue(getPseudo(5,2));
		device->AddValue(getDevice(5,2));
		if((i+1)%100==0){
			printf("i=%i\n",i+1);
			printf("average: %f, %f\n",pseudo->getAverage(),device->getAverage());
			printf("Sigma: %f, %f\n",pseudo->getSigma(),device->getSigma());
		}
	}
	delete pseudo;
	delete device;
	return 0;
}
