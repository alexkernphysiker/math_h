#ifndef EMEFWAYNIGJGENCP
#define EMEFWAYNIGJGENCP
#include <math.h>
template<class numt=double>
double KExpLX(double k, double l, double x){
	return k*exp(l*x);
}
template<class numt=double>
numt Gaussian(numt x, numt X_max, numt sigma){
	numt koef= 1/(sigma*sqrt(2*3.1415926));
	return koef*exp(-(pow(X_max-x,2))/(2*sigma*sigma));
}
template<class numt=double>
numt BreitWigner(numt x, numt ampl, numt gamma, numt pos){
	return (ampl*pow(gamma/2,2))/(pow(x-pos,2)+pow(gamma/2,2));
}
template<class numt=double>
numt FermiFunc(numt x, numt X_border, numt diffuse){
	return 1.0/(1.0+exp((x-X_border)/diffuse));
}
template<class numt, class indexer>
numt Polynom(numt x,indexer  p,unsigned int P, int index_offset=0){
	numt res=0;numt c=1;
	for(unsigned int i=0; i<=P;i++){res+=c*p[index_offset+i];if(i<P)c*=x;}
	return res;
}
template<unsigned int P,class numt, class indexer, int index_offset=0>
numt Polynom(numt x,indexer  p){
	return Polynom<numt,indexer>(x,p,P,index_offset);
}
#endif
