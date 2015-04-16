#ifndef EMEFWAYNIGJGENCP
#define EMEFWAYNIGJGENCP
#include <math.h>

template<class numt=double>
numt Gaussian(numt x, numt X_max, numt sigma){
	numt koef= 1/(sigma*sqrt(2*3.1415926));
	return koef*exp(-(pow(X_max-x,2))/(2*sigma*sigma));
}
template<class numt=double>
numt BreitWigner(numt x, numt ampl, numt pos, numt gamma){
	return (ampl*pow(gamma/2,2))/(pow(x-pos,2)+pow(gamma/2,2));
}
template<class numt=double>
numt BukinFunction(numt x,numt ampl,numt pos,numt sigma,numt assymetryparam,numt rho1,numt rho2){
	numt sq2ln2=sqrt(2.0*log(2));
	numt x1=pos+sigma*sq2ln2*(assymetryparam/sqrt(assymetryparam+1) -1);
	numt x2=pos+sigma*sq2ln2*(assymetryparam/sqrt(assymetryparam+1) +1);
	numt in_exp=assymetryparam*sqrt(assymetryparam*assymetryparam+1)*(x-x1)*sq2ln2;
	in_exp/=sigma*pow(sqrt(assymetryparam*assymetryparam+1)-assymetryparam,2)*log(sqrt(assymetryparam*assymetryparam+1)+assymetryparam);
	in_exp-=log(2);
	if(x<x1)
		in_exp+=rho1*pow((x-x1)/(pos-x1),2);
	if(x>=x2)
		in_exp+=rho2*pow((x-x2)/(pos-x2),2);
	return ampl*exp(in_exp);
}
template<class numt=double>
numt Novosibirsk(numt x,numt pos,numt sigma,numt asym){
	numt sqlog4=sqrt(log(4));
	numt lambda = sinh(sigma*sqlog4)*pow(asym*sqlog4,-1);
	if(1+lambda*(x-pos) >0) 
		return exp(-0.5*pow(log(1+lambda*(x-pos)),2)*pow(sigma,-2)+pow(sigma,2));
	return 0;
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
