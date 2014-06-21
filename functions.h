#ifndef FUNCTIONS_H
#	define FUNCTIONS_H
#include <math.h>
template<class numt=double>//cannot use integer types
numt Gaussian(numt x, numt X_max, numt sigma){
	numt koef= 1/(sigma*::sqrt(2*3.1415926));
	return koef*exp(-(pow(X_max-x,2))/(2*sigma*sigma));
}
template<class numt=double, class indexer=numt*>
numt Gauss(indexer x, indexer p){return Gaussian(x[0],p[0],p[1])*p[2];}
template<class numt=double, class indexer=numt*>
numt gauss(numt x, indexer p){return Gaussian(x,p[0],p[1])*p[2];}

template<class numt=double>//cannot use integer types
numt FermiFunc(numt x, numt X_border, numt diffuse){
	return 1.0/(1.0+exp((x-X_border)/diffuse));
}
template<class numt=double, class indexer=numt*>
numt Fermi(indexer x, indexer p){return FermiFunc(x[0],p[0],p[1])*p[2];}
template<class numt=double, class indexer=numt*>
numt fermi(numt x, indexer p){return FermiFunc(x,p[0],p[1])*p[2];}

template<class numt=double>//cannot use integer types
numt KExpLX(numt K, numt l, numt x){return K*exp(l*x);}
template<class numt=double, class indexer=numt*>
numt Exp(indexer x, indexer p){return KExpLX(p[0],p[1],x[0]);}
template<class numt=double, class indexer=numt*>
numt _exp(numt x, indexer p){return KExpLX(p[0],p[1],x);}

template<class numt=double, class indexer=numt*>
numt polynom(numt x,indexer  p,unsigned int P){
	numt res=0;numt c=1;
	for(unsigned int i=0; i<=P;i++){res+=c*p[i];if(i<P)c*=x;}
	return res;
}
template<unsigned int P,class numt=double, class indexer=numt*>
numt polynom(numt x,indexer  p){return polynom<numt,indexer>(x,p,P);}
template<unsigned int P,class numt=double, class indexer=numt*>
numt Polynom(indexer x, indexer p){return polynom<P,numt,indexer>(x[0],p);}

#endif // FUNCTIONS_H