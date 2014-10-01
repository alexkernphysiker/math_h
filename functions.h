// https://github.com/alexkernphysiker/math_h
#ifndef FUNCTIONS_H
#	define FUNCTIONS_H
#include <math.h>
/////// FUNCTIONS
template<class numt=double>//cannot use integer types
numt Gaussian(numt x, numt X_max, numt sigma){
	numt koef= 1/(sigma*::sqrt(2*3.1415926));
	return koef*exp(-(pow(X_max-x,2))/(2*sigma*sigma));
}
template<class numt=double>//cannot use integer types
numt FermiFunc(numt x, numt X_border, numt diffuse){
	return 1.0/(1.0+exp((x-X_border)/diffuse));
}
template<class numt=double>//cannot use integer types
numt KExpLX(numt K, numt l, numt x){return K*exp(l*x);}
template<class numt=double, class indexer=numt*>
numt polynom(numt x,indexer  p,unsigned int P, int index_offset=0){
	numt res=0;numt c=1;
	for(unsigned int i=0; i<=P;i++){res+=c*p[index_offset+i];if(i<P)c*=x;}
	return res;
}
////////// WRAPPERS
template<class numt=double, class indexer=numt*, int offset_x=0, int offset_p=0>
numt Gauss(indexer x, indexer p){return Gaussian(x[offset_x],p[offset_p+0],p[offset_p+1])*p[offset_p+2];}
template<class numt=double, class indexer=numt*>
numt Gauss(indexer x, indexer p, int offset_x=0, int offset_p=0){return Gaussian(x[offset_x],p[offset_p+0],p[offset_p+1])*p[offset_p+2];}
template<class numt=double, class indexer=numt*,int offset_p=0>
numt gauss(numt x, indexer p){return Gaussian(x,p[offset_p+0],p[offset_p+1])*p[offset_p+2];}
template<class numt=double, class indexer=numt*>
numt gauss(numt x, indexer p,int offset_p=0){return Gaussian(x,p[offset_p+0],p[offset_p+1])*p[offset_p+2];}
template<class numt=double, class indexer=numt*, int offset_x=0, int offset_p=0>
numt Fermi(indexer x, indexer p){return FermiFunc(x[offset_x],p[offset_p],p[offset_p+1])*p[offset_p+2];}
template<class numt=double, class indexer=numt*>
numt Fermi(indexer x, indexer p, int offset_x=0, int offset_p=0){return FermiFunc(x[offset_x],p[offset_p],p[offset_p+1])*p[offset_p+2];}
template<class numt=double, class indexer=numt*,int offset_p=0>
numt fermi(numt x, indexer p){return FermiFunc(x,p[offset_p],p[offset_p+1])*p[offset_p+2];}
template<class numt=double, class indexer=numt*>
numt fermi(numt x, indexer p,int offset_p=0){return FermiFunc(x,p[offset_p],p[offset_p+1])*p[offset_p+2];}
template<class numt=double, class indexer=numt*, int offset_x=0, int offset_p=0>
numt Exp(indexer x, indexer p){return KExpLX(p[offset_p],p[1+offset_p],x[offset_x]);}
template<class numt=double, class indexer=numt*>
numt Exp(indexer x, indexer p, int offset_x=0, int offset_p=0){return KExpLX(p[offset_p],p[1+offset_p],x[offset_x]);}
template<class numt=double, class indexer=numt*,int offset_p=0>
numt _exp(numt x, indexer p){return KExpLX(p[offset_p],p[1+offset_p],x);}
template<class numt=double, class indexer=numt*>
numt _exp(numt x, indexer p,int offset_p=0){return KExpLX(p[offset_p],p[1+offset_p],x);}
template<unsigned int P,class numt=double, class indexer=numt*, int index_offset=0>
numt polynom(numt x,indexer  p){return polynom<numt,indexer>(x,p,P,index_offset);}
template<unsigned int P,class numt=double, class indexer=numt*, int offset_x=0, int offset_p=0>
numt Polynom(indexer x, indexer p){return polynom<P,numt,indexer,offset_p>(x[offset_x],p);}
template<unsigned int P,class numt=double, class indexer=numt*>
numt Polynom(indexer x, indexer p, int offset_x=0, int offset_p=0){return polynom<P,numt,indexer,offset_p>(x[offset_x],p);}
#endif // FUNCTIONS_H
