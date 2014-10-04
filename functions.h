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
template<class numt, class indexer>
numt Polynom(numt x,indexer  p,unsigned int P, int index_offset=0){
	numt res=0;numt c=1;
	for(unsigned int i=0; i<=P;i++){res+=c*p[index_offset+i];if(i<P)c*=x;}
	return res;
}
template<unsigned int P,class numt, class indexer, int index_offset=0>
numt Polynom(numt x,indexer  p){return Polynom<numt,indexer>(x,p,P,index_offset);}
namespace FuncWrappers{
template<class numt,class indexer,int i>
numt par(indexer P){return P[i];}
template<class numt,class indexer,numt (f)(numt),numt(F)(indexer)>
numt func(indexer P){return f(F(P));}
template<class numt,class indexer,numt (f)(numt,numt),numt(F1)(indexer),numt(F2)(indexer)>
numt func2(indexer P){return f(F1(P),F2(P));}
template<class numt,class indexer,numt (f)(numt,numt,numt),numt(F1)(indexer),numt(F2)(indexer),numt(F3)(indexer)>
numt func3(indexer P){return f(F1(P),F2(P),F3(P));}
template<class numt,class indexer,numt(F1)(indexer),numt(F2)(indexer)>
numt add(indexer P){return F1(P)+F2(P);}
template<class numt,class indexer,numt(F1)(indexer),numt(F2)(indexer)>
numt sub(indexer P){return F1(P)-F2(P);}
template<class numt,class indexer,numt(F1)(indexer),numt(F2)(indexer)>
numt mul(indexer P){return F1(P)*F2(P);}
template<class numt,class indexer,numt(F1)(indexer),numt(F2)(indexer)>
numt div(indexer P){return F1(P)/F2(P);}
template<class numt,class indexer,numt(F1)(indexer),numt(F2)(indexer)>
numt power(indexer P){return pow(F1(P),F2(P));}
}
#endif // FUNCTIONS_H
