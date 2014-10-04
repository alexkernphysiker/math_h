// https://github.com/alexkernphysiker/math_h
#ifndef FUNCTIONS_H
#	define FUNCTIONS_H
#include <math.h>
template<class numt=double>double KExpLX(double k, double l, double x){return k*exp(l*x);}
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
namespace FuncWrappers2{
	template<class numt,class indexer,int i>
	numt arg(indexer X,indexer){return X[i];}
	template<class numt,class indexer,int i>
	numt par(indexer,indexer P){return P[i];}
	template<class numt,class indexer,numt (f)(numt),numt(F)(indexer,indexer)>
	numt func(indexer X,indexer P){return f(F(X,P));}
	template<class numt,class indexer,numt (f)(numt,numt),numt(F1)(indexer,indexer),numt(F2)(indexer,indexer)>
	numt func2(indexer X,indexer P){return f(F1(X,P),F2(X,P));}
	template<class numt,class indexer,numt (f)(numt,numt,numt),numt(F1)(indexer,indexer),numt(F2)(indexer,indexer),numt(F3)(indexer,indexer)>
	numt func3(indexer X,indexer P){return f(F1(X,P),F2(X,P),F3(X,P));}
	template<class numt,class indexer,numt(F1)(indexer,indexer),numt(F2)(indexer,indexer)>
	numt add(indexer X,indexer P){return F1(X,P)+F2(X,P);}
	template<class numt,class indexer,numt(F1)(indexer,indexer),numt(F2)(indexer,indexer)>
	numt sub(indexer X,indexer P){return F1(X,P)-F2(X,P);}
	template<class numt,class indexer,numt(F1)(indexer,indexer),numt(F2)(indexer,indexer)>
	numt mul(indexer X,indexer P){return F1(X,P)*F2(X,P);}
	template<class numt,class indexer,numt(F1)(indexer,indexer),numt(F2)(indexer,indexer)>
	numt div(indexer X,indexer P){return F1(X,P)/F2(X,P);}
	template<class numt,class indexer,numt(F1)(indexer,indexer),numt(F2)(indexer,indexer)>
	numt power(indexer X,indexer P){return pow(F1(X,P),F2(X,P));}
}
#endif // FUNCTIONS_H
