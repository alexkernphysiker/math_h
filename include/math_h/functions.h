// this file is distributed under 
// MIT license
#ifndef EMEFWAYNIGJGENCP
#define EMEFWAYNIGJGENCP
#include <math.h>

template<class numt=double>
numt Gaussian(numt x, numt X_max, numt sigma){
	return exp(-pow((X_max-x)/sigma,2)/2)/(sigma*sqrt(2*3.1415926));
}
template<class numt=double>
numt BreitWigner(numt x, numt pos, numt gamma){
	return gamma/(2*3.1415926*(pow(x-pos,2)+pow(gamma/2,2)));
}
template<class numt=double>
numt Novosibirsk(numt x,numt pos,numt sigma,numt asym){
	if(pow(asym/sigma,2)<=0.000001)
		return Gaussian<numt>(x,pos,sigma);
	numt Lsqlog4=asym*sqrt(log(4));
	numt qb=sinh(Lsqlog4)/Lsqlog4;
	numt q4=qb*asym*(x-pos)/sigma + 1;
	if(q4<=0)return 0;
	numt slq4=pow(log(q4),2);
	numt k=sigma*sqrt(2*3.1415926*(pow(asym,2)+1.0));
	return exp(-slq4/(pow(asym,2)*2)+pow(asym,2))/k;
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
inline numt Polynom(numt x,indexer  p){
	return Polynom<numt,indexer>(x,p,P,index_offset);
}
#endif
