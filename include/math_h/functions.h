// this file is distributed under 
// MIT license
#ifndef EMEFWAYNIGJGENCP
#define EMEFWAYNIGJGENCP
#include <math.h>
#include <functional>
namespace MathTemplates{
	template<class numX, class numY=numX>
	class IFunction{
	public:
		virtual numY operator()(const numX x)const=0;
		virtual const std::function<numY(numX)> func()const=0;		
	};
	
	
	template<class numt=double>
	numt Gaussian(const numt x,const numt X_max,const numt sigma){
		return exp(-pow((X_max-x)/sigma,2)/2)/(sigma*sqrt(2*3.1415926));
	}
	template<class numt=double>
	numt BreitWigner(const numt x,const numt pos,const numt gamma){
		return gamma/(2*3.1415926*(pow(x-pos,2)+pow(gamma/2,2)));
	}
	template<class numt=double>
	numt Novosibirsk(const numt x,const numt pos,const numt sigma,const numt asym){
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
	numt FermiFunc(const numt x,const numt X_border,const numt diffuse){
		return 1.0/(1.0+exp((x-X_border)/diffuse));
	}
	template<class numt, class indexer>
	numt Polynom(const numt x,const indexer  p,const unsigned int P,const int index_offset=0){
		numt res=0;numt c=1;
		for(unsigned int i=0; i<=P;i++){res+=c*p[index_offset+i];if(i<P)c*=x;}
		return res;
	}
	template<unsigned int P,class numt, class indexer, int index_offset=0>
	inline numt Polynom(const numt x,const indexer  p){
		return Polynom<numt,indexer>(x,p,P,index_offset);
	}
};
#endif
