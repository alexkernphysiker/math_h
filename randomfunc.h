/* The newest version of this file can be got by this link
https://github.com/alexkernphysiker/MathLibs/blob/master/functions
There are also some examples how to use these templates in this repository
author: alex_kernphysiker@privatdemail.net */
#ifndef RANDOMFUNC_H
#define RANDOMFUNC_H
#include "interpolate.h"
#include "sympson.h"
#include <random>
template<class numt,class func>//Generates random values distributed by given formula
// func cannot be a lambda-expression
class RandomValueGenerator{
private:
	func m_distr;	int N;
	numt *values;numt *distrib_func;
public:
	RandomValueGenerator(){}
	RandomValueGenerator(RandomValueGenerator &R){
		m_distr=R.m_distr;N=R.N;
		if(N>0){
			values=new numt[N];distrib_func=new numt[N];
			for(int i=0; i<N;i++){
				values[i]=R.values[i];distrib_func[i]=R.distrib_func[i];
			}
		}
	}
	RandomValueGenerator(func distribution_density){
		m_distr=distribution_density;N=0;values=nullptr;distrib_func=nullptr;
	}
	virtual ~RandomValueGenerator(){
		if(values!=nullptr){delete[] values;delete[] distrib_func;}
	}
	void Init(numt x1, numt x2, numt step){
		if(step<=0)throw; //wrong parameter value
		if(values!=nullptr){delete[] values;delete[] distrib_func;}
		N=int((x2-x1)/step)+1;
		if(N<0)N=-N;
		values=new numt[N];
		for(int i=0; i<N; i++)values[i]=x1+(step*i);
		distrib_func=SympsonTable<numt,numt*,func>(N,values,m_distr);
	}
	numt operator ()(){
		return Interpolate_Linear(0,N-1,
					distrib_func,values,
					distrib_func[N-1]*numt(rand()) / numt(RAND_MAX)
				);
	}
	numt operator ()(numt){return (*this)();}
};
template<class numt>
double RandomUniformly(numt x1, numt x2){
	return ((x2-x1)*numt(rand()) / numt(RAND_MAX))+x1;
}
template<class numt>// cannot use integer types
numt RandomGauss(numt sigma, numt average=0, unsigned int precision=12){
	if(sigma==0)return average;
	numt res=0.0;
	numt coeff=1.0/::sqrt(12.0);
	for(unsigned int i=0;i<precision;i++)
		res+=RandomUniformly(-0.5,0.5);
	res*=(sigma/(coeff*precision))*::sqrt(numt(precision));
	res+=average;
	return res;
}
#endif // RANDOMFUNC_H
