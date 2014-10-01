// https://github.com/alexkernphysiker/math_h
#ifndef RANDOMFUNC_H
#	define RANDOMFUNC_H
#include "interpolate.h"
#include "sympson.h"
#include <random>
template<class numt>
double RandomUniformly(numt x1, numt x2){
#ifdef USE_RANDOM_DEVICE
	std::random_device random;
	auto val=random()-random.min();
	auto max=random.max()-random.min();
#else
	auto val=rand();
	auto max=RAND_MAX;
#endif
	return ((x2-x1)*numt(val) / numt(max))+x1;
}
template<class numt>// cannot use integer types
numt RandomGauss(numt sigma, numt average=0, unsigned int precision=12){
	if(sigma==0)return average;
	numt res=0.0;
	numt coeff=1.0/::sqrt(12.0);
	for(unsigned int i=0;i<precision;i++)
		res+=RandomUniformly<numt>(-0.5,0.5);
	res*=(sigma/(coeff*precision))*::sqrt(numt(precision));
	res+=average;
	return res;
}
template<class numt,class func>//Generates random values distributed by given formula
// func cannot be a lambda-expression
class RandomValueGenerator{
private:
	func m_distr;int N;
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
					RandomUniformly<numt>(0,distrib_func[N-1])
				);
	}
	numt operator ()(double){return (*this)();}
	numt operator |(numt x){return m_distr(x)/distrib_func[N-1];}
};
#endif // RANDOMFUNC_H
