// https://github.com/alexkernphysiker/math_h
#ifndef RANDOMFUNC_H
#	define RANDOMFUNC_H
#include "interpolate.h"
#include "sympson.h"
#include <random>
#include <vector>
template<class numt>
numt RandomUniformly(numt x1, numt x2){
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
	func m_distr;
	std::vector<numt> values,distrib_func;
public:
	RandomValueGenerator(){}
	RandomValueGenerator(RandomValueGenerator &R){
		m_distr=R.m_distr;
		for(auto v:R.values)
			values.push_back(v);
		for(auto v:R.distrib_func)
			distrib_func.push_back(v);
	}
	RandomValueGenerator(func distribution_density,numt x1, numt x2, numt step){
		m_distr=distribution_density;
		for(numt x=x1;x<=x2;x+=step)
			values.push_back(x);
		auto N=values.size();
		numt* tbl=SympsonTable<numt,std::vector<numt>,func>(N,values,m_distr);
		for(int i=0;i<N;i++)
			distrib_func.push_back(tbl[i]);
	}
	virtual ~RandomValueGenerator(){}
	numt operator ()(){
		auto N=values.size();
		return Interpolate_Linear<numt,std::vector<numt>>(0,N-1,
					distrib_func,values,
					RandomUniformly<numt>(0,distrib_func[N-1])
				);
	}
	numt operator ()(double){return (*this)();}
	numt operator |(numt x){
		auto N=values.size();
		return m_distr(x)/distrib_func[N-1];
	}
};
#endif // RANDOMFUNC_H
