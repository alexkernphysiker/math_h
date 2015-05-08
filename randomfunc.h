#ifndef XRYPRAVJWTJCYPQI
#define XRYPRAVJWTJCYPQI
#ifdef USE_RANDOM_DEVICE
# include <random>
#endif
#include <vector>
#include <functional>
#include <math.h>
#include "interpolate.h"
#include "sympson.h"
template<class inumt>
inumt RandomUniformlyI(inumt x1,inumt x2){
	if(x1>x2)throw std::exception();
	if(x1==x2)return x1;
	#ifdef USE_RANDOM_DEVICE
	std::random_device random;
	inumt val=random();
	if(val<0)
		val=-val;
	#else
	inumt val=rand();
	#endif
	return val%(x2+1-x1)+x1;
}
template<class numt>
numt RandomUniformlyR(numt x1, numt x2){
	if(x1>x2)throw std::exception();
	if(x1==x2)return x1;
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
template<class numt>
numt RandomGauss(numt sigma, numt average=0, unsigned int precision=12){
	if(sigma<0)throw std::exception();
	if(sigma==0)return average;
	numt res=0.0;
	numt coeff=1.0/sqrt(12.0);
	for(unsigned int i=0;i<precision;i++)
		res+=RandomUniformlyR<numt>(-0.5,0.5);
	res*=(sigma/(coeff*precision))*sqrt(numt(precision));
	if(0==average)
		return res;
	res+=average;
	return res;
}
template<class numt>
class RandomValueGenerator{
private:
	LinearInterpolation<numt> distrib;
	double C;
public:
	RandomValueGenerator():C(0){}
	RandomValueGenerator(RandomValueGenerator &R):C(R.C){
		for(auto p:R.distrib)
			distrib<<p;
	}
	RandomValueGenerator(std::function<numt(numt)> distribution_density,numt x1, numt x2, numt step){
		using namespace std;
		vector<numt> X;
		for(numt x=x1;x<=x2;x+=step)
			X.push_back(x);
		numt* Y=SympsonTable<numt,vector<numt>>(distribution_density,X,X.size());
		for(int i=0;i<X.size();i++)
			distrib<<make_pair(Y[i],X[i]);
		C=Y[X.size()-1];
		delete[] Y;
	}
	virtual ~RandomValueGenerator(){}
	numt operator ()(){
		return distrib(RandomUniformlyR<numt>(0,C));
	}
	numt operator ()(double){
		return operator()();
	}
};
#endif
