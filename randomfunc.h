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
	LinearInterpolation_fixedsize<numt,numt> distrib;
	double C;
public:
	RandomValueGenerator():C(1),distrib(0,1,2){}
	RandomValueGenerator(RandomValueGenerator &R):C(R.C),distrib(R.distrib.min(),R.distrib.max(),R.distrib.size()){
		for(int i=0,n=distrib.size();i<n;i++)
			distrib.setY(i,R.distrib.getY(i));
	}
	RandomValueGenerator(std::function<numt(numt)> distribution_density,numt x1, numt x2, int bins):distrib(x1,x2,bins){
		using namespace std;
		if(bins<=1)throw exception();
		if(x2<=x1)throw exception();
		numt X[bins];
		for(int i=0;i<bins;i++)
			X[i]=distrib.getX(i);
		numt* Y=SympsonTable<numt,numt*>(distribution_density,X,bins);
		for(int i=0;i<bins;i++)
			distrib.setY(i,Y[i]);
		C=Y[bins-1];
		delete[] Y;
		if(C<=0)throw exception();
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
