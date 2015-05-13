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
numt RandomGauss(numt sigma, numt average=0){
	if(sigma<0)throw std::exception();
	if(sigma==0)return average;
	numt phi=RandomUniformlyR<numt>(0,2.0*3.1415926);
	numt r=0;while(r==0)r=RandomUniformlyR<numt>(0,1);
	numt Z=sigma*cos(phi)*sqrt(-2.0*log(r));
	if(average==0)return Z;
	return average+Z;
}
template<class numt>
class RandomValueGenerator{
private:
	LinearInterpolation_fixedsize<numt,numt> distrib;
public:
	RandomValueGenerator(const RandomValueGenerator &R):distrib(R.distrib){}
	RandomValueGenerator(std::function<numt(numt)> distribution_density,numt x1, numt x2, int bins):distrib(x1,x2,bins){
		using namespace std;
		if(bins<=1)throw exception();
		if(x2<=x1)throw exception();
		numt X[bins];
		for(int i=0;i<bins;i++)
			X[i]=distrib.getX(i);
		numt* Y=SympsonTable<numt,numt*>(distribution_density,X,bins);
		for(int i=0;i<bins;i++){
			distrib.setX(i,Y[i]);
			distrib.setY(i,X[i]);
		}
		delete[] Y;
		for(int i=1;i<bins;i++){
			if(distrib.getX(i)<=distrib.getX(i-1))
				throw exception();
		}
	}
	virtual ~RandomValueGenerator(){}
	numt operator ()(){
		return distrib(RandomUniformlyR<numt>(distrib.min(),distrib.max()));
	}
	numt operator ()(double){
		return operator()();
	}
};
#endif
