#ifndef XRYPRAVJWTJCYPQI
#define XRYPRAVJWTJCYPQI
#include <random>
#include <functional>
#include <math.h>
#include "interpolate.h"
#include "sympson.h"
template<class numt,class RG=std::default_random_engine,typename... Args>
class RandomValueGenerator:protected LinearInterpolation_fixedsize<numt,numt>{
private:
	std::uniform_real_distribution<numt> distr;
	RG generator;
public:
	RandomValueGenerator(const RandomValueGenerator &R)
		:LinearInterpolation_fixedsize<numt,numt>(R),distr(R.distr),generator(R.generator){}
	RandomValueGenerator(std::function<numt(numt)> distribution_density,numt x1, numt x2, int bins,Args... args)
		:LinearInterpolation_fixedsize<numt,numt>(x1,x2,bins),generator(args...){
		using namespace std;
		if(bins<=1)throw exception();
		if(x2<=x1)throw exception();
		numt X[bins];
		for(int i=0;i<bins;i++)
			X[i]=LinearInterpolation_fixedsize<numt,numt>::getX(i);
		numt* Y=SympsonTable<numt,numt*>(distribution_density,X,bins);
		for(int i=0;i<bins;i++){
			LinearInterpolation_fixedsize<numt,numt>::setX(i,Y[i]);
			LinearInterpolation_fixedsize<numt,numt>::setY(i,X[i]);
		}
		delete[] Y;
		for(int i=1;i<bins;i++){
			if(LinearInterpolation_fixedsize<numt,numt>::getX(i)<LinearInterpolation_fixedsize<numt,numt>::getX(i-1))
				throw exception();
		}
		if(LinearInterpolation_fixedsize<numt,numt>::min()>=LinearInterpolation_fixedsize<numt,numt>::max())
			throw exception();
		distr=uniform_real_distribution<numt>(LinearInterpolation_fixedsize<numt,numt>::min(),LinearInterpolation_fixedsize<numt,numt>::max());
	}
	RandomValueGenerator(numt x1, numt x2,Args... args):RandomValueGenerator([](double){return 1.0;},x1,x2,2,args...){}
	virtual ~RandomValueGenerator(){}
	numt operator ()(){
		return LinearInterpolation_fixedsize<numt,numt>::operator()(distr(generator));
	}
	numt operator ()(double){
		return operator()();
	}
};
#endif
