#ifndef XRYPRAVJWTJCYPQI
#define XRYPRAVJWTJCYPQI
#include <random>
#include <functional>
#include <math.h>
#include "interpolate.h"
#include "sympson.h"
#include "exception_math_h.h"
template<class numt,class RG=std::mt19937>
class RandomValueGenerator:protected LinearInterpolation_fixedsize<numt,numt>{
private:
	std::uniform_real_distribution<numt> distr;
public:
	RandomValueGenerator(const RandomValueGenerator &R)
		:LinearInterpolation_fixedsize<numt,numt>(R),distr(R.distr){}
	RandomValueGenerator(std::function<numt(numt)> distribution_density,numt x1, numt x2, int bins)
		:LinearInterpolation_fixedsize<numt,numt>(x1,x2,bins){
		using namespace std;
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
				throw math_h_error<RandomValueGenerator>("Probability density function has points below zero");
		}
		if(LinearInterpolation_fixedsize<numt,numt>::min()>=LinearInterpolation_fixedsize<numt,numt>::max())
			throw math_h_error<RandomValueGenerator>("Probability density function has points below zero");
		distr=uniform_real_distribution<numt>(LinearInterpolation_fixedsize<numt,numt>::min(),LinearInterpolation_fixedsize<numt,numt>::max());
	}
	RandomValueGenerator(numt x1, numt x2):RandomValueGenerator([](double){return 1.0;},x1,x2,2){}
	virtual ~RandomValueGenerator(){}
	numt operator ()(RG&generator){
		return LinearInterpolation_fixedsize<numt,numt>::operator()(distr(generator));
	}
};
#endif
