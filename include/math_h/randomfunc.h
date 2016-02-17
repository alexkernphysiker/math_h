// this file is distributed under 
// MIT license
#ifndef XRYPRAVJWTJCYPQI
#define XRYPRAVJWTJCYPQI
#include <random>
#include <functional>
#include <math.h>
#include "interpolate.h"
#include "integrate.h"
#include "error.h"
namespace MathTemplates{
	using namespace std;
	template<class numt,class RG=mt19937>
	class RandomValueGenerator{
	private:
		LinearInterpolation<numt,numt> reverse_distr_func;
		uniform_real_distribution<numt> uniform;
	public:
		RandomValueGenerator(const RandomValueGenerator &R):reverse_distr_func(R.reverse_distr_func),uniform(R.uniform){}
		RandomValueGenerator(const LinearInterpolation<numt,numt>&distribution_density)
			:reverse_distr_func(Int_Trapez_Table_PositiveStrict(distribution_density).Transponate())
			,uniform(reverse_distr_func.min(),reverse_distr_func.max()){}
		RandomValueGenerator(function<numt(numt)> distribution_density,const vector<numt>&chain)
			:RandomValueGenerator(LinearInterpolation<numt>(distribution_density,chain)){}
		RandomValueGenerator(function<numt(numt)> distribution_density,vector<numt>&&chain)
			:RandomValueGenerator(LinearInterpolation<numt>(distribution_density,chain)){}
		RandomValueGenerator(function<numt(numt)> distribution_density,initializer_list<numt>&&chain)
			:RandomValueGenerator(LinearInterpolation<numt>(distribution_density,chain)){}
		RandomValueGenerator(numt x1, numt x2):RandomValueGenerator({make_pair(x1,numt(1)),make_pair(x2,numt(1))}){}
		virtual ~RandomValueGenerator(){}
		numt operator ()(RG&generator){
			return reverse_distr_func(uniform(generator));
		}
		function<numt()>func(RG&generator){
			return [&generator,this]()->numt{return operator()(generator);};
		}
	};
};
#endif
