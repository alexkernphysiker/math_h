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
		std::pair<numt,numt> range;
	public:
		RandomValueGenerator(const RandomValueGenerator &R):reverse_distr_func(R.reverse_distr_func),range(R.range){}
		RandomValueGenerator(const SortedPoints<numt,numt>&distribution_density)
			:reverse_distr_func(Int_Trapez_Table_PositiveStrict(distribution_density).Transponate())
			,range(reverse_distr_func.min(),reverse_distr_func.max()){}
		RandomValueGenerator(const function<numt(numt)> distribution_density,const vector<numt>&chain)
		:RandomValueGenerator(SortedPoints<numt>(distribution_density,chain)){}
		RandomValueGenerator(const function<numt(numt)> distribution_density,const vector<numt>&&chain)
		:RandomValueGenerator(SortedPoints<numt>(distribution_density,chain)){}
		RandomValueGenerator(const function<numt(numt)> distribution_density,const initializer_list<numt>&&chain)
		:RandomValueGenerator(SortedPoints<numt>(distribution_density,chain)){}
		RandomValueGenerator(const numt x1,const numt x2):RandomValueGenerator({point<numt>(x1,numt(1)),point<numt>(x2,numt(1))}){}
		virtual ~RandomValueGenerator(){}
		numt operator ()(RG&generator)const{
			uniform_real_distribution<numt> uniform(range.first,range.second);
			return reverse_distr_func(uniform(generator));
		}
		const function<numt()>func(RG&generator)const{
			return [&generator,this]()->numt{return operator()(generator);};
		}
	};
};
#endif
