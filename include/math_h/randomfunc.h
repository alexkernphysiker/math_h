// this file is distributed under 
// MIT license
#ifndef XRYPRAVJWTJCYPQI
#define XRYPRAVJWTJCYPQI
#include <random>
#include <functional>
#include <math.h>
#include "tabledata.h"
#include "interpolate.h"
#include "integrate.h"
namespace MathTemplates{
	template<class numt,class RG=std::mt19937>
	class RandomValueGenerator:public IFunction<numt,RG&>{
	private:
		LinearInterpolation<numt,numt> reverse_distr_func;
		std::pair<numt,numt> range;
	public:
		RandomValueGenerator(const RandomValueGenerator &R):reverse_distr_func(R.reverse_distr_func),range(R.range){}
		RandomValueGenerator(const SortedPoints<numt>&distribution_density)
		:reverse_distr_func(Int_Trapez_Table_PositiveStrict(distribution_density).TransponateAndSort())
		,range(reverse_distr_func.left().X(),reverse_distr_func.right().X()){}
		RandomValueGenerator(const std::function<numt(numt)> distribution_density,const SortedChain<numt>&chain)
		:RandomValueGenerator(SortedPoints<numt>(distribution_density,chain)){}
		RandomValueGenerator(const std::function<numt(numt)> distribution_density,const SortedChain<numt>&&chain)
		:RandomValueGenerator(SortedPoints<numt>(distribution_density,chain)){}
		RandomValueGenerator(const std::function<numt(numt)> distribution_density,const std::initializer_list<numt>&&chain)
		:RandomValueGenerator(SortedPoints<numt>(distribution_density,chain)){}
		RandomValueGenerator(const numt x1,const numt x2):RandomValueGenerator({point<numt>(x1,numt(1)),point<numt>(x2,numt(1))}){}
		virtual ~RandomValueGenerator(){}
		virtual numt operator ()(RG&generator)const override{
			std::uniform_real_distribution<numt> uniform(range.first,range.second);
			return reverse_distr_func(uniform(generator));
		}
	};
};
#endif
