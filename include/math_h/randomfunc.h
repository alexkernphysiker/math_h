// this file is distributed under 
// MIT license
#ifndef XRYPRAVJWTJCYPQI
#define XRYPRAVJWTJCYPQI
#include <random>
#include <memory>
#include <functional>
#include <math.h>
#include "tabledata.h"
#include "interpolate.h"
#include "integrate.h"
namespace MathTemplates{
    template<class numt=double,class RG=std::mt19937>
    class RandomValueTableDistr:public IFunction<numt,RG&>{
    private:
	LinearInterpolation<numt,numt> reverse_distr_func;
	std::pair<numt,numt> range;
    public:
	RandomValueTableDistr(const RandomValueTableDistr &R):reverse_distr_func(R.reverse_distr_func),range(R.range){}
	RandomValueTableDistr(const SortedPoints<numt>&distribution_density)
	    :reverse_distr_func(Int_Trapez_Table_PositiveStrict(distribution_density).TransponateAndSort())
	    ,range(reverse_distr_func.left().X(),reverse_distr_func.right().X()){}
	RandomValueTableDistr(const std::function<numt(numt)> distribution_density,const SortedChain<numt>&chain)
	    :RandomValueTableDistr(SortedPoints<numt>(distribution_density,chain)){}
	RandomValueTableDistr(const std::function<numt(numt)> distribution_density,const SortedChain<numt>&&chain)
	    :RandomValueTableDistr(SortedPoints<numt>(distribution_density,chain)){}
	RandomValueTableDistr(const std::function<numt(numt)> distribution_density,const std::initializer_list<numt>&&chain)
	    :RandomValueTableDistr(SortedPoints<numt>(distribution_density,chain)){}
	RandomValueTableDistr(const numt x1,const numt x2)
	    :RandomValueTableDistr({point<numt>(x1,numt(1)),point<numt>(x2,numt(1))}){}
	virtual ~RandomValueTableDistr(){}
	virtual numt operator ()(RG&generator)const override{
	    std::uniform_real_distribution<numt> uniform(range.first,range.second);
	    return reverse_distr_func(uniform(generator));
	}
    };
    template<class numt=double,class RG=std::mt19937>
    class RandomUniform:public IFunction<numt,RG&>{
    private:
	std::shared_ptr<std::uniform_real_distribution<numt>> f_distr;
    public:
	RandomUniform(const numt x1,const numt x2):f_distr(std::make_shared<std::uniform_real_distribution<numt>>(x1,x2)){}
	RandomUniform(const RandomUniform&source):f_distr(source.f_distr){}
	virtual ~RandomUniform(){}
	virtual numt operator ()(RG&generator)const override{
	    return f_distr->operator()(generator);
	}
    };
    template<class numt=double,class RG=std::mt19937>
    class RandomGauss:public IFunction<numt,RG&>{
    private:
	std::shared_ptr<std::normal_distribution<numt>> f_distr;
    public:
	RandomGauss(const numt x1,const numt x2):f_distr(std::make_shared<std::normal_distribution<numt>>(x1,x2)){}
	RandomGauss(const RandomGauss&source):f_distr(source.f_distr){}
	virtual ~RandomGauss(){}
	virtual numt operator ()(RG&generator)const override{
	    return f_distr->operator()(generator);
	}
    };
};
#endif
