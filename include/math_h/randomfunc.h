// this file is distributed under
// LGPLv3 license
#ifndef XRYPRAVJWTJCYPQI
#define XRYPRAVJWTJCYPQI
#if __cplusplus<201100L
#error c++>=11 is needed for using math_h headers
#endif
#include <random>
#include <memory>
#include "tabledata.h"
#include "interpolate.h"
namespace MathTemplates
{
template<typename Engine=std::mt19937>
class RandomEngine{
public:
    static Engine&Instance(){
	static Engine res;
	return res;
    }
};
template<class numt = double>
using RandomValueGenerator = IFunction<numt>;
template<class numt = double, class RG = RandomEngine<>>
class RandomValueTableDistr: public RandomValueGenerator<numt>
{
private:
    const ReverseIntegratedLinearInterpolation<numt> reverse_distr_func;
    const std::shared_ptr<std::uniform_real_distribution<numt>> f_distr;
public:
    RandomValueTableDistr(const RandomValueTableDistr &R):
        reverse_distr_func(R.reverse_distr_func), f_distr(R.f_distr) {}
    RandomValueTableDistr(const ReverseIntegratedLinearInterpolation<numt> &distribution_density):
        reverse_distr_func(distribution_density),
        f_distr(
            std::make_shared<std::uniform_real_distribution<numt>>(
                reverse_distr_func.left().X(), reverse_distr_func.right().X()
            )
        ) {}
    RandomValueTableDistr(const LinearInterpolation<numt> &source):
        RandomValueTableDistr(ReverseIntegratedLinearInterpolation<numt>(source)) {}
    RandomValueTableDistr(const Points<numt>&source):
        RandomValueTableDistr(LinearInterpolation<numt>(source)) {}
    RandomValueTableDistr(const std::function<numt(numt)> distribution_density, const SortedChain<numt> &chain):
        RandomValueTableDistr(LinearInterpolation<numt>(distribution_density, chain)) {}
    virtual ~RandomValueTableDistr() {}
    virtual numt operator()()const override
    {
        return reverse_distr_func(f_distr->operator()(RG::Instance()));
    }
};
template<class numt = double, class RG = RandomEngine<>>
class RandomUniform: public RandomValueGenerator<numt>
{
private:
    std::shared_ptr<std::uniform_real_distribution<numt>> f_distr;
public:
    RandomUniform(const numt x1, const numt x2): f_distr(std::make_shared<std::uniform_real_distribution<numt>>(x1, x2)) {}
    RandomUniform(const RandomUniform &source): f_distr(source.f_distr) {}
    virtual ~RandomUniform() {}
    virtual numt operator()()const override
    {
        return f_distr->operator()(RG::Instance());
    }
};
template<class numt = double, class RG = RandomEngine<>>
class RandomGauss: public RandomValueGenerator<numt>
{
private:
    std::shared_ptr<std::normal_distribution<numt>> f_distr;
public:
    RandomGauss(const numt x1, const numt x2): f_distr(std::make_shared<std::normal_distribution<numt>>(x1, x2)) {}
    RandomGauss(const RandomGauss &source): f_distr(source.f_distr) {}
    virtual ~RandomGauss() {}
    virtual numt operator()()const override
    {
        return f_distr->operator()(RG::Instance());
    }
};
};
#endif
