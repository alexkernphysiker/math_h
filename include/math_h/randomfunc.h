// this file is distributed under
// LGPLv3 license
#ifndef XRYPRAVJWTJCYPQI
#define XRYPRAVJWTJCYPQI
#include <random>
#include <memory>
#include "error.h"
#include "functions.h"
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
    interpolation_details::ReverseIntegratedLinearInterpolation<numt> reverse_distr_func;
    mutable std::uniform_real_distribution<numt> f_distr;
public:
    typedef numt NumberType;
    RandomValueTableDistr(const RandomValueTableDistr &R) = delete;
    RandomValueTableDistr& operator=(const RandomValueTableDistr &R) = delete;
    RandomValueTableDistr(RandomValueTableDistr &&R) = default;
    RandomValueTableDistr(interpolation_details::ReverseIntegratedLinearInterpolation<numt> &&distribution_density):
        reverse_distr_func(std::move(distribution_density)),
        f_distr(reverse_distr_func.left().X(), reverse_distr_func.right().X())
        {}
    RandomValueTableDistr(LinearInterpolation<numt> &&source):
        RandomValueTableDistr(interpolation_details::ReverseIntegratedLinearInterpolation<numt>(std::move(source))) {}
    RandomValueTableDistr(Points<numt>&&source):
        RandomValueTableDistr(LinearInterpolation<numt>(std::move(source))) {}
    RandomValueTableDistr(const IFunction<numt,const numt&>& distribution_density, const SortedChain<numt> &chain):
        RandomValueTableDistr(LinearInterpolation<numt>(distribution_density, chain)) {}
    template<class FUNC>
    RandomValueTableDistr(FUNC F,const SortedChain<numt> &c):
	RandomValueTableDistr(static_cast<const IFunction<numt,const numt&>&>(FunctionWrap<numt,const numt&>(F)),c){}

    virtual ~RandomValueTableDistr() {}
    virtual numt operator()()const override
    {
        return reverse_distr_func(f_distr(RG::Instance()));
    }
};
template<class numt = double, class RG = RandomEngine<>>
class RandomUniform: public RandomValueGenerator<numt>
{
private:
    mutable std::uniform_real_distribution<numt> f_distr;
public:
    typedef numt NumberType;
    RandomUniform(const numt& x1, const numt& x2): f_distr(x1, x2) {}
    RandomUniform(const RandomUniform &source) = default;
    virtual ~RandomUniform() {}
    virtual numt operator()()const override
    {
        return f_distr(RG::Instance());
    }
};
template<class numt = double, class RG = RandomEngine<>>
class RandomGauss: public RandomValueGenerator<numt>
{
private:
    mutable std::normal_distribution<numt> f_distr;
public:
    typedef numt NumberType;
    RandomGauss(const numt& x1, const numt& x2): f_distr(x1, x2) {}
    RandomGauss(const RandomGauss &source) = default;
    virtual ~RandomGauss() {}
    virtual numt operator()()const override
    {
        return f_distr(RG::Instance());
    }
};
template<class numt = double, class RG = RandomEngine<>>
auto Poisson(const numt&avr){
    const auto chain=ChainWithCount(size_t(avr*3+1),numt(0),avr*3);
    return RandomValueTableDistr<numt,RG>([&avr](const numt&k){
	    numt kf=1;
	    for(numt i=2;i<=k;i+=1) {
            kf *= i;
        }
	    return pow(avr,k) / kf;
    },chain);
}
};
#endif
