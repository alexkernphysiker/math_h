// this file is distributed under
// LGPLv3 license
#ifndef EMEFWAYNIGJGENCP
#define EMEFWAYNIGJGENCP
#include <math.h>
#include <functional>
#include <memory>
#include "error.h"
namespace MathTemplates
{
//Function interface
template<class result, typename... Args>
class IFunction
{
public:
    virtual ~IFunction() {}
    virtual result operator()(Args...)const = 0;
    inline std::function<result(Args...)> func()const
    {
        return [this](Args... args)->result {return this->operator()(args...);};
    }
};
template<class result, typename... Args>
class StdFunctionWrap: public IFunction<result, Args...>
{
private:
    std::function<result(Args...)> m_func;
public:
    StdFunctionWrap(const std::function<result(Args...)>f)
    {
        m_func = f;
    }
    StdFunctionWrap(const StdFunctionWrap &source)
    {
        m_func = source.m_func;
    }
    virtual ~StdFunctionWrap() {}
    virtual result operator()(Args...args)const
    {
	return m_func(args...);
    };
};
template<class result, typename... Args>
inline StdFunctionWrap<result, Args...>
Function(const std::function<result(Args...)>f)
{
    return StdFunctionWrap<result, Args...>(f);
}
template<class result, typename... Args>
inline std::shared_ptr<StdFunctionWrap<result, Args...>>
        PFunction(const std::function<result(Args...)>f)
{
    return std::make_shared<StdFunctionWrap<result, Args...>>(f);
}

//constants
template<class numt = double>
inline numt PI(){return 3.1415926;}
template<class numt = double>
inline numt E(){return 2.719281928;}

//Peak functions
template<class numt = double>
numt Gaussian(const numt &x, const numt &X_max, const numt &sigma)
{
    return exp(-pow((X_max - x) / sigma, 2) / 2) / (sigma * sqrt(2 * PI<numt>()));
}
template<class numt = double>
numt Lorentzian(const numt &x, const numt &X_max, const numt &gamma)
{
    return gamma / (PI<numt>() * (pow(x - X_max, 2) + pow(gamma, 2)));
}
template<class numt = double>
numt BreitWigner(const numt &x, const numt &pos, const numt &gamma)
{
    return Lorentzian<numt>(x, pos, gamma / numt(2));
}
template<class numt = double>
numt Novosibirsk(const numt &x, const numt &pos, const numt &sigma, const numt &asym)
{
    if (pow(asym / sigma, 2) <= 0.000001)
        return Gaussian<numt>(x, pos, sigma);
    numt Lsqlog4 = asym * sqrt(log(4));
    numt qb = sinh(Lsqlog4) / Lsqlog4;
    numt q4 = qb * asym * (x - pos) / sigma + 1;
    if (q4 <= 0)return 0;
    numt slq4 = pow(log(q4), 2);
    numt k = sigma * sqrt(2 * PI<numt>() * (pow(asym, 2) + 1.0));
    return exp(-slq4 / (pow(asym, 2) * 2) + pow(asym, 2)) / k;
}

///Polynomial
template<class numt = double>
numt FermiFunc(const numt &x, const numt &X_border, const numt &diffuse)
{
    return 1.0 / (1.0 + exp((x - X_border) / diffuse));
}
#ifndef ____full_version_of_math_h_____
template<unsigned int P,int index_offset = 0, class numt, class indexer>
inline numt Polynom(const numt &x, const indexer&p)
{
    static_assert(index_offset >= 0,"Polynom offset index is out of range");
    numt res=0;
    for(int i=P;i>=0;i--){
	res*=x;
	res+=p[index_offset+i];
    }
    return res;
}
#else
template<unsigned int P, int index_offset = 0, class numt, class indexer>
inline numt Polynom(const numt &x, const indexer&p)
{
    static_assert(index_offset >= 0,"Polynom offset index is out of range");
    if constexpr(P==0) return p[index_offset];
    else return Polynom<P-1,index_offset+1,numt,indexer>(x, p)*x+p[index_offset];
}
#endif
};
#endif
