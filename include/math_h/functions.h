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
class FunctionWrap:public IFunction<result,Args...>
{
private:
    std::function<result(Args...)> func;
public:
    template<class FUNC>
    FunctionWrap(FUNC F){func=[F](Args... args)->result{return F(args...);};}
    virtual ~FunctionWrap() {}
    virtual result operator()(Args...args)const override{return func(args...);}
};

template<class numtx=double,class numty=numtx>
class Opr{
private:
    FunctionWrap<numty,const IFunction<numty,const numtx&>&> m_func;
public:
    template<class OPR>Opr(OPR O):m_func(O){}
    virtual ~Opr(){}
    inline numty operator*(const IFunction<numty,const numtx&>& f)const{return m_func(f);}
    template<class FUNC>
    inline numty operator*(FUNC F)const{
	FunctionWrap<numty,const numtx&> f(F);
	return operator*(static_cast<const IFunction<numty,const numtx&>&>(f));
    }
};

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
};
#endif
