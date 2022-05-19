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
    FunctionWrap(const FunctionWrap&source):func(source.func){}
    FunctionWrap(std::function<result(Args...)>source):func(source){}
    FunctionWrap(const result&v):FunctionWrap([v](Args...){return v;}){}
    FunctionWrap(int v):FunctionWrap([v](Args...){return result(v);}){}
    template<class FUNC> 
    FunctionWrap(FUNC F):func([F](Args... args)->result{return F(args...);}){}
    virtual ~FunctionWrap() {}
    virtual result operator()(Args...args)const override{return func(args...);}
    template<class FUNC> 
    FunctionWrap operator+(const FUNC& second)const{
	    FunctionWrap first(*this);
	    return FunctionWrap([first,second](Args... args){return first(args...)+second(args...);});
    }
    template<class FUNC> 
    FunctionWrap operator-(const FUNC& second)const{
	    FunctionWrap first(*this);
	    return FunctionWrap([first,second](Args... args){return first(args...)-second(args...);});
    }
};
template<class result, typename... Args>
inline FunctionWrap<result,Args...> operator-(const FunctionWrap<result,Args...>&source){
    return FunctionWrap<result,Args...>([source](Args...args){return -source(args...);});
}
template<class result, typename... Args>
inline const FunctionWrap<result,Args...>& operator+(const FunctionWrap<result,Args...>&source){
    return source;
}
template<class numtz, class numty, typename... Args>
inline auto operator*(const numtz&alpha, const FunctionWrap<numty,Args...>&source)->FunctionWrap<decltype(alpha*numty(0)),Args...>{
    return FunctionWrap<decltype(alpha*numty(0)),Args...>([alpha,source](Args...args){return alpha*source(args...);});
}

template<class numty=double,class numtx=numty>
class Opr{
private:
    FunctionWrap<numty,const numtx&> m_func;
public:
    Opr(const FunctionWrap<numty,const numtx&>&source):m_func(source){}
    Opr(const Opr&source):m_func(source.m_func){}
    Opr(const numty&v):m_func([v](const numtx&){return v;}){}
    Opr(int v):m_func([v](const numtx&){return numty(v);}){}
    template<class OPR>Opr(OPR O):m_func(O){}
    virtual ~Opr(){}
    inline numty operator*(const numtx& f)const{return m_func(f);}
    inline Opr operator+(const Opr&second)const{
	    return Opr([*this,second](const numtx&F){return (operator*(F))+(second*F);});
    }
    inline Opr operator-(const Opr&second)const{
	    return Opr([*this,second](const numtx&F){return (operator*(F))-(second*F);});
    }
};
template<class numty,class numtx>
inline Opr<numty,numtx> operator-(const Opr<numty,numtx>&source){
    return Opr<numty,numtx>([source](const numtx&F){return -(source*F);});
}
template<class numty,class numtx>
inline const Opr<numty,numtx>& operator+(const Opr<numty,numtx>&source){
    return source;
}
template<class numtz,class numty,class numtx>
inline auto operator*(const numtz&alpha, const Opr<numty,numtx>&source){
    return Opr<decltype(alpha*numty(0)),numtx>([alpha,source](const numtx&F){return alpha*(source*F);});
}
//constants
template<class numt = double>
inline numt PI(){return 3.14159265358;}
template<class numt = double>
inline numt E(){return 2.718281828459;}


template<class numt = double>
numt FermiFunc(const numt &x, const numt &X_border, const numt &diffuse)
{
    return 1.0 / (1.0 + exp((x - X_border) / diffuse));
}

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
template<unsigned int Power,int index_offset = 0, class numt, class indexer>
inline numt Polynom(const numt &x, const indexer&coefficients)
{
    static_assert(index_offset >= 0,"Polynom offset index is out of range");
    numt res=0;
    for(int i=Power;i>=0;i--){
        res*=x;
        res+=coefficients[index_offset+i];
    }
    return res;
}
};
#endif
