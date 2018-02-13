// this file is distributed under
// LGPLv3 license
#ifndef _____EXTENDED_UNCERTAINTIES_ESTIMATION_____________
#define _____EXTENDED_UNCERTAINTIES_ESTIMATION_____________
#include <memory>
#include <math.h>
#include "error.h"
#include "sigma.h"
#include "randomfunc.h"
namespace MathTemplates
{
template<size_t d,class n=double>class value_f;
template<class n=double>class value_f_chain;
template<typename numt = double>
class abstract_value_with_uncertainty_numeric:public abstract_value_with_uncertainty<numt>
{
public:
    typedef numt NumType;
    template<size_t d,class n>friend class value_f;
    template<class n>friend class value_f_chain;
protected:
    virtual numt mmeasure()const=0;
private:
    std::shared_ptr<StandardDeviation<numt>> m_cache;
    static size_t mm_number;
    const StandardDeviation<numt> &VAL()const
    {
        using namespace std;
        if (!m_cache) {
	    const auto tmp=make_shared<StandardDeviation<numt>>();
            const_cast<abstract_value_with_uncertainty_numeric&>(*this).m_cache =tmp;
	    for (size_t i=0;i<mm_number;i++)
		tmp->operator<<(mmeasure());
        }
        return *m_cache;
    }
public:
    abstract_value_with_uncertainty_numeric():m_cache(nullptr){}
    static inline size_t micro_measurements_count(){return abstract_value_with_uncertainty_numeric::mm_number;}
    static inline void set_micro_measurements_count(size_t v){abstract_value_with_uncertainty_numeric::mm_number=v;}
    virtual ~abstract_value_with_uncertainty_numeric(){}
    virtual const numt&val()const override
    {
        return VAL().val();
    }
    virtual const numt&uncertainty()const override
    {
        return VAL().uncertainty();
    }
};
template<typename numt>
size_t abstract_value_with_uncertainty_numeric<numt>::mm_number=1000;
//measured values
template<class numt=double>
class value_numeric_const:public abstract_value_with_uncertainty_numeric<numt>{
private:
    numt m_val;
protected:
    virtual numt mmeasure()const override{return m_val;}
public:
    value_numeric_const(const value_numeric_const&source):abstract_value_with_uncertainty_numeric<numt>(),m_val(source.m_val){}
    value_numeric_const(const numt&v):abstract_value_with_uncertainty_numeric<numt>(),m_val(v){}
    virtual ~value_numeric_const(){}
};
template<class numt=double,class DISTR=RandomGauss<numt>>
class value_numeric_distr:public abstract_value_with_uncertainty_numeric<numt>{
private:
    DISTR m_distr;
protected:
    virtual numt mmeasure()const override{
	return m_distr();
    }
public:
    value_numeric_distr(const value_numeric_distr&S):abstract_value_with_uncertainty_numeric<numt>(),m_distr(S.m_distr){}
    value_numeric_distr(const DISTR&D):abstract_value_with_uncertainty_numeric<numt>(),m_distr(D){}
    template<typename... Args>
    value_numeric_distr(Args... args):abstract_value_with_uncertainty_numeric<numt>(),m_distr(args...){}
    virtual ~value_numeric_distr(){}
};
template<class numt>
inline value_numeric_distr<numt,RandomGauss<numt>> value_numeric(const abstract_value_with_uncertainty<numt>&a){
    return value_numeric_distr<numt,RandomGauss<numt>>(a.val(),a.uncertainty());
}
//calculated values
template<class numt>
class value_f_chain:public abstract_value_with_uncertainty_numeric<numt>{
private:
    std::function<numt(const Chain<numt>&)> m_func;
    Chain<const abstract_value_with_uncertainty_numeric<numt>*> m_chain;
protected:
    virtual numt mmeasure()const override{
	Chain<numt> P;
	for(const auto A:m_chain)P.push_back(A->mmeasure());
	return m_func(P);
    }
public:
    value_f_chain(const value_f_chain&source)
	:abstract_value_with_uncertainty_numeric<numt>()
	    ,m_func(source.m_func),m_chain(source.m_chain){}
    template<class VType>
    value_f_chain(std::function<numt(const Chain<numt>&)> F,const Chain<VType>&chain)
	:abstract_value_with_uncertainty_numeric<numt>(),m_func(F){
	    for(const auto&a:chain)m_chain.push_back(&a);
	}
    virtual ~value_f_chain(){}
};
template<class numt,class VType>
inline value_f_chain<numt> FUNC(std::function<numt(const Chain<numt>&)> F,const Chain<VType>&chain){
    return value_f_chain<numt>(F,chain);
}
template<class numt>
class value_f<0,numt>:public abstract_value_with_uncertainty_numeric<numt>{
private:
    std::function<numt()> m_func;
protected:
    virtual numt mmeasure()const override{return m_func();}
public:
    value_f(const value_f&s)
	:abstract_value_with_uncertainty_numeric<numt>(),m_func(s.m_func){}
    inline value_f(std::function<numt()> F)
	:abstract_value_with_uncertainty_numeric<numt>(),m_func(F){}
    virtual ~value_f(){}
};


template<class numt>
class value_f<1,numt>:public value_f<0,numt>{
public:
    value_f(const value_f&s):value_f<0,numt>(s){}
    template<class ArgType1>
    inline value_f(std::function<numt(const numt&)> F,ArgType1 a)
	:value_f<0,numt>([a,F](){return F(dynamic_cast<const abstract_value_with_uncertainty_numeric<numt>*>(&a)->mmeasure());}){}
    virtual ~value_f(){}
};
template<class numt>
class value_f<2,numt>:public value_f<1,numt>{
public:
    value_f(const value_f&s):value_f<1,numt>(s){}
    template<class ArgType1,class ArgType2>
    inline value_f(std::function<numt(const numt&,const numt&)> F,ArgType1 a,ArgType2 b)
	:value_f<1,numt>([a,F](const numt&y){
	    return F(dynamic_cast<const abstract_value_with_uncertainty_numeric<numt>*>(&a)->mmeasure(),y);
	},b){}
    virtual ~value_f(){}
};


template<class ArgType1>
inline value_f<1,typename ArgType1::NumType> SQRT(ArgType1 a){
    return value_f<1,typename ArgType1::NumType>([](const typename ArgType1::NumType&x){return sqrt(x);},a);
}
template<class ArgType1>
inline value_f<1,typename ArgType1::NumType> SQR(ArgType1 a){
    return value_f<1,typename ArgType1::NumType>([](const typename ArgType1::NumType&x){return x*x;},a);
}
template<class ArgType1>
inline value_f<1,typename ArgType1::NumType> EXP(ArgType1 a){
    return value_f<1,typename ArgType1::NumType>([](const typename ArgType1::NumType&x){return exp(x);},a);
}
template<class ArgType1>
inline value_f<1,typename ArgType1::NumType> LOG(ArgType1 a){
    return value_f<1,typename ArgType1::NumType>([](const typename ArgType1::NumType&x){return log(x);},a);
}
template<class ArgType1>
inline value_f<1,typename ArgType1::NumType> SIN(ArgType1 a){
    return value_f<1,typename ArgType1::NumType>([](const typename ArgType1::NumType&x){return sin(x);},a);
}
template<class ArgType1>
inline value_f<1,typename ArgType1::NumType> COS(ArgType1 a){
    return value_f<1,typename ArgType1::NumType>([](const typename ArgType1::NumType&x){return cos(x);},a);
}
template<class ArgType1>
inline value_f<1,typename ArgType1::NumType> TAN(ArgType1 a){
    return value_f<1,typename ArgType1::NumType>([](const typename ArgType1::NumType&x){return tan(x);},a);
}
template<class ArgType1>
inline value_f<1,typename ArgType1::NumType> ATAN(ArgType1 a){
    return value_f<1,typename ArgType1::NumType>([](const typename ArgType1::NumType&x){return atan(x);},a);
}
template<class ArgType1,class ArgType2>
inline value_f<2,typename ArgType1::NumType> ATAN2(ArgType1 a,ArgType2 b){
    return value_f<2,typename ArgType1::NumType>([](
	const typename ArgType1::NumType&x,const typename ArgType1::NumType&y
    ){return atan2(x,y);},a,b);
}
template<class ArgType1,class ArgType2>
inline value_f<2,typename ArgType1::NumType> POW(ArgType1 a,ArgType2 b
){
    return value_f<2,typename ArgType1::NumType>([](
	const typename ArgType1::NumType&x,const typename ArgType1::NumType&y
    ){return pow(x,y);},a,b);
}
template<class ArgType1>
inline value_f<1,typename ArgType1::NumType> POW(ArgType1 a,const typename ArgType1::NumType&b){
    return value_f<1,typename ArgType1::NumType>([&b](const typename ArgType1::NumType&x){return pow(x,b);},a);
}
template<class ArgType1,class ArgType2>
inline value_f<2,typename ArgType1::NumType> operator+(ArgType1 a,ArgType2 b){
    return value_f<2,typename ArgType1::NumType>([](
	const typename ArgType1::NumType&x,const typename ArgType1::NumType&y
    ){return x+y;},a,b);
}
template<class ArgType1,class ArgType2>
inline value_f<2,typename ArgType1::NumType> operator-(ArgType1 a,ArgType2 b){
    return value_f<2,typename ArgType1::NumType>([](
	const typename ArgType1::NumType&x,const typename ArgType1::NumType&y
    ){return x-y;},a,b);
}
template<class ArgType1,class ArgType2>
inline value_f<2,typename ArgType1::NumType> operator*(ArgType1 a,ArgType2 b){
    return value_f<2,typename ArgType1::NumType>([](
	const typename ArgType1::NumType&x,const typename ArgType1::NumType&y
    ){return x*y;},a,b);
}
template<class ArgType1,class ArgType2>
inline value_f<2,typename ArgType1::NumType> operator/(ArgType1 a,ArgType2 b){
    return value_f<2,typename ArgType1::NumType>([](
	const typename ArgType1::NumType&x,const typename ArgType1::NumType&y
    ){return x/y;},a,b);
}

#ifdef ____middle_version_of_math_h_____
template<size_t dim,class numt>
class value_f:public value_f<dim-1,numt>{
public:
    value_f(const value_f&s):value_f<dim-1,numt>(s){}
    template<class FUNC,class ArgType,typename... Args>
    inline value_f(FUNC F,ArgType a,Args... args)
	:value_f<dim-1,numt>([a,F](auto...z){
	    return F(dynamic_cast<const abstract_value_with_uncertainty_numeric<numt>*>(&a)->mmeasure(),z...);
	},args...){}
    virtual ~value_f(){}
};
#endif
template<class FF, class ArgType,typename... Args>
inline value_f<sizeof...(Args)+1,typename ArgType::NumType> FUNC(FF F,ArgType a,Args... args){
    return value_f<sizeof...(Args)+1,typename ArgType::NumType>(F,a,args...);
}

};
#endif
