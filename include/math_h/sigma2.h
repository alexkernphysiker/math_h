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
template<typename numt = double>
class abstract_value_with_extended_uncertainty:public abstract_value_with_uncertainty<numt>
{
    template<size_t d,class n>friend class value_f;
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
            const_cast<abstract_value_with_extended_uncertainty&>(*this).m_cache =tmp;
	    for (size_t i=0;i<mm_number;i++)
		tmp->operator<<(mmeasure());
        }
        return *m_cache;
    }
public:
    abstract_value_with_extended_uncertainty():m_cache(nullptr){}
    static inline size_t micro_measurements_count(){return abstract_value_with_extended_uncertainty::mm_number;}
    static inline void set_micro_measurements_count(size_t v){abstract_value_with_extended_uncertainty::mm_number=v;}
    virtual ~abstract_value_with_extended_uncertainty(){}
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
size_t abstract_value_with_extended_uncertainty<numt>::mm_number=1000;
//measured values
template<class numt=double>
class value_const:public abstract_value_with_extended_uncertainty<numt>{
private:
    numt m_val;
protected:
    virtual numt mmeasure()const override{return m_val;}
public:
    value_const(const numt&v):abstract_value_with_extended_uncertainty<numt>(),m_val(v){}
    virtual ~value_const(){}
};
template<class numt=double,class DISTR=RandomGauss<numt>>
class value_ext:public abstract_value_with_extended_uncertainty<numt>{
private:
    DISTR m_distr;
protected:
    virtual numt mmeasure()const override{
	return m_distr();
    }
public:
    value_ext(const DISTR&D):abstract_value_with_extended_uncertainty<numt>(),m_distr(D){}
    template<typename... Args>
    value_ext(Args... args):abstract_value_with_extended_uncertainty<numt>(),m_distr(args...){}
    virtual ~value_ext(){}
};
template<class numt>
inline value_ext<numt,RandomGauss<numt>> ext_value(const abstract_value_with_uncertainty<numt>&a){
    return value_ext<numt,RandomGauss<numt>>(a.val(),a.uncertainty());
}
//calculated values
template<class numt>
class value_f<0,numt>:public abstract_value_with_extended_uncertainty<numt>{
private:
    std::function<numt()> m_func;
protected:
    virtual numt mmeasure()const override{return m_func();}
public:
    value_f(const std::function<numt()> F)
	:abstract_value_with_extended_uncertainty<numt>(),m_func(F){}
    virtual ~value_f(){}
};
template<class numt>
class value_f<1,numt>:public value_f<0,numt>{
public:
    value_f(
	const std::function<numt(const numt&)> F,
	const abstract_value_with_extended_uncertainty<numt>&a
    )
	:value_f<0,numt>([&a,F](){return F(a.mmeasure());}){}
    virtual ~value_f(){}
};
template<class numt>
class value_f<2,numt>:public value_f<0,numt>{
public:
    value_f(
	const std::function<numt(const numt&,const numt&)> F,
	const abstract_value_with_extended_uncertainty<numt>&a,
	const abstract_value_with_extended_uncertainty<numt>&b
    )
	:value_f<0,numt>([&a,&b,F](){return F(a.mmeasure(),b.mmeasure());}){}
    virtual ~value_f(){}
};
template<class numt>
inline value_f<1,numt> SQRT(
    	const abstract_value_with_extended_uncertainty<numt>&a
){
    return value_f<1,numt>([](const numt&x){return sqrt(x);},a);
}
template<class numt>
inline value_f<1,numt> SQR(
    	const abstract_value_with_extended_uncertainty<numt>&a
){
    return value_f<1,numt>([](const numt&x){return x*x;},a);
}
template<class numt>
inline value_f<1,numt> EXP(
    	const abstract_value_with_extended_uncertainty<numt>&a
){
    return value_f<1,numt>([](const numt&x){return exp(x);},a);
}
template<class numt>
inline value_f<1,numt> LOG(
    	const abstract_value_with_extended_uncertainty<numt>&a
){
    return value_f<1,numt>([](const numt&x){return log(x);},a);
}
template<class numt>
inline value_f<1,numt> SIN(
    	const abstract_value_with_extended_uncertainty<numt>&a
){
    return value_f<1,numt>([](const numt&x){return sin(x);},a);
}
template<class numt>
inline value_f<1,numt> COS(
    	const abstract_value_with_extended_uncertainty<numt>&a
){
    return value_f<1,numt>([](const numt&x){return cos(x);},a);
}
template<class numt>
inline value_f<1,numt> TAN(
    	const abstract_value_with_extended_uncertainty<numt>&a
){
    return value_f<1,numt>([](const numt&x){return tan(x);},a);
}
template<class numt>
inline value_f<1,numt> ATAN(
    	const abstract_value_with_extended_uncertainty<numt>&a
){
    return value_f<1,numt>([](const numt&x){return atan(x);},a);
}
template<class numt>
inline value_f<2,numt> ATAN2(
    	const abstract_value_with_extended_uncertainty<numt>&a,
	const abstract_value_with_extended_uncertainty<numt>&b
){
    return value_f<2,numt>([](const numt&x,const numt&y){return atan2(x,y);},a,b);
}
template<class numt>
inline value_f<2,numt> POW(
    	const abstract_value_with_extended_uncertainty<numt>&a,
	const abstract_value_with_extended_uncertainty<numt>&b
){
    return value_f<2,numt>([](const numt&x,const numt&y){return pow(x,y);},a,b);
}
template<class numt>
inline value_f<1,numt> POW(
    	const abstract_value_with_extended_uncertainty<numt>&a,
	const numt&b
){
    return value_f<1,numt>([&b](const numt&x){return pow(x,b);},a);
}
template<class numt>
inline value_f<2,numt> operator+(
    	const abstract_value_with_extended_uncertainty<numt>&a,
	const abstract_value_with_extended_uncertainty<numt>&b
){
    return value_f<2,numt>([](const numt&x,const numt&y){return x+y;},a,b);
}
template<class numt>
inline value_f<2,numt> operator-(
    	const abstract_value_with_extended_uncertainty<numt>&a,
	const abstract_value_with_extended_uncertainty<numt>&b
){
    return value_f<2,numt>([](const numt&x,const numt&y){return x-y;},a,b);
}
template<class numt>
inline value_f<2,numt> operator*(
    	const abstract_value_with_extended_uncertainty<numt>&a,
	const abstract_value_with_extended_uncertainty<numt>&b
){
    return value_f<2,numt>([](const numt&x,const numt&y){return x*y;},a,b);
}
template<class numt>
inline value_f<2,numt> operator/(
    	const abstract_value_with_extended_uncertainty<numt>&a,
	const abstract_value_with_extended_uncertainty<numt>&b
){
    return value_f<2,numt>([](const numt&x,const numt&y){return x/y;},a,b);
}

}
#endif
