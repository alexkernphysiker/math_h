// this file is distributed under
// LGPLv3 license
#ifndef YJIOGPKSIIYJVKND
#define YJIOGPKSIIYJVKND
#include <iostream>
#include <memory>
#include <vector>
#include <math.h>
#include "error.h"
#include "statistics.h"
namespace MathTemplates
{
template<typename numt = double>
class abstract_value_with_uncertainty
{
public:
    typedef numt NumberType;
    virtual ~abstract_value_with_uncertainty(){}
    virtual const numt&val()const=0;
    virtual const numt&uncertainty()const=0;
    inline numt epsilon()const
    {
        return uncertainty()/val();
    }
    inline numt min()const
    {
        return val()-uncertainty();
    }
    inline numt max()const
    {
        return val()+uncertainty();
    }

    //Physical comparing of magnitudes with uncertainties
    inline bool Contains(const numt &x)const
    {
        return (x >= min()) && (x <= max());
    }
    inline bool Contains(const abstract_value_with_uncertainty &x)const
    {
        return (x.max() >= min()) && (x.min() <= max());
    }
    inline bool NotEqual(const numt &x)const
    {
        return (x < min()) || (x > max());
    }
    inline bool NotEqual(const abstract_value_with_uncertainty &x)const
    {
        return (x.max() < min()) || (x.min() > max());
    }
    inline bool Below(const numt &x)const
    {
        return max() < x;
    }
    inline bool Below(const abstract_value_with_uncertainty &x)const
    {
        return max() < x.min();
    }
    inline bool Above(const numt &x)const
    {
        return min() > x;
    }
    inline bool Above(const abstract_value_with_uncertainty &x)const
    {
        return min() > x.max();
    }
    //chi-square-like numeric comparing of magnitudes
    numt NumCompare(const numt &x)const
    {
        return pow((val() - x) / uncertainty(), 2);
    }
    numt NumCompare(const abstract_value_with_uncertainty &x)const
    {
        return pow((val() - x.val()) / (uncertainty() + x.uncertainty()), 2);
    }
    //Inheriting number-like comparing
    inline bool operator<(const abstract_value_with_uncertainty &other)const
    {
        return val() < other.val();
    }
    inline bool operator>(const abstract_value_with_uncertainty &other)const
    {
        return val() > other.val();
    }
    inline bool operator<(const numt &other)const
    {
        return val() < other;
    }
    inline bool operator>(const numt &other)const
    {
        return val() > other;
    }
    inline bool operator==(const abstract_value_with_uncertainty&other)const
    {
        return val() == other.val();
    }
    inline bool operator>=(const abstract_value_with_uncertainty&other)const
    {
        return val() >= other.val();
    }
    inline bool operator<=(const abstract_value_with_uncertainty &other)const
    {
        return val() <= other.val();
    }
    inline bool operator==(const numt&other)const
    {
        return val() == other;
    }
    inline bool operator>=(const numt&other)const
    {
        return val() >= other;
    }
    inline bool operator<=(const numt &other)const
    {
        return val() <= other;
    }
};
template<typename numt>
inline std::ostream &operator<<(std::ostream &str, const abstract_value_with_uncertainty<numt> &P)
{
    return str << P.val() << " " << P.uncertainty();
}

template<typename numt = double>
class value:public abstract_value_with_uncertainty<numt>
{
private:
    numt m_val, m_uncertainty;
public:
    virtual ~value() {}
    virtual const numt&val()const override
    {
        return m_val;
    }
    virtual const numt&uncertainty()const override
    {
        return m_uncertainty;
    }
    value(): m_val(0), m_uncertainty(0){}
    value(const numt&v): m_val(v), m_uncertainty(0){}
    value(const numt&v, const numt&err): m_val(v), m_uncertainty(err){
        if (m_uncertainty < 0)
            m_uncertainty = INFINITY;
    }
    template<class numt2>
    value(const abstract_value_with_uncertainty<numt2> &source): value(source.val(),source.uncertainty()){}

    value make_wider(const numt &scale)const
    {
        return value(m_val, scale * m_uncertainty);
    }
    //arithmetic actions
    value &operator+=(const abstract_value_with_uncertainty<numt> &other)
    {
        m_uncertainty = sqrt(pow(m_uncertainty, 2) + pow(other.uncertainty(), 2));
        m_val += other.val();
        return *this;
    }
    value &operator-=(const abstract_value_with_uncertainty<numt> &other)
    {
        m_uncertainty = sqrt(pow(m_uncertainty, 2) + pow(other.uncertainty(), 2));
        m_val -= other.val();
        return *this;
    }
    value &operator*=(const abstract_value_with_uncertainty<numt> &other)
    {
        m_uncertainty = sqrt(pow(m_uncertainty * other.val(), 2) + pow(other.uncertainty() * m_val, 2));
        m_val *= other.val();
        return *this;
    }
    value &operator/=(const abstract_value_with_uncertainty<numt>&other)
    {
        m_uncertainty = sqrt(pow(m_uncertainty / other.val(), 2)
                     + pow(other.uncertainty() * m_val / pow(other.val(), 2), 2));
        m_val /= other.val();
        return *this;
    }
    //extending arithmetic actions
    inline value &operator+=(const numt&v){return operator+=(value(v));}
    inline value &operator-=(const numt&v){return operator-=(value(v));}
    inline value &operator*=(const numt&v){return operator*=(value(v));}
    inline value &operator/=(const numt&v){return operator/=(value(v));}
    inline value operator+(const abstract_value_with_uncertainty<numt>&other)const{return value(*this) += other;}
    inline value operator-(const abstract_value_with_uncertainty<numt>&other)const{return value(*this) -= other;}
    inline value operator*(const abstract_value_with_uncertainty<numt>&other)const{return value(*this) *= other;}
    inline value operator/(const abstract_value_with_uncertainty<numt>&other)const{return value(*this) /= other;}
    inline value operator+(const numt&other)const{return value(*this) += other;}
    inline value operator-(const numt&other)const{return value(*this) -= other;}
    inline value operator*(const numt&other)const{return value(*this) *= other;}
    inline value operator/(const numt&other)const{return value(*this) /= other;}
};

template<class numt>
inline value<numt> std_error(const numt &v)
{
        if (v < 0)throw Exception<value<numt>>("Cannot calculate std error for negative value");
        if (v == 0)return value<numt>(0,1);
        return value<numt>(v, sqrt(v));
}
template<class numt>
inline value<numt> value_in_range(const numt &a,const numt &b)
{
        if (b < a)throw Exception<value<numt>>("Cannot use empty range for value with uncertainty");
        return value<numt>((a+b)/numt(2),(b-a)/numt(2));
}
namespace sigma_details{
    template<class numt,class Func>
    inline numt get_val(Func F,const abstract_value_with_uncertainty<numt>&v)
    {
	return F(v.val());
    }
    template<class numt,class Func>
    inline numt get_val_u(Func F,const abstract_value_with_uncertainty<numt>&v)
    {
	return F(v.max());
    }
    template<class numt,class Func>
    inline numt get_val_d(Func F,const abstract_value_with_uncertainty<numt>&v)
    {
	return F(v.min());
    }
    template<class numt,class Func>
    inline numt uncertainty_sqr(Func F,const abstract_value_with_uncertainty<numt>&v)
    {
	return (
	    pow(get_val_u(F,v)-get_val(F,v),2)+
	    pow(get_val_d(F,v)-get_val(F,v),2)
	)/numt(2);
    }
#ifdef ____middle_version_of_math_h_____
    template<class numt,class Func,typename... Args>
    inline numt get_val(Func F,const abstract_value_with_uncertainty<numt> &v,Args... args)
    {
	return get_val([&v,&F](auto... a){return F(v.val(),a...);},args...);
    }
    template<class numt,class Func,typename... Args>
    inline numt get_val_u(Func F,const abstract_value_with_uncertainty<numt> &v,Args... args)
    {
	return get_val([&v,&F](auto... a){return F(v.max(),a...);},args...);
    }
    template<class numt,class Func,typename... Args>
    inline numt get_val_d(Func F,const abstract_value_with_uncertainty<numt> &v,Args... args)
    {
	return get_val([&v,&F](auto... a){return F(v.min(),a...);},args...);
    }
    template<class numt,class Func,typename... Args>
    inline numt uncertainty_sqr(Func F,const abstract_value_with_uncertainty<numt>&v,Args... args)
    {
	return (
	    pow(get_val_u(F,v,args...)-get_val(F,v,args...),2)+
	    pow(get_val_d(F,v,args...)-get_val(F,v,args...),2)
	)/numt(2)
	    +uncertainty_sqr([&v,&F](auto... a){return F(v.val(),a...);},args...);
    }
#endif
};
#ifdef ____middle_version_of_math_h_____
template<class numt,class Func,typename... Args>
inline value<numt> func_with_uncertainty(Func F,const abstract_value_with_uncertainty<numt>&v,Args... args)
{
    return value<numt>(
	sigma_details::get_val([F](auto...a){return F(a...);},v,args...),
	sqrt(sigma_details::uncertainty_sqr([F](auto...a){return F(a...);},v,args...))
    );
}
#else
template<class numt,class Func>
inline value<numt> func_with_uncertainty(Func F,const abstract_value_with_uncertainty<numt>&v)
{
    return value<numt>(
	sigma_details::get_val([F](const numt&a){return F(a);},v),
	sqrt(sigma_details::uncertainty_sqr([F](const numt&a){return F(a);},v))
    );
}
#endif
template<typename numt>
inline std::istream &operator>>(std::istream &str, value<numt> &P)
{
    numt v, u;
    str >> v >> u;
    P = value<numt>(v, u);
    return str;
}

template<typename numt = double>
class StandardDeviation:public abstract_value_with_uncertainty<numt>
{
private:
    Sampling<1,numt> m_sample;
    std::pair<size_t,numt> m_cache;
public:
    StandardDeviation():m_cache(0,0){}
    virtual ~StandardDeviation() {}
    inline StandardDeviation(const numt &x):StandardDeviation(){m_sample.Fill(x);}
    inline StandardDeviation&operator<<(const numt &x)
    {
        m_sample.Fill(x);
        return *this;
    }
    inline const Sampling<1,numt>&Sample()const{return m_sample;}
    virtual const numt&val()const override
    {
        return m_sample.Average().x();
    }
    virtual const numt&uncertainty()const override
    {
	const size_t sz = m_sample.count();
	if(sz<2)throw Exception<StandardDeviation>("Cannot obtain standard deviation for sample less than 2 elements");
        if (m_cache.first!=sz) {
	    const_cast<StandardDeviation&>(*this).m_cache =
                std::make_pair(sz,sqrt(m_sample.Cov().template element<1,1>()));
	}
        return m_cache.second;
    }
    //extending arithmetic actions
    inline value<numt> operator+(const abstract_value_with_uncertainty<numt>&other)const{return value<numt>(*this) += other;}
    inline value<numt> operator-(const abstract_value_with_uncertainty<numt>&other)const{return value<numt>(*this) -= other;}
    inline value<numt> operator*(const abstract_value_with_uncertainty<numt>&other)const{return value<numt>(*this) *= other;}
    inline value<numt> operator/(const abstract_value_with_uncertainty<numt>&other)const{return value<numt>(*this) /= other;}
    inline value<numt> operator+(const numt&other)const{return value<numt>(*this) += other;}
    inline value<numt> operator-(const numt&other)const{return value<numt>(*this) -= other;}
    inline value<numt> operator*(const numt&other)const{return value<numt>(*this) *= other;}
    inline value<numt> operator/(const numt&other)const{return value<numt>(*this) /= other;}
};

template<typename numt = double>
class WeightedAverage:public abstract_value_with_uncertainty<numt>
{
private:
    numt Sum;
    numt Wnorm;
    std::shared_ptr<value<numt>>m_cache;
public:
    WeightedAverage()
    {
        Sum = 0;
        Wnorm = 0;
    }
    virtual ~WeightedAverage() {}
    WeightedAverage &operator<<(const abstract_value_with_uncertainty<numt> &X)
    {
        if (X.uncertainty() == 0)
            throw Exception<WeightedAverage>("Cannot add value with zero error");
        numt w = 1.0 / pow(X.uncertainty(), 2);
        Sum += w * X.val();
        Wnorm += w;
        m_cache = nullptr;
        return *this;
    }
    inline WeightedAverage(const abstract_value_with_uncertainty<numt> &X):WeightedAverage()
    {
	operator<<(X);
    }
private:
    const value<numt> &VAL()const
    {
        using namespace std;
        if (!m_cache) {
            if (Wnorm <= 0)throw Exception<WeightedAverage>("Attempt to check empty data");
            const_cast<WeightedAverage &>(*this).m_cache = make_shared<value<numt>>(Sum / Wnorm, 1.0 / sqrt(Wnorm));
        }
        return *m_cache;
    }
public:
    virtual const numt&val()const override
    {
        return VAL().val();
    }
    virtual const numt&uncertainty()const override
    {
        return VAL().uncertainty();
    }
    //extending arithmetic actions
    inline value<numt> operator+(const abstract_value_with_uncertainty<numt>&other)const{return value<numt>(*this) += other;}
    inline value<numt> operator-(const abstract_value_with_uncertainty<numt>&other)const{return value<numt>(*this) -= other;}
    inline value<numt> operator*(const abstract_value_with_uncertainty<numt>&other)const{return value<numt>(*this) *= other;}
    inline value<numt> operator/(const abstract_value_with_uncertainty<numt>&other)const{return value<numt>(*this) /= other;}
    inline value<numt> operator+(const numt&other)const{return value<numt>(*this) += other;}
    inline value<numt> operator-(const numt&other)const{return value<numt>(*this) -= other;}
    inline value<numt> operator*(const numt&other)const{return value<numt>(*this) *= other;}
    inline value<numt> operator/(const numt&other)const{return value<numt>(*this) /= other;}
};
};
#endif
