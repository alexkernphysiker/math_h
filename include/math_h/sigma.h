// this file is distributed under
// LGPLv3 license
#ifndef YJIOGPKSIIYJVKND
#define YJIOGPKSIIYJVKND
#include <iostream>
#include <memory>
#include <vector>
#include <math.h>
#include "error.h"
namespace MathTemplates
{
template<class comparable = double>
using Chain=std::vector<comparable>;
template<typename numt = double>
class abstract_value_with_uncertainty
{
public:
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
namespace details{
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
	details::get_val([F](auto...a){return F(a...);},v,args...),
	sqrt(details::uncertainty_sqr([F](auto...a){return F(a...);},v,args...))
    );
}
#else
template<class numt,class Func>
inline value<numt> func_with_uncertainty(Func F,const abstract_value_with_uncertainty<numt>&v)
{
    return value<numt>(
	details::get_val([F](const numt&a){return F(a);},v),
	sqrt(details::uncertainty_sqr([F](const numt&a){return F(a);},v))
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
    Chain<numt> m_list;
    numt m_sum;
    std::shared_ptr<value<numt>>m_cache;
public:
    typedef typename Chain<numt>::const_iterator const_iterator;
    const_iterator begin()const{return m_list.begin();}
    const_iterator end() const{return m_list.end();}
    const numt&operator[](const size_t index){
	if(index>=m_list.size())throw Exception<StandardDeviation>("range check error");
	return m_list[index];
    }
    inline size_t size()const
    {
        return m_list.size();
    }

    StandardDeviation():m_sum(0){}
    virtual ~StandardDeviation() {}
    StandardDeviation &operator<<(const numt &x)
    {
        m_list.push_back(x);
        m_sum += x;
        m_cache = nullptr;
        return *this;
    }
    inline StandardDeviation(const numt &x):StandardDeviation()
    {
        operator<<(x);
    }
    inline size_t count()const
    {
        return m_list.size();
    }
private:
    const value<numt> &VAL()const
    {
        using namespace std;
        if (!m_cache) {
            size_t sz = m_list.size();
            if (sz <= 1)
                throw Exception<StandardDeviation>("No data to check. for sigma needed at least two elements.");
            numt average = m_sum / sz;
            numt m_sigsqr = 0;
            for (auto value : m_list)
                m_sigsqr += pow(value - average, 2);
            m_sigsqr /= sz - 1;
            const_cast<StandardDeviation &>(*this).m_cache =
                make_shared<value<numt>>(average, sqrt(m_sigsqr));
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

template<typename numt = double>
class CorrelationLinear
{
private:
    StandardDeviation<numt> _X, _Y, _XY;
public:
    CorrelationLinear(const numt &scale = 1):
        _X(scale), _Y(scale), _XY(scale) {}
    virtual ~CorrelationLinear() {}
    CorrelationLinear &operator<<(const std::pair<numt, numt> &P)
    {
        _X << P.first;
        _Y << P.second;
        _XY << P.first *P.second;
        return *this;
    }
    inline numt Covariance()const
    {
        return _XY.val() - _X.val() * _Y.val();
    }
    inline numt R()const
    {
        return Covariance() / (_X.uncertainty() * _Y.uncertainty());
    }
    const StandardDeviation<numt>&X()const{return _X;}
    const StandardDeviation<numt>&Y()const{return _Y;}
};
};
#endif
