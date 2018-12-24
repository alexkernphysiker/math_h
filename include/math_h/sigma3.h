// this file is distributed under
// LGPLv3 license
#ifndef _____EXTENDED_UNCERTAINTIES_ESTIMATION_MORE_THAN_ONE_UNCERTAINTY_____________
#define _____EXTENDED_UNCERTAINTIES_ESTIMATION_MORE_THAN_ONE_UNCERTAINTY_____________
#include <memory>
#include <math.h>
#include "error.h"
#include "sigma.h"
#include "tabledata.h"
namespace MathTemplates
{
template<size_t size, class numt = double>class Uncertainties;
template<class numt>
class Uncertainties<1, numt>
{
    template<size_t sizef, class n>friend class Uncertainties;
public:
    enum {UncertaintiesCount = 1};
    typedef numt NumberType;
private:
    numt m_value,m_uncertainty;
    bool m_use_maximum_estimation;
public:
    virtual ~Uncertainties() {}
    inline std::tuple<numt> to_tuple()const
    {
        return std::make_tuple(m_value,m_uncertainty);
    }
    inline Uncertainties():m_value(0),m_uncertainty(0),m_use_maximum_estimation(false){}
    inline Uncertainties(const numt&val):m_value(val),m_uncertainty(0),m_use_maximum_estimation(false){}
    template<class... Args>
    inline Uncertainties(const std::tuple<Args...> &v): m_value(std::get<0>(v)),m_uncertainty(std::get<1>(v)),m_use_maximum_estimation(false){}
    template<class Arg>
    inline Uncertainties(const Arg&arg):m_value(arg),m_uncertainty(0),m_use_maximum_estimation(false){}
    template<class Arg,class Arg2>
    inline Uncertainties(const Arg&arg,const Arg2&arg2):m_value(arg),m_uncertainty(arg2),m_use_maximum_estimation(false){}
    template<class numt2>
    inline Uncertainties(const Uncertainties<UncertaintiesCount,numt2>&V):m_value(V.m_value)
	    ,m_uncertainty(V.m_uncertainty),m_use_maximum_estimation(V.m_use_maximum_estimation){}
    template<size_t index>
    inline const NumberType&uncertainty()const
    {
	static_assert(index == 1,"uncertainty index is out of range");
        return m_uncertainty;
    }
    template<size_t index>
    inline bool using_maximum_estimation()const
    {
	static_assert(index == 1,"uncertainty index is out of range");
        return m_use_maximum_estimation;
    }
    template<size_t index>
    inline void use_maximum_estimation(bool v=true)
    {
	static_assert(index == 1,"uncertainty index is out of range");
        m_use_maximum_estimation=v;
    }
    inline const NumberType&val()const{return m_value;}
    template<size_t index>
    inline value<NumberType> take_uncertainty_component()const
    {
       return value<NumberType>(val(),uncertainty<index>());
    }
    inline value<NumberType> wrap()const{return value<NumberType>(m_value,m_uncertainty);}
    //number-like comparing
    inline bool operator<(const Uncertainties &other)const{return val() < other.val();}
    inline bool operator>(const Uncertainties &other)const{return val() > other.val();}
    inline bool operator==(const Uncertainties&other)const{return val() == other.val();}
    inline bool operator>=(const Uncertainties&other)const{return val() >= other.val();}
    inline bool operator<=(const Uncertainties&other)const{return val() <= other.val();}
    inline bool operator<(const NumberType&other)const{return val() < other;}
    inline bool operator>(const NumberType&other)const{return val() > other;}
    inline bool operator==(const NumberType&other)const{return val() == other;}
    inline bool operator>=(const NumberType&other)const{return val() >= other;}
    inline bool operator<=(const NumberType&other)const{return val() <= other;}
    //arithmetic actions
    Uncertainties &operator+=(const Uncertainties&other)
    {
	if(m_use_maximum_estimation)m_uncertainty += other.m_uncertainty;
        else m_uncertainty = sqrt(pow(m_uncertainty, 2) + pow(other.m_uncertainty, 2));
        m_value += other.m_value;
        return *this;
    }
    Uncertainties&operator-=(const Uncertainties&other)
    {
	if(m_use_maximum_estimation)m_uncertainty += other.m_uncertainty;
        else m_uncertainty = sqrt(pow(m_uncertainty, 2) + pow(other.m_uncertainty, 2));
        m_value -= other.m_value;
        return *this;
    }
    Uncertainties&operator*=(const Uncertainties&other)
    {
	if(m_use_maximum_estimation)m_uncertainty = ((m_uncertainty * other.val()) + (other.m_uncertainty * val()));
        else m_uncertainty = sqrt(pow(m_uncertainty * other.val(), 2) + pow(other.m_uncertainty * val(), 2));
        m_value *= other.m_value;
        return *this;
    }
    Uncertainties&operator/=(const Uncertainties&other)
    {
	if(m_use_maximum_estimation)m_uncertainty = ((m_uncertainty / other.val()) + (other.m_uncertainty * val() / pow(other.val(), 2)));
        else m_uncertainty = sqrt(pow(m_uncertainty / other.val(), 2) + pow(other.m_uncertainty * val() / pow(other.val(), 2), 2));
        m_value /= other.m_value;
        return *this;
    }
    //extending arithmetic actions
    inline Uncertainties&operator+=(const NumberType&v){return operator+=(Uncertainties(v));}
    inline Uncertainties&operator-=(const NumberType&v){return operator-=(Uncertainties(v));}
    inline Uncertainties&operator*=(const NumberType&v){return operator*=(Uncertainties(v));}
    inline Uncertainties&operator/=(const NumberType&v){return operator/=(Uncertainties(v));}
    inline Uncertainties operator+(const Uncertainties&other)const{return Uncertainties(*this) += other;}
    inline Uncertainties operator-(const Uncertainties&other)const{return Uncertainties(*this) -= other;}
    inline Uncertainties operator*(const Uncertainties&other)const{return Uncertainties(*this) *= other;}
    inline Uncertainties operator/(const Uncertainties&other)const{return Uncertainties(*this) /= other;}
    inline Uncertainties operator+(const NumberType&other)const{return Uncertainties(*this) += other;}
    inline Uncertainties operator-(const NumberType&other)const{return Uncertainties(*this) -= other;}
    inline Uncertainties operator*(const NumberType&other)const{return Uncertainties(*this) *= other;}
    inline Uncertainties operator/(const NumberType&other)const{return Uncertainties(*this) /= other;}
    //IO
    inline void output(std::ostream&str)const{
	str<<m_value<<" "<<m_uncertainty;
    }
    inline void input(std::istream&str){
	str>>m_value>>m_uncertainty;
    }
};
template<size_t size,class numt>
class Uncertainties
{
    template<size_t sizef, class n>friend class Uncertainties;
public:
    enum {UncertaintiesCount = size};
    typedef numt NumberType;
private:
    Uncertainties<UncertaintiesCount-1,NumberType> m_lesser;
    numt m_uncertainty;
    bool m_use_maximum_estimation;
public:
    virtual ~Uncertainties() {}
    inline auto to_tuple()const->decltype(std::tuple_cat(m_lesser.to_tuple(), std::make_tuple(m_uncertainty)))
    {
        return std::tuple_cat(m_lesser.to_tuple(), std::make_tuple(m_uncertainty));
    }
    inline Uncertainties():m_lesser(),m_uncertainty(0),m_use_maximum_estimation(false){}
    inline Uncertainties(const numt&val):m_lesser(val),m_uncertainty(0),m_use_maximum_estimation(false){}
    template<class... Args>
    inline Uncertainties(const std::tuple<Args...> &v): m_lesser(v),m_uncertainty(std::get<UncertaintiesCount>(v)),m_use_maximum_estimation(false) {}
    template<class Arg>
    inline Uncertainties(const Arg&arg):m_lesser(arg),m_uncertainty(0),m_use_maximum_estimation(false){}
    template<class Arg,class Arg2>
    inline Uncertainties(const Arg&arg,const Arg2&arg2):m_lesser(arg),m_uncertainty(arg2),m_use_maximum_estimation(false){}
    template<class numt2>
    inline Uncertainties(const Uncertainties<UncertaintiesCount,numt2>&V):m_lesser(V.m_lesser),
	    m_uncertainty(V.m_uncertainty),m_use_maximum_estimation(V.m_use_maximum_estimation){}
    template<size_t index>
    inline const NumberType&uncertainty()const
    {
	static_assert(index > 0,"uncertainty index is out of range");
	static_assert(index<=UncertaintiesCount,"uncertainty index is out of range");
        if constexpr(index == UncertaintiesCount)return m_uncertainty;
	else return m_lesser.template uncertainty<index>();
    }
    template<size_t index>
    inline bool using_maximum_estimation()const
    {
	static_assert(index > 0,"uncertainty index is out of range");
	static_assert(index<=UncertaintiesCount,"uncertainty index is out of range");
        if constexpr(index == UncertaintiesCount)return m_use_maximum_estimation;
	else return m_lesser.template using_maximum_estimation<index>();
    }
    template<size_t index>
    inline void use_maximum_estimation(bool v=true)
    {
	static_assert(index > 0,"uncertainty index is out of range");
	static_assert(index<=UncertaintiesCount,"uncertainty index is out of range");
        if constexpr(index == UncertaintiesCount)m_use_maximum_estimation=v;
	else m_lesser.template use_maximum_estimation<index>(v);
    }
    inline const NumberType&val()const{return m_lesser.val();}
    template<size_t index>
    inline value<NumberType> take_uncertainty_component()const
    {
       return value<NumberType>(val(),uncertainty<index>());
    }
    inline value<NumberType> wrap()const{
	return value<NumberType>(val(),sqrt(pow(m_uncertainty,2)+pow(m_lesser.wrap().uncertainty(),2)));
    }
    //number-like comparing
    inline bool operator<(const Uncertainties &other)const{return val() < other.val();}
    inline bool operator>(const Uncertainties &other)const{return val() > other.val();}
    inline bool operator==(const Uncertainties&other)const{return val() == other.val();}
    inline bool operator>=(const Uncertainties&other)const{return val() >= other.val();}
    inline bool operator<=(const Uncertainties&other)const{return val() <= other.val();}
    inline bool operator<(const NumberType&other)const{return val() < other;}
    inline bool operator>(const NumberType&other)const{return val() > other;}
    inline bool operator==(const NumberType&other)const{return val() == other;}
    inline bool operator>=(const NumberType&other)const{return val() >= other;}
    inline bool operator<=(const NumberType&other)const{return val() <= other;}
    //arithmetic actions
    Uncertainties &operator+=(const Uncertainties&other)
    {
	if(m_use_maximum_estimation)m_uncertainty += other.m_uncertainty;
        else m_uncertainty = sqrt(pow(m_uncertainty, 2) + pow(other.m_uncertainty, 2));
        m_lesser += other.m_lesser;
        return *this;
    }
    Uncertainties&operator-=(const Uncertainties&other)
    {
	if(m_use_maximum_estimation)m_uncertainty += other.m_uncertainty;
        else m_uncertainty = sqrt(pow(m_uncertainty, 2) + pow(other.m_uncertainty, 2));
        m_lesser -= other.m_lesser;
        return *this;
    }
    Uncertainties&operator*=(const Uncertainties&other)
    {
	if(m_use_maximum_estimation)m_uncertainty = ((m_uncertainty * other.val()) + (other.m_uncertainty * val()));
        else m_uncertainty = sqrt(pow(m_uncertainty * other.val(), 2) + pow(other.m_uncertainty * val(), 2));
        m_lesser *= other.m_lesser;
        return *this;
    }
    Uncertainties&operator/=(const Uncertainties&other)
    {
	if(m_use_maximum_estimation)m_uncertainty = ((m_uncertainty / other.val()) + (other.m_uncertainty * val() / pow(other.val(), 2)));
        else m_uncertainty = sqrt(pow(m_uncertainty / other.val(), 2) + pow(other.m_uncertainty * val() / pow(other.val(), 2), 2));
        m_lesser /= other.m_lesser;
        return *this;
    }
    //extending arithmetic actions
    inline Uncertainties&operator+=(const NumberType&v){return operator+=(Uncertainties(v));}
    inline Uncertainties&operator-=(const NumberType&v){return operator-=(Uncertainties(v));}
    inline Uncertainties&operator*=(const NumberType&v){return operator*=(Uncertainties(v));}
    inline Uncertainties&operator/=(const NumberType&v){return operator/=(Uncertainties(v));}
    inline Uncertainties operator+(const Uncertainties&other)const{return Uncertainties(*this) += other;}
    inline Uncertainties operator-(const Uncertainties&other)const{return Uncertainties(*this) -= other;}
    inline Uncertainties operator*(const Uncertainties&other)const{return Uncertainties(*this) *= other;}
    inline Uncertainties operator/(const Uncertainties&other)const{return Uncertainties(*this) /= other;}
    inline Uncertainties operator+(const NumberType&other)const{return Uncertainties(*this) += other;}
    inline Uncertainties operator-(const NumberType&other)const{return Uncertainties(*this) -= other;}
    inline Uncertainties operator*(const NumberType&other)const{return Uncertainties(*this) *= other;}
    inline Uncertainties operator/(const NumberType&other)const{return Uncertainties(*this) /= other;}
    //IO
    inline void output(std::ostream&str)const{
	m_lesser.output(str);
	str<<" "<<m_uncertainty;
    }
    inline void input(std::istream&str){
	m_lesser.input(str);
	str>>m_uncertainty;
    }

};
template<class numt, class... Args>
inline Uncertainties < sizeof...(Args), numt > uncertainties(const numt &x, Args... args)
{
    return Uncertainties< sizeof...(Args), numt > (std::make_tuple(x, args...));
}
template<size_t index,size_t sz=index,class numt=double>
inline Uncertainties<sz,numt> extend_value(const value<numt>&source)
{
    static_assert(sz>=index,"uncertainty indexing error");
    return Uncertainties<sz,numt>(Uncertainties<index,numt>(source.val(),source.uncertainty()));
}

template<size_t i,class numt>
inline std::ostream& operator<<(std::ostream&str,const Uncertainties<i,numt>&X)
{
    X.output(str);
    return str;
}
template<size_t i,class numt>
inline std::istream& operator>>(std::istream&str,Uncertainties<i,numt>&X)
{
    X.input(str);
    return str;
}

template<size_t sz,class numtX = double, class numtY = numtX>
using ext_hist=SortedPoints<value<numtX>,Uncertainties<sz,numtY>>;

template<size_t index,size_t sz,class numtX, class numtY>
ext_hist<sz,numtX,numtY> extend_hist(const hist<numtX,numtY>&source){
    static_assert(sz>=index,"uncertainty indexing error");
    ext_hist<sz,numtX,numtY> res;
    for(const auto&p:source)res<<make_point(p.X(),extend_value<index,sz,numtY>(p.Y()));
    return res;
}
template<size_t sz,class numtX, class numtY>
ext_hist<sz+1,numtX,numtY> add_one_uncertainty(const ext_hist<sz,numtX,numtY>&source){
    ext_hist<sz+1,numtX,numtY> res;
    for(const auto&p:source)res<<make_point(p.X(),Uncertainties<sz+1>(p.Y(),numtY(0)));
    return res;
}

template<size_t sz,class numtX, class numtY>
hist<numtX,numtY> wrap_hist(const SortedPoints<value<numtX>,Uncertainties<sz,numtY>>&source){
    hist<numtX,numtY> res;
    for(const auto&p:source)res<<make_point(p.X(),p.Y().wrap());
    return res;
}
template<size_t sz,class numt>
inline value<numt> wrap_value(const Uncertainties<sz,numt>&source){
    return source.wrap();
}
template<size_t index,size_t sz,class numt>
inline value<numt> take_uncertainty_component(const Uncertainties<sz,numt>&source){
    return source.template take_uncertainty_component<index>();
}
template<size_t index,size_t sz,class numtX, class numtY>
hist<numtX,numtY> take_uncertainty_component(const SortedPoints<value<numtX>,Uncertainties<sz,numtY>>&source){
    hist<numtX,numtY> res;
    for(const auto&p:source)res<<make_point(p.X(),p.Y().template take_uncertainty_component<index>());
    return res;
}

};
#endif
