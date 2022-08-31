// this file is distributed under
// LGPLv3 license
#ifndef ___________VECTORS_H_____
#	define ___________VECTORS_H_____
#include <tuple>
#include <iostream>
#include <math.h>
#include <type_traits>
#include "error.h"
namespace MathTemplates
{

template<size_t size = 3, class numt = double>class Vector;
template<size_t size = 3, class linetype = Vector<3, double>>class Matrix;
template<size_t size = 3, class numt = double>class Direction;
template<class numt, class... Args>
inline Vector < sizeof...(Args) + 1, numt > vec(const numt &x, Args... args);

template<class numt>
class Vector<1, numt>
{
    template<size_t sizef, class n>friend class Vector;
    template<size_t sizef, class n>friend class Direction;
    template<size_t sizef, class n>friend class Matrix;

public:
    enum {Dimensions = 1};

    typedef numt NumberType;
    typedef Direction<1, numt> DType;
    typedef Vector<2,numt> PlusOneComponent;

private:
    numt m_x;

protected:
    inline const numt &___last_component()const
    {
        return m_x;
    }

public:
    inline Vector(const numt &x)
        : m_x(x) {}
    virtual ~Vector() {}

    inline std::tuple<numt> to_tuple()const
    {
        return std::make_tuple(m_x);
    }

    template<class numt2>
    inline Vector(const Vector<Dimensions, numt2> &source)
        : m_x(source.___last_component()) {}

    template<class... Args>
    inline Vector(const std::tuple<Args...> &v)
        : m_x(std::get<0>(v)) {}

    inline static Vector zero()
    {
        return Vector(std::make_tuple(numt(0)));
    }

    template<size_t index>
    inline auto InsertComponent(const NumberType&c)const
    {
	    static_assert(index>0,"Range check error for insertion position");
	    static_assert(index<=(Dimensions+1),"Range check error for insertion position");
	    if constexpr(index==(Dimensions+1))return vec(m_x,c);
	    if constexpr(index==Dimensions) return vec(c,m_x);
    }

    template<size_t index>
    inline static Vector basis_vector()
    {
	    static_assert(index ==1,"dimension index is out of range");
        return Vector(std::make_tuple(numt(1)));
    }

    template<size_t index>
    inline const auto &component()const
    {
	    static_assert(index == 1,"dimension index is out of range");
        return m_x;
    }

    inline const numt &x()const
    {
        return m_x;
    }
    inline operator numt()const
    {
        return m_x;
    }

    template<size_t index, typename... Args>
    inline auto component(Args... args)const
    {
	    static_assert(index == 1,"dimension index is out of range");
        return m_x(args...);
    }

    template<typename... Args>
    inline auto x(Args... args)const
    {
        return m_x(args...);
    }
    
    Vector &operator=(const Vector &source) = default;
    
    Vector &operator+=(const Vector &second)
    {
        m_x += second.m_x;
        return *this;
    }
    
    Vector operator+(const Vector&second)const
    {
        return Vector(std::make_tuple(m_x + second.___last_component()));
    }
    
    Vector &operator-=(const Vector &second)
    {
        m_x -= second.m_x;
        return *this;
    }
    
    Vector operator-(const Vector&second)const
    {
        return Vector(std::make_tuple(m_x - second.___last_component()));
    }
    
    Vector operator-()const
    {
        return Vector(std::make_tuple(-m_x));
    }
    
    inline const Vector& operator+()const
    {
        return *this;
    }

    Vector &operator*=(const numt &second)
    {
        m_x *= second;
        return *this;
    }
    
    Vector &operator/=(const numt &second)
    {
        m_x /= second;
        return *this;
    }
    
    Vector operator/(const numt &second)const
    {
        return Vector(std::make_tuple(m_x / second));
    }
    
    template<class numt2>
    inline auto operator*(const Vector<1,numt2> &second)const
    {
        return m_x * second.___last_component();
    }
    
    template<class numt2>
    inline auto operator*(const numt2&second)const
    {
        return Vector<1,decltype(x()*second)>(std::make_tuple(x()*second));
    }
    
    template<class numt2>
    inline auto mul_inv(const numt2&second)const
    {
        return Vector<1,decltype(second*x())>(std::make_tuple(second*x()));
    }
    
    template<typename... Args>
    inline auto operator()(Args... args)const{
	    return vec(m_x(args...));
    }
    
    template<typename... Args>
    inline operator numt()const{
	    return m_x;
    }
    
    inline numt length_sqr()const
    {
        return operator*(*this);
    }
    
    inline numt length()const
    {
        return sqrt(length_sqr());
    }
    
    bool operator==(const Vector &second)const
    {
        return (m_x == second.m_x);
    }
    
    inline bool CloseTo(const Vector &second, const numt &epsilon)const
    {
        return operator-(second).length() < epsilon;
    }
    
    inline void output(std::ostream&str)const
    {
        str<<m_x;
    }
};

template<size_t size, class numt>
class Vector
{
    template<size_t sizef, class n>friend class Vector;
    template<size_t sizef, class n>friend class Direction;
    template<size_t sizef, class n>friend class Matrix;

public:
    enum {Dimensions = size};

    typedef numt NumberType;
    typedef Vector < size - 1, numt > MinusOneComponent;
    typedef Direction<size, numt> DType;
    typedef Vector<size+1,numt> PlusOneComponent;

private:
    MinusOneComponent m_other;
    numt m_x;

    inline Vector(const MinusOneComponent &other, const numt &x)
        : m_other(other), m_x(x) {}

protected:
    inline const auto &___minus_one_component()const
    {
        return m_other;
    }

    inline const auto &___last_component()const
    {
        return m_x;
    }

public:
    virtual ~Vector() {}

    inline auto to_tuple()const
    {
        return std::tuple_cat(m_other.to_tuple(), std::make_tuple(m_x));
    }

    template<class numt2>
    inline Vector(const Vector<Dimensions, numt2> &source)
        : m_other(source.___minus_one_component()), m_x(source.___last_component()) {}

    template<class... Args>
    inline Vector(const std::tuple<Args...> &v)
        : m_other(v), m_x(std::get < Dimensions - 1 > (v)) {}
    
    inline static auto zero()
    {
        return Vector(MinusOneComponent::zero(), numt(0));
    }
    
    template<size_t index>
    inline static auto basis_vector()
    {
	    static_assert(index > 0,"dimension index is out of range");
	    static_assert(index <= Dimensions,"dimension index is out of range");
        if constexpr(index == Dimensions) return Vector(MinusOneComponent::zero(), numt(1));
	    else return Vector(MinusOneComponent::template basis_vector<index>(), numt(0));
    }
    
    template<size_t index>
    inline const auto &component()const
    {
	    static_assert(index > 0,"dimension index is out of range");
	    static_assert(index <= Dimensions,"dimension index is out of range");
        if constexpr(index == Dimensions)return m_x;
	    else return m_other.template component<index>();
    }
    
    template<size_t index, typename... Args>
    inline auto component(Args... args)const
    {
	    static_assert(index > 0,"dimension index is out of range");
	    static_assert(index<=Dimensions,"dimension index is out of range");
        if constexpr(index == Dimensions)return m_x(args...);
	    else return m_other.template component<index>(args...);
    }
    
    template<size_t index>
    inline auto RemoveComponent()const
    {
	    static_assert(index > 0,"dimension index is out of range");
	    static_assert(index<=Dimensions,"dimension index is out of range");
        if constexpr(index == Dimensions)return m_other;
        else if constexpr(Dimensions ==2)return MinusOneComponent(m_x);
	    else return MinusOneComponent(m_other.template RemoveComponent<index>(),m_x);
    }

    template<size_t index>
    inline auto InsertComponent(const NumberType&c)const
    {
    	static_assert(index>0,"Range check error for insertion position");
	    static_assert(index<=(Dimensions+1),"Range check error for insertion position");
	    if constexpr(index==(Dimensions+1))return PlusOneComponent(*this,c);
	    if constexpr(index==Dimensions)return PlusOneComponent(Vector(m_other,c),m_x);
	    if constexpr(index<Dimensions)return PlusOneComponent(m_other.template InsertComponent<index>(c),m_x);
    }

    inline const auto &x()const
    {
        return component<1>();
    }

    inline const auto &y()const
    {
        return component<2>();
    }

    inline const auto &z()const
    {
        return component<3>();
    
    }
    template<typename... Args>
    inline auto x(Args... args)const
    {
        return component<1>(args...);
    }

    template<typename... Args>
    inline auto y(Args... args)const
    {
        return component<2>(args...);
    }

    template<typename... Args>
    inline auto z(Args... args)const
    {
        return component<3>(args...);
    }

    Vector &operator=(const Vector &source) = default;

    auto &operator+=(const Vector &second)
    {
        m_x += second.m_x;
        m_other += second.m_other;
        return *this;
    }

    auto operator+(const Vector&second)const
    {
        return Vector(m_other + second.___minus_one_component(), m_x + second.___last_component());
    }

    auto &operator-=(const Vector &second)
    {
        m_x -= second.m_x;
        m_other -= second.m_other;
        return *this;
    }

    auto operator-(const Vector&second)const
    {
        return Vector(m_other - second.___minus_one_component(), m_x - second.___last_component());
    }

    auto operator-()const
    {
        return Vector(-m_other, -m_x);
    }

    inline const auto& operator+()const
    {
        return *this;
    }

    auto &operator*=(const numt &second)
    {
        m_x *= second;
        m_other *= second;
        return *this;
    }

    auto &operator/=(const numt &second)
    {
        m_x /= second;
        m_other /= second;
        return *this;
    }

    auto operator/(const numt &second)const
    {
        return Vector(m_other / second, m_x / second);
    }

    template<class numt2>
    inline auto operator*(const Vector<size,numt2> &second)const
    {
        return (m_other * second.___minus_one_component()) + (m_x * second.___last_component());
    }

    template<class numt2>
    inline auto operator*(const numt2 &second)const
    {
        return Vector<size,decltype(x()*second)>(m_other * second, m_x * second);
    }

    template<class numt2>
    inline auto mul_inv(const numt2 &second)const
    {
        return Vector<size,decltype(second*x())>(m_other.mul_inv(second), second*m_x);
    }

    template<typename... Args>
    inline auto operator()(Args... args)const{
	    return Vector<Dimensions,decltype(m_x(args...))>(m_other(args...),m_x(args...));
    }

    inline auto length_sqr()const
    {
        return operator*(*this);
    }

    inline auto length()const
    {
        return sqrt(length_sqr());
    }

    bool operator==(const Vector &second)const
    {
        return (m_x == second.m_x) && (m_other == second.m_other);
    }

    inline bool CloseTo(const Vector &second, const numt &epsilon)const
    {
        return operator-(second).length() < epsilon;
    }

    inline void output(std::ostream&str)const{
	    m_other.output(str);
	    str<<" "<<m_x;
    }
};

template<size_t size = 3,class numta, class numtb>
inline auto operator*(const numta&a,const Vector<size,numtb>&B)
{
    return B.mul_inv(a);
}

template<class numt, class... Args>
inline Vector < sizeof...(Args) + 1, numt > vec(const numt &x, Args... args)
{
    return Vector < sizeof...(Args) + 1, numt > (std::make_tuple(x, args...));
}

template<size_t i, size_t d, class numt = double>
inline auto axis()
{
    return Vector<d, numt>::template basis_vector<i>();
}

template<class numt = double>
inline auto x()
{
    return Vector<2, numt>::template basis_vector<1>();
}

template<class numt = double>
inline auto y()
{
    return Vector<2, numt>::template basis_vector<2>();
}

template<class numt = double>
inline auto zero()
{
    return Vector<2, numt>::zero();
}

template<class numt = double>
inline auto X()
{
    return Vector<3, numt>::template basis_vector<1>();
}

template<class numt = double>
inline auto Y()
{
    return Vector<3, numt>::template basis_vector<2>();
}

template<class numt = double>
inline auto Z()
{
    return Vector<3, numt>::template basis_vector<3>();
}

template<class numt = double>
inline auto Zero()
{
    return Vector<3, numt>::zero();
}

template<size_t i,class numt>
inline auto& operator<<(std::ostream&str,const Vector<i,numt>&X)
{
    X.output(str);
    return str;
}

};
#endif
