// this file is distributed under
// LGPLv3 license
#ifndef ___________VECTORS_H_____
#	define ___________VECTORS_H_____
#include <tuple>
#include <math.h>
#include <type_traits>
#include "error.h"
#include "randomfunc.h"
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
    inline Vector(const numt &x): m_x(x) {}
public:
    virtual ~Vector() {}
    inline std::tuple<numt> to_tuple()const
    {
        return std::make_tuple(m_x);
    }
    template<class numt2>
    inline Vector(const Vector<Dimensions, numt2> &source): m_x(source.___last_component()) {}
    template<class... Args>
    inline Vector(const std::tuple<Args...> &v): m_x(std::get<0>(v)) {}
    inline static Vector zero()
    {
        return Vector(std::make_tuple(numt(0)));
    }
#ifdef ____full_version_of_math_h_____
    template<size_t index>
    inline PlusOneComponent InsertComponent(const NumberType&c)const
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
    inline const numt &component()const
    {
	static_assert(index == 1,"dimension index is out of range");
        return m_x;
    }
#else
    template<size_t index>
    static Vector basis_vector()
    {
	static_assert(index>0,"dimension index is out of range");
	if(index>1)throw Exception<Vector>("dimension index is out of range");
        return Vector(std::make_tuple(numt(1)));
    }
    template<size_t index>
    const numt &component()const
    {
	static_assert(index>0,"dimension index is out of range");
	if(index>1)throw Exception<Vector>("dimension index is out of range");
        return m_x;
    }
#endif
    inline const numt &x()const
    {
        return m_x;
    }
    Vector &operator=(const Vector &source)
    {
        m_x = source.m_x;
        return *this;
    }
    Vector &operator+=(const Vector &second)
    {
        m_x += second.m_x;
        return *this;
    }
    Vector operator+(const Vector &second)const
    {
        return Vector(std::make_tuple(m_x + second.m_x));
    }
    Vector &operator-=(const Vector &second)
    {
        m_x -= second.m_x;
        return *this;
    }
    Vector operator-(const Vector &second)const
    {
        return Vector(std::make_tuple(m_x - second.m_x));
    }

    Vector &operator*=(const numt &second)
    {
        m_x *= second;
        return *this;
    }
    Vector operator*(const numt &second)const
    {
        return Vector(std::make_tuple(m_x * second));
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
    numt operator*(const Vector &second)const
    {
        return m_x * second.m_x;
    }
    inline numt M_sqr()const
    {
        return operator*(*this);
    }
    inline numt M()const
    {
        return sqrt(M_sqr());
    }
    bool operator==(const Vector &second)const
    {
        return (m_x == second.m_x);
    }
    bool CloseTo(const Vector &second, const numt &epsilon)const
    {
        return operator-(second).M() < epsilon;
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
    inline Vector(const MinusOneComponent &other, const numt &x): m_other(other), m_x(x) {}
protected:
    inline const MinusOneComponent &___minus_one_component()const
    {
        return m_other;
    }
    inline const numt &___last_component()const
    {
        return m_x;
    }
public:
    virtual ~Vector() {}
    inline auto to_tuple()const->decltype(std::tuple_cat(m_other.to_tuple(), std::make_tuple(m_x)))
    {
        return std::tuple_cat(m_other.to_tuple(), std::make_tuple(m_x));
    }
    template<class numt2>
    inline Vector(const Vector<Dimensions, numt2> &source): m_other(source.___minus_one_component()), m_x(source.___last_component()) {}
    template<class... Args>
    inline Vector(const std::tuple<Args...> &v): m_other(v), m_x(std::get < Dimensions - 1 > (v)) {}
    inline static Vector zero()
    {
        return Vector(MinusOneComponent::zero(), numt(0));
    }
#ifdef ____full_version_of_math_h_____
    template<size_t index>
    inline static Vector basis_vector()
    {
	static_assert(index > 0,"dimension index is out of range");
	static_assert(index<=Dimensions,"dimension index is out of range");
        if constexpr(index == Dimensions) return Vector(MinusOneComponent::zero(), numt(1));
	else return Vector(MinusOneComponent::template basis_vector<index>(), numt(0));
    }
    template<size_t index>
    inline const numt &component()const
    {
	static_assert(index > 0,"dimension index is out of range");
	static_assert(index<=Dimensions,"dimension index is out of range");
        if constexpr(index == Dimensions)return m_x;
	else return m_other.template component<index>();
    }
    template<size_t index>
    inline MinusOneComponent RemoveComponent()const
    {
	static_assert(index > 0,"dimension index is out of range");
	static_assert(index<=Dimensions,"dimension index is out of range");
        if constexpr(index == Dimensions)return m_other;
        else if constexpr(Dimensions ==2)return MinusOneComponent(m_x);
	else return MinusOneComponent(m_other.template RemoveComponent<index>(),m_x);
    }
    template<size_t index>
    inline PlusOneComponent InsertComponent(const NumberType&c)const
    {
	static_assert(index>0,"Range check error for insertion position");
	static_assert(index<=(Dimensions+1),"Range check error for insertion position");
	if constexpr(index==(Dimensions+1))return PlusOneComponent(*this,c);
	if constexpr(index==Dimensions)return PlusOneComponent(Vector(m_other,c),m_x);
	if constexpr(index<Dimensions)return PlusOneComponent(m_other.template InsertComponent<index>(c),m_x);
    }
#else
    template<size_t index>
    static Vector basis_vector()
    {
	static_assert(index > 0,"dimension index is out of range");
	if(index>Dimensions)throw Exception<Vector>("dimension index is out of range");
        if(index == Dimensions) return Vector(MinusOneComponent::zero(), numt(1));
	else return Vector(MinusOneComponent::template basis_vector<index>(), numt(0));
    }
    template<size_t index>
    const numt &component()const
    {
	static_assert(index > 0,"dimension index is out of range");
	if(index>Dimensions)throw Exception<Vector>("dimension index is out of range");
        if(index == Dimensions)return m_x;
	else return m_other.template component<index>();
    }
#endif
    inline const numt &x()const
    {
        return component<1>();
    }
    inline const numt &y()const
    {
        return component<2>();
    }
    inline const numt &z()const
    {
        return component<3>();
    }
    Vector &operator+=(const Vector &second)
    {
        m_x += second.m_x;
        m_other += second.m_other;
        return *this;
    }
    Vector operator+(const Vector &second)const
    {
        return Vector(m_other + second.m_other, m_x + second.m_x);
    }
    Vector &operator-=(const Vector &second)
    {
        m_x -= second.m_x;
        m_other -= second.m_other;
        return *this;
    }
    Vector operator-(const Vector &second)const
    {
        return Vector(m_other - second.m_other, m_x - second.m_x);
    }

    Vector &operator*=(const numt &second)
    {
        m_x *= second;
        m_other *= second.m_other;
        return *this;
    }
    Vector operator*(const numt &second)const
    {
        return Vector(m_other * second, m_x * second);
    }
    Vector &operator/=(const numt &second)
    {
        m_x /= second;
        m_other /= second.m_other;
        return *this;
    }
    Vector operator/(const numt &second)const
    {
        return Vector(m_other / second, m_x / second);
    }
    numt operator*(const Vector &second)const
    {
        return (m_other * second.m_other) + (m_x * second.m_x);
    }
    inline numt M_sqr()const
    {
        return operator*(*this);
    }
    inline numt M()const
    {
        return sqrt(M_sqr());
    }
    bool operator==(const Vector &second)const
    {
        return (m_x == second.m_x) && (m_other == second.m_other);
    }
    bool CloseTo(const Vector &second, const numt &epsilon)const
    {
        return operator-(second).M() < epsilon;
    }
};
template<class numt, class... Args>
inline Vector < sizeof...(Args) + 1, numt > vec(const numt &x, Args... args)
{
    return Vector < sizeof...(Args) + 1, numt > (std::make_tuple(x, args...));
}
template<size_t i, size_t d, class numt = double>
inline Vector<d, numt> axis()
{
    return Vector<d, numt>::template basis_vector<i>();
}
template<class numt = double>
inline Vector<2, numt> x()
{
    return Vector<2, numt>::template basis_vector<1>();
}
template<class numt = double>
inline Vector<2, numt> y()
{
    return Vector<2, numt>::template basis_vector<2>();
}
template<class numt = double>
inline Vector<2, numt> zero()
{
    return Vector<2, numt>::zero();
}
template<class numt = double>
inline Vector<3, numt> X()
{
    return Vector<3, numt>::template basis_vector<1>();
}
template<class numt = double>
inline Vector<3, numt> Y()
{
    return Vector<3, numt>::template basis_vector<2>();
}
template<class numt = double>
inline Vector<3, numt> Z()
{
    return Vector<3, numt>::template basis_vector<3>();
}
template<class numt = double>
inline Vector<3, numt> Zero()
{
    return Vector<3, numt>::zero();
}
template<size_t i, class numt = double>
inline Vector<i, numt> operator-(const Vector<i, numt> &V)
{
    return V * numt(-1);
}
};
#endif
