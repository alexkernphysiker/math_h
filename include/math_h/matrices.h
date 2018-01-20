// this file is distributed under
// LGPLv3 license
#ifndef ______MATRICES_H_______
#	define ______MATRICES_H_______
#if __cplusplus<201100L
#error c++>=11 is needed for using math_h headers
#endif
#include <iostream>
#include <math.h>
#include "error.h"
#include "vectors.h"
#if __cplusplus>201700L
#define ____optimized_version_of_matrices_h_____
#else
#warning compiler does not support "if constexpr(...)". c++>=17 is needed. classes from vectors.h will work slower
#endif
namespace MathTemplates
{
template<class linetype>
class Matrix<1, linetype>
{
    template<size_t sizef, class n>friend class Vector;
    template<size_t sizef, class n>friend class Direction;
    template<size_t sizef, class n>friend class Matrix;
public:
    enum {DimensionsFinal = 1};
    enum {DimensionsInitial = linetype::Dimensions};
    typedef typename linetype::NumberType NumberType;
    typedef linetype VIType;
    typedef Vector<1, NumberType> VFType;
    typedef Vector < DimensionsInitial + 1, NumberType > VIWType;
    typedef Matrix<DimensionsFinal, VIWType> PlusOneColumn;
    typedef Vector < DimensionsInitial - 1, NumberType > VINType;
    typedef Matrix<DimensionsFinal, VINType> MinusOneColumn;
private:
    VIType m_line;
protected:
    inline Matrix(const VIType &line): m_line(line) {}
    inline MinusOneColumn ___minus_one_column()const
    {
        return MinusOneColumn(m_line.___recursive());
    }
    inline VFType ___last_column()const
    {
        return desCartes(m_line.___last_component());
    }
    inline PlusOneColumn ___add_column(const VFType &col)const
    {
        return VIWType(m_line, col.___last_component());
    }
    inline const VIType &___last_line()const
    {
        return m_line;
    }
    inline const Matrix &AddColumns()const
    {
        return *this;
    }
public:
    virtual ~Matrix() {}
#ifdef ____optimized_version_of_matrices_h_____
    template<size_t index, size_t jindex>
    inline const NumberType &element()const
    {
	static_assert(index ==1 ,"dimension index is out of range");
        return m_line.template component<jindex>();
    }
#else
    template<size_t index, size_t jindex>
    const NumberType&element()const
    {
	static_assert(index>0, "dimension index is out of range");
	if(index>1)throw Exception<Matrix>("dimension index out of range");
	return m_line.template component<jindex>();
    }
#endif
    VFType operator*(const VIType &v)const
    {
        return desCartes(m_line * v);
    }
    Matrix<DimensionsFinal, Vector<1, NumberType>>
            operator*(const Matrix<DimensionsInitial, Vector<1, NumberType>> &B)const
    {
        return desCartes(m_line * B.___last_column());
    }
    template<size_t third_size>
    Matrix<DimensionsFinal, Vector<third_size, NumberType>>
            operator*(const Matrix<DimensionsInitial, Vector<third_size, NumberType>> &B)const
    {
        return operator*(B.___minus_one_column()).___add_column(operator*(B.___last_column()));
    }
    template<class otherlinetype>
    inline Matrix(const Matrix<DimensionsFinal, otherlinetype> &source): m_line(source.___last_line()) {}
    template<class... Args>
    inline Matrix(const std::tuple<Args...> &v): m_line(std::get < DimensionsFinal - 1 > (v)) {}
    inline Matrix(const VFType &A, const VIType &B): m_line(B *A.___last_component()) {}
    bool operator==(const Matrix &B)const
    {
        return m_line == B.m_line;
    }
    Matrix operator*(const NumberType &v)const
    {
        return Matrix(m_line * v);
    }
    Matrix operator/(const NumberType &v)const
    {
        return Matrix(m_line / v);
    }
    Matrix operator+(const Matrix &B)const
    {
        return Matrix(m_line + B.m_line);
    }
    Matrix operator-(const Matrix &B)const
    {
        return Matrix(m_line - B.m_line);
    }
    static inline Matrix zero()
    {
        return Matrix(VIType::zero());
    }
    static inline Matrix one()
    {
        return Matrix(VIType::template basis_vector<DimensionsFinal>());
    }
    template<size_t x, size_t y>
    static inline Matrix RotationInPlane(const NumberType &angle)
    {
        return Matrix(
                   (x == DimensionsFinal) ? ((VIType::template basis_vector<x>() * cos(angle)) - (VIType::template basis_vector<y>() * sin(angle))) :
                   (y == DimensionsFinal) ? ((VIType::template basis_vector<y>() * cos(angle)) + (VIType::template basis_vector<x>() * sin(angle))) :
                   VIType::template basis_vector<DimensionsFinal>()
                                                );
    }
    template<class... Args>
    inline Matrix < DimensionsFinal, Vector < DimensionsInitial + 1 + sizeof...(Args), NumberType >> AddColumns(const VFType &col, Args...args)const
    {
        return ___add_column(col).AddColumns(args...);
    }
};
template<size_t sizef, class linetype>
class Matrix
{
    template<size_t sizeff, class n>friend class Vector;
    template<size_t sizeff, class n>friend class Direction;
    template<size_t sizeff, class n>friend class Matrix;
public:
    enum {DimensionsFinal = sizef};
    enum {DimensionsInitial = linetype::Dimensions};
    typedef typename linetype::NumberType NumberType;
    typedef linetype VIType;
    typedef Vector<sizef, NumberType> VFType;
    typedef Matrix < sizef - 1, linetype > MinorTransformation;
    typedef Vector < DimensionsInitial + 1, NumberType > VIWType;
    typedef Matrix<DimensionsFinal, VIWType> PlusOneColumn;
    typedef Vector < DimensionsInitial - 1, NumberType > VINType;
    typedef Matrix<DimensionsFinal, VINType> MinusOneColumn;
private:
    MinorTransformation m_minor;
    VIType m_line;
protected:
    inline Matrix(const MinorTransformation &minor, const VIType &line): m_minor(minor), m_line(line) {}
    inline MinusOneColumn ___minus_one_column()const
    {
        const auto new_minor = m_minor.___minus_one_column();
        const auto new_line = m_line.___recursive();
        return MinusOneColumn(new_minor, new_line);
    }
    inline VFType ___last_column()const
    {
        return VFType(m_minor.___last_column(), m_line.___last_component());
    }
    inline PlusOneColumn ___add_column(const VFType &col)const
    {
        const auto new_minor = m_minor.___add_column(col.___recursive());
        const auto new_line = VIWType(m_line, col.___last_component());
        return PlusOneColumn(new_minor, new_line);
    }
    inline const VIType &___last_line()const
    {
        return m_line;
    }
    inline const MinorTransformation &___recursive()const
    {
        return m_minor;
    }
    inline const Matrix &AddColumns()const
    {
        return *this;
    }
public:
    virtual ~Matrix() {}
#ifdef ____optimized_version_of_matrices_h_____
    template<size_t index, size_t jindex>
    inline const NumberType &element()const
    {
	static_assert(index > 0,"dimension index is out of range");
	static_assert(index<=DimensionsFinal,"dimension index is out of range");
	if constexpr(index==DimensionsFinal) return m_line.template component<jindex>();
	else return m_minor.template element<index, jindex>();
    }
#else
    template<size_t index, size_t jindex>
    const NumberType &element()const
    {
	static_assert(index > 0,"dimension index is out of range");
	if(index>DimensionsFinal)throw Exception<Matrix>("dimension index is out of range");
	if(index==DimensionsFinal) return m_line.template component<jindex>();
	else return m_minor.template element<index, jindex>();
    }
#endif
    VFType operator*(const VIType &v)const
    {
        return VFType(m_minor * v, m_line * v);
    }
    Matrix<DimensionsFinal, Vector<1, NumberType>>
            operator*(const Matrix<DimensionsInitial, Vector<1, NumberType>> &B)const
    {
        const auto P = m_minor * B;
        const auto C = desCartes(m_line * B.___last_column());
        return Matrix<DimensionsFinal, Vector<1, NumberType>>(P, C);
    }
    template<size_t third_size>
    Matrix<DimensionsFinal, Vector<third_size, NumberType>>
            operator*(const Matrix<DimensionsInitial, Vector<third_size, NumberType>> &B)const
    {
        const auto P = operator*(B.___minus_one_column());
        const VFType C = operator*(B.___last_column());
        return P.___add_column(C);
    }
    template<class otherlinetype>
    inline Matrix(const Matrix<DimensionsFinal, otherlinetype> &source): m_minor(source.___recursive()), m_line(source.___last_line()) {}
    template<class... Args>
    inline Matrix(const std::tuple<Args...> &v): m_minor(v), m_line(std::get < DimensionsFinal - 1 > (v)) {}
    inline Matrix(const VFType &A, const VIType &B): m_minor(A.___recursive(), B), m_line(B *A.___last_component()) {}
    bool operator==(const Matrix &B)const
    {
        return (m_line == B.m_line) && (m_minor == B.m_minor);
    }
    Matrix operator*(const NumberType &v)const
    {
        return Matrix(m_minor * v, m_line * v);
    }
    Matrix operator/(const NumberType &v)const
    {
        return Matrix(m_minor / v, m_line / v);
    }
    Matrix operator+(const Matrix &B)const
    {
        return Matrix(m_minor + B.m_minor, m_line + B.m_line);
    }
    Matrix operator-(const Matrix &B)const
    {
        return Matrix(m_minor - B.m_minor, m_line - B.m_line);
    }
    static inline Matrix zero()
    {
        return Matrix(MinorTransformation::zero(), VIType::zero());
    }
    static inline Matrix one()
    {
        return Matrix(MinorTransformation::one(), VIType::template basis_vector<DimensionsFinal>());
    }
    template<size_t x, size_t y>
    static inline Matrix RotationInPlane(const NumberType &angle)
    {
        return Matrix(
                   MinorTransformation::template RotationInPlane<x, y>(angle),
                   (x == DimensionsFinal) ? ((VIType::template basis_vector<x>() * cos(angle)) - (VIType::template basis_vector<y>() * sin(angle))) :
                   (y == DimensionsFinal) ? ((VIType::template basis_vector<y>() * cos(angle)) + (VIType::template basis_vector<x>() * sin(angle))) :
                   VIType::template basis_vector<DimensionsFinal>()
                                                );
    }
    template<class... Args>
    inline Matrix < DimensionsFinal, Vector < DimensionsInitial + 1 + sizeof...(Args), NumberType >> AddColumns(const VFType &col, Args...args)const
    {
        return ___add_column(col).AddColumns(args...);
    }
};
template<class linetype, class... Args>
inline Matrix < sizeof...(Args) + 1, linetype > lines(const linetype &x, Args... args)
{
    return Matrix < sizeof...(Args) + 1, linetype > (std::make_tuple(x, args...));
}
template<class numt, class... Args>
inline Matrix < 1, Vector < sizeof...(Args) + 1, numt >> line(const numt &x, Args... args)
{
    return lines(desCartes(x, args...));
}
template<class numt, class...Args>
inline Matrix < 1 + sizeof...(Args), Vector<1, numt >> column(const numt &x, Args...args)
{
    return Matrix < 1 + sizeof...(Args), Vector<1, numt >> (std::make_tuple(x, args...));
}
template<class VecT, class...Args>
inline Matrix < VecT::Dimensions, Vector < 1 + sizeof...(Args), typename VecT::NumberType >> columns(const VecT &x, Args...args)
{
    return Matrix<VecT::Dimensions, Vector<1, typename VecT::NumberType>>(x.to_tuple()).AddColumns(args...);
}
template<size_t s, class numt = double>
inline Matrix<s, Vector<s, numt>> ZERO()
{
    return Matrix<s, Vector<s, numt>>::zero();
}
template<size_t s, class numt = double>
inline Matrix<s, Vector<s, numt>> ONE()
{
    return Matrix<s, Vector<s, numt>>::one();
}
template<size_t size, class VIType>
inline Matrix<size, VIType> TensorProduct(const Vector<size, typename VIType::NumberType> &A, const VIType &B)
{
    return Matrix<size, VIType>(A, B);
}

};
#endif
