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
#warning compiler does not support "if constexpr(...)". c++>=17 is needed. Some features of matrices.h will be absent
#endif
namespace MathTemplates
{
template<class linetype, class... Args>
inline Matrix < sizeof...(Args) + 1, linetype > lines(const linetype &x, Args... args);
template<class linetype>
class Matrix<1, linetype>
{
    template<size_t sizef, class n>friend class Vector;
    template<size_t sizef, class n>friend class Direction;
    template<size_t sizef, class n>friend class Matrix;
public:
    enum {RowsCount = 1};
    enum {ColumnsCount = linetype::Dimensions};
    typedef typename linetype::NumberType NumberType;
    typedef linetype LineType;
    typedef Vector<1, NumberType> ColumnType;
    typedef Vector < ColumnsCount + 1, NumberType > LongerLine;
    typedef Matrix<RowsCount, LongerLine> PlusOneColumn;
    typedef Vector < ColumnsCount - 1, NumberType > ShorterLine;
    typedef Matrix<RowsCount, ShorterLine> MinusOneColumn;
    typedef Matrix<RowsCount+1, LineType> PlusOneRow;
private:
    LineType m_line;
protected:
    inline Matrix(const LineType &line): m_line(line) {}
    inline MinusOneColumn ___minus_one_column()const
    {
        return MinusOneColumn(m_line.___recursive());
    }
    inline ColumnType ___last_column()const
    {
        return desCartes(m_line.___last_component());
    }
    inline PlusOneColumn ___add_column(const ColumnType &col)const
    {
        return LongerLine(m_line, col.___last_component());
    }
    inline const LineType &___last_line()const
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
    template<size_t index>
    inline MinusOneColumn RemoveColumn()const
    {
	static_assert(index > 0,"dimension index is out of range");
	static_assert(index<=ColumnsCount,"dimension index is out of range");
	return MinusOneColumn(m_line.template RemoveComponent<index>());
    }
    template<size_t index>
    inline PlusOneColumn InsertColumn(const ColumnType&C)const
    {
	static_assert(index > 0,"dimension index is out of range");
	static_assert(index<=(ColumnsCount+1),"dimension index is out of range");
	return PlusOneColumn(m_line.template InsertComponent<index>(C.x()));
    }
    template<size_t index>
    inline PlusOneRow InsertRow(const LineType&L)const
    {
	static_assert(index > 0,"dimension index is out of range");
	static_assert(index<=(RowsCount+1),"dimension index is out of range");
	if constexpr(index==(RowsCount+1)) return lines(m_line,L);
	if constexpr(index==RowsCount) return lines(L,m_line);
    }
    template<size_t index, size_t jindex>
    inline const NumberType &element()const
    {
	static_assert(index == 1 ,"dimension index is out of range");
        return m_line.template component<jindex>();
    }
    inline NumberType Determinant()const
    {
	static_assert(size_t(ColumnsCount)==size_t(RowsCount),"cannot calculate a non-squared matrix determinant");
	return m_line.x();
    }
    inline LineType Cramer(const LineType&X)const
    {
	static_assert(size_t(ColumnsCount)==size_t(RowsCount),"cannot calculate a non-squared matrix determinant");
	const NumberType D=Determinant();
	if(D==0)throw Exception<Matrix>("System of equations cannot be solved");
	return desCartes(X.x()/D);
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
    ColumnType operator*(const LineType &v)const
    {
        return m_line * v;
    }
    Matrix<RowsCount, Vector<1, NumberType>>
            operator*(const Matrix<ColumnsCount, Vector<1, NumberType>> &B)const
    {
        return desCartes(m_line * B.___last_column());
    }
    template<size_t third_size>
    Matrix<RowsCount, Vector<third_size, NumberType>>
	operator*(const Matrix<ColumnsCount, Vector<third_size, NumberType>> &B)const
    {
        return operator*(B.___minus_one_column()).___add_column(operator*(B.___last_column()));
    }
    template<class otherlinetype>
    inline Matrix(const Matrix<RowsCount, otherlinetype> &source): m_line(source.___last_line()) {}
    template<class... Args>
    inline Matrix(const std::tuple<Args...> &v): m_line(std::get < RowsCount - 1 > (v)) {}
    inline Matrix(const ColumnType &A, const LineType &B): m_line(B *A.___last_component()) {}
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
        return Matrix(LineType::zero());
    }
    static inline Matrix one()
    {
        return Matrix(LineType::template basis_vector<RowsCount>());
    }
    template<size_t x, size_t y>
    static inline Matrix RotationInPlane(const NumberType &angle)
    {
        return Matrix(
	    (x == RowsCount) ? ((LineType::template basis_vector<x>() * cos(angle)) - (LineType::template basis_vector<y>() * sin(angle))) :
	    (y == RowsCount) ? ((LineType::template basis_vector<y>() * cos(angle)) + (LineType::template basis_vector<x>() * sin(angle))) :
	    LineType::template basis_vector<RowsCount>()
	);
    }
    template<class... Args>
    inline Matrix < RowsCount, Vector < ColumnsCount + 1 + sizeof...(Args), NumberType >> AddColumns(const ColumnType &col, Args...args)const
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
    enum {RowsCount = sizef};
    enum {ColumnsCount = linetype::Dimensions};
    typedef typename linetype::NumberType NumberType;
    typedef linetype LineType;
    typedef Vector<sizef, NumberType> ColumnType;
    typedef Matrix < sizef - 1, linetype > MinusOneRow;
    typedef Matrix<RowsCount+1, LineType> PlusOneRow;
    typedef Vector < ColumnsCount + 1, NumberType > LongerLine;
    typedef Matrix<RowsCount, LongerLine> PlusOneColumn;
    typedef Vector < ColumnsCount - 1, NumberType > ShorterLine;
    typedef Matrix<RowsCount, ShorterLine> MinusOneColumn;
    typedef typename MinusOneRow::MinusOneColumn Minor;
private:
    MinusOneRow m_other_lines;
    LineType m_line;
protected:
    inline Matrix(const MinusOneRow &minor, const LineType &line): m_other_lines(minor), m_line(line) {}
    inline MinusOneColumn ___minus_one_column()const
    {
        const auto new_minor = m_other_lines.___minus_one_column();
        const auto new_line = m_line.___recursive();
        return MinusOneColumn(new_minor, new_line);
    }
    inline ColumnType ___last_column()const
    {
        return ColumnType(m_other_lines.___last_column(), m_line.___last_component());
    }
    inline PlusOneColumn ___add_column(const ColumnType &col)const
    {
        const auto new_minor = m_other_lines.___add_column(col.___recursive());
        const auto new_line = LongerLine(m_line, col.___last_component());
        return PlusOneColumn(new_minor, new_line);
    }
    inline const LineType &___last_line()const
    {
        return m_line;
    }
    inline const MinusOneRow &___recursive()const
    {
        return m_other_lines;
    }
    inline const Matrix &AddColumns()const
    {
        return *this;
    }
public:
    virtual ~Matrix() {}
#ifdef ____optimized_version_of_matrices_h_____
    template<size_t index>
    inline MinusOneColumn RemoveColumn()const
    {
	static_assert(index > 0,"dimension index is out of range");
	static_assert(index<=ColumnsCount,"dimension index is out of range");
        const auto new_minor = m_other_lines.template RemoveColumn<index>();
        const auto new_line = m_line.template RemoveComponent<index>();
        return MinusOneColumn(new_minor, new_line);
    }
    template<size_t index>
    inline PlusOneColumn InsertColumn(const ColumnType&C)const
    {
	static_assert(index > 0,"dimension index is out of range");
	static_assert(index<=(ColumnsCount+1),"dimension index is out of range");
	return PlusOneColumn(
		m_other_lines.template InsertColumn<index>(C.___recursive()),
		m_line.template InsertComponent<index>(C.___last_component())
	);
    }
    template<size_t index>
    inline MinusOneRow RemoveRow()const
    {
	static_assert(index > 0,"dimension index is out of range");
	static_assert(index<=RowsCount,"dimension index is out of range");
        if constexpr(index == RowsCount)return m_other_lines;
	if constexpr(RowsCount==2) return MinusOneRow(m_line);
	else return MinusOneRow(m_other_lines.template RemoveRow<index>(),m_line);
    }
    template<size_t index>
    inline PlusOneRow InsertRow(const LineType&L)const
    {
	static_assert(index > 0,"dimension index is out of range");
	static_assert(index<=(RowsCount+1),"dimension index is out of range");
	if constexpr(index==(RowsCount+1)) return PlusOneRow(*this,L);
	if constexpr(index==RowsCount) return PlusOneRow(Matrix(m_other_lines,L),m_line);
	if constexpr(index<RowsCount) return PlusOneRow(m_other_lines.template InsertRow<index>(L),m_line);
    }
    template<size_t index, size_t jindex>
    inline const NumberType &element()const
    {
	static_assert(index > 0,"dimension index is out of range");
	static_assert(index<=RowsCount,"dimension index is out of range");
	if constexpr(index==RowsCount) return m_line.template component<jindex>();
	if constexpr(index<RowsCount) return m_other_lines.template element<index, jindex>();
    }
    template<size_t row,size_t col>
    inline Minor GetMinor()const{return RemoveColumn<col>().template RemoveRow<row>();}
private:
    template<size_t index>
    inline NumberType __det_until()const
    {
	static_assert(index <= RowsCount,"dimension index is out of range");
	static_assert(index > 0,"dimension index is out of range");
	const NumberType res=element<1,index>()*GetMinor<1,index>().Determinant();
	if constexpr(index==RowsCount) return res;
	if constexpr(index<RowsCount) return res-__det_until<index+1>();
    }
public:
    inline NumberType Determinant()const
    {
	static_assert(size_t(ColumnsCount)==size_t(RowsCount),"cannot calculate a non-squared matrix determinant");
	return __det_until<1>();
    }
private:
    template<size_t index>
    inline Vector<index,NumberType> ___cramer(const LineType&X,const NumberType&D)const
    {
	static_assert(index <= ColumnsCount,"dimension index is out of range");
	static_assert(index > 0,"dimension index is out of range");
	const NumberType d=RemoveColumn<index>().template InsertColumn<index>(X).Determinant();
	if constexpr(index==1) return desCartes(d/D);
	if constexpr(index>1) return Vector<index,NumberType>(___cramer<index-1>(X,D),d/D);
    }
public:
    inline LineType Cramer(const LineType&X)const
    {
	static_assert(size_t(ColumnsCount)==size_t(RowsCount),"cannot calculate a non-squared matrix determinant");
	const NumberType D=Determinant();
	if(D==0)throw Exception<Matrix>("System of equations cannot be solved");
	return ___cramer<ColumnsCount>(X,D);
    }
#else
    template<size_t index, size_t jindex>
    const NumberType &element()const
    {
	static_assert(index > 0,"dimension index is out of range");
	if(index>RowsCount)throw Exception<Matrix>("dimension index is out of range");
	if(index==RowsCount) return m_line.template component<jindex>();
	else return m_other_lines.template element<index, jindex>();
    }
#endif
    ColumnType operator*(const LineType &v)const
    {
        return ColumnType(m_other_lines * v, m_line * v);
    }
    Matrix<RowsCount, Vector<1, NumberType>>
	operator*(const Matrix<ColumnsCount, Vector<1, NumberType>> &B)const
    {
        const auto P = m_other_lines * B;
        const auto C = desCartes(m_line * B.___last_column());
        return Matrix<RowsCount, Vector<1, NumberType>>(P, C);
    }
    template<size_t third_size>
    Matrix<RowsCount, Vector<third_size, NumberType>>
	operator*(const Matrix<ColumnsCount, Vector<third_size, NumberType>> &B)const
    {
        const auto P = operator*(B.___minus_one_column());
        const ColumnType C = operator*(B.___last_column());
        return P.___add_column(C);
    }
    template<class otherlinetype>
    inline Matrix(const Matrix<RowsCount, otherlinetype> &source): m_other_lines(source.___recursive()), m_line(source.___last_line()) {}
    template<class... Args>
    inline Matrix(const std::tuple<Args...> &v): m_other_lines(v), m_line(std::get < RowsCount - 1 > (v)) {}
    inline Matrix(const ColumnType &A, const LineType &B): m_other_lines(A.___recursive(), B), m_line(B *A.___last_component()) {}
    bool operator==(const Matrix &B)const
    {
        return (m_line == B.m_line) && (m_other_lines == B.m_other_lines);
    }
    Matrix operator*(const NumberType &v)const
    {
        return Matrix(m_other_lines * v, m_line * v);
    }
    Matrix operator/(const NumberType &v)const
    {
        return Matrix(m_other_lines / v, m_line / v);
    }
    Matrix operator+(const Matrix &B)const
    {
        return Matrix(m_other_lines + B.m_other_lines, m_line + B.m_line);
    }
    Matrix operator-(const Matrix &B)const
    {
        return Matrix(m_other_lines - B.m_other_lines, m_line - B.m_line);
    }
    static inline Matrix zero()
    {
        return Matrix(MinusOneRow::zero(), LineType::zero());
    }
    static inline Matrix one()
    {
        return Matrix(MinusOneRow::one(), LineType::template basis_vector<RowsCount>());
    }
    template<size_t x, size_t y>
    static inline Matrix RotationInPlane(const NumberType &angle)
    {
        return Matrix(
                   MinusOneRow::template RotationInPlane<x, y>(angle),
                   (x == RowsCount) ? ((LineType::template basis_vector<x>() * cos(angle)) - (LineType::template basis_vector<y>() * sin(angle))) :
                   (y == RowsCount) ? ((LineType::template basis_vector<y>() * cos(angle)) + (LineType::template basis_vector<x>() * sin(angle))) :
                   LineType::template basis_vector<RowsCount>()
                                                );
    }
    template<class... Args>
    inline Matrix < RowsCount, Vector < ColumnsCount + 1 + sizeof...(Args), NumberType >> AddColumns(const ColumnType &col, Args...args)const
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
template<size_t size, class LineType>
inline Matrix<size, LineType> TensorProduct(const Vector<size, typename LineType::NumberType> &A, const LineType &B)
{
    return Matrix<size, LineType>(A, B);
}

};
#endif
