// this file is distributed under
// LGPLv3 license
#ifndef ______MATRICES_H_______
#	define ______MATRICES_H_______
#include <iostream>
#include <math.h>
#include "error.h"
#include "vectors.h"
namespace MathTemplates
{
template<class linetype, class... Args>
inline Matrix < sizeof...(Args) + 1, linetype > rows(const linetype &x, Args... args);
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
    typedef Vector < ColumnsCount + 1, NumberType > LongerLine;
    typedef Matrix<RowsCount, LongerLine> PlusOneColumn;
    typedef Vector < ColumnsCount - 1, NumberType > ShorterLine;
    typedef Matrix<RowsCount, ShorterLine> MinusOneColumn;
    typedef Matrix<RowsCount+1,linetype> PlusOneRow;
    typedef Vector<1, NumberType> DefaultColumnType;
    typedef linetype DefaultRowType;
    template<class numt2=NumberType> using ColumnType=Vector<1, numt2>;
    template<class numt2=NumberType> using RowType=Vector<linetype::Dimensions,numt2>;
private:
    DefaultRowType m_row;
protected:
    inline Matrix(const DefaultRowType &line): m_row(line) {}
    inline MinusOneColumn ___minus_one_column()const
    {
        return MinusOneColumn(m_row.___minus_one_component());
    }
    inline DefaultColumnType ___last_column()const
    {
        return vec(m_row.___last_component());
    }
    inline const DefaultRowType &___last_row()const
    {
        return m_row;
    }
public:
    virtual ~Matrix() {}
    inline PlusOneColumn AddColumn(const DefaultColumnType &col)const
    {
        return LongerLine(m_row, col.___last_component());
    }
    inline Matrix<RowsCount,Vector<ColumnsCount+1,NumberType>> AddColumns(const Matrix<RowsCount,Vector<1,NumberType>>&M)const
    {
	return AddColumn(M.___last_column());
    }
    template<size_t third_size>
    inline Matrix<RowsCount,Vector<ColumnsCount+third_size,NumberType>> AddColumns(const Matrix<RowsCount,Vector<third_size,NumberType>>&M)const
    {
	return AddColumns(M.___minus_one_column()).AddColumn(M.___last_column());
    }
    inline PlusOneRow AddRow(const DefaultRowType &row)const
    {
        return rows(m_row,row);
    }
    inline Matrix<RowsCount+1,Vector<ColumnsCount,NumberType>> AddRows(const Matrix<1,Vector<ColumnsCount,NumberType>>&M)const
    {
	return AddRow(M.___last_row());
    }
    template<size_t third_size>
    inline Matrix<RowsCount+third_size,Vector<ColumnsCount,NumberType>> AddRows(const Matrix<third_size,Vector<ColumnsCount,NumberType>>&M)const
    {
	return AddRows(M.___minus_one_component()).AddRow(M.___last_row());
    }
    template<size_t index, size_t jindex>
    inline const NumberType &element()const
    {
	static_assert(index == 1 ,"dimension index is out of range");
        return m_row.template component<jindex>();
    }
    template<size_t index>
    inline const DefaultRowType &row()const
    {
	static_assert(index == 1 ,"dimension index is out of range");
        return m_row;
    }

    template<size_t index>
    inline MinusOneColumn RemoveColumn()const
    {
	static_assert(index > 0,"dimension index is out of range");
	static_assert(index<=ColumnsCount,"dimension index is out of range");
	return MinusOneColumn(m_row.template RemoveComponent<index>());
    }
    template<size_t index,size_t count>
    inline auto RemoveColumns()const
    {
	static_assert(index > 0,"dimension index is out of range");
	static_assert(index<=ColumnsCount,"dimension index is out of range");
	static_assert(count>0,"cannot remove less that 1 column");
	static_assert(count<ColumnsCount,"too many columns to remove");
	if constexpr(count==1)return RemoveColumn<index>();
	if constexpr(count>1)return RemoveColumns<index,count-1>().template RemoveColumn<index>();
    }
    template<size_t index>
    inline PlusOneColumn InsertColumn(const DefaultColumnType&C)const
    {
	static_assert(index > 0,"dimension index is out of range");
	static_assert(index<=(ColumnsCount+1),"dimension index is out of range");
	return PlusOneColumn(m_row.template InsertComponent<index>(C.x()));
    }
    template<size_t index,size_t third_size>
    inline auto InsertColumns(const Matrix<RowsCount,Vector<third_size,NumberType>>&C)const
    {
	static_assert(index > 0,"dimension index is out of range");
	static_assert(index<=(ColumnsCount+1),"dimension index is out of range");
	if constexpr(third_size==1)return InsertColumn<index>(C.___last_column());
	if constexpr(third_size>1)return InsertColumn<index>(C.___last_column())
	    .template InsertColumns<index,third_size-1>(C.___minus_one_column());
    }

    template<size_t index>
    inline PlusOneRow InsertRow(const DefaultRowType&L)const
    {
	static_assert(index > 0,"dimension index is out of range");
	static_assert(index<=(RowsCount+1),"dimension index is out of range");
	if constexpr(index==(RowsCount+1)) return rows(m_row,L);
	if constexpr(index==RowsCount) return rows(L,m_row);
    }
    template<size_t index,size_t third_size>
    inline auto InsertRows(const Matrix<third_size,Vector<ColumnsCount,NumberType>>&C)const
    {
	static_assert(index > 0,"dimension index is out of range");
	static_assert(index<=(RowsCount+1),"dimension index is out of range");
	if constexpr(third_size==1)return InsertRow<index>(C.___last_row());
	if constexpr(third_size>1)return InsertRow<index>(C.___last_row())
	    .template InsertRows<index,third_size-1>(C.___minus_one_component());
    }
    inline NumberType Determinant()const
    {
	static_assert(size_t(ColumnsCount)==size_t(RowsCount),"cannot calculate a non-squared matrix determinant");
	return m_row.x();
    }
    inline DefaultRowType Cramer(const DefaultRowType&X)const
    {
	static_assert(size_t(ColumnsCount)==size_t(RowsCount),"cannot calculate a non-squared matrix determinant");
	const NumberType D=Determinant();
	if(D==0)throw Exception<Matrix>("System of equations cannot be solved");
	return vec(X.x()/D);
    }
    inline Matrix<ColumnsCount,Vector<1,NumberType>> transponate()const{
	return Matrix <ColumnsCount,Vector<1, NumberType>> (m_row.to_tuple());
    }
    inline Vector<RowsCount,NumberType> diagonal()const
    {
	return vec(m_row.x());
    }
    inline NumberType trace()const
    {
	return m_row.x();
    }
    template<class numt2>
    inline DefaultColumnType operator*(const RowType<numt2> &v)const
    {
        return m_row * v;
    }
    template<class numt2>
    inline Matrix<RowsCount, Vector<1, NumberType>>
            operator*(const Matrix<ColumnsCount, Vector<1, numt2>> &B)const
    {
        return vec(m_row * B.___last_column());
    }
    template<size_t third_size,class numt2>
    inline Matrix<RowsCount, Vector<third_size, NumberType>>
	operator*(const Matrix<ColumnsCount, Vector<third_size, numt2>> &B)const
    {
        return operator*(B.___minus_one_column()).___add_column(operator*(B.___last_column()));
    }
    template<class otherlinetype>
    inline Matrix(const Matrix<RowsCount, otherlinetype> &source): m_row(source.___last_row()) {}
    template<class... Args>
    inline Matrix(const std::tuple<Args...> &v): m_row(std::get < RowsCount - 1 > (v)) {}
    inline Matrix(const DefaultColumnType &A, const DefaultRowType &B): m_row(B *A.___last_component()) {}
    inline bool operator==(const Matrix &B)const
    {
        return m_row == B.m_row;
    }
    template<class numt2>
    inline auto operator*(const numt2 &v)const->Matrix<RowsCount,RowType<decltype(element<1,1>()*v)>>
    {
        return Matrix<RowsCount,RowType<decltype(element<1,1>()*v)>>(m_row * v);
    }
    inline Matrix operator/(const NumberType &v)const
    {
        return Matrix(m_row / v);
    }
    inline Matrix operator+(const Matrix&B)const
    {
        return Matrix(m_row + B.___last_row());
    }
    inline Matrix operator-(const Matrix&B)const
    {
        return Matrix(m_row - B.___last_row());
    }
    static inline Matrix zero()
    {
        return Matrix(DefaultRowType::zero());
    }
    static inline Matrix one()
    {
        return Matrix(DefaultRowType::template basis_vector<RowsCount>());
    }
    template<size_t x, size_t y>
    static inline Matrix RotationInPlane(const NumberType &angle)
    {
        return Matrix(
	    (x == RowsCount) ? ((DefaultRowType::template basis_vector<x>() * cos(angle)) - (DefaultRowType::template basis_vector<y>() * sin(angle))) :
	    (y == RowsCount) ? ((DefaultRowType::template basis_vector<y>() * cos(angle)) + (DefaultRowType::template basis_vector<x>() * sin(angle))) :
	    DefaultRowType::template basis_vector<RowsCount>()
	);
    }
    inline void output(std::ostream&str)const{str<<m_row;}
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
    typedef Matrix < sizef - 1, linetype > MinusOneRow;
    typedef Matrix<RowsCount+1, linetype> PlusOneRow;
    typedef Vector < ColumnsCount + 1, NumberType > LongerLine;
    typedef Matrix<RowsCount, LongerLine> PlusOneColumn;
    typedef Vector < ColumnsCount - 1, NumberType > ShorterLine;
    typedef Matrix<RowsCount, ShorterLine> MinusOneColumn;
    typedef typename MinusOneRow::MinusOneColumn Minor;
    typedef Vector<sizef, NumberType> DefaultColumnType;
    typedef linetype DefaultRowType;
    template<class numt2=NumberType> using ColumnType=Vector<1, numt2>;
    template<class numt2=NumberType> using RowType=Vector<linetype::Dimensions,numt2>;
private:
    MinusOneRow m_other_rows;
    DefaultRowType m_row;
protected:
    inline Matrix(const MinusOneRow &minor, const DefaultRowType &line): m_other_rows(minor), m_row(line) {}
    inline MinusOneColumn ___minus_one_column()const
    {
        return MinusOneColumn(m_other_rows.___minus_one_column(),m_row.___minus_one_component());
    }
    inline DefaultColumnType ___last_column()const
    {
        return DefaultColumnType(m_other_rows.___last_column(), m_row.___last_component());
    }
    inline const DefaultRowType &___last_row()const
    {
        return m_row;
    }
    inline const MinusOneRow &___minus_one_component()const
    {
        return m_other_rows;
    }
public:
    virtual ~Matrix() {}
    inline PlusOneColumn AddColumn(const DefaultColumnType &col)const
    {
        return PlusOneColumn(m_other_rows.AddColumn(col.___minus_one_component()),LongerLine(m_row, col.___last_component()));
    }
    inline Matrix<RowsCount,Vector<ColumnsCount+1,NumberType>> AddColumns(const Matrix<RowsCount,Vector<1,NumberType>>&M)const
    {
	return AddColumn(M.___last_column());
    }
    template<size_t third_size>
    inline Matrix<RowsCount,Vector<ColumnsCount+third_size,NumberType>> AddColumns(const Matrix<RowsCount,Vector<third_size,NumberType>>&M)const
    {
	return AddColumns(M.___minus_one_column()).AddColumn(M.___last_column());
    }
    inline PlusOneRow AddRow(const DefaultRowType &row)const
    {
        return PlusOneRow(*this,row);
    }
    inline Matrix<RowsCount+1,Vector<ColumnsCount,NumberType>> AddRows(const Matrix<1,Vector<ColumnsCount,NumberType>>&M)const
    {
	return AddRow(M.___last_row());
    }
    template<size_t third_size>
    inline Matrix<RowsCount+third_size,Vector<ColumnsCount,NumberType>> AddRows(const Matrix<third_size,Vector<ColumnsCount,NumberType>>&M)const
    {
	return AddRows(M.___minus_one_component()).AddRow(M.___last_row());
    }
    template<size_t index, size_t jindex>
    inline const NumberType &element()const
    {
	static_assert(index > 0,"dimension index is out of range");
	static_assert(index<=RowsCount,"dimension index is out of range");
	if constexpr(index==RowsCount) return m_row.template component<jindex>();
	if constexpr(index<RowsCount) return m_other_rows.template element<index, jindex>();
    }
    template<size_t index>
    inline const DefaultRowType &row()const
    {
	static_assert(index > 0,"dimension index is out of range");
	static_assert(index<=RowsCount,"dimension index is out of range");
	if constexpr(index==RowsCount) return m_row;
	if constexpr(index<RowsCount) return m_other_rows.template row<index>();
    }
    template<size_t index>
    inline MinusOneColumn RemoveColumn()const
    {
	static_assert(index > 0,"dimension index is out of range");
	static_assert(index<=ColumnsCount,"dimension index is out of range");
        return MinusOneColumn(m_other_rows.template RemoveColumn<index>(),m_row.template RemoveComponent<index>());
    }
    template<size_t index,size_t count>
    inline auto RemoveColumns()const
    {
	static_assert(index > 0,"dimension index is out of range");
	static_assert(index<=ColumnsCount,"dimension index is out of range");
	static_assert(count>0,"cannot remove less that 1 column");
	static_assert(count<ColumnsCount,"too many columns to remove");
	if constexpr(count==1)return RemoveColumn<index>();
	if constexpr(count>1)return RemoveColumns<index,count-1>().template RemoveColumn<index>();
    }

    template<size_t index>
    inline MinusOneRow RemoveRow()const
    {
	static_assert(index > 0,"dimension index is out of range");
	static_assert(index<=RowsCount,"dimension index is out of range");
        if constexpr(index == RowsCount)return m_other_rows;
	if constexpr(RowsCount==2) return MinusOneRow(m_row);
	else return MinusOneRow(m_other_rows.template RemoveRow<index>(),m_row);
    }
    template<size_t index,size_t count>
    inline auto RemoveRows()const
    {
	static_assert(index > 0,"dimension index is out of range");
	static_assert(index<=RowsCount,"dimension index is out of range");
	static_assert(count>0,"cannot remove less that 1 column");
	static_assert(count<RowsCount,"too many columns to remove");
	if constexpr(count==1)return RemoveRow<index>();
	if constexpr(count>1)return RemoveRows<index,count-1>().template RemoveRow<index>();
    }

    template<size_t index>
    inline PlusOneColumn InsertColumn(const DefaultColumnType&C)const
    {
	static_assert(index > 0,"dimension index is out of range");
	static_assert(index<=(ColumnsCount+1),"dimension index is out of range");
	return PlusOneColumn(
		m_other_rows.template InsertColumn<index>(C.___minus_one_component()),
		m_row.template InsertComponent<index>(C.___last_component())
	);
    }
    template<size_t index,size_t third_size>
    inline auto InsertColumns(const Matrix<RowsCount,Vector<third_size,NumberType>>&C)const
    {
	static_assert(index > 0,"dimension index is out of range");
	static_assert(index<=(ColumnsCount+1),"dimension index is out of range");
	if constexpr(third_size==1)return InsertColumn<index>(C.___last_column());
	if constexpr(third_size>1)return InsertColumn<index>(C.___last_column())
	    .template InsertColumns<index,third_size-1>(C.___minus_one_column());
    }
    template<size_t index>
    inline PlusOneRow InsertRow(const DefaultRowType&L)const
    {
	static_assert(index > 0,"dimension index is out of range");
	static_assert(index<=(RowsCount+1),"dimension index is out of range");
	if constexpr(index==(RowsCount+1)) return PlusOneRow(*this,L);
	if constexpr(index==RowsCount) return PlusOneRow(Matrix(m_other_rows,L),m_row);
	if constexpr(index<RowsCount) return PlusOneRow(m_other_rows.template InsertRow<index>(L),m_row);
    }
    template<size_t index,size_t third_size>
    inline auto InsertRows(const Matrix<third_size,Vector<ColumnsCount,NumberType>>&C)const
    {
	static_assert(index > 0,"dimension index is out of range");
	static_assert(index<=(RowsCount+1),"dimension index is out of range");
	if constexpr(third_size==1)return InsertRow<index>(C.___last_row());
	if constexpr(third_size>1)return InsertRow<index>(C.___last_row())
	    .template InsertRowss<index,third_size-1>(C.___minus_one_component());
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
    template<size_t index>
    inline Vector<index,NumberType> ___cramer(const DefaultRowType&X,const NumberType&D)const
    {
	static_assert(index <= ColumnsCount,"dimension index is out of range");
	static_assert(index > 0,"dimension index is out of range");
	const NumberType d=RemoveColumn<index>().template InsertColumn<index>(X).Determinant();
	if constexpr(index==1) return vec(d/D);
	if constexpr(index>1) return Vector<index,NumberType>(___cramer<index-1>(X,D),d/D);
    }
public:
    inline NumberType Determinant()const
    {
	static_assert(size_t(ColumnsCount)==size_t(RowsCount),"cannot calculate a non-squared matrix determinant");
	return __det_until<1>();
    }
    inline DefaultRowType Cramer(const DefaultRowType&X)const
    {
	static_assert(size_t(ColumnsCount)==size_t(RowsCount),"cannot calculate a non-squared matrix determinant");
	const NumberType D=Determinant();
	if(D==0)throw Exception<Matrix>("System of equations cannot be solved");
	return ___cramer<ColumnsCount>(X,D);
    }
    inline Matrix<ColumnsCount,Vector<RowsCount,NumberType>> transponate()const{
	return m_other_rows.transponate().template InsertColumn<RowsCount>(m_row);
    }
    inline Vector<RowsCount,NumberType> diagonal()const
    {
	return Vector<RowsCount,NumberType>(m_other_rows.diagonal(),m_row.template component<RowsCount>());
    }
    inline NumberType trace()const
    {
	return m_row.template component<RowsCount>()+m_other_rows.trace();
    }
    template<class numt2>
    inline DefaultColumnType operator*(const RowType<numt2> &v)const
    {
        return DefaultColumnType(typename MinusOneRow::DefaultColumnType(m_other_rows * v), NumberType(m_row * v));
    }
    template<class numt2>
    Matrix<RowsCount, Vector<1, NumberType>>
	operator*(const Matrix<ColumnsCount, Vector<1, numt2>> &B)const
    {
        const auto P = m_other_rows * B;
        const auto C = vec(m_row * B.___last_column());
        return Matrix<RowsCount, Vector<1, NumberType>>(P, C);
    }
    template<size_t third_size, class numt2>
    Matrix<RowsCount, Vector<third_size, NumberType>>
	operator*(const Matrix<ColumnsCount, Vector<third_size, numt2>> &B)const
    {
        const auto P = operator*(B.___minus_one_column());
        const DefaultColumnType C = operator*(B.___last_column());
        return P.AddColumn(C);
    }
    template<class otherlinetype>
    inline Matrix(const Matrix<RowsCount, otherlinetype> &source): m_other_rows(source.___minus_one_component()), m_row(source.___last_row()) {}
    template<class... Args>
    inline Matrix(const std::tuple<Args...> &v): m_other_rows(v), m_row(std::get < RowsCount - 1 > (v)) {}
    inline Matrix(const DefaultColumnType &A, const DefaultRowType &B): m_other_rows(A.___minus_one_component(), B), m_row(B *A.___last_component()) {}
    inline bool operator==(const Matrix &B)const
    {
        return (m_row == B.m_row) && (m_other_rows == B.m_other_rows);
    }
    template<class numt2>
    inline auto operator*(const numt2 &v)const->Matrix<RowsCount,RowType<decltype(element<1,1>()*v)>>
    {
        return Matrix<RowsCount,RowType<decltype(element<1,1>()*v)>>(m_other_rows * v, m_row * v);
    }
    inline Matrix operator/(const NumberType &v)const
    {
        return Matrix(m_other_rows / v, m_row / v);
    }
    template<class numt2>
    inline Matrix operator+(const Matrix<RowsCount,RowType<numt2>> &B)const
    {
        return Matrix(m_other_rows + B.___minus_one_component(), m_row + B.___last_row());
    }
    template<class numt2>
    inline Matrix operator-(const Matrix<RowsCount,RowType<numt2>> &B)const
    {
        return Matrix(m_other_rows - B.___minus_one_component(), m_row - B.___last_row());
    }
    static inline Matrix zero()
    {
        return Matrix(MinusOneRow::zero(), DefaultRowType::zero());
    }
    static inline Matrix one()
    {
        return Matrix(MinusOneRow::one(), DefaultRowType::template basis_vector<RowsCount>());
    }
    template<size_t x, size_t y>
    static inline Matrix RotationInPlane(const NumberType &angle)
    {
        return Matrix(
	    MinusOneRow::template RotationInPlane<x, y>(angle),
	    (x == RowsCount) ? ((DefaultRowType::template basis_vector<x>() * cos(angle)) - (DefaultRowType::template basis_vector<y>() * sin(angle))) :
		(y == RowsCount) ? ((DefaultRowType::template basis_vector<y>() * cos(angle)) + (DefaultRowType::template basis_vector<x>() * sin(angle))) :
		DefaultRowType::template basis_vector<RowsCount>()
	);
    }
    inline void output(std::ostream&str)const{
	m_other_rows.output(str);
	str<<std::endl<<m_row;
    }
};
template<class linetype, class... Args>
inline Matrix < sizeof...(Args) + 1, linetype > rows(const linetype &x, Args... args)
{
    return Matrix < sizeof...(Args) + 1, linetype > (std::make_tuple(x, args...));
}
template<class numt, class... Args>
inline Matrix < 1, Vector < sizeof...(Args) + 1, numt >> row(const numt &x, Args... args)
{
    return rows(vec(x, args...));
}
template<class numt, class...Args>
inline Matrix < 1 + sizeof...(Args), Vector<1, numt >> column(const numt &x, Args...args)
{
    return Matrix < 1 + sizeof...(Args), Vector<1, numt >> (std::make_tuple(x, args...));
}
template<class VecT>
inline Matrix < VecT::Dimensions, Vector < 1, typename VecT::NumberType >> columns(const VecT &x)
{
    return Matrix<VecT::Dimensions, Vector<1, typename VecT::NumberType>>(x.to_tuple());
}
template<class VecT, class...Args>
inline Matrix < VecT::Dimensions, Vector < 1 + sizeof...(Args), typename VecT::NumberType >> columns(const VecT &x, Args...args)
{
    return Matrix<VecT::Dimensions, Vector<1, typename VecT::NumberType>>(x.to_tuple()).AddColumns(columns(args...));
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
template<size_t size, class RowType>
inline Matrix<size, RowType> TensorProduct(const Vector<size, typename RowType::NumberType> &A, const RowType &B)
{
    return Matrix<size, RowType>(A, B);
}
template<size_t i,class VT>
inline std::ostream& operator<<(std::ostream&str,const Matrix<i,VT>&X)
{
    X.output(str);
    return str;
}
};
#endif
