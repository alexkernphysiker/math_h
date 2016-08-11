// this file is distributed under 
// MIT license
#ifndef ______MATRICES_H_______
#	define ______MATRICES_H_______
#include <list>
#include <vector>
#include <functional>
#include "error.h"
namespace MathTemplates{
	template<class numt>
	class Matrix{
	public:
		typedef std::function<const numt(const size_t,const size_t)> function;
		virtual ~Matrix(){}
		virtual const size_t height()const=0;
		virtual const size_t width()const=0;
	protected:
		virtual const numt get_element(const size_t i,const size_t j)const=0;
	public:
		const numt operator()(const size_t i,const size_t j)const{
			if(i>=height())throw Exception<Matrix>("Range check error");
			if(j>=width())throw Exception<Matrix>("Range check error");
			return get_element(i,j);
		}
		inline const bool HasSizeAs(const Matrix&other)const{
			return (height()==other.height())&&(width()==other.width());
		}
		const bool operator==(const Matrix&other)const{
			if(!HasSizeAs(other))return false;
			for(size_t i=0;i<height();i++)
				for(size_t j=0;j<width();j++)
					if(operator()(i,j)!=other(i,j))
						return false;
					return true;
		}
		inline const bool operator!=(const Matrix&other)const{return !operator==(other);}
		inline const bool operator==(const numt&c)const{
			return (width()==1)&&(height()==1)&&(operator()(0,0)==c);
		}
		inline const bool operator!=(const numt&c)const{return !operator==(c);}
	};
	template<class numt>
	class MatrixData:public Matrix<numt>{
	public:
		typedef std::vector<std::vector<numt>> container;
		typedef std::function<const numt(const size_t,const size_t,const numt&)> transform_function;
	private:
		container f_data;
	protected:
		virtual const numt get_element(const size_t i,const size_t j)const override{return f_data[i][j];}
	public:
		virtual const size_t height()const override{return f_data.size();}
		virtual const size_t width()const override{return f_data[0].size();}
		virtual ~MatrixData(){}
		MatrixData(const numt&A){
			f_data.push_back({A});
		}
		MatrixData(const container&A){
				if(A.size()==0)
				throw Exception<MatrixData,0>("invalid matrix size");
			if((A.size()==1)&&(A[0].size()==0))
				throw Exception<MatrixData,0>("invalid matrix size");
			const size_t first_row_size=A[0].size();
			if(first_row_size==0)
				throw Exception<MatrixData,0>("invalid matrix size");
			for(const auto&row:A){
				if(row.size()!=first_row_size)
					throw Exception<MatrixData,0>("invalid matrix size");
				f_data.push_back(std::vector<numt>());
				const auto row_index=f_data.size()-1;
				for(const auto&item:row)
					f_data[row_index].push_back(item);
			}
		}
		MatrixData(const Matrix<numt>&source){
			for(size_t i=0;i<source.height();i++){
				f_data.push_back(std::vector<numt>());
				for(size_t j=0;j<source.width();j++)
					f_data[i].push_back(source(i,j));
			}
		}
		MatrixData&Transform(const transform_function F){
			for(size_t i=0;i<height();i++)
				for(size_t j=0;j<width();j++)
					f_data[i][j]=F(i,j,f_data[i][j]);
			return *this;
		}
	};
	template<class numt>
	class MatrixByFormula:public Matrix<numt>{
	private:
		size_t f_N,f_M;
		typename Matrix<numt>::function f_func;
	public:
		MatrixByFormula(const size_t N,const size_t M,const typename Matrix<numt>::function F):f_N(N),f_M(M),f_func(F){}
		virtual ~MatrixByFormula(){}
		virtual const size_t height()const override{return f_N;}
		virtual const size_t width()const override{return f_M;}
	protected:
		virtual const numt get_element(const size_t i,const size_t j)const{return f_func(i,j);}
	};
	
	template<class numt>
	const MatrixByFormula<numt> Unitary(const size_t N){
		return MatrixByFormula<numt>(N,N,[](size_t i,size_t j)->numt{return (i==j)?1:0;});
	}
	template<class numt>
	const MatrixByFormula<numt> Zeros(const size_t N,const size_t M){
		return MatrixByFormula<numt>(N,M,[](size_t,size_t)->numt{return 0;});
	}
	template<class numt>
	const MatrixByFormula<numt> RVec(const size_t N,const size_t i){
		if(i>=N)throw Exception<Matrix<numt>>("Invalid reference vector");
		return MatrixByFormula<numt>(N,1,[i](size_t ii,size_t)->numt{return (ii==i)?1:0;});
	}
	template<class numt>
	const MatrixByFormula<numt> Diagonal(const std::vector<numt>&V){
		return MatrixByFormula<numt>(V.size(),V.size(),[&V](size_t i,size_t j)->numt{return (i==j)?V[i]:0;});
	}
	template<class numt>
	const MatrixByFormula<numt> Permutation(const std::vector<size_t>&V){
		for(const size_t i:V)if(i>=V.size())
			throw Exception<MatrixByFormula<numt>>("invalid permutation matrix");
		return MatrixByFormula<numt>(V.size(),V.size(),[&V](size_t i,size_t j)->numt{return j==V[i]?1:0;});
	}
	
	
	template<class numt>
	const MatrixByFormula<numt> Transponate(const Matrix<numt>&source){
		return MatrixByFormula<numt>(source.width(),source.height(),[&source](size_t i,size_t j)->numt{return source(j,i);});
	}
	template<class numt>
	const MatrixByFormula<numt> Minor(const Matrix<numt>&source,const size_t i,const size_t j){
		if((source.height()<2)||(source.width()<2)||(i>=source.height())||(j>=source.width()))
			throw Exception<MatrixByFormula<numt>>("Invalid minor");
		return MatrixByFormula<numt>(source.height()-1,source.width()-1,[&source,i,j](size_t ii,size_t jj)->numt{
			auto new_ii=ii,new_jj=jj;
			if(ii>=i)new_ii++;if(jj>=j)new_jj++;
			return source(new_ii,new_jj);
		});
	}
	template<class numt>
	const numt Determinant(const Matrix<numt>&source){
		if((source.height()!=source.width())||(source.height()==0))
			throw Exception<MatrixByFormula<numt>>("Cannot calculate the determinant");
		if(source.height()==1)return source(0,0);
		numt result=0,k=1;
		for(size_t i=0;i<source.width();i++){
			result+=k*source(0,i)*Determinant(Minor(source,0,i));
			k=-k;
		}
		return result;
	}
	template<class numt>
	const MatrixByFormula<numt> Multiply(const Matrix<numt>&A,const Matrix<numt>&B){
		if(A.width()!=B.height())
			throw Exception<Matrix<numt>>("Matrix Multiplication: size mismatch");
		size_t cycle_length=A.width();
		return MatrixByFormula<numt>(A.height(),B.width(),[&A,&B,cycle_length](size_t i,size_t j)->numt{
			numt result=0;
			for(size_t k=0;k<cycle_length;k++)
				result+=A(i,k)*B(k,j);
			return result;
		});
	}

}
#endif