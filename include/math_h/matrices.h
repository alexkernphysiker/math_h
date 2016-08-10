// this file is distributed under 
// MIT license
#ifndef ______MATRICES_H_______
#	define ______MATRICES_H_______
#include <vector>
#include <functional>
#include "error.h"
namespace MathTemplates{
	template<class numt>
	const bool MatrixValid(const std::vector<std::vector<numt>>&A){
		if(A.size()==0)return false;
		if(A.size()==1)return A[0].size()>0;
		const size_t first_row_size=A[0].size();
		if(first_row_size==0)return false;
		for(const auto&row:A)
			if(row.size()!=first_row_size)
				return false;
			return true;
	}
	template<class numt>
	const bool MatrixSquare(const std::vector<std::vector<numt>>&A){
		if(A.size()==0)return false;
		if(A.size()==1)return A[0].size()==1;
		const size_t first_row_size=A[0].size();
		if(first_row_size!=A.size())return false;
		for(const auto&row:A)
			if(row.size()!=first_row_size)
				return false;
			return true;
	}
	template<class numt>
	inline const bool MatrixSizeEqual(
		const std::vector<std::vector<numt>>&A,
		const std::vector<std::vector<numt>>&B
	){
		if(!MatrixValid<numt>(A))throw Exception<std::vector<std::vector<numt>>,0>("Multiplication: first operand is not valid");
		if(!MatrixValid<numt>(B))throw Exception<std::vector<std::vector<numt>>,0>("Multiplication: second operand is not valid");
		if(A.size()!=B.size())return false;
		return (A[0].size()==B[0].size());
	}
	template<class numt>
	inline const bool MatricesEqual(
		const std::vector<std::vector<numt>>&A,
		const std::vector<std::vector<numt>>&B
	){
		if(!MatrixSizeEqual<numt>(A,B))
			return false;
		for(size_t r=0;r<A.size();r++)
			for(size_t c=0;c<A[0].size();c++)
				if(A[r][c]!=B[r][c])
					return false;
		return true;
	}
	template<class numt>
	const std::vector<std::vector<numt>>MulMatrices(
		const std::vector<std::vector<numt>>&A,
		const std::vector<std::vector<numt>>&B
	){
		if(!MatrixValid<numt>(A))throw Exception<std::vector<std::vector<numt>>,0>("Multiplication: first operand is not valid");
		if(!MatrixValid<numt>(B))throw Exception<std::vector<std::vector<numt>>,0>("Multiplication: second operand is not valid");
		if(A[0].size()!=B.size())throw Exception<std::vector<std::vector<numt>>,1>("Multiplication: matrices size mismatch");
		std::vector<std::vector<numt>> result;
		const auto mul_cycle_length=B.size();
		const auto res_height=A.size();
		const auto res_width=B[0].size();
		for(size_t r=0;r<res_height;r++){
			result.push_back(std::vector<numt>());
			for(size_t c=0;c<res_width;c++){
				result[r].push_back(0);
				for(size_t m=0;m<mul_cycle_length;m++)
					result[r][c]+=A[r][m]*B[m][c];
			}
		}
		return result;
	}
}
#endif