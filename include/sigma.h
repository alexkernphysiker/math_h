// this file is distributed under 
// MIT license
#ifndef YJIOGPKSIIYJVKND
#define YJIOGPKSIIYJVKND
#include <list>
#include <vector>
#include <utility>
#include <math.h>
#include "exception_math_h.h"
template<typename numt>
class Sigma{
private:
	std::list<numt> m_list;
	numt m_sum;
public:
	Sigma(){m_sum=0;}
	virtual ~Sigma(){}
	Sigma &AddValue(numt value){
		m_list.push_back(value);
		m_sum+=value;
		return *this;
	}
	numt getAverage()const{
		int sz=m_list.size();
		if(sz<=0)
			throw math_h_error<Sigma>("No data to check. For average needed at least one element.");
		return m_sum/sz;
	}
	int count()const{return m_list.size();}
	numt getSigmaSqr()const{
		int sz=m_list.size();
		if(sz<=1)
			throw math_h_error<Sigma>("No data to check. for sigma needed at least two elements.");
		numt average=getAverage();
		numt m_sigsqr=0;
		for(auto value:m_list)
			m_sigsqr+=pow(value-average,2);
		m_sigsqr/=sz-1;
		return m_sigsqr;
	}
	numt getSigma()const{return sqrt(getSigmaSqr());}
};
template<typename numt>
class WeightedAverageCalculator{
private:
	numt Sum;
	numt Wnorm;
public:
	WeightedAverageCalculator(){Sum=0;Wnorm=0;}
	virtual ~WeightedAverageCalculator(){}
	WeightedAverageCalculator &AddValue(numt X,numt deltaX){
		if(deltaX<=0)
			throw math_h_error<WeightedAverageCalculator>("Attempt to add point with zero error.");
		numt w=1.0/pow(deltaX,2);
		Sum+=w*X;
		Wnorm+=w;
		return *this;
	}
	numt Average()const{
		if(Wnorm<=0)throw math_h_error<WeightedAverageCalculator>("Attempt to check empty data");
		return Sum/Wnorm;
	}
	numt Sigma()const{
		if(Wnorm<=0)throw math_h_error<WeightedAverageCalculator>("Attempt to check empty data");
		return 1.0/sqrt(Wnorm);
	}
};
#endif
