// this file is distributed under 
// MIT license
#ifndef YJIOGPKSIIYJVKND
#define YJIOGPKSIIYJVKND
#include <list>
#include <vector>
#include <utility>
#include <math.h>
#include "error.h"
namespace MathTemplates{
	using namespace std;
	template<typename numt>
	class value{
	private:
		numt Value,Error;
	public:
		value():Value(0),Error(0){}
		value(double v):Value(v),Error(sqrt(v)){if(Error<1)Error=1;}
		value(double v,double err):Value(v),Error(err){if(Error<0)throw Exception<value>("Error cannot be negative");}
		value(const value&source):Value(source.Value),Error(source.Error){}
		
		double val()const{return Value;}
		double delta()const{return Error;}
		double epsilon()const{return Error/Value;}
		double min()const{return Value-Error;}
		double max()const{return Value+Error;}
		bool contains(double x)const{return (x>=min())&&(x<=max());}
		bool contains(const value&x)const{return (x.max()>=min())&&(x.min()<=max());}
		
		value&operator+=(const value&other){
			Error=sqrt(pow(Error,2)+pow(other.Error,2));
			Value+=other.Value;
		}
		value&operator-=(const value&other){
			Error=sqrt(pow(Error,2)+pow(other.Error,2));
			Value-=other.Value;
		}
		value&operator*=(const value&other){
			Error=sqrt(pow(Error*other.Value,2)+pow(other.Error*Value,2));
			Value*=other.Value;
		}
		value&operator/=(const value&other){
			Error=sqrt(pow(Error/other.Value,2)+pow(other.Error*Value/pow(other.Value,2),2));
			Value/=other.Value;
		}
	};
	
	template<typename numt>
	value<numt> operator+(const value<numt>&a,const value<numt>&b){auto res=a;res+=b;return res;}
	template<typename numt>
	value<numt> operator+(value<numt>&&a,const value<numt>&b){return a+b;}
	template<typename numt>
	value<numt> operator+(const value<numt>&a,value<numt>&&b){return a+b;}
	template<typename numt>
	value<numt> operator+(value<numt>&&a,value<numt>&&b){return a+b;}
	
	template<typename numt>
	value<numt> operator-(const value<numt>&a,const value<numt>&b){auto res=a;res-=b;return res;}
	template<typename numt>
	value<numt> operator-(value<numt>&&a,const value<numt>&b){return a-b;}
	template<typename numt>
	value<numt> operator-(const value<numt>&a,value<numt>&&b){return a-b;}
	template<typename numt>
	value<numt> operator-(value<numt>&&a,value<numt>&&b){return a-b;}
	
	template<typename numt>
	value<numt> operator*(const value<numt>&a,const value<numt>&b){auto res=a;res*=b;return res;}
	template<typename numt>
	value<numt> operator*(value<numt>&&a,const value<numt>&b){return a*b;}
	template<typename numt>
	value<numt> operator*(const value<numt>&a,value<numt>&&b){return a*b;}
	template<typename numt>
	value<numt> operator*(value<numt>&&a,value<numt>&&b){return a*b;}

	template<typename numt>
	value<numt> operator/(const value<numt>&a,const value<numt>&b){auto res=a;res/=b;return res;}
	template<typename numt>
	value<numt> operator/(value<numt>&&a,const value<numt>&b){return a/b;}
	template<typename numt>
	value<numt> operator/(const value<numt>&a,value<numt>&&b){return a/b;}
	template<typename numt>
	value<numt> operator/(value<numt>&&a,value<numt>&&b){return a/b;}

	template<typename numt>
	bool operator<(const value<numt>&a,const value<numt>&b){return a.val()<b.val();}
	template<typename numt>
	bool operator<(value<numt>&&a,const value<numt>&b){return a.val()<b.val();}
	template<typename numt>
	bool operator<(const value<numt>&a,value<numt>&&b){return a.val()<b.val();}
	template<typename numt>
	bool operator<(value<numt>&&a,value<numt>&&b){return a.val()<b.val();}

	template<typename numt>
	bool operator>(const value<numt>&a,const value<numt>&b){return a.val()>b.val();}
	template<typename numt>
	bool operator>(value<numt>&&a,const value<numt>&b){return a.val()>b.val();}
	template<typename numt>
	bool operator>(const value<numt>&a,value<numt>&&b){return a.val()>b.val();}
	template<typename numt>
	bool operator>(value<numt>&&a,value<numt>&&b){return a.val()>b.val();}
	
	template<typename numt>
	class Sigma{
	private:
		list<numt> m_list;
		numt m_sum;
		value<numt>*m_cache;
	public:
		Sigma(){
			m_sum=0;
			m_cache=new value<numt>(INFINITY,INFINITY);
		}
		virtual ~Sigma(){
			delete m_cache;
		}
		Sigma &operator<<(numt x){
			m_list.push_back(x);
			m_sum+=x;
			(*m_cache)=value<numt>(INFINITY,INFINITY);
			return *this;
		}
		int count()const{return m_list.size();}
		value<numt>&get()const{
			if(isfinite(m_cache->val())){
				return const_cast<value<numt>&>(*m_cache);
			}else{
				int sz=m_list.size();
				if(sz<=1)
					throw Exception<Sigma>("No data to check. for sigma needed at least two elements.");
				numt average=m_sum/sz;
				numt m_sigsqr=0;
				for(auto value:m_list)
					m_sigsqr+=pow(value-average,2);
				m_sigsqr/=sz-1;
				(*m_cache)=value<numt>(average,sqrt(m_sigsqr));
				return const_cast<value<numt>&>(*m_cache);
			}
		}
	};
	template<typename numt>
	class WeightedAverageCalculator{
	private:
		numt Sum;
		numt Wnorm;
		value<numt>*m_cache;
	public:
		WeightedAverageCalculator(){
			Sum=0;Wnorm=0;
			m_cache=new value<numt>(INFINITY,INFINITY);
		}
		virtual ~WeightedAverageCalculator(){delete m_cache;}
		WeightedAverageCalculator &operator<<(const value<numt>&X){
			if(X.delta()<=0)
				throw Exception<WeightedAverageCalculator>("Cannot add value with zero error");
			numt w=1.0/pow(X.delta(),2);
			Sum+=w*X.val();
			Wnorm+=w;
			(*m_cache)=value<numt>(INFINITY,INFINITY);
			return *this;
		}
		WeightedAverageCalculator &operator<<(value<numt>&&X){
			return operator<<(X);
		}
		value<numt>&get()const{
			if(isfinite(m_cache->val())){
				return const_cast<value<numt>&>(*m_cache);
			}else{
				if(Wnorm<=0)throw Exception<WeightedAverageCalculator>("Attempt to check empty data");
				(*m_cache)=value<numt>(Sum/Wnorm,1.0/sqrt(Wnorm));
				return const_cast<value<numt>&>(*m_cache);
			}
		}
	};
};
#endif
