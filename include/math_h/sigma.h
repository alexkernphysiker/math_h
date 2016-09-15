// this file is distributed under 
// MIT license
#ifndef YJIOGPKSIIYJVKND
#define YJIOGPKSIIYJVKND
#include <iostream>
#include <list>
#include <vector>
#include <functional>
#include <memory>
#include <math.h>
#include "error.h"
namespace MathTemplates{
	template<typename numt>
	class value{
	private:
		numt Value,Error;
		struct cache{
			numt epsilon,min,max;
			cache(const numt&&e,const numt&&b,const numt&&a)
			:epsilon(e),min(b),max(a){}
		};
		std::shared_ptr<cache> f_cache;
		inline void invalidate(){
			if(Error<0)throw Exception<value>("Error cannot be negative");
			f_cache=nullptr;
		}
		inline void calc()const{
			if(!f_cache){
				const_cast<value&>(*this).f_cache=std::make_shared<cache>(
					Error/Value,Value-Error,Value+Error
				);
			}
		}
	public:
		virtual ~value(){}
		value():Value(numt(0)),Error(numt(0)){invalidate();}
		value(const numt v):Value(v),Error(numt(0)){invalidate();}
		value(const numt v,const numt err):Value(v),Error(err){invalidate();}
		value(const std::initializer_list<numt>&source){
			if(source.size()==0)
				throw Exception<value>("wrong initialization of value from emply list");
			if(source.size()>2)
				throw Exception<value>("wrong initialization of value from list with more than two numbers");
			std::vector<numt> v;for(const numt&x:source)v.push_back(x);
			Value=v[0];
			if(v.size()==1)Error=numt(0);
			else Error=v[1];
			invalidate();
		}
		static const value std_error(const numt&v){
			if(v<0)throw Exception<value>("Cannot calculate std error for negative value");
			auto res=value(v,sqrt(v));
			if(res.Error<numt(1))res.Error=numt(1);
			return res;
		}
		inline static const value std_error(const numt&&v){return std_error(v);}
		value(const value&source):Value(source.Value),Error(source.Error){invalidate();}
		value&operator=(const value&source){
			Value=source.Value;
			Error=source.Error;
			invalidate();
			return *this;
		}
		
		const numt&val()const{return Value;}
		const numt&uncertainty()const{return Error;}
		const numt&epsilon()const{calc();return f_cache->epsilon;}
		const numt&min()const{calc();return f_cache->min;}
		const numt&max()const{calc();return f_cache->max;}
		//Physical comparing of magnitudes with uncertainties
		const bool Contains(const numt& x)const{return (x>=min())&&(x<=max());}
		const bool Contains(const value&x)const{return (x.max()>=min())&&(x.min()<=max());}
		inline const bool Contains(const numt&&x)const{return Contains(x);}
		inline const bool Contains(const value&&x)const{return Contains(x);}
		const bool NotEqual(const numt& x)const{return (x<min())||(x>max());}
		const bool NotEqual(const value&x)const{return (x.max()<min())||(x.min()>max());}
		inline const bool NotEqual(const numt&&x)const{return NotEqual(x);}
		inline const bool NotEqual(const value&&x)const{return NotEqual(x);}
		const bool Below(const numt& x)const{return max()<x;}
		const bool Below(const value&x)const{return max()<x.min();}
		inline const bool Below(const numt&&x)const{return Below(x);}
		inline const bool Below(const value&&x)const{return Below(x);}
		const bool Above(const numt& x)const{return min()>x;}
		const bool Above(const value&x)const{return min()>x.max();}
		inline const bool Above(const numt&&x)const{return Above(x);}
		inline const bool Above(const value&&x)const{return Above(x);}
		
		//Inheriting number-like comparing
		const bool operator<(const value&other)const{return Value<other.Value;}
		const bool operator<(const value&&other)const{return Value<other.Value;}
		const bool operator>(const value&other)const{return Value>other.Value;}
		const bool operator>(const value&&other)const{return Value>other.Value;}
		const bool operator==(const value&other)const{return Value==other.Value;}
		const bool operator==(const value&&other)const{return Value==other.Value;}
		const bool operator>=(const value&other)const{return Value>=other.Value;}
		const bool operator>=(const value&&other)const{return Value>=other.Value;}
		const bool operator<=(const value&other)const{return Value<=other.Value;}
		const bool operator<=(const value&&other)const{return Value<=other.Value;}
		//arithmetic actions 
		value&operator+=(const value&other){
			Error=sqrt(pow(Error,2)+pow(other.Error,2));
			Value+=other.Value;
			invalidate();
			return *this;
		}
		inline value&operator+=(const value&&other){return operator+=(other);}
		value&operator-=(const value&other){
			Error=sqrt(pow(Error,2)+pow(other.Error,2));
			Value-=other.Value;
			invalidate();
			return *this;
		}
		inline value&operator-=(const value&&other){return operator-=(other);}
		value&operator*=(const value&other){
			Error=sqrt(pow(Error*other.Value,2)+pow(other.Error*Value,2));
			Value*=other.Value;
			invalidate();
			return *this;
		}
		inline value&operator*=(const value&&other){return operator*=(other);}
		value&operator/=(const value&other){
			Error=sqrt(pow(Error/other.Value,2)+pow(other.Error*Value/pow(other.Value,2),2));
			Value/=other.Value;
			invalidate();
			return *this;
		}
		inline value&operator/=(const value&&other){return operator/=(other);}
		
		const value operator+(const value&other)const{return value(*this)+=other;}
		const value operator+(const value&&other)const{return value(*this)+=other;}
		const value operator-(const value&other)const{return value(*this)-=other;}
		const value operator-(const value&&other)const{return value(*this)-=other;}
		const value operator*(const value&other)const{return value(*this)*=other;}
		const value operator*(const value&&other)const{return value(*this)*=other;}
		const value operator/(const value&other)const{return value(*this)/=other;}
		const value operator/(const value&&other)const{return value(*this)/=other;}
		const value Func(const std::function<numt(const numt&)>F)const{
			numt V=F(val());
			return value(V,sqrt( (pow(F(min())-V,2)+pow(F(max())-V,2)) / numt(2) ));
		}
	};
	template<typename numt>
	inline std::istream&operator>>(std::istream&str,value<numt>&P){
		numt v,u;
		str>>v>>u;
		P={v,u};
		return str;
	}
	template<typename numt>
	inline std::ostream&operator<<(std::ostream&str,const value<numt>&P){
		return str<<P.val()<<" "<<P.uncertainty();
	}
	template<typename numt>
	inline std::ostream&operator<<(std::ostream&str,const value<numt>&&P){return str<<P;}
	template<typename numt>
	const value<numt> func_value(const std::function<numt(const numt&)>F,const value<numt>&X){return X.Func(F);}
	template<typename numt>
	inline const value<numt> func_value(const std::function<numt(const numt&)> F,const value<numt>&&X){return X.Func(F);}
	template<typename numt>
	const std::function<value<numt>(const value<numt>&)> value_func(const std::function<numt(const numt&)>F){
		return [F](const value<numt>&X)->value<numt>{return X.Func(F);};
	}
	
	template<typename numt>
	class StandardDeviation{
	private:
		std::list<numt> m_list;
		numt m_sum;
		std::shared_ptr<value<numt>>m_cache;
		numt m_scale;
	public:
		StandardDeviation(const numt scale=1){
			if(scale<=0)throw Exception<StandardDeviation>("Uncertainty scaling factor must be greater than zero");
			m_scale=scale;
			m_sum=0;
		}
		virtual ~StandardDeviation(){}
		StandardDeviation &operator<<(const numt&x){
			m_list.push_back(x);
			m_sum+=x;
			m_cache=nullptr;
			return *this;
		}
		inline StandardDeviation&operator<<(const numt&&x){return operator<<(x);}
		const size_t count()const{return m_list.size();}
		const numt&scaling_factor()const{return m_scale;}
		const value<numt>&operator()()const{
			using namespace std;
			if(!m_cache){
				size_t sz=m_list.size();
				if(sz<=1)
					throw Exception<StandardDeviation>("No data to check. for sigma needed at least two elements.");
				numt average=m_sum/sz;
				numt m_sigsqr=0;
				for(auto value:m_list)
					m_sigsqr+=pow(value-average,2);
				m_sigsqr/=sz-1;
				const_cast<StandardDeviation&>(*this).m_cache=make_shared<value<numt>>(average,sqrt(m_sigsqr)*m_scale);
			}
			return *m_cache;
		}
	};
	template<typename numt>
	class WeightedAverage{
	private:
		numt Sum;
		numt Wnorm;
		std::shared_ptr<value<numt>>m_cache;
	public:
		WeightedAverage(){
			Sum=0;Wnorm=0;
		}
		virtual ~WeightedAverage(){}
		WeightedAverage&operator<<(const value<numt>&X){
			if(X.uncertainty()==0)
				throw Exception<WeightedAverage>("Cannot add value with zero error");
			numt w=1.0/pow(X.uncertainty(),2);
			Sum+=w*X.val();
			Wnorm+=w;
			m_cache=nullptr;
			return *this;
		}
		inline WeightedAverage&operator<<(const value<numt>&&X){return operator<<(X);}
		const value<numt>&operator()()const{
			using namespace std;
			if(!m_cache){
				if(Wnorm<=0)throw Exception<WeightedAverage>("Attempt to check empty data");
				const_cast<WeightedAverage&>(*this).m_cache=make_shared<value<numt>>(Sum/Wnorm,1.0/sqrt(Wnorm));
			}
			return *m_cache;
		}
	};
};
#endif
