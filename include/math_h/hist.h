// this file is distributed under 
// MIT license
#ifndef ___HISTS_HEADER______
#	define ___HISTS_HEADER______
#include <list>
#include <vector>
#include <utility>
#include <functional>
#include "error.h"
#include "sigma.h"
namespace MathTemplates{
	using namespace std;
	template<class numtX,class numtY=numtX>class point{
	private:
		value<numtX> x;
		value<numtY> y;
	public:
		point(const value<numtX>&pos):x(pos){}
		point(value<numtX>&&pos):x(pos){}
		point(const value<numtX>&pos,const value<numtY>&val):x(pos),y(val){}
		point(value<numtX>&&pos,const value<numtY>&val):x(pos),y(val){}
		point(const value<numtX>&pos,value<numtY>&&val):x(pos),y(val){}
		point(value<numtX>&&pos,value<numtY>&&val):x(pos),y(val){}
		point(const point&source):x(source.x),y(source.y){}
		value<numtX>&X()const{return const_cast<value<numtX>&>(x);}
		value<numtY>&Y()const{return const_cast<value<numtY>&>(y);}
	protected:
		value<numtX>&__X(){return x;}
		value<numtY>&__Y(){return y;}
	};
	template<class numtX,class numtY=numtX> class point_editable_x:public point<numtX,numtY>{
	public:
		point_editable_x(const point<numtX,numtY>&source):point<numtX,numtY>(source){}
		point_editable_x(point<numtX,numtY>&&source):point<numtX,numtY>(source){}
		value<numtX>&varX(){return point<numtX,numtY>::__X();}
	};
	template<class numtX,class numtY=numtX> class point_editable_y:public point<numtX,numtY>{
	public:
		point_editable_y(const point<numtX,numtY>&source):point<numtX,numtY>(source){}
		point_editable_y(point<numtX,numtY>&&source):point<numtX,numtY>(source){}
		value<numtY>&varY(){return point<numtX,numtY>::__Y();}
	};
	template<class numtX,class numtY=numtX> class point_editable:public point<numtX,numtY>{
	public:
		point_editable(const point<numtX,numtY>&source):point<numtX,numtY>(source){}
		point_editable(point<numtX,numtY>&&source):point<numtX,numtY>(source){}
		value<numtX>&varX(){return point<numtX,numtY>::__X();}
		value<numtY>&varY(){return point<numtX,numtY>::__Y();}
	};
	template<class numtX,class numtY=numtX>class hist{
	public:
		typedef point_editable_y<numtX,numtY> Point;
	private:
		vector<Point> m_data;
	public:
		hist(const initializer_list<value<numtX>>&data){
			for(const auto& v:data)m_data.push_back(Point(v));
		}
		hist(initializer_list<value<numtX>>&&data):hist(data){}
		hist(const vector<value<numtX>>&data){
			for(const auto& v:data)m_data.push_back(Point(v));
		}
		hist(vector<value<numtX>>&&data):hist(data){}
		hist(const initializer_list<point<numtX,numtY>>&data){
			for(const auto& P:data)m_data.push_back(Point(P));
		}
		hist(initializer_list<point<numtX,numtY>>&&data):hist(data){}
		hist(const vector<point<numtX,numtY>>&data){
			for(const auto& P:data)m_data.push_back(Point(P));
		}
		hist(vector<point<numtX,numtY>>&&data):hist(data){}
		hist(const initializer_list<Point>&data){
			for(const auto& P:data)m_data.push_back(P);
		}
		hist(initializer_list<Point>&&data):hist(data){}
		hist(const vector<Point>&data){
			for(const auto& P:data)m_data.push_back(P);
		}
		hist(vector<Point>&&data):hist(data){}
		hist(const hist&source){
			for(const Point& P:source.m_data)m_data.push_back(P);
		}
		virtual ~hist(){}
		hist&operator=(const hist& source){
			m_data.clear();
			for(const Point& P:source.m_data)m_data.push_back(P);
			return *this;
		}
		typedef typename vector<Point>::const_iterator const_iterator;
		const_iterator begin()const{return m_data.cbegin();}
		const_iterator cbegin()const{return m_data.cbegin();}
		const_iterator end() const{return m_data.cend();}
		const_iterator cend() const{return m_data.cend();}
		size_t size()const{return m_data.size();}
		Point&operator[](size_t i)const{
			if(m_data.size()<=i)
				throw Exception<hist>("range check error");
			return const_cast<Point&>(m_data[i]);
		}
		hist CloneEmptyBins()const{
			vector<Point> initer;
			for(const Point&P:m_data)initer.push_back(P);
			return hist(initer);
		}
		numtY Total()const{
			numtY res=0;
			for(const Point&P:m_data)res+=P.Y().val();
			return res;
		}
		value<numtY> TotalSum()const{
			value<numtY> res(0,0);
			for(const Point&P:m_data)res+=P.Y();
			return res;
		}
		numtY ChiSq_only_y_error(function<numtY(numtX)>f,size_t paramcount){
			numtY res=0,k=numtY(size())-numtY(paramcount);
			if(k<=0)throw Exception<hist>("ChiSq error: too few points or too many parameters");
			for(const Point&p:m_data)
				res+=pow((p.Y().val()-f(p.X().val()))/p.Y().delta(),2);
			return res/k;
		}
		numtY ChiSq(function<numtY(numtX)>f,size_t paramcount){
			numtY res=0,k=numtY(size())-numtY(paramcount);
			if(k<=0)throw Exception<hist>("ChiSq error: too few points or too many parameters");
			for(const Point&p:m_data){
				numtY w=pow(p.Y().delta(),2)+pow((f(p.X().min())-f(p.X().max()))/numtY(2),2);
				res+=pow(p.Y().val()-f(p.X().val()),2)/w;
			}
			return res/k;
		}
		
		hist&operator+=(const hist& second){
			for(size_t i=0,n=size();i<n;i++){
				if(m_data[i].X().val()==second[i].X().val()){
					m_data[i].varY()+=second[i].Y();
				}else
					throw Exception<hist>("Cannot add histogram. bins differ");
			}
			return *this;
		}
		hist&operator+=(function<numtY(numtX)>f){
			for(size_t i=0,n=size();i<n;i++)
				m_data[i].varY()+=value<numtY>(f(m_data[i].X().val()),0);
			return *this;
		}
		hist&operator+=(const value<numtY>&c){
			for(size_t i=0,n=size();i<n;i++)m_data[i].varY()+=c;
			return *this;
		}
		hist&operator+=(value<numtY>&&c){return operator+=(c);}
		
		hist&operator-=(const hist& second){
			for(size_t i=0,n=size();i<n;i++){
				if(m_data[i].X().val()==second[i].X().val()){
					m_data[i].varY()-=second[i].Y();
				}else
					throw Exception<hist>("Cannot add histogram. bins differ");
			}
			return *this;
		}
		hist&operator-=(function<numtY(numtX)>f){
			for(size_t i=0,n=size();i<n;i++)
				m_data[i].varY()-=value<numtY>(f(m_data[i].X().val()),0);
			return *this;
		}
		hist&operator-=(const value<numtY>&c){
			for(size_t i=0,n=size();i<n;i++)m_data[i].varY()-=c;
			return *this;
		}
		hist&operator-=(value<numtY>&&c){return operator-=(c);}
		
		hist&operator*=(const hist& second){
			for(size_t i=0,n=size();i<n;i++){
				if(m_data[i].X().val()==second[i].X().val()){
					m_data[i].varY()*=second[i].Y();
				}else
					throw Exception<hist>("Cannot add histogram. bins differ");
			}
			return *this;
		}
		hist&operator*=(function<numtY(numtX)>f){
			for(size_t i=0,n=size();i<n;i++)
				m_data[i].varY()*=value<numtY>(f(m_data[i].X().val()),0);
			return *this;
		}
		hist&operator*=(const value<numtY>&c){
			for(size_t i=0,n=size();i<n;i++)m_data[i].varY()*=c;
			return *this;
		}
		hist&operator*=(value<numtY>&&c){return operator*=(c);}
		
		hist&operator/=(const hist& second){
			for(size_t i=0,n=size();i<n;i++){
				if(m_data[i].X().val()==second[i].X().val()){
					m_data[i].varY()*=second[i].Y();
				}else
					throw Exception<hist>("Cannot add histogram. bins differ");
			}
			return *this;
		}
		hist&operator/=(function<numtY(numtX)>f){
			for(size_t i=0,n=size();i<n;i++)
				m_data[i].varY()/=value<numtY>(f(m_data[i].X().val()),0);
			return *this;
		}
		hist&operator/=(const value<numtY>&c){
			for(size_t i=0,n=size();i<n;i++)m_data[i].varY()/=c;
			return *this;
		}
		hist&operator/=(value<numtY>&&c){return operator/=(c);}
	protected:
		hist&FillWithValues(const value<numtY>&v){
			for(Point&P:m_data)P.varY()=v;
			return *this;
		}
		hist&FillWithValues(value<numtY>&&v){
			return FillWithValues(v);
		}
		hist&imbibe(const hist& second){//the uncertanties are set standard way (sqrt)
			for(int i=0,n=size();i<n;i++){
				if(m_data[i].X().val()==second[i].X().val()){
					m_data[i].varY()=value<numtY>(m_data[i].Y().val()+second[i].Y().val());
				}else
					throw Exception<hist>("Cannot imbibe histogram. bins differ");
			}
			return *this;
		}
		hist&fill(numtX v){
			for(Point&P:m_data)
				if(P.X().contains(v))
					P.varY()=value<numtY>(P.Y().val()+numtY(1));
			return *this;
		}
		
	};
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator+(const hist<numtX,numtY>&a,const hist<numtX,numtY>&b){auto res=a;res+=b;return res;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator+(hist<numtX,numtY>&&a,const hist<numtX,numtY>&b){auto res=a;res+=b;return res;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator+(const hist<numtX,numtY>&a,hist<numtX,numtY>&&b){auto res=a;res+=b;return res;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator+(hist<numtX,numtY>&&a,hist<numtX,numtY>&&b){auto res=a;res+=b;return res;}
	
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator+(const hist<numtX,numtY>&a,function<numtY(numtX)>b){auto res=a;res+=b;return res;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator+(hist<numtX,numtY>&&a,function<numtY(numtX)>b){auto res=a;res+=b;return res;}
	
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator+(const hist<numtX,numtY>&a,const value<numtY>&b){auto res=a;res+=b;return res;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator+(hist<numtX,numtY>&&a,const value<numtY>&b){return a+b;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator+(const hist<numtX,numtY>&a,value<double>&&b){return a+b;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator+(hist<numtX,numtY>&&a,value<double>&&b){return a+b;}
	
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator-(const hist<numtX,numtY>&a,const hist<numtX,numtY>&b){auto res=a;res-=b;return res;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator-(hist<numtX,numtY>&&a,const hist<numtX,numtY>&b){auto res=a;res-=b;return res;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator-(const hist<numtX,numtY>&a,hist<numtX,numtY>&&b){auto res=a;res-=b;return res;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator-(hist<numtX,numtY>&&a,hist<numtX,numtY>&&b){auto res=a;res-=b;return res;}
	
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator-(const hist<numtX,numtY>&a,function<numtY(numtX)>b){auto res=a;res-=b;return res;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator-(hist<numtX,numtY>&&a,function<numtY(numtX)>b){auto res=a;res-=b;return res;}
	
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator-(const hist<numtX,numtY>&a,const value<numtY>&b){auto res=a;res-=b;return res;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator-(hist<numtX,numtY>&&a,const value<numtY>&b){return a-b;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator-(const hist<numtX,numtY>&a,value<numtY>&&b){return a-b;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator-(hist<numtX,numtY>&&a,value<numtY>&&b){return a-b;}
	
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator*(const hist<numtX,numtY>&a,const hist<numtX,numtY>&b){auto res=a;res*=b;return res;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator*(hist<numtX,numtY>&&a,const hist<numtX,numtY>&b){auto res=a;res*=b;return res;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator*(const hist<numtX,numtY>&a,hist<numtX,numtY>&&b){auto res=a;res*=b;return res;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator*(hist<numtX,numtY>&&a,hist<numtX,numtY>&&b){auto res=a;res*=b;return res;}
	
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator*(const hist<numtX,numtY>&a,function<numtY(numtX)>b){auto res=a;res*=b;return res;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator*(hist<numtX,numtY>&&a,function<numtY(numtX)>b){auto res=a;res*=b;return res;}
	
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator*(const hist<numtX,numtY>&a,const value<numtY>&b){auto res=a;res*=b;return res;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator*(hist<numtX,numtY>&&a,const value<numtY>&b){return a*b;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator*(const hist<numtX,numtY>&a,value<numtY>&&b){return a*b;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator*(hist<numtX,numtY>&&a,value<numtY>&&b){return a*b;}
	
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator/(const hist<numtX,numtY>&a,const hist<numtX,numtY>&b){auto res=a;res/=b;return res;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator/(hist<numtX,numtY>&&a,const hist<numtX,numtY>&b){auto res=a;res/=b;return res;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator/(const hist<numtX,numtY>&a,hist<numtX,numtY>&&b){auto res=a;res/=b;return res;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator/(hist<numtX,numtY>&&a,hist<numtX,numtY>&&b){auto res=a;res/=b;return res;}
	
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator/(const hist<numtX,numtY>&a,function<numtY(numtX)>b){auto res=a;res/=b;return res;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator/(hist<numtX,numtY>&&a,function<numtY(numtX)>b){auto res=a;res/=b;return res;}
	
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator/(const hist<numtX,numtY>&a,const value<numtY>&b){auto res=a;res/=b;return res;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator/(hist<numtX,numtY>&&a,const value<numtY>&b){return a/b;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator/(const hist<numtX,numtY>&a,value<numtY>&&b){return a/b;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator/(hist<numtX,numtY>&&a,value<numtY>&&b){return a/b;}
	
	template<class numtX,class numtY=numtX>
	class Distribution1D:public hist<numtX,numtY>{
	public:
		Distribution1D(const initializer_list<value<numtX>>&data):hist<numtX,numtY>(data){hist<numtX,numtY>::FillWithValues(value<numtY>(0));}
		Distribution1D(initializer_list<value<numtX>>&&data):hist<numtX,numtY>(data){hist<numtX,numtY>::FillWithValues(value<numtY>(0));}
		Distribution1D(const vector<value<numtX>>&data):hist<numtX,numtY>(data){hist<numtX,numtY>::FillWithValues(value<numtY>(0));}
		Distribution1D(vector<value<numtX>>&&data):hist<numtX,numtY>(data){hist<numtX,numtY>::FillWithValues(value<numtY>(0));}
		Distribution1D&operator<<(numtX v){
			hist<numtX,numtY>::fill(v);
			return *this;
		}
	};
	
	template<class numtX,class numtY=numtX,class numtZ=numtY>
	class Distribution2D{
	private:
		vector<value<numtX>> m_x_axis;
		vector<value<numtY>> m_y_axis;
		vector<vector<value<numtZ>>> m_data;
		void init(){
			m_data.clear();
			for(size_t i=0,I=m_x_axis.size();i<I;i++){
				m_data.push_back(vector<value<numtZ>>());
				for(size_t j=0,J=m_y_axis.size();j<J;j++)
					m_data[m_data.size()-1].push_back(value<numtZ>(0));
			}
		}
	public:
		Distribution2D(const initializer_list<value<numtX>>&X,const initializer_list<value<numtY>>&Y){
			for(const auto&x:X)m_x_axis.push_back(x);
			for(const auto&y:Y)m_y_axis.push_back(y);
			init();
		}
		Distribution2D(initializer_list<value<numtX>>&&X,initializer_list<value<numtY>>&&Y):Distribution2D(X,Y){}
		Distribution2D(const vector<value<numtX>>&X,const vector<value<numtY>>&Y){
			for(const auto&x:X)m_x_axis.push_back(x);
			for(const auto&y:Y)m_y_axis.push_back(y);
			init();
		}
		Distribution2D(vector<value<numtX>>&&X,vector<value<numtY>>&&Y):Distribution2D(X,Y){}
		virtual ~Distribution2D(){}
		typedef typename vector<vector<value<numtZ>>>::const_iterator const_iterator;
		const_iterator begin()const{return m_data.cbegin();}
		const_iterator cbegin()const{return m_data.cbegin();}
		const_iterator end() const{return m_data.cend();}
		const_iterator cend() const{return m_data.cend();}
		size_t size()const{return m_data.size();}
		vector<value<numtZ>>&operator[](size_t i)const{
			if(size()<=i)throw Exception<Distribution2D>("range check error");
			return const_cast<vector<value<numtZ>>&>(m_data[i]);
		}
		vector<value<numtX>>&X()const{return const_cast<vector<value<numtX>>&>(m_x_axis);}
		vector<value<numtY>>&Y()const{return const_cast<vector<value<numtY>>&>(m_y_axis);}
		class Point{
			friend class Distribution2D;
		private:
			value<numtX> x;
			value<numtY> y;
			value<numtZ> z;
		protected:
			Point(const value<numtX>&_x,const value<numtY>&_y,const value<numtZ>&_z)
				:x(_x),y(_y),z(_z){}
		public:
			virtual ~Point(){}
			value<numtX>&X()const{return const_cast<value<numtX>&>(x);}
			value<numtY>&Y()const{return const_cast<value<numtY>&>(y);}
			value<numtZ>&Z()const{return const_cast<value<numtZ>&>(z);}
		};
		Point operator()(size_t i,size_t j)const{
			if(size()<=i)throw Exception<Distribution2D>("range check error");
			if(m_y_axis.size()<=j)throw Exception<Distribution2D>("range check error");
			return Point(m_x_axis[i],m_y_axis[j],m_data[i][j]);
		}
		void FullCycle(function<void(Point&&)>f)const{
			for(size_t i=0,I=m_x_axis.size();i<I;i++)
				for(size_t j=0,J=m_y_axis.size();j<J;j++)
					f(Point(m_x_axis[i],m_y_axis[j],m_data[i][j]));
		}
		Distribution2D&operator<<(const pair<numtX,numtY>&p){
			for(size_t i=0,I=m_x_axis.size();i<I;i++)if(m_x_axis[i].contains(p.first))
				for(size_t j=0,J=m_y_axis.size();j<J;j++)if(m_y_axis[j].contains(p.second))
					m_data[i][j]=value<numtZ>(m_data[i][j].val()+numtZ(1));
			return *this;
		}
		Distribution2D&operator<<(pair<numtX,numtY>&&p){
			return operator<<(p);
		}
	};
}
#endif