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
#include "interpolate.h"
namespace MathTemplates{
	using namespace std;
	template<class numtX,class numtY=numtX>class point{
	private:
		value<numtX> x;
		value<numtY> y;
	public:
		point(const value<numtX>&pos):x(pos){}
		point(const value<numtX>&&pos):x(pos){}
		point(const value<numtX>&pos,const value<numtY>&val):x(pos),y(val){}
		point(const value<numtX>&&pos,const value<numtY>&val):x(pos),y(val){}
		point(const value<numtX>&pos,const value<numtY>&&val):x(pos),y(val){}
		point(const value<numtX>&&pos,const value<numtY>&&val):x(pos),y(val){}
		point(const point&source):x(source.x),y(source.y){}
		const value<numtX>&X()const{return x;}
		const value<numtY>&Y()const{return y;}
	protected:
		value<numtX>&__X(){return x;}
		value<numtY>&__Y(){return y;}
	};
	template<class numtX,class numtY=numtX> class point_editable_x:public point<numtX,numtY>{
	public:
		point_editable_x(const point<numtX,numtY>&source):point<numtX,numtY>(source){}
		point_editable_x(const point<numtX,numtY>&&source):point<numtX,numtY>(source){}
		value<numtX>&varX(){return point<numtX,numtY>::__X();}
	};
	template<class numtX,class numtY=numtX> class point_editable_y:public point<numtX,numtY>{
	public:
		point_editable_y(const point<numtX,numtY>&source):point<numtX,numtY>(source){}
		point_editable_y(const point<numtX,numtY>&&source):point<numtX,numtY>(source){}
		value<numtY>&varY(){return point<numtX,numtY>::__Y();}
	};
	template<class numtX,class numtY=numtX> class point_editable:public point<numtX,numtY>{
	public:
		point_editable(const point<numtX,numtY>&source):point<numtX,numtY>(source){}
		point_editable(const point<numtX,numtY>&&source):point<numtX,numtY>(source){}
		value<numtX>&varX(){return point<numtX,numtY>::__X();}
		value<numtY>&varY(){return point<numtX,numtY>::__Y();}
	};
	template<class numtX,class numtY=numtX,class numtZ=numtY>
	class point3d{
	private:
		value<numtX> x;
		value<numtY> y;
		value<numtZ> z;
	public: 
		point3d(const value<numtX>&_x,const value<numtY>&_y,const value<numtZ>&_z)
		:x(_x),y(_y),z(_z){}
		point3d(const value<numtX>&&_x,const value<numtY>&&_y,const value<numtZ>&&_z)
		:x(_x),y(_y),z(_z){}
		point3d(const numtX _x,const numtY _y,const numtZ _z)
		:x(_x,0),y(_y,0),z(_z,0){}
		virtual ~point3d(){}
		const value<numtX>&X()const{return x;}
		const value<numtY>&Y()const{return y;}
		const value<numtZ>&Z()const{return z;}
	};
	
	template<class numtX,class numtY=numtX>class hist{
	public:
		typedef point_editable_y<numtX,numtY> Point;
	private:
		vector<Point> m_data;
	public:
		//Copying
		hist(){}
		hist(const initializer_list<value<numtX>>&data){
			for(const auto& v:data)m_data.push_back(Point(v));
		}
		hist(const initializer_list<value<numtX>>&&data):hist(data){}
		hist(const vector<value<numtX>>&data){
			for(const auto& v:data)m_data.push_back(Point(v));
		}
		hist(const vector<value<numtX>>&&data):hist(data){}
		hist(const initializer_list<point<numtX,numtY>>&data){
			for(const auto& P:data)m_data.push_back(Point(P));
		}
		hist(const initializer_list<point<numtX,numtY>>&&data):hist(data){}
		hist(const vector<point<numtX,numtY>>&data){
			for(const auto& P:data)m_data.push_back(Point(P));
		}
		hist(const vector<point<numtX,numtY>>&&data):hist(data){}
		hist(const initializer_list<Point>&data){
			for(const auto& P:data)m_data.push_back(P);
		}
		hist(const initializer_list<Point>&&data):hist(data){}
		hist(const vector<Point>&data){
			for(const auto& P:data)m_data.push_back(P);
		}
		hist(const vector<Point>&&data):hist(data){}
		hist(const hist&source){
			for(const Point& P:source.m_data)m_data.push_back(P);
		}
		virtual ~hist(){}
		hist&operator=(const hist& source){
			m_data.clear();
			for(const Point& P:source.m_data)m_data.push_back(P);
			return *this;
		}
		//Iterating
		typedef typename vector<Point>::const_iterator const_iterator;
		const_iterator begin()const{return m_data.cbegin();}
		const_iterator cbegin()const{return m_data.cbegin();}
		const_iterator end() const{return m_data.cend();}
		const_iterator cend() const{return m_data.cend();}
		size_t size()const{return m_data.size();}
		const Point&operator[](const size_t i)const{
			if(m_data.size()<=i)
				throw Exception<hist>("range check error");
			return m_data[i];
		}
		Point&Bin(const size_t i){
			if(m_data.size()<=i)
				throw Exception<hist>("range check error");
			return m_data[i];
		}
		//Simple transform
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
		LinearInterpolation<numtX,numtY> Line()const{
			LinearInterpolation<numtX,numtY> res;
			for(const Point&P:m_data)
				res<<make_pair(P.X().val(),P.Y().val());
			return res;
		}
		
		//Comparing
		numtY ChiSq_only_y_error(function<numtY(numtX)>f,size_t paramcount)const{
			numtY res=0,k=numtY(size())-numtY(paramcount);
			if(k<=0)throw Exception<hist>("ChiSq error: too few points or too many parameters");
			for(const Point&p:m_data)
				res+=pow((p.Y().val()-f(p.X().val()))/p.Y().delta(),2);
			return res/k;
		}
		numtY ChiSq_only_y_error(const hist<numtX,numtY>&H,size_t paramcount)const{
			if(H.size()!=size())throw Exception<hist>("ChiSq: histograms size mismatch");
			numtY res=0,k=numtY(size())-numtY(paramcount);
			if(k<=0)throw Exception<hist>("ChiSq error: too few points or too many parameters");
			for(size_t i=0,n=size();i<n;i++){
				if(
					(m_data[i].X().val()!=H[i].X().val())&&
					(m_data[i].X().delta()!=H[i].X().delta())
				)throw Exception<hist>("ChiSq: histograms bins mismatch");
				res+=pow(m_data[i].Y().val()-H[i].Y().val(),2)/(pow(m_data[i].Y().delta(),2)+pow(H[i].Y().delta(),2));
			}
			return res/k;
		}
		numtY ChiSq_only_y_error(hist<numtX,numtY>&&H,size_t paramcount)const{
			ChiSq_only_y_error(H,paramcount);
		}
		numtY ChiSq(function<numtY(numtX)>f,size_t paramcount)const{
			numtY res=0,k=numtY(size())-numtY(paramcount);
			if(k<=0)throw Exception<hist>("ChiSq error: too few points or too many parameters");
			for(const Point&p:m_data){
				numtY w=pow(p.Y().delta(),2)+pow((f(p.X().min())-f(p.X().max()))/numtY(2),2);
				res+=pow(p.Y().val()-f(p.X().val()),2)/w;
			}
			return res/k;
		}
		//Arithmetic operations
		hist&operator+=(const hist& second){
			for(size_t i=0,n=size();i<n;i++){
				if(m_data[i].X().val()==second[i].X().val()){
					m_data[i].varY()+=second[i].Y();
				}else
					throw Exception<hist>("Cannot add histogram. bins differ");
			}
			return *this;
		}
		hist&operator+=(const function<numtY(numtX)>f){
			for(size_t i=0,n=size();i<n;i++)
				m_data[i].varY()+=value<numtY>(f(m_data[i].X().val()),0);
			return *this;
		}
		hist&operator+=(const value<numtY>&c){
			for(size_t i=0,n=size();i<n;i++)m_data[i].varY()+=c;
			return *this;
		}
		hist&operator+=(const value<numtY>&&c){return operator+=(c);}
		
		hist&operator-=(const hist& second){
			for(size_t i=0,n=size();i<n;i++){
				if(m_data[i].X().val()==second[i].X().val()){
					m_data[i].varY()-=second[i].Y();
				}else
					throw Exception<hist>("Cannot add histogram. bins differ");
			}
			return *this;
		}
		hist&operator-=(const function<numtY(numtX)>f){
			for(size_t i=0,n=size();i<n;i++)
				m_data[i].varY()-=value<numtY>(f(m_data[i].X().val()),0);
			return *this;
		}
		hist&operator-=(const value<numtY>&c){
			for(size_t i=0,n=size();i<n;i++)m_data[i].varY()-=c;
			return *this;
		}
		hist&operator-=(const value<numtY>&&c){return operator-=(c);}
		
		hist&operator*=(const hist& second){
			for(size_t i=0,n=size();i<n;i++){
				if(m_data[i].X().val()==second[i].X().val()){
					m_data[i].varY()*=second[i].Y();
				}else
					throw Exception<hist>("Cannot add histogram. bins differ");
			}
			return *this;
		}
		hist&operator*=(const function<numtY(numtX)>f){
			for(size_t i=0,n=size();i<n;i++)
				m_data[i].varY()*=value<numtY>(f(m_data[i].X().val()),0);
			return *this;
		}
		hist&operator*=(const value<numtY>&c){
			for(size_t i=0,n=size();i<n;i++)m_data[i].varY()*=c;
			return *this;
		}
		hist&operator*=(const value<numtY>&&c){return operator*=(c);}
		
		hist&operator/=(const hist& second){
			for(size_t i=0,n=size();i<n;i++){
				if(m_data[i].X().val()==second[i].X().val()){
					m_data[i].varY()*=second[i].Y();
				}else
					throw Exception<hist>("Cannot add histogram. bins differ");
			}
			return *this;
		}
		hist&operator/=(const function<numtY(numtX)>f){
			for(size_t i=0,n=size();i<n;i++)
				m_data[i].varY()/=value<numtY>(f(m_data[i].X().val()),0);
			return *this;
		}
		hist&operator/=(const value<numtY>&c){
			for(size_t i=0,n=size();i<n;i++)m_data[i].varY()/=c;
			return *this;
		}
		hist&operator/=(const value<numtY>&&c){return operator/=(c);}
	
		//Advanced transrormations
		hist&FillWithValues(const value<numtY>&v){
			for(Point&P:m_data)P.varY()=v;
			return *this;
		}
		hist&FillWithValues(const value<numtY>&&v){
			return FillWithValues(v);
		}
		hist Scale(const size_t sc_x)const{
			//uncertanties are set to standard sqrt
			vector<value<numtX>> new_x,sorted_x;
			for(const auto&item:m_data)
				InsertSorted(item.X(),sorted_x,std_size(sorted_x),std_insert(sorted_x,value<numtX>));
			for(size_t i=sc_x-1,n=sorted_x.size();i<n;i+=sc_x){
				auto min=sorted_x[i+1-sc_x].min();
				auto max=sorted_x[i].max();
				new_x.push_back(value<numtX>((max+min)/numtX(2),(max-min)/numtX(2)));
			}
			hist res(new_x);
			for(size_t i=0;i<new_x.size();i++){
				numtY v=0;
				for(size_t ii=0;ii<sc_x;ii++)
					v+=m_data[i*sc_x+ii].Y().val();
				res.Bin(i).varY()=value<numtY>(v,sqrt(v));
			}
			return res;
		}
		hist&imbibe(const hist& second){
			//Sum of histograms. the uncertanties are set to standard sqrt
			for(int i=0,n=size();i<n;i++){
				if(m_data[i].X().val()==second[i].X().val()){
					m_data[i].varY()=value<numtY>(m_data[i].Y().val()+second[i].Y().val());
				}else
					throw Exception<hist>("Cannot imbibe histogram. bins differ");
			}
			return *this;
		}
	};
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator+(const hist<numtX,numtY>&a,const hist<numtX,numtY>&b){auto res=a;res+=b;return res;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator+(const hist<numtX,numtY>&&a,const hist<numtX,numtY>&b){auto res=a;res+=b;return res;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator+(const hist<numtX,numtY>&a,const hist<numtX,numtY>&&b){auto res=a;res+=b;return res;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator+(const hist<numtX,numtY>&&a,const hist<numtX,numtY>&&b){auto res=a;res+=b;return res;}
	
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator+(const hist<numtX,numtY>&a,const function<numtY(numtX)>b){auto res=a;res+=b;return res;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator+(const hist<numtX,numtY>&&a,const function<numtY(numtX)>b){auto res=a;res+=b;return res;}
	
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator+(const hist<numtX,numtY>&a,const value<numtY>&b){auto res=a;res+=b;return res;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator+(const hist<numtX,numtY>&&a,const value<numtY>&b){return a+b;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator+(const hist<numtX,numtY>&a,const value<double>&&b){return a+b;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator+(const hist<numtX,numtY>&&a,const value<double>&&b){return a+b;}
	
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator-(const hist<numtX,numtY>&a,const hist<numtX,numtY>&b){auto res=a;res-=b;return res;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator-(const hist<numtX,numtY>&&a,const hist<numtX,numtY>&b){auto res=a;res-=b;return res;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator-(const hist<numtX,numtY>&a,const hist<numtX,numtY>&&b){auto res=a;res-=b;return res;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator-(const hist<numtX,numtY>&&a,const hist<numtX,numtY>&&b){auto res=a;res-=b;return res;}
	
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator-(const hist<numtX,numtY>&a,const function<numtY(numtX)>b){auto res=a;res-=b;return res;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator-(const hist<numtX,numtY>&&a,const function<numtY(numtX)>b){auto res=a;res-=b;return res;}
	
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator-(const hist<numtX,numtY>&a,const value<numtY>&b){auto res=a;res-=b;return res;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator-(const hist<numtX,numtY>&&a,const value<numtY>&b){return a-b;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator-(const hist<numtX,numtY>&a,const value<numtY>&&b){return a-b;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator-(const hist<numtX,numtY>&&a,const value<numtY>&&b){return a-b;}
	
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator*(const hist<numtX,numtY>&a,const hist<numtX,numtY>&b){auto res=a;res*=b;return res;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator*(const hist<numtX,numtY>&&a,const hist<numtX,numtY>&b){auto res=a;res*=b;return res;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator*(const hist<numtX,numtY>&a,const hist<numtX,numtY>&&b){auto res=a;res*=b;return res;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator*(const hist<numtX,numtY>&&a,const hist<numtX,numtY>&&b){auto res=a;res*=b;return res;}
	
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator*(const hist<numtX,numtY>&a,const function<numtY(numtX)>b){auto res=a;res*=b;return res;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator*(const hist<numtX,numtY>&&a,const function<numtY(numtX)>b){auto res=a;res*=b;return res;}
	
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator*(const hist<numtX,numtY>&a,const value<numtY>&b){auto res=a;res*=b;return res;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator*(const hist<numtX,numtY>&&a,const value<numtY>&b){return a*b;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator*(const hist<numtX,numtY>&a,const value<numtY>&&b){return a*b;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator*(const hist<numtX,numtY>&&a,const value<numtY>&&b){return a*b;}
	
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator/(const hist<numtX,numtY>&a,const hist<numtX,numtY>&b){auto res=a;res/=b;return res;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator/(const hist<numtX,numtY>&&a,const hist<numtX,numtY>&b){auto res=a;res/=b;return res;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator/(const hist<numtX,numtY>&a,const hist<numtX,numtY>&&b){auto res=a;res/=b;return res;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator/(const hist<numtX,numtY>&&a,const hist<numtX,numtY>&&b){auto res=a;res/=b;return res;}
	
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator/(const hist<numtX,numtY>&a,const function<numtY(numtX)>b){auto res=a;res/=b;return res;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator/(const hist<numtX,numtY>&&a,const function<numtY(numtX)>b){auto res=a;res/=b;return res;}
	
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator/(const hist<numtX,numtY>&a,const value<numtY>&b){auto res=a;res/=b;return res;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator/(const hist<numtX,numtY>&&a,const value<numtY>&b){return a/b;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator/(const hist<numtX,numtY>&a,const value<numtY>&&b){return a/b;}
	template<class numtX,class numtY>
	inline hist<numtX,numtY> operator/(const hist<numtX,numtY>&&a,const value<numtY>&&b){return a/b;}

	template<class numtX,class numtY=numtX,class numtZ=numtY>
	class hist2d{
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
		hist2d(const initializer_list<value<numtX>>&X,const initializer_list<value<numtY>>&Y){
			for(const auto&x:X)m_x_axis.push_back(x);
			for(const auto&y:Y)m_y_axis.push_back(y);
			init();
		}
		hist2d(const initializer_list<value<numtX>>&&X,const initializer_list<value<numtY>>&&Y):hist2d(X,Y){}
		hist2d(const vector<value<numtX>>&X,const vector<value<numtY>>&Y){
			for(const auto&x:X)m_x_axis.push_back(x);
			for(const auto&y:Y)m_y_axis.push_back(y);
			init();
		}
		hist2d(const vector<value<numtX>>&&X,const vector<value<numtY>>&&Y):hist2d(X,Y){}
		hist2d():hist2d({},{}){}
		hist2d(const hist2d&source){
			for(const auto&item:source.m_x_axis)m_x_axis.push_back(item);
			for(const auto&item:source.m_y_axis)m_y_axis.push_back(item);
			for(size_t i=0,I=source.m_data.size();i<I;i++){
				m_data.push_back(vector<value<numtZ>>());
				for(const auto&item:source.m_data[i])
					m_data[i].push_back(item);
			}
		}
		virtual ~hist2d(){}
		typedef typename vector<vector<value<numtZ>>>::const_iterator const_iterator;
		const_iterator begin()const{return m_data.cbegin();}
		const_iterator cbegin()const{return m_data.cbegin();}
		const_iterator end() const{return m_data.cend();}
		const_iterator cend() const{return m_data.cend();}
		size_t size()const{return m_data.size();}
		const vector<value<numtZ>>&operator[](const size_t i)const{
			if(size()<=i)throw Exception<hist2d>("range check error");
			return m_data[i];
		}
		value<numtZ>&Bin(const size_t i,const size_t j){
			if(size()<=i)throw Exception<hist2d>("range check error");
			if(m_data[i].size()<=j)throw Exception<hist2d>("range check error");
			return m_data[i][j];
		}
		const vector<value<numtX>>&X()const{return m_x_axis;}
		const vector<value<numtY>>&Y()const{return m_y_axis;}
		
		hist2d Scale(const size_t sc_x,const size_t sc_y)const{
			//uncertanties are set to standard sqrt
			vector<value<numtX>> new_x,sorted_x;
			for(const auto&item:X())
				InsertSorted(item,sorted_x,std_size(sorted_x),std_insert(sorted_x,value<numtX>));
			for(size_t i=sc_x-1,n=sorted_x.size();i<n;i+=sc_x){
				auto min=sorted_x[i+1-sc_x].min();
				auto max=sorted_x[i].max();
				new_x.push_back(value<numtX>((max+min)/numtX(2),(max-min)/numtX(2)));
			}
			vector<value<numtY>> new_y,sorted_y;
			for(const auto&item:Y())
				InsertSorted(item,sorted_y,std_size(sorted_y),std_insert(sorted_y,value<numtY>));
			for(size_t i=sc_y-1,n=sorted_y.size();i<n;i+=sc_y){
				auto min=sorted_y[i+1-sc_y].min();
				auto max=sorted_y[i].max();
				new_y.push_back(value<numtY>((max+min)/numtY(2),(max-min)/numtY(2)));
			}
			hist2d res(new_x,new_y);
			for(size_t i=0;i<new_x.size();i++)for(size_t j=0;j<new_y.size();j++){
				numtZ v=0;
				for(size_t ii=0;ii<sc_x;ii++)
					for(size_t jj=0;jj<sc_y;jj++)
						v+=m_data[i*sc_x+ii][j*sc_y+jj].val();
				res.Bin(i,j)=value<numtZ>(v,sqrt(v));
			}
			return res;
		}
		hist2d&imbibe(const hist2d& second){
			//Sum of histograms. the uncertanties are set to standard sqrt
			if((X().size()!=second.X().size())||(Y().size()!=second.Y().size()))
				throw Exception<hist2d>("cannot imbibe second histogram: bins differ");
			for(int i=0,n=size();i<n;i++)for(int j=0,m=m_data[i].size();j<m;j++)
				m_data[i][j]=value<numtY>(m_data[i][j].val()+second[i][j].val());
			return *this;
		}
		
		const point3d<numtX,numtY,numtZ> operator()(const size_t i,const size_t j)const{
			if(size()<=i)throw Exception<hist2d>("range check error");
			if(m_y_axis.size()<=j)throw Exception<hist2d>("range check error");
			return point3d<numtX,numtY,numtZ>(m_x_axis[i],m_y_axis[j],m_data[i][j]);
		}
		void FullCycle(const function<void(const point3d<numtX,numtY,numtZ>&)>f)const{
			for(size_t i=0,I=m_x_axis.size();i<I;i++)
				for(size_t j=0,J=m_y_axis.size();j<J;j++){
					point3d<numtX,numtY,numtZ> P(m_x_axis[i],m_y_axis[j],m_data[i][j]);
					f(P);
				}
		}
		void FullCycle(const function<void(const value<numtX>&,const value<numtY>&,const value<numtZ>&)>f)const{
			for(size_t i=0,I=m_x_axis.size();i<I;i++)
				for(size_t j=0,J=m_y_axis.size();j<J;j++)
					f(m_x_axis[i],m_y_axis[j],m_data[i][j]);
		}
		void FullCycleVar(const function<void(const value<numtX>&,const value<numtY>&,value<numtZ>&)>f){
			for(size_t i=0,I=m_x_axis.size();i<I;i++)
				for(size_t j=0,J=m_y_axis.size();j<J;j++)
					f(m_x_axis[i],m_y_axis[j],m_data[i][j]);
		}
	};

	
	template<class numtX,class numtY=numtX>
	class Distribution1D:public hist<numtX,numtY>{
	public:
		Distribution1D(const initializer_list<value<numtX>>&data):hist<numtX,numtY>(data){hist<numtX,numtY>::FillWithValues(value<numtY>(0));}
		Distribution1D(const initializer_list<value<numtX>>&&data):hist<numtX,numtY>(data){hist<numtX,numtY>::FillWithValues(value<numtY>(0));}
		Distribution1D(const vector<value<numtX>>&data):hist<numtX,numtY>(data){hist<numtX,numtY>::FillWithValues(value<numtY>(0));}
		Distribution1D(const vector<value<numtX>>&&data):hist<numtX,numtY>(data){hist<numtX,numtY>::FillWithValues(value<numtY>(0));}
		Distribution1D&operator<<(const numtX v){
			for(size_t i=0,n=hist<double>::size();i<n;i++)
				if(hist<double>::Bin(i).X().contains(v))
					hist<double>::Bin(i).varY()=value<numtY>(hist<double>::Bin(i).Y().val()+numtY(1));
				return *this;
		}
	};
	template<class numtX,class numtY=numtX,class numtZ=numtY>
	class Distribution2D:public hist2d<numtX,numtY,numtZ>{
	public:
		Distribution2D(const initializer_list<value<numtX>>&X,const initializer_list<value<numtY>>&Y):hist2d<numtX,numtY,numtZ>(X,Y){}
		Distribution2D(const initializer_list<value<numtX>>&&X,initializer_list<value<numtY>>&&Y):hist2d<numtX,numtY,numtZ>(X,Y){}
		Distribution2D(const vector<value<numtX>>&X,const vector<value<numtY>>&Y):hist2d<numtX,numtY,numtZ>(X,Y){}
		Distribution2D(const vector<value<numtX>>&&X,vector<value<numtY>>&&Y):hist2d<numtX,numtY,numtZ>(X,Y){}
		Distribution2D&operator<<(const pair<numtX,numtY>&p){
			for(size_t i=0,I=hist2d<numtX,numtY,numtZ>::X().size();i<I;i++)if(hist2d<numtX,numtY,numtZ>::X()[i].contains(p.first))
				for(size_t j=0,J=hist2d<numtX,numtY,numtZ>::Y().size();j<J;j++)if(hist2d<numtX,numtY,numtZ>::Y()[j].contains(p.second))
					hist2d<numtX,numtY,numtZ>::Bin(i,j)=value<numtZ>(hist2d<numtX,numtY,numtZ>::operator[](i)[j].val()+numtZ(1));
				return *this;
		}
		Distribution2D&operator<<(const pair<numtX,numtY>&&p){
			return operator<<(p);
		}
	};
}
#endif