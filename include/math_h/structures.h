// this file is distributed under 
// MIT license
#ifndef ____STRUCTURES_H_____
#	define ____STRUCTURES_H_____
#include <vector>
#include <functional>
#include "sigma.h"
#include "error.h"
namespace MathTemplates{
	using namespace std;
	template<class numX>
	const vector<numX> ChainWithStep(const numX from,const numX step,const numX to){
		if(from>=to)throw Exception<vector<numX>>("wrong binning ranges");
		if(step<=0)throw Exception<vector<numX>>("wrong binning step");
		vector<numX> res;
		for(numX x=from;x<=to;x+=step)res.push_back(x);
		return res;
	}
	template<class numX>
	const vector<numX> ChainWithCount(const size_t cont,const numX from,const numX to){
		if(from>=to)throw Exception<vector<numX>>("wrong binning ranges");
		if(0==cont)throw Exception<vector<numX>>("wrong bins count");
		numX step=(to-from)/numX(cont);
		vector<numX> res;
		for(numX x=from;x<=to;x+=step)res.push_back(x);
		return res;
	}
	template<class numt>
	const vector<value<numt>> BinsByStep(const numt from,const numt step,const numt to){
		if(0>=step)throw Exception<vector<value<numt>>>("wrong bin width");
		if(to<=from)throw Exception<vector<value<numt>>>("wrong range");
		numt delta=step/numt(2);
		vector<value<numt>> res;
		for(numt x=from+delta;x<to;x+=step)
			res.push_back(value<numt>(x,delta));
		return res;
	}
	template<class numt>
	const vector<value<numt>> BinsByCount(const size_t count,const numt from,const numt to){
		if(0==count)throw Exception<vector<value<numt>>>("wrong bins count");
		return BinsByStep(from,(to-from)/numt(count),to);
	}
	
	namespace details{
		template<class comparable, class indexer=vector<comparable>>
		int  WhereToInsert(const int from,const int to,const indexer&X,const comparable&x){
			if(from>to) return from;
			int beg=from,end=to;
			if(x>X[end]) return end+1;
			if(x<X[beg]) return beg;
			while(1<(end-beg)){
				int mid=(beg+end)/2;
				if(x<X[mid]) end=mid;
				else
					if(x>X[mid]) beg=mid;
					else return mid;
			}
			return end;
		}
		template<class comparable, class indexer=vector<comparable>>
		int ClosestPosition(const int from,const int to,const indexer&X,const comparable&x){
			if(from>to) return from-1;
			int beg=from,end=to;
			if(x>X[end]) return from-1;
			if(x<X[beg]) return from-1;
			while(1<(end-beg)){
				int mid=(beg+end)/2;
				if(x<X[mid]) end=mid;
				else
					if(x>X[mid]) beg=mid;
					else return mid;
			}
			if(X[beg]==x)return beg;
			if(X[end]==x)return end;
			return from-1;
		}
	}
	template<class comparable,class indexer, class Size, class Insert>
	void InsertSorted(const comparable&x,indexer&X,const Size size,const Insert insert){
		insert(details::WhereToInsert(0,size()-1,X,x),x);
	}
	template<class comparable,class indexer, class Size, class Insert>
	void InsertSorted(const comparable&&x,indexer&X,const Size size,const Insert insert){
		InsertSorted(x,X,size,insert);
	}
	#define std_size(vector) [&vector](){return (vector).size();}
	#define std_insert(vector,type) [&vector](int pos,type x){(vector).insert((vector).begin()+pos,x);}
	#define field_size(vector)  [this](){return (vector).size();}
	#define field_insert(vector,type)  [this](int pos,type x){(vector).insert((vector).begin()+pos,x);}
	
	template<class comparable>
	class SortedChain{
	private:
		vector<comparable> data;
	public:
		SortedChain(){}
		SortedChain&operator<<(const comparable&p){
			InsertSorted(p,data,field_size(data),field_insert(data,comparable));
			return *this;
		}
		SortedChain&operator<<(const comparable&&p){
			return operator<<(p);
		}
		SortedChain(const initializer_list<comparable>&points){
			for(const auto&p:points)
				operator<<(p);
		}
		SortedChain(const initializer_list<comparable>&&points):SortedChain(points){}
		SortedChain(const vector<comparable>&points){
			for(const auto&p:points)
				operator<<(p);
		}
		SortedChain(const vector<comparable>&&points):SortedChain(points){}
		SortedChain(const SortedChain&points){
			for(const auto&p:points.data)
				operator<<(p);
		}
		SortedChain(const SortedChain&&points):SortedChain(points){}
		virtual ~SortedChain(){}
		SortedChain& operator=(const SortedChain&points){
			data.clear();
			for(const auto&p:points.data)
				operator<<(p);
			return *this;
		}
		
		const size_t size()const{return data.size();}
		const comparable&operator[](const int i)const{
			if(size()<=i)
				throw Exception<SortedChain>("Range check error");
			return data[i];
		}
		const comparable&left()const{
			if(size()<1)
				throw Exception<SortedChain>("Attempt to obtain empty properties.");
			return data[0];
		}
		const comparable&right()const{
			if(size()<1)
				throw Exception<SortedChain>("Attempt to obtain empty properties.");
			return data[size()-1];
		}
		typedef typename vector<comparable>::const_iterator const_iterator;
		const_iterator begin()const{return data.begin();}
		const_iterator cbegin()const{return data.cbegin();}
		const_iterator end() const{return data.end();}
		const_iterator cend() const{return data.cend();}
	protected:
		comparable&accessBin(const size_t i){
			if(data.size()<=i)
				throw Exception<SortedChain>("range check error");
			return data[i];
		}
	};
	
	template<class numtX,class numtY=numtX>class point{
	private:
		numtX x;
		numtY y;
	public:
		point(const numtX&pos):x(pos),y(numtY(0)){}
		point(const numtX&&pos):x(pos),y(numtY(0)){}
		point(const numtX&pos,const numtY&val):x(pos),y(val){}
		point(const numtX&&pos,const numtY&val):x(pos),y(val){}
		point(const numtX&pos,const numtY&&val):x(pos),y(val){}
		point(const numtX&&pos,const numtY&&val):x(pos),y(val){}
		point(const point&source):x(source.x),y(source.y){}
		const numtX&X()const{return x;}
		const numtY&Y()const{return y;}
	protected:
		numtX&__X(){return x;}
		numtY&__Y(){return y;}
	};
	template<class numtX,class numtY=numtX> class point_editable_x:public point<numtX,numtY>{
	public:
		point_editable_x(const point<numtX,numtY>&source):point<numtX,numtY>(source){}
		point_editable_x(const point<numtX,numtY>&&source):point<numtX,numtY>(source){}
		numtX&varX(){return point<numtX,numtY>::__X();}
	};
	template<class numtX,class numtY=numtX> class point_editable_y:public point<numtX,numtY>{
	public:
		point_editable_y(const point<numtX,numtY>&source):point<numtX,numtY>(source){}
		point_editable_y(const point<numtX,numtY>&&source):point<numtX,numtY>(source){}
		numtY&varY(){return point<numtX,numtY>::__Y();}
	};
	template<class numtX,class numtY=numtX> class point_editable:public point<numtX,numtY>{
	public:
		point_editable(const point<numtX,numtY>&source):point<numtX,numtY>(source){}
		point_editable(const point<numtX,numtY>&&source):point<numtX,numtY>(source){}
		numtX&varX(){return point<numtX,numtY>::__X();}
		numtY&varY(){return point<numtX,numtY>::__Y();}
	};
	template<class numX, class numY=numX>
	bool operator<(const point<numX,numY>&a,const point<numX,numY>&b){
		return a.X()<b.X();
	}
	template<class numX, class numY=numX>
	bool operator>(const point<numX,numY>&a,const point<numX,numY>&b){
		return a.X()>b.X();
	}
	template<class numX, class numY=numX>
	class SortedPoints:public SortedChain<point_editable_y<numX,numY>>{
	public:
		typedef function<numY(const numX&)> Func;
		SortedPoints(){}
		SortedPoints&operator<<(const point<numX,numY>&p){
			SortedChain<point_editable_y<numX,numY>>::operator<<(point_editable_y<numX,numY>(p));
			return *this;
		}
		SortedPoints&operator<<(const point<numX,numY>&&p){return operator<<(p);}
		SortedPoints(const initializer_list<point<numX,numY>>&points){
			for(const auto&p:points)operator<<(p);
		}
		SortedPoints(const initializer_list<point<numX,numY>>&&points):SortedPoints(points){}
		
		SortedPoints(const vector<point<numX,numY>>&points){
			for(const auto&p:points)operator<<(p);
		}
		SortedPoints(const vector<point<numX,numY>>&&points):SortedPoints(points){}
		
		SortedPoints(const Func f,const vector<numX>&chain){
			for(numX x:chain)operator<<(point<numX,numY>(x,f(x)));
		}
		SortedPoints(const Func f,const vector<numX>&&chain):SortedPoints(f,chain){}
		SortedPoints(const Func f,const initializer_list<numX>&&chain){
			for(numX x:chain)operator<<(point<numX,numY>(x,f(x)));
		}
		
		virtual ~SortedPoints(){}
		
		SortedPoints(const SortedChain<point_editable_y<numX,numY>>&points):SortedChain<point_editable_y<numX,numY>>(points){}
		SortedPoints(const SortedChain<point_editable_y<numX,numY>>&&points):SortedPoints(points){}
		SortedPoints& operator=(const SortedChain<point_editable_y<numX,numY>>&points){
			SortedChain<point_editable_y<numX,numY>>::operator=(points);
			return *this;
		}
		
		point_editable_y<numX,numY>&Bin(const size_t i){return SortedChain<point_editable_y<numX,numY>>::accessBin(i);}

		const  vector<point<numY,numX>> Transponate()const{
			vector<point<numY,numX>> res;
			for(const auto&p:*this)
				res.push_back(point<numX,numY>(p.Y(),p.X()));
			return res;
		}
		const SortedPoints XRange(const numX from,const numX to)const{
			SortedPoints res;
			for(const auto&P:*this)
				if((P.X()>=from)&&(P.X()<=to))
					res<<P;
				return res;
		}
		const SortedPoints YRange(const numY from,const numY to)const{
			SortedPoints res;
			for(const auto&P:*this)
				if((P.Y()>=from)&&(P.Y()<=to))
					res<<P;
				return res;
		}
		SortedPoints&FillWithValues(const numY&v){
			for(size_t i=0,n=this->size();i<n;i++)
				Bin(i).varY()=v;
			return *this;
		}
		SortedPoints&FillWithValues(const numY&&v){
			return FillWithValues(v);
		}
		SortedPoints&Transform(const function<numY(const numX&,const numY&)>&F){
			for(size_t i=0,n=this->size();i<n;i++)
				Bin(i).varY()=F(Bin(i).X(),Bin(i).Y());
			return *this;
		}
		SortedPoints&operator+=(const SortedPoints& second){
			for(size_t i=0,n=this->size();i<n;i++){
				if(Bin(i).X()==second[i].X())
					Bin(i).varY()+=second[i].Y();
				else
					throw Exception<SortedPoints>("Cannot perform arithmetic actions two chains. X coordinates differ");
			}
			return *this;
		}
		SortedPoints&operator+=(const Func f){
			for(size_t i=0,n=this->size();i<n;i++)
				Bin(i).varY()+=f(Bin(i).X());
			return *this;
		}
		SortedPoints&operator+=(const numY&c){
			for(size_t i=0,n=this->size();i<n;i++)
				Bin(i).varY()+=c;
			return *this;
		}
		SortedPoints&operator+=(const numY&&c){
			return operator+=(c);
		}
		
		SortedPoints&operator-=(const SortedPoints& second){
			for(size_t i=0,n=this->size();i<n;i++){
				if(Bin(i).X()==second[i].X())
					Bin(i).varY()-=second[i].Y();
				else
					throw Exception<SortedPoints>("Cannot perform arithmetic actions two chains. X coordinates differ");
			}
			return *this;
		}
		SortedPoints&operator-=(const Func f){
			for(size_t i=0,n=this->size();i<n;i++)
				Bin(i).varY()-=f(Bin(i).X());
			return *this;
		}
		SortedPoints&operator-=(const numY&c){
			for(size_t i=0,n=this->size();i<n;i++)
				Bin(i).varY()-=c;
			return *this;
		}
		SortedPoints&operator-=(const numY&&c){
			return operator-=(c);
		}
		
		SortedPoints&operator*=(const SortedPoints& second){
			for(size_t i=0,n=this->size();i<n;i++){
				if(Bin(i).X()==second[i].X())
					Bin(i).varY()*=second[i].Y();
				else
					throw Exception<SortedPoints>("Cannot perform arithmetic actions two chains. X coordinates differ");
			}
			return *this;
		}
		SortedPoints&operator*=(const Func f){
			for(size_t i=0,n=this->size();i<n;i++)
				Bin(i).varY()*=f(Bin(i).X());
			return *this;
		}
		SortedPoints&operator*=(const numY&c){
			for(size_t i=0,n=this->size();i<n;i++)
				Bin(i).varY()*=c;
			return *this;
		}
		SortedPoints&operator*=(const numY&&c){
			return operator*=(c);
		}
		
		SortedPoints&operator/=(const SortedPoints& second){
			for(size_t i=0,n=this->size();i<n;i++){
				if(Bin(i).X()==second[i].X())
					Bin(i).varY()*=second[i].Y();
				else
					throw Exception<SortedPoints>("Cannot perform arithmetic actions two chains. X coordinates differ");
			}
			return *this;
		}
		SortedPoints&operator/=(const Func f){
			for(size_t i=0,n=this->size();i<n;i++)
				Bin(i).varY()/=f(Bin(i).X());
			return *this;
		}
		SortedPoints&operator/=(const numY&c){
			for(size_t i=0,n=this->size();i<n;i++)
				Bin(i).varY()/=c;
			return *this;
		}
		SortedPoints&operator/=(const numY&&c){
			return operator/=(c);
		}
	};
	template<class numtX,class numtY>
	inline const SortedPoints<numtX,numtY> operator+(const SortedPoints<numtX,numtY>&a,const SortedPoints<numtX,numtY>&b){auto res=a;res+=b;return res;}
	template<class numtX,class numtY>
	inline const SortedPoints<numtX,numtY> operator+(const SortedPoints<numtX,numtY>&&a,const SortedPoints<numtX,numtY>&b){auto res=a;res+=b;return res;}
	template<class numtX,class numtY>
	inline const SortedPoints<numtX,numtY> operator+(const SortedPoints<numtX,numtY>&a,const SortedPoints<numtX,numtY>&&b){auto res=a;res+=b;return res;}
	template<class numtX,class numtY>
	inline const SortedPoints<numtX,numtY> operator+(const SortedPoints<numtX,numtY>&&a,const SortedPoints<numtX,numtY>&&b){auto res=a;res+=b;return res;}
	
	template<class numtX,class numtY>
	inline const SortedPoints<numtX,numtY> operator+(const SortedPoints<numtX,numtY>&a,const typename SortedPoints<numtX,numtY>::Func b){auto res=a;res+=b;return res;}
	template<class numtX,class numtY>
	inline const SortedPoints<numtX,numtY> operator+(const SortedPoints<numtX,numtY>&&a,const typename SortedPoints<numtX,numtY>::Func b){auto res=a;res+=b;return res;}
	
	template<class numtX,class numtY>
	inline const SortedPoints<numtX,numtY> operator+(const SortedPoints<numtX,numtY>&a,const numtY&b){auto res=a;res+=b;return res;}
	template<class numtX,class numtY>
	inline const SortedPoints<numtX,numtY> operator+(const SortedPoints<numtX,numtY>&&a,const numtY&b){return a+b;}
	template<class numtX,class numtY>
	inline const SortedPoints<numtX,numtY> operator+(const SortedPoints<numtX,numtY>&a,const numtY&&b){return a+b;}
	template<class numtX,class numtY>
	inline const SortedPoints<numtX,numtY> operator+(const SortedPoints<numtX,numtY>&&a,const numtY&&b){return a+b;}
	
	template<class numtX,class numtY>
	inline const SortedPoints<numtX,numtY> operator-(const SortedPoints<numtX,numtY>&a,const SortedPoints<numtX,numtY>&b){auto res=a;res-=b;return res;}
	template<class numtX,class numtY>
	inline const SortedPoints<numtX,numtY> operator-(const SortedPoints<numtX,numtY>&&a,const SortedPoints<numtX,numtY>&b){auto res=a;res-=b;return res;}
	template<class numtX,class numtY>
	inline const SortedPoints<numtX,numtY> operator-(const SortedPoints<numtX,numtY>&a,const SortedPoints<numtX,numtY>&&b){auto res=a;res-=b;return res;}
	template<class numtX,class numtY>
	inline const SortedPoints<numtX,numtY> operator-(const SortedPoints<numtX,numtY>&&a,const SortedPoints<numtX,numtY>&&b){auto res=a;res-=b;return res;}
	
	template<class numtX,class numtY>
	inline const SortedPoints<numtX,numtY> operator-(const SortedPoints<numtX,numtY>&a,const typename SortedPoints<numtX,numtY>::Func b){auto res=a;res-=b;return res;}
	template<class numtX,class numtY>
	inline const SortedPoints<numtX,numtY> operator-(const SortedPoints<numtX,numtY>&&a,const typename SortedPoints<numtX,numtY>::Func b){auto res=a;res-=b;return res;}
	
	template<class numtX,class numtY>
	inline const SortedPoints<numtX,numtY> operator-(const SortedPoints<numtX,numtY>&a,const numtY&b){auto res=a;res-=b;return res;}
	template<class numtX,class numtY>
	inline const SortedPoints<numtX,numtY> operator-(const SortedPoints<numtX,numtY>&&a,const numtY&b){return a-b;}
	template<class numtX,class numtY>
	inline const SortedPoints<numtX,numtY> operator-(const SortedPoints<numtX,numtY>&a,const numtY&&b){return a-b;}
	template<class numtX,class numtY>
	inline const SortedPoints<numtX,numtY> operator-(const SortedPoints<numtX,numtY>&&a,const numtY&&b){return a-b;}
	
	template<class numtX,class numtY>
	inline const SortedPoints<numtX,numtY> operator*(const SortedPoints<numtX,numtY>&a,const SortedPoints<numtX,numtY>&b){auto res=a;res*=b;return res;}
	template<class numtX,class numtY>
	inline const SortedPoints<numtX,numtY> operator*(const SortedPoints<numtX,numtY>&&a,const SortedPoints<numtX,numtY>&b){auto res=a;res*=b;return res;}
	template<class numtX,class numtY>
	inline const SortedPoints<numtX,numtY> operator*(const SortedPoints<numtX,numtY>&a,const SortedPoints<numtX,numtY>&&b){auto res=a;res*=b;return res;}
	template<class numtX,class numtY>
	inline const SortedPoints<numtX,numtY> operator*(const SortedPoints<numtX,numtY>&&a,const SortedPoints<numtX,numtY>&&b){auto res=a;res*=b;return res;}
	
	template<class numtX,class numtY>
	inline const SortedPoints<numtX,numtY> operator*(const SortedPoints<numtX,numtY>&a,const typename SortedPoints<numtX,numtY>::Func b){auto res=a;res*=b;return res;}
	template<class numtX,class numtY>
	inline const SortedPoints<numtX,numtY> operator*(const SortedPoints<numtX,numtY>&&a,const typename SortedPoints<numtX,numtY>::Func b){auto res=a;res*=b;return res;}
	
	template<class numtX,class numtY>
	inline const SortedPoints<numtX,numtY> operator*(const SortedPoints<numtX,numtY>&a,const numtY&b){auto res=a;res*=b;return res;}
	template<class numtX,class numtY>
	inline const SortedPoints<numtX,numtY> operator*(const SortedPoints<numtX,numtY>&&a,const numtY&b){return a*b;}
	template<class numtX,class numtY>
	inline const SortedPoints<numtX,numtY> operator*(const SortedPoints<numtX,numtY>&a,const numtY&&b){return a*b;}
	template<class numtX,class numtY>
	inline const SortedPoints<numtX,numtY> operator*(const SortedPoints<numtX,numtY>&&a,const numtY&&b){return a*b;}
	
	template<class numtX,class numtY>
	inline const SortedPoints<numtX,numtY> operator/(const SortedPoints<numtX,numtY>&a,const SortedPoints<numtX,numtY>&b){auto res=a;res/=b;return res;}
	template<class numtX,class numtY>
	inline const SortedPoints<numtX,numtY> operator/(const SortedPoints<numtX,numtY>&&a,const SortedPoints<numtX,numtY>&b){auto res=a;res/=b;return res;}
	template<class numtX,class numtY>
	inline const SortedPoints<numtX,numtY> operator/(const SortedPoints<numtX,numtY>&a,const SortedPoints<numtX,numtY>&&b){auto res=a;res/=b;return res;}
	template<class numtX,class numtY>
	inline const SortedPoints<numtX,numtY> operator/(const SortedPoints<numtX,numtY>&&a,const SortedPoints<numtX,numtY>&&b){auto res=a;res/=b;return res;}
	
	template<class numtX,class numtY>
	inline const SortedPoints<numtX,numtY> operator/(const SortedPoints<numtX,numtY>&a,const typename SortedPoints<numtX,numtY>::Func b){auto res=a;res/=b;return res;}
	template<class numtX,class numtY>
	inline const SortedPoints<numtX,numtY> operator/(const SortedPoints<numtX,numtY>&&a,const typename SortedPoints<numtX,numtY>::Func b){auto res=a;res/=b;return res;}
	
	template<class numtX,class numtY>
	inline const SortedPoints<numtX,numtY> operator/(const SortedPoints<numtX,numtY>&a,const numtY&b){auto res=a;res/=b;return res;}
	template<class numtX,class numtY>
	inline const SortedPoints<numtX,numtY> operator/(const SortedPoints<numtX,numtY>&&a,const numtY&b){return a/b;}
	template<class numtX,class numtY>
	inline const SortedPoints<numtX,numtY> operator/(const SortedPoints<numtX,numtY>&a,const numtY&&b){return a/b;}
	template<class numtX,class numtY>
	inline const SortedPoints<numtX,numtY> operator/(const SortedPoints<numtX,numtY>&&a,const numtY&&b){return a/b;}
	
	template<class numtX,class numtY=numtX>
	class hist:public SortedPoints<value<numtX>,value<numtY>>{
	public:
		typedef point<value<numtX>,value<numtY>> Point;
		hist(){}
		hist(const initializer_list<value<numtX>>&data){
			for(const auto& v:data)this->operator<<(Point(v));
		}
		hist(const initializer_list<value<numtX>>&&data):hist(data){}
		hist(const vector<value<numtX>>&data){
			for(const auto& v:data)this->operator<<(Point(v));
		}
		hist(const vector<value<numtX>>&&data):hist(data){}
		hist(const initializer_list<point<numtX,numtY>>&data){
			for(const auto& P:data)this->operator<<(Point(P));
		}
		hist(const initializer_list<point<numtX,numtY>>&&data):hist(data){}
		hist(const vector<point<numtX,numtY>>&data){
			for(const auto& P:data)this->operator<<(Point(P));
		}
		hist(const vector<point<numtX,numtY>>&&data):hist(data){}
		hist(const initializer_list<Point>&data){
			for(const auto& P:data)this->operator<<(P);
		}
		hist(const initializer_list<Point>&&data):hist(data){}
		hist(const vector<Point>&data){
			for(const auto& P:data)this->operator<<(P);
		}
		hist(const vector<Point>&&data):hist(data){}
		hist(const SortedPoints<value<numtX>,value<numtY>>&source):SortedPoints<value<numtX>,value<numtY>>(source){}
		virtual ~hist(){}
		hist& operator=(const SortedPoints<value<numtX>,value<numtY>>&points){
			SortedPoints<value<numtX>,value<numtY>>::operator=(points);
			return *this;
		}
		//Simple transform
		hist CloneEmptyBins()const{
			vector<Point> initer;
			for(const Point&P:*this)initer.push_back(P);
			return hist(initer);
		}
		numtY Total()const{
			numtY res=0;
			for(const Point&P:*this)res+=P.Y().val();
			return res;
		}
		value<numtY> TotalSum()const{
			value<numtY> res=0;
			for(const Point&P:*this)res+=P.Y();
			return res;
		}
		const SortedPoints<numtX,numtY> Line()const{
			SortedPoints<numtX,numtY> res;
			for(const Point&P:*this)
				res<<point<numtX,numtY>(P.X().val(),P.Y().val());
			return res;
		}
		//Advanced transrormations
		const hist Scale(const size_t sc_x)const{
			//uncertanties are set to standard sqrt
			vector<value<numtX>> new_x,sorted_x;
			for(const auto&item:*this)
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
					v+=this->operator[](i*sc_x+ii).Y().val();
				res.Bin(i).varY()=value<numtY>::std_error(v);
			}
			return res;
		}
		hist&imbibe(const hist& second){
			//Sum of histograms. the uncertanties are set to standard sqrt
			for(int i=0,n=this->size();i<n;i++){
				if(this->operator[](i).X()==second[i].X()){
					this->Bin(i).varY()=value<numtY>::std_error(this->operator[](i).Y().val()+second[i].Y().val());
				}else
					throw Exception<hist>("Cannot imbibe histogram. bins differ");
			}
			return *this;
		}
	};
	
	
	
	
	template<class numtX,class numtY=numtX,class numtZ=numtY>
	class point3d{
	private:
		numtX x;
		numtY y;
		numtZ z;
	public: 
		point3d(const numtX&_x,const numtY&_y,const numtZ&_z)
		:x(_x),y(_y),z(_z){}
		point3d(const numtX&&_x,const numtY&&_y,const numtZ&&_z)
		:x(_x),y(_y),z(_z){}
		virtual ~point3d(){}
		const numtX&X()const{return x;}
		const numtY&Y()const{return y;}
		const numtZ&Z()const{return z;}
	};
	template<class numtX,class numtY=numtX,class numtZ=numtY>
	class BiSortedPoints{
	private:
		SortedChain<numtX> m_x_axis;
		SortedChain<numtY> m_y_axis;
		vector<vector<numtZ>> m_data;
		void init(){
			m_data.clear();
			for(size_t i=0,I=m_x_axis.size();i<I;i++){
				m_data.push_back(vector<numtZ>());
				for(size_t j=0,J=m_y_axis.size();j<J;j++)
					m_data[m_data.size()-1].push_back(numtZ(0));
			}
		}
	public:
		const SortedChain<numtX>&X()const{return m_x_axis;}
		const SortedChain<numtY>&Y()const{return m_y_axis;}
		BiSortedPoints(const initializer_list<numtX>&X,const initializer_list<numtY>&Y):m_x_axis(X),m_y_axis(Y){init();}
		BiSortedPoints(const initializer_list<numtX>&&X,const initializer_list<numtY>&&Y):m_x_axis(X),m_y_axis(Y){init();}
		BiSortedPoints(const vector<numtX>&X,const vector<numtY>&Y):m_x_axis(X),m_y_axis(Y){init();}
		BiSortedPoints(const vector<numtX>&&X,const vector<numtY>&&Y):m_x_axis(X),m_y_axis(Y){init();}
		BiSortedPoints():BiSortedPoints({},{}){}
		BiSortedPoints(const BiSortedPoints&source):m_x_axis(source.X()),m_y_axis(source.Y()){
			for(size_t i=0,I=source.m_data.size();i<I;i++){
				m_data.push_back(vector<numtZ>());
				for(const auto&item:source.m_data[i])
					m_data[i].push_back(item);
			}
		}
		virtual ~BiSortedPoints(){}
		typedef typename vector<vector<numtZ>>::const_iterator const_iterator;
		const_iterator begin()const{return m_data.cbegin();}
		const_iterator cbegin()const{return m_data.cbegin();}
		const_iterator end() const{return m_data.cend();}
		const_iterator cend() const{return m_data.cend();}
		size_t size()const{return m_data.size();}
		const vector<numtZ>&operator[](const size_t i)const{
			if(size()<=i)throw Exception<BiSortedPoints>("range check error");
			return m_data[i];
		}
		numtZ&Bin(const size_t i,const size_t j){
			if(size()<=i)throw Exception<BiSortedPoints>("range check error");
			if(m_data[i].size()<=j)throw Exception<BiSortedPoints>("range check error");
			return m_data[i][j];
		}
		const point3d<numtX,numtY,numtZ> operator()(const size_t i,const size_t j)const{
			if(size()<=i)throw Exception<BiSortedPoints>("range check error");
			if(m_y_axis.size()<=j)throw Exception<BiSortedPoints>("range check error");
			return point3d<numtX,numtY,numtZ>(m_x_axis[i],m_y_axis[j],m_data[i][j]);
		}
		void FullCycle(const function<void(const point3d<numtX,numtY,numtZ>&)>f)const{
			for(size_t i=0,I=m_x_axis.size();i<I;i++)
				for(size_t j=0,J=m_y_axis.size();j<J;j++){
					point3d<numtX,numtY,numtZ> P(m_x_axis[i],m_y_axis[j],m_data[i][j]);
					f(P);
				}
		}
		void FullCycle(const function<void(const numtX&,const numtY&,const numtZ&)>f)const{
			for(size_t i=0,I=m_x_axis.size();i<I;i++)
				for(size_t j=0,J=m_y_axis.size();j<J;j++)
					f(m_x_axis[i],m_y_axis[j],m_data[i][j]);
		}
		void FullCycleVar(const function<void(const numtX&,const numtY&,numtZ&)>f){
			for(size_t i=0,I=m_x_axis.size();i<I;i++)
				for(size_t j=0,J=m_y_axis.size();j<J;j++)
					f(m_x_axis[i],m_y_axis[j],m_data[i][j]);
		}
		
	};
	template<class numtX,class numtY=numtX,class numtZ=numtY>
	class hist2d:public BiSortedPoints<value<numtX>,value<numtY>,value<numtZ>>{
	public:
		hist2d(const initializer_list<value<numtX>>&X,const initializer_list<value<numtY>>&Y):BiSortedPoints<value<numtX>,value<numtY>,value<numtZ>>(X,Y){}
		hist2d(const initializer_list<value<numtX>>&&X,const initializer_list<value<numtY>>&&Y):hist2d(X,Y){}
		hist2d(const vector<value<numtX>>&X,const vector<value<numtY>>&Y):BiSortedPoints<value<numtX>,value<numtY>,value<numtZ>>(X,Y){}
		hist2d(const vector<value<numtX>>&&X,const vector<value<numtY>>&&Y):hist2d(X,Y){}
		hist2d():hist2d({},{}){}
		hist2d(const BiSortedPoints<value<numtX>,value<numtY>,value<numtZ>>&source):BiSortedPoints<value<numtX>,value<numtY>,value<numtZ>>(source){}
		virtual ~hist2d(){}
		
		const hist2d Scale(const size_t sc_x,const size_t sc_y)const{
			//uncertanties are set to standard sqrt
			vector<value<numtX>> new_x,sorted_x;
			for(const auto&item:this->X())
				InsertSorted(item,sorted_x,std_size(sorted_x),std_insert(sorted_x,value<numtX>));
			for(size_t i=sc_x-1,n=sorted_x.size();i<n;i+=sc_x){
				auto min=sorted_x[i+1-sc_x].min();
				auto max=sorted_x[i].max();
				new_x.push_back(value<numtX>((max+min)/numtX(2),(max-min)/numtX(2)));
			}
			vector<value<numtY>> new_y,sorted_y;
			for(const auto&item:this->Y())
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
						v+=this->operator[](i*sc_x+ii)[j*sc_y+jj].val();
					res.Bin(i,j)=value<numtZ>(v,sqrt(v));
			}
			return res;
		}
		hist2d&imbibe(const hist2d& second){
			//Sum of histograms. the uncertanties are set to standard sqrt
			if((this->X().size()!=second.X().size())||(this->Y().size()!=second.Y().size()))
				throw Exception<hist2d>("cannot imbibe second histogram: bins differ");
			for(int i=0,n=this->size();i<n;i++)for(int j=0,m=this->operator[](i).size();j<m;j++)
				this->Bin(i,j)=value<numtY>::std_error(this->operator[](i)[j].val()+second[i][j].val());
			return *this;
		}
	};
	
	
	
	
	
	
	template<class numtX,class numtY=numtX>
	class Distribution1D:public hist<numtX,numtY>{
	private:
		unsigned long long counter;
	public:
		Distribution1D(const initializer_list<value<numtX>>&data):hist<numtX,numtY>(data){this->FillWithValues(value<numtY>(0));counter=0;}
		Distribution1D(const initializer_list<value<numtX>>&&data):hist<numtX,numtY>(data){this->FillWithValues(value<numtY>(0));counter=0;}
		Distribution1D(const vector<value<numtX>>&data):hist<numtX,numtY>(data){this->FillWithValues(value<numtY>(0));counter=0;}
		Distribution1D(const vector<value<numtX>>&&data):hist<numtX,numtY>(data){this->FillWithValues(value<numtY>(0));counter=0;}
		Distribution1D&Fill(const numtX& v){
			counter++;
			for(size_t i=0,n=this->size();i<n;i++)
				if(this->Bin(i).X().contains(v))
					this->Bin(i).varY()=value<numtY>::std_error(this->Bin(i).Y().val()+numtY(1));
			return *this;
		}
		Distribution1D&Fill(const numtX&&v){
			return Fill(v);
		}
		const unsigned long long Entries()const{
			return counter;
		}
	};
	template<class numtX,class numtY=numtX,class numtZ=numtY>
	class Distribution2D:public hist2d<numtX,numtY,numtZ>{
	private:
		unsigned long long counter;
	public:
		Distribution2D(const initializer_list<value<numtX>>&X,const initializer_list<value<numtY>>&Y):hist2d<numtX,numtY,numtZ>(X,Y){counter=0;}
		Distribution2D(const initializer_list<value<numtX>>&&X,initializer_list<value<numtY>>&&Y):hist2d<numtX,numtY,numtZ>(X,Y){counter=0;}
		Distribution2D(const vector<value<numtX>>&X,const vector<value<numtY>>&Y):hist2d<numtX,numtY,numtZ>(X,Y){counter=0;}
		Distribution2D(const vector<value<numtX>>&&X,vector<value<numtY>>&&Y):hist2d<numtX,numtY,numtZ>(X,Y){counter=0;}
		Distribution2D&Fill(const pair<numtX,numtY>&p){
			counter++;
			for(size_t i=0,I=this->X().size();i<I;i++)if(this->X()[i].contains(p.first))
				for(size_t j=0,J=this->Y().size();j<J;j++)if(this->Y()[j].contains(p.second))
					this->Bin(i,j)=value<numtZ>::std_error(this->operator[](i)[j].val()+numtZ(1));
				return *this;
		}
		Distribution2D&Fill(const pair<numtX,numtY>&&p){
			return Fill(p);
		}
		const unsigned long long Entries()const{
			return counter;
		}
	};
};

#endif