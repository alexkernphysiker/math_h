// this file is distributed under 
// MIT license
#ifndef PPRNSANGJGXVGERD
#define PPRNSANGJGXVGERD
#include <vector>
#include <utility>
#include <functional>
#include "functions.h"
#include "error.h"
namespace MathTemplates{
	using namespace std;
	template<class comparable, class indexer=std::vector<comparable>>
	int  WhereToInsert(const int from,const int to,const indexer X,const comparable x){
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
	template<class comparable, class indexer=std::vector<comparable>>
	int Search(const int from,const int to,const indexer X,const comparable x){
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
	template<class comparable,class indexer, class Size, class Insert>
	void InsertSorted(const comparable x,indexer&X,const Size size,const Insert insert){
		insert(WhereToInsert(0,size()-1,X,x),x);
	}
	#define std_size(vector) [&vector](){return (vector).size();}
	#define std_insert(vector,type) [&vector](int pos,type x){(vector).insert((vector).begin()+pos,x);}
	#define field_size(vector)  [this](){return (vector).size();}
	#define field_insert(vector,type)  [this](int pos,type x){(vector).insert((vector).begin()+pos,x);}
	
	namespace details{
		template<class numX, class numY=numX, class Pair=std::pair<numX,numY>>
		inline bool operator<(const Pair&a,const Pair&b){
			return a.first<b.first;
		}
		template<class numX, class numY=numX, class Pair=std::pair<numX,numY>>
		inline bool operator<(const Pair&a,const Pair&&b){return a<b;}
		template<class numX, class numY=numX, class Pair=std::pair<numX,numY>>
		inline bool operator<(const Pair&&a,const Pair&b){return a<b;}
		template<class numX, class numY=numX, class Pair=std::pair<numX,numY>>
		inline bool operator<(const Pair&&a,const Pair&&b){return a<b;}
		
		template<class numX, class numY=numX, class Pair=std::pair<numX,numY>>
		inline bool operator>(const Pair&a,const Pair&b){
			return a.first>b.first;
		}
		template<class numX, class numY=numX, class Pair=std::pair<numX,numY>>
		inline bool operator>(const Pair&a,const Pair&&b){return a>b;}
		template<class numX, class numY=numX, class Pair=std::pair<numX,numY>>
		inline bool operator>(const Pair&&a,const Pair&b){return a>b;}
		template<class numX, class numY=numX, class Pair=std::pair<numX,numY>>
		inline bool operator>(const Pair&&a,const Pair&&b){return a>b;}
	}
	template<class numX, class numY=numX, class PairIndexer=std::vector<std::pair<numX,numY>>>
	numY  Interpolate_Linear(const int from,const int to,const PairIndexer tbl,const numX x){
		using namespace std;
		using namespace details;
		if(x==tbl[from].first)
			return tbl[from].second;
		if(x==tbl[to].first)
			return tbl[to].second;
		pair<numX,numY> p=make_pair(x,numY(0));
		int i=WhereToInsert<pair<numX,numY>>(from,to,tbl,p);
		if((i<=from)||(i>to))
			throw Exception<std::pair<numX,numY>>("Attempt to interpolate outside the given region.");
		numX k=(x-tbl[i-1].first)/(tbl[i].first-tbl[i-1].first);
		return tbl[i-1].second+(tbl[i].second-tbl[i-1].second)*numY(k);
	}
	template<class numX, class numY=numX, class PairIndexer=std::vector<std::pair<numX,numY>>,class Size=std::function<int()>>
	numY InterpolateLinear(const numX x,const PairIndexer tbl,const Size size){
		return Interpolate_Linear(0,size()-1,tbl,x);
	}
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
	template<class numX, class numY=numX>
	class SortedPoints{
	public:
		typedef pair<numX,numY> Point;
	private:
		vector<Point> data;
	public:
		SortedPoints(){}
		SortedPoints &operator<<(const Point&p){
			InsertSorted(p,data,field_size(data),field_insert(data,Point));
			return *this;
		}
		SortedPoints &operator<<(const Point&&p){
			return operator<<(p);
		}
		SortedPoints(const initializer_list<Point>&points){
			for(const Point&p:points)operator<<(p);
		}
		SortedPoints(const initializer_list<Point>&&points):SortedPoints(points){}
		SortedPoints(const SortedPoints&points){
			for(const Point&p:points.data)operator<<(p);
		}
		SortedPoints(const SortedPoints&&points):SortedPoints(points){}
		SortedPoints(const function<numY(numX)> f,const vector<numX>&chain){
			for(numX x:chain)operator<<(make_pair(x,f(x)));
		}
		SortedPoints(const function<numY(numX)> f,const vector<numX>&&chain):SortedPoints(f,chain){}
		SortedPoints(const function<numY(numX)> f,const initializer_list<numX>&&chain){
			for(numX x:chain)operator<<(make_pair(x,f(x)));
		}
		virtual ~SortedPoints(){}
		//Points access
		int size()const{return data.size();}
		const Point&operator[](const int i)const{
			if(size()<=i)
				throw Exception<SortedPoints>("Range check error");
			return data[i];
		}
		//typedef typename vector<Point>::iterator iterator;
		typedef typename vector<Point>::const_iterator const_iterator;
		//iterator begin(){return data.begin();}
		const_iterator begin()const{return data.begin();}
		const_iterator cbegin()const{return data.cbegin();}
		//iterator end(){return data.end();}
		const_iterator end() const{return data.end();}
		const_iterator cend() const{return data.cend();}
		const Point&left()const{
			if(size()<1)
				throw Exception<SortedPoints>("Attempt to obtain empty properties.");
			return data[0];
		}
		const Point&right()const{
			if(size()<1)
				throw Exception<SortedPoints>("Attempt to obtain empty properties.");
			return data[size()-1];
		}
		numX min()const{return left().first;}
		numX max()const{return right().first;}
		//Arithmetic actions
		const SortedPoints<numY,numX> Transponate()const{
			SortedPoints<numY,numX> res;
			for(const Point&p:data)
				res<<make_pair(p.second,p.first);
			return res;
		}
		SortedPoints &operator+=(const numY value){
			for(Point&point:data)
				point.second+=value;
			return *this;
		}
		SortedPoints &operator-=(const numY value){
			for(Point&point:data)
				point.second-=value;
			return *this;
		}
		SortedPoints &operator*=(const numY value){
			for(Point&point:data)
				point.second*=value;
			return *this;
		}
		SortedPoints &operator/=(const numY value){
			for(Point&point:data)
				point.second/=value;
			return *this;
		}
		SortedPoints &transform(const std::function<numY(numY)>F){
			for(Point&point:data)
				point.second=F(point.second);
			return *this;
		}
		SortedPoints &transform(const std::function<numY(numX,numY)>F){
			for(Point&point:data)
				point.second=F(point.first,point.second);
			return *this;
		}
	};
	template<class numX, class numY=numX>
	class LinearInterpolation:public SortedPoints<numX,numY>,public IFunction<numX,numY>{
	public:
		typedef pair<numX,numY> Point;
	private:
		vector<Point> data;
	public:
		LinearInterpolation(){}
		LinearInterpolation(const SortedPoints<numX,numY>&source):SortedPoints<numX,numY>(source){}
		LinearInterpolation(const SortedPoints<numX,numY>&&source):SortedPoints<numX,numY>(source){}
		virtual ~LinearInterpolation(){}
		virtual numY operator()(const numX x)const override{
			using namespace details;
			return InterpolateLinear<numX,numY>(x,*this,field_size(*this));
		}
		virtual const function<numY(numX)> func()const override{
			return [this](double x){return operator()(x);};
		}
	};
};
#endif
