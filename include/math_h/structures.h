// this file is distributed under 
// MIT license
#ifndef ____STRUCTURES_H_____
#	define ____STRUCTURES_H_____
#include <vector>
#include "error.h"
namespace MathTemplates{
	using namespace std;
	template<class numtX,class numtY=numtX>class point{
	private:
		numtX x;
		numtY y;
	public:
		point(const numtX&pos):x(pos){}
		point(const numtX&&pos):x(pos){}
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
		int Search(const int from,const int to,const indexer&X,const comparable&x){
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
	private:
		vector<point_editable_y<numX,numY>> data;
	public:
		SortedPoints(){}
		SortedPoints&operator<<(const point<numX,numY>&p){
			typedef point_editable_y<numX,numY> P;
			InsertSorted(P(p),data,field_size(data),field_insert(data,P));
			return *this;
		}
		SortedPoints&operator<<(const point<numX,numY>&&p){
			return operator<<(p);
		}
		SortedPoints(const initializer_list<point<numX,numY>>&points){
			for(const auto&p:points)operator<<(p);
		}
		SortedPoints(const initializer_list<point<numX,numY>>&&points):SortedPoints(points){}
		SortedPoints(const SortedPoints&points){
			for(const auto&p:points.data)operator<<(p);
		}
		SortedPoints(const SortedPoints&&points):SortedPoints(points){}
		SortedPoints(const function<numY(numX)> f,const vector<numX>&chain){
			for(numX x:chain)operator<<(point<numX,numY>(x,f(x)));
		}
		SortedPoints(const function<numY(numX)> f,const vector<numX>&&chain):SortedPoints(f,chain){}
		SortedPoints(const function<numY(numX)> f,const initializer_list<numX>&&chain){
			for(numX x:chain)operator<<(point<numX,numY>(x,f(x)));
		}
		virtual ~SortedPoints(){}

		int size()const{return data.size();}
		const point<numX,numY>&operator[](const int i)const{
			if(size()<=i)
				throw Exception<SortedPoints>("Range check error");
			return data[i];
		}
		typedef typename vector<point_editable_y<numX,numY>>::const_iterator const_iterator;
		const_iterator begin()const{return data.begin();}
		const_iterator cbegin()const{return data.cbegin();}
		const_iterator end() const{return data.end();}
		const_iterator cend() const{return data.cend();}
		const point<numX,numY>&left()const{
			if(size()<1)
				throw Exception<SortedPoints>("Attempt to obtain empty properties.");
			return data[0];
		}
		const point<numX,numY>&right()const{
			if(size()<1)
				throw Exception<SortedPoints>("Attempt to obtain empty properties.");
			return data[size()-1];
		}
		numX min()const{return left().X();}
		numX max()const{return right().X();}
		//Arithmetic actions
		const SortedPoints<numY,numX> Transponate()const{
			SortedPoints<numY,numX> res;
			for(const point_editable_y<numX,numY>&p:data)
				res<<point<numX,numY>(p.Y(),p.X());
			return res;
		}
		SortedPoints &operator+=(const numY value){
			for(point_editable_y<numX,numY>&point:data)
				point.varY()+=value;
			return *this;
		}
		SortedPoints &operator-=(const numY value){
			for(point_editable_y<numX,numY>&point:data)
				point.varY()-=value;
			return *this;
		}
		SortedPoints &operator*=(const numY value){
			for(point_editable_y<numX,numY>&point:data)
				point.varY()*=value;
			return *this;
		}
		SortedPoints &operator/=(const numY value){
			for(point_editable_y<numX,numY>&point:data)
				point.varY()/=value;
			return *this;
		}
		SortedPoints &transform(const std::function<numY(numY)>F){
			for(point_editable_y<numX,numY>&point:data)
				point.varY()=F(point.Y());
			return *this;
		}
		SortedPoints &transform(const std::function<numY(numX,numY)>F){
			for(point_editable_y<numX,numY>&point:data)
				point.varY()=F(point.X(),point.Y());
			return *this;
		}
	};
	
};

#endif