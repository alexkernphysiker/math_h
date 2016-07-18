// this file is distributed under 
// MIT license
#ifndef ____CHAINS_H_____
#	define ____CHAINS_H_____
#include <vector>
#include <functional>
#include "error.h"
namespace MathTemplates{
	namespace details{
		template<class comparable, class indexer=std::vector<comparable>>
		size_t  WhereToInsert(const size_t from,const size_t to,const indexer&X,const comparable&x){
			if(from>to) return from;
			size_t beg=from,end=to;
			if(x>X[end]) return end+1;
			if(x<X[beg]) return beg;
			while(1<(end-beg)){
				size_t mid=(beg+end)/2;
				if(x<X[mid]) end=mid;
				else
					if(x>X[mid]) beg=mid;
					else return mid;
			}
			return end;
		}
		template<class comparable,class indexer, class Size, class Insert>
		void InsertSorted(const comparable&x,indexer&X,const Size size,const Insert insert){
			if(size()==0) insert(0,x);
			else insert(WhereToInsert(0,size()-1,X,x),x);
		}
		template<class comparable,class indexer, class Size, class Insert>
		void InsertSorted(const comparable&&x,indexer&X,const Size size,const Insert insert){
			InsertSorted(x,X,size,insert);
		}
	}
	#define std_size(vector) [&vector](){return (vector).size();}
	#define std_insert(vector,type) [&vector](int pos,type x){(vector).insert((vector).begin()+pos,x);}
	#define field_size(vector)  [this](){return (vector).size();}
	#define field_insert(vector,type)  [this](int pos,type x){(vector).insert((vector).begin()+pos,x);}
	
	template<class comparable>
	class SortedChain{
	private:
		std::vector<comparable> data;
	public:
		SortedChain(){}
		SortedChain(const SortedChain&points){
			for(const auto&p:points.data)
				data.push_back(p);
		}
		SortedChain(const SortedChain&&points):SortedChain(points){}
		SortedChain&operator<<(const comparable&p){
			details::InsertSorted(p,data,field_size(data),field_insert(data,comparable));
			return *this;
		}
		SortedChain&operator<<(const comparable&&p){
			return operator<<(p);
		}
		SortedChain(const std::initializer_list<comparable>&points){
			for(const auto&p:points)
				operator<<(p);
		}
		virtual ~SortedChain(){}
		SortedChain& operator=(const SortedChain&points){
			data.clear();
			for(const auto&p:points.data)
				data.push_back(p);
			return *this;
		}
		void clear(){data.clear();}
		
		const size_t size()const{return data.size();}
		const comparable&operator[](const size_t i)const{
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
		typedef typename std::vector<comparable>::const_iterator const_iterator;
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
		void append_item_from_sorted(const comparable&c){data.push_back(c);}
	};
	template<class numX>
	const SortedChain<numX> ChainWithStep(const numX from,const numX step,const numX to){
		if(from>=to)throw Exception<SortedChain<numX>>("wrong binning ranges");
		if(step<=0)throw Exception<SortedChain<numX>>("wrong binning step");
		SortedChain<numX> res;
		for(numX x=from;x<=to;x+=step)res<<x;
		return res;
	}
	template<class numX>
	const SortedChain<numX> ChainWithCount(const size_t cont,const numX from,const numX to){
		if(from>=to)throw Exception<SortedChain<numX>>("wrong binning ranges");
		if(0==cont)throw Exception<SortedChain<numX>>("wrong bins count");
		numX step=(to-from)/numX(cont);
		SortedChain<numX> res;
		for(numX x=from;x<=to;x+=step)res<<x;
		return res;
	}
};
#endif