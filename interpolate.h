#ifndef PPRNSANGJGXVGERD
#define PPRNSANGJGXVGERD
#include <vector>
#include <utility>
#include <functional>
template<class comparable, class indexer=std::vector<comparable>>
int  WhereToInsert(int from, int to, indexer X, comparable x){
	int beg=from;
	int end=to;
	if(beg>end)
		return beg;
	if(x>X[end])
		return end+1;
	if(x<X[beg])
		return beg;
	while(1<(end-beg)){
		int mid=(beg+end)/2;
		if(x<X[mid])
			end=mid;
		else
			if(x>X[mid])
				beg=mid;
			else
				return mid;
	}
	return end;
}
template<class comparable,class indexer, class Size, class Insert>
void InsertSorted(comparable x,indexer X,Size size,Insert insert){
	insert(WhereToInsert(0,size()-1,X,x),x);
}
#define std_size(vector) [&vector](){return vector.size();}
#define std_insert(vector,type) [&vector](int pos,type x){vector.insert(vector.begin()+pos,x);}
#define field_size(vector)  [this](){return vector.size();}
#define field_insert(vector,type)  [this](int pos,type x){vector.insert(vector.begin()+pos,x);}
template<class numX, class indexerX, class numY=numX, class indexerY=indexerX>
numY  Interpolate_Linear(int from, int to, indexerX X, indexerY Y, numX x){
	int i=WhereToInsert(from,to,X,x);
	if(i<=from)
		throw;
	if(i>to)
		throw;
	numX k=(x-X[i-1])/(X[i]-X[i-1]);
	return Y[i-1]+(Y[i]-Y[i-1])*numY(k);
}
template<class numX, class Size,class indexerX, class numY, class indexerY=indexerX>
numY InterpolateLinear(numX x, indexerX X, indexerY Y,Size size){
	return Interpolate_Linear(0,size()-1,X,Y,x);
}
namespace details{
	template<class numX, class numY=numX, class PairIndexer=std::vector<std::pair<numX,numY>>>
	bool operator<(PairIndexer a,PairIndexer b){
		return a.first<b.first;
	}
	template<class numX, class numY=numX, class PairIndexer=std::vector<std::pair<numX,numY>>>
	bool operator>(PairIndexer a,PairIndexer b){
		return a.first>b.first;
	}
}
template<class numX, class numY=numX, class PairIndexer=std::vector<std::pair<numX,numY>>>
numY  Interpolate_Linear2(int from, int to, PairIndexer tbl, numX x){
	using namespace std;
	using namespace details;
	if(x==tbl[from].first)
		return tbl[from].second;
	if(x==tbl[to].first)
		return tbl[to].second;
	pair<numX,numY> p=make_pair(x,numY(0));
	int i=WhereToInsert<pair<numX,numY>>(from,to,tbl,p);
	if(i<=from)
		throw;
	if(i>to)
		throw;
	numX k=(x-tbl[i-1].first)/(tbl[i].first-tbl[i-1].first);
	return tbl[i-1].second+(tbl[i].second-tbl[i-1].second)*numY(k);
}
template<class numX, class numY=numX, class PairIndexer=std::vector<std::pair<numX,numY>>,class Size=std::function<int()>>
numY InterpolateLinear2(numX x, PairIndexer tbl,Size size){
	return Interpolate_Linear2(0,size()-1,tbl,x);
}
template<class numX, class numY=numX>
class LinearInterpolation{
public:
	typedef std::pair<numX,numY> Point;
private:
	std::vector<Point> data;
public:
	LinearInterpolation(){}
	LinearInterpolation &operator<<(Point p){
		InsertSorted(p,data,field_size(data),field_insert(data,Point));
		return *this;
	}
	numY operator()(numX x){
		using namespace details;
		return InterpolateLinear2<numX,numY>(x,data,field_size(data));
	}
	Point& operator[](int i){
		return data[i];
	}
	int size(){
		return data.size();
	}
	typedef typename std::vector<Point>::iterator iterator;
	typedef typename std::vector<Point>::const_iterator const_iterator;
	iterator begin(){
		return data.begin();
	}
	const_iterator cbegin()const{
		return data.cbegin();
	}
	iterator end(){
		return data.end();
	}
	const_iterator cend() const{
		return data.cend();
	}
};
#endif
