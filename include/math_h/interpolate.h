// this file is distributed under 
// MIT license
#ifndef PPRNSANGJGXVGERD
#define PPRNSANGJGXVGERD
#include <vector>
#include <utility>
#include <functional>
#include "error.h"
template<class comparable, class indexer=std::vector<comparable>>
int  WhereToInsert(int from, int to, indexer X, comparable x){
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
int Search(int from, int to, indexer X, comparable x){
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
void InsertSorted(comparable x,indexer X,Size size,Insert insert){
	insert(WhereToInsert(0,size()-1,X,x),x);
}
#define std_size(vector) [&vector](){return vector.size();}
#define std_insert(vector,type) [&vector](int pos,type x){vector.insert(vector.begin()+pos,x);}
#define field_size(vector)  [this](){return vector.size();}
#define field_insert(vector,type)  [this](int pos,type x){vector.insert(vector.begin()+pos,x);}
template<class numX, class indexerX, class numY=numX, class indexerY=indexerX>
numY  Interpolate_Linear(int from, int to, indexerX X, indexerY Y, numX x){
	if(x==X[from])return Y[from];
	if(x==X[to])return Y[to];
	int i=WhereToInsert(from,to,X,x);
	if((i<=from)||(i>to))
		throw Error<std::pair<numX,numY>>("Attempt to interpolate outside given region.");
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
	if((i<=from)||(i>to))
		throw Error<std::pair<numX,numY>>("Attempt to interpolate outside the given region.");
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
	LinearInterpolation &operator<<(Point&&p){
		InsertSorted(p,data,field_size(data),field_insert(data,Point));
		return *this;
	}
	LinearInterpolation &operator+=(numY value){
		for(Point&point:data)
			point.second+=value;
		return *this;
	}
	LinearInterpolation &operator-=(numY value){
		for(Point&point:data)
			point.second-=value;
		return *this;
	}
	LinearInterpolation &operator*=(numY value){
		for(Point&point:data)
			point.second*=value;
		return *this;
	}
	LinearInterpolation &operator/=(numY value){
		for(Point&point:data)
			point.second/=value;
		return *this;
	}
	LinearInterpolation &transform(std::function<numY(numY)>F){
		for(Point&point:data)
			point.second=F(point.second);
		return *this;
	}
	LinearInterpolation &transform(std::function<numY(numX,numY)>F){
		for(Point&point:data)
			point.second=F(point.first,point.second);
		return *this;
	}
	LinearInterpolation(std::vector<Point>&&points){
		for(Point&p:points)
			operator<<(static_cast<Point&&>(p));
	}
	virtual ~LinearInterpolation(){}
	int size()const{return data.size();}
	numX min()const{
		if(size()<1)
			throw Error<LinearInterpolation>("Attempt to obtain empty properties.");
		return data[0].first;
	}
	numX max()const{
		if(size()<1)
			throw Error<LinearInterpolation>("Attempt to obtain empty properties.");
		return data[size()-1].first;
	}
	numY operator()(numX x)const{
		using namespace details;
		return InterpolateLinear2<numX,numY>(x,data,field_size(data));
	}
	std::function<numY(numX)> func()const{
		return [this](double x){return operator()(x);};
	}
	Point&operator[](int i)const{return data[i];}
	typedef typename std::vector<Point>::iterator iterator;
	typedef typename std::vector<Point>::const_iterator const_iterator;
	iterator begin(){return data.begin();}
	const_iterator begin()const{return data.begin();}
	const_iterator cbegin()const{return data.cbegin();}
	iterator end(){return data.end();}
	const_iterator end() const{return data.end();}
	const_iterator cend() const{return data.cend();}
};
template<class numX, class numY=numX>
class LinearInterpolation_fixedsize{
private:
	numX *X;
	numY *Y;
	int n;
public:
	LinearInterpolation_fixedsize(const LinearInterpolation_fixedsize &source){
		n=source.n;
		X=new numX[n];
		Y=new numY[n];
		for(int i=0;i<n;i++){
			X[i]=source.X[i];
			Y[i]=source.Y[i];
		}
	}
	LinearInterpolation_fixedsize(std::function<numY(numX)> F, numX from,numX to,int points){
		if(from>=to)
			throw Error<LinearInterpolation_fixedsize>("Wrong ranges given in constructor");
		if(points<2)
			throw Error<LinearInterpolation_fixedsize>("Too few bins given in constructor");
		numX step= (to-from)/numX(points-1);
		n=points;
		X=new numX[n];
		Y=new numY[n];
		for(int i=0;i<n;i++){
			numX x=from+step*numX(i);
			X[i]=x;
			Y[i]=F(x);
		}
	}
	LinearInterpolation_fixedsize(numX from,numX to,int points):LinearInterpolation_fixedsize([](numX){return numY(0);},from,to,points){}
	virtual ~LinearInterpolation_fixedsize(){
		delete[] X;
		delete[] Y;
	}
	int size()const{return n;}
	numX min()const{
		if(size()<1)
			throw Error<LinearInterpolation_fixedsize>("Attempt to obtain empty proterties");
		return X[0];
	}
	numX max()const{
		if(size()<1)
			throw Error<LinearInterpolation_fixedsize>("Attempt to obtain empty proterties");
		return X[size()-1];
	}
	numY operator()(numX x)const{
		using namespace details;
		return Interpolate_Linear(0,n-1,X,Y,x);
	}
	std::function<numY(numX)> func()const{
		return [this](double x){return operator()(x);};
	}
	numX getX(int i)const{
		if((i<0)||(i>=size()))
			throw Error<LinearInterpolation_fixedsize>("Range check error");
		return X[i];
	}
	numY getY(int i)const{
		if((i<0)||(i>=size()))
			throw Error<LinearInterpolation_fixedsize<numX,numY>>("Range check error");
		return Y[i];
	}
	void setY(int i,numY v){
		if((i<0)||(i>=size()))
			throw Error<LinearInterpolation_fixedsize<numX,numY>>("Range check error");
		Y[i]=v;
	}
protected:
	void setX(int i,numX v){
		if((i<0)||(i>=size()))
			throw Error<LinearInterpolation_fixedsize<numX,numY>>("Range check error");
		X[i]=v;
	}
};
template<class numX, class numY=numX>
class Distribution:public LinearInterpolation_fixedsize<numX,numY>{
private:
	numX bindelta;
public:
	Distribution(numX from,numX to,int bincount):LinearInterpolation_fixedsize<numX,numY>(
		from+(to-from)/numX(bincount*2),
		to-(to-from)/numX(bincount*2),
		bincount
	){
		bindelta=(to-from)/numX(bincount*2);
	}
	virtual ~Distribution(){}
	Distribution &AddValue(numX v){
		for(int i=0,n=LinearInterpolation_fixedsize<numX,numY>::size();i<n;i++){
			numX x=LinearInterpolation_fixedsize<numX,numY>::getX(i);
			if(((x-bindelta)<=v)&&(v<(x+bindelta)))
				LinearInterpolation_fixedsize<numX,numY>::setY(i,LinearInterpolation_fixedsize<numX,numY>::getY(i)+1);
		}
		return *this;
	}
	numX BinWidth()const{return bindelta*numX(2);}
};
#endif
