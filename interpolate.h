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
	numY k=numY(x-X[i-1])/numY(X[i]-X[i-1]);
	return Y[i-1]+k*(Y[i]-Y[i-1]);
}
template<class numt, class Size,class indexerX, class indexerY=indexerX>
numt InterpolateLinear(numt x, indexerX X, indexerY Y,Size size){
	return Interpolate_Linear(0,size()-1,X,Y,x);
}
template<class numX, class numY=numX>
class FuncTable{
private:
	numX* X;
	numY* Y;
	int cnt;
public:
	FuncTable(){cnt=2;X=new numX[cnt];Y=new numY[cnt];}
	FuncTable(int sz){
		if(sz<=1)throw;
		cnt=sz;
		X=new numX[cnt];
		Y=new numY[cnt];
	}
	FuncTable(FuncTable &f){
		cnt=f.cnt;
		X=new numX[cnt];
		Y=new numY[cnt];
		for(int i=0; i<cnt;i++){
			X[i]=(f.X[i]);
			Y[i]=(f.Y[i]);
		}
	}
	virtual ~FuncTable(){delete[] X; delete[] Y;}
	void set(int i,numX x, numY y){X[i]=(x);Y[i]=(y);}
	void sety(int i,numY y){Y[i]=(y);}
	numX getx(int i){return X[i];}
	numY gety(int i){return Y[i];}
	int size(){return cnt;}
	numY operator()(numX x){
		return Interpolate_Linear<numX,numX*,numY,numY*>(0,cnt-1,X,Y,x);
	}
};
template<class numX, class func, class numY=numX>
void fillFuncTable(FuncTable<numX,numY> &tbl,numX from, numX to,func y){
	numX step=(to-from)/(tbl.size()-1);
	for(int i=0;i<tbl.size();i++){
		numX x=from+(step*numX(i));
		tbl.set(i,x,y(x));
	}
}
template<class numX, class numY=numX>
void fillFuncTableWithZeros(FuncTable<numX,numY> &tbl,numX from, numX to){
	numX step=(to-from)/(tbl.size()-1);
	for(int i=0;i<tbl.size();i++){
		numX x=from+(step*numX(i));
		tbl.set(i,x,numY(0));
	}
}
template<class numt>
void add2hist(FuncTable<numt> &tbl, numt value){
	int index=tbl.size()-1;int beg=0;
	if(value<tbl.getx(index)){
		if(value>tbl.getx(beg)){
			while(1<(index-beg)){
				int mid=(beg+index)/2;
				if(value<tbl.getx(mid))index=mid;else
					if(value>=tbl.getx(mid))beg=mid;
			}
		}else index=0;
	}else index=tbl.size();
	if((tbl.size()-1)<=index){
		int i=tbl.size()-1;
		tbl.sety(i,tbl.gety(i)+1);
	}
	else if(0==index){
		tbl.sety(0,tbl.gety(0)+1);
	}else {
	  double x1=tbl.getx(index-1);
	  double x2=tbl.getx(index);
	  if((value-x1)<(x2-value))index--;
	  tbl.sety(index,tbl.gety(index)+1);
	}
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
	pair<numX,numY> p=make_pair(x,numY(0));
	int i=WhereToInsert<pair<numX,numY>>(from,to,tbl,p);
	if(i<=from)
		throw;
	if(i>to)
		throw;
	numY k=numY(x-tbl[i-1].first)/numY(tbl[i].first-tbl[i-1].first);
	return tbl[i-1].second+k*(tbl[i].second-tbl[i-1].second);
}
template<class numX, class numY, class PairIndexer,class Size>
numY InterpolateLinear2(numX x, PairIndexer tbl,Size size){
	return Interpolate_Linear2(0,size()-1,tbl,x);
}
template<class numX, class numY=numX>
class LinearInterpolation{
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
		using namespace std;
		using namespace details;
		return InterpolateLinear2<numX,numY,vector<Point>,function<int()>>(x,data,field_size(data));
	}
	Point operator[](int i){
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
