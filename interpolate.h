// https://github.com/alexkernphysiker/math_h
#ifndef ___INTERPOLATE_H
#	define ___INTERPOLATE_H
//X must be sorted (lesser index - lesser value)
// Looks for a position where next value x should be inserted into X for X stayed sorted
template<class comparable, class indexer>
int  WhereToInsert(int from, int to, indexer X, comparable x){
	int beg=from;int end=to;if(beg>end)return beg;//the array seems to be empty
	if(x>X[end])return end+1;// x should be inserted after all elements
	if(x<X[beg])return beg;// x should be inserted before all elements
	while(1<(end-beg)){// looking for place between beg and end
		int mid=(beg+end)/2;
		if(x<X[mid])end=mid;else
			if(x>X[mid])beg=mid;else
				return mid;
	}
	return end;//beg+1, new element x should be inserted between beg and end
}
//Linear interpolation based on previous algorithm
template<class numX, class indexerX, class numY, class indexerY>
numY  Interpolate_Linear(int from, int to, indexerX X, indexerY Y, numX x){
	int i=WhereToInsert(from,to,X,x);
	if(i<=from)throw;//x out of border; leftside
	if(i>to)throw;//x out of border; rightside
	numY k=numY(x-X[i-1])/numY(X[i]-X[i-1]);return Y[i-1]+k*(Y[i]-Y[i-1]);
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
		return Interpolate_Linear(0,cnt-1,X,Y,x);
	}
};
template<class numX, class func, class numY=numX>
void fillFuncTable(FuncTable<numX,numY> &tbl,numX from, numX to,func y){
	numX step=(to-from)/(tbl.size()-1);
	for(int i=0;i<tbl.size();i++){
		numX x=from+(step*numt(i));
		tbl.set(i,x,y(x));
	}
}
template<class numX, class numY=numX>
void fillFuncTableWithZeros(FuncTable<numX,numY> &tbl,numt from, numt to){
	numX step=(to-from)/(tbl.size()-1);
	for(int i=0;i<tbl.size();i++){
		numX x=from+(step*numt(i));
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
#endif
