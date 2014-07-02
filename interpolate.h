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
template<class comparable, class indexer>
int  InsertIndex(indexer X, comparable x){
	return WhereToInsert<comparable,indexer>(0,X.Count()-1,X,x);
}
template<class comparable, class indexer>
int  insertIndex(indexer X, comparable x){
	return WhereToInsert<comparable,indexer>(0,X.count()-1,X,x);
}

//Linear interpolation based on previous algorithm
template<class numt, class indexer, class indexer2>
numt  Interpolate_Linear(int from, int to, indexer X, indexer2 Y, numt x){
	int i=WhereToInsert(from,to,X,x);
	if(i<=from)throw;//x out of border; leftside
	if(i>to)throw;//x out of border; rightside
	numt k=(x-X[i-1])/(X[i]-X[i-1]);return Y[i-1]+k*(Y[i]-Y[i-1]);
}
template<class numt, class indexer, class indexer2>
numt  InterpolateLinear(indexer X, indexer2 Y, numt x){
	int i=InsertIndex(X,x);if(i<=0)throw;if(i>=X.Count())throw;
	numt k=(x-X[i-1])/(X[i]-X[i-1]);return Y[i-1]+k*(Y[i]-Y[i-1]);
}
template<class numt, class indexer, class indexer2>
numt  interpolateLinear(indexer X, indexer2 Y, numt x){
	int i=insertIndex(X,x);if(i<=0)throw;if(i>=X.count())throw;
	numt k=(x-X[i-1])/(X[i]-X[i-1]);return Y[i-1]+k*(Y[i]-Y[i-1]);
}
#endif