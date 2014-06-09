/* The newest version of this file can be got by this link
https://github.com/alexkernphysiker/MathLibs/blob/master/functions
There are also some examples how to use these templates in this repository
author: alex_kernphysiker@privatdemail.net */
#ifndef MATH_TEMPLATES_H
#define MATH_TEMPLATES_H
#include <utility>
namespace Math_{
namespace detail{// implementation details
	template <
		typename ret,typename F, typename Tuple,
		bool Done,
		 int Total, int... N
	>
	struct call_impl{
		static ret call(F f, Tuple && t){
			return call_impl<
					ret, F, Tuple,
					Total == 1 + sizeof...(N),
					Total, N..., sizeof...(N)
			>::call(f, std::forward<Tuple>(t));
		}
	};
	template <
			typename ret,typename F, typename Tuple,
			int Total, int... N
	>
	struct call_impl<ret, F, Tuple, true, Total, N...>{
		static ret call(F f, Tuple && t){
			return f(std::get<N>(std::forward<Tuple>(t))...);
		}
	};
}
//allows to call multi-parameter function as a single-parameter one
template<typename restype,int parn, typename... Args>
struct SingleParam{
private:
	typedef restype (*F)(Args...);
	typedef std::tuple<Args...> Tuple;
	typedef typename std::decay<Tuple>::type ttype;
	typedef typename std::tuple_element<parn,Tuple>::type partype;
	F m_func;
	Tuple m_args;
public:
	SingleParam(){}
	SingleParam(SingleParam &S){m_args=S.m_args;m_func=S.m_func;}
	SingleParam(F func,Args... args):m_args(args...){m_func=func;}
	restype operator()(partype x){
		std::get<parn>(m_args)=x;
		return detail::call_impl<
				restype,F, Tuple,
				0 == std::tuple_size<ttype>::value,
				std::tuple_size<ttype>::value
		>::call(m_func, std::forward<Tuple>(m_args));
	}
};

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

// Sympson integral
template<class numt,class functype>
numt Sympson(functype y,numt a, numt b, numt step){
	numt res=0;
	if(step<=0)return 0;
	numt A=a;numt B=b;
	if(A>B){numt C=A;A=B;B=C;}
	numt lastfunc=y(A);
	numt halfstep=step/2;
	for(numt x=A+step;x<B;x+=step){
		numt midfunc=y(x-halfstep);
		numt nextfunc=y(x);
		res+=(lastfunc+4*midfunc+nextfunc)*step/6;
		lastfunc=nextfunc;
	}
	return res;
}
// Sympson integral (output table)
template<class numt,class indexer,class functype>
numt* SympsonTable(int N,indexer x,functype y){
	numt res=0;
	if(N<=0)return nullptr;
	numt *table=new numt[N];
	table[0]=0;
	numt lastx=x[0];numt lastfunc=y(lastx);
	for(int i=1;i<N;i++){
		numt thisarg=x[i];
		numt midfunc=y((thisarg+lastx)/2);
		numt nextfunc=y(thisarg);
		res+=(lastfunc+4*midfunc+nextfunc)*(thisarg-lastx)/6;
		lastfunc=nextfunc;
		lastx=thisarg;
		table[i]=res;
	}
	return table;
}

// Convolution integral
// func1 and func2 cannot be lambda-expressions
// Use SingleParam template class instead
template<class numt,class func1, class func2>
class Convolution{
private:
	func1 A;func2 B;numt Ksi1;numt Ksi2;numt Step;
	class ConvUInt{
	private:	numt X;Convolution *master;
	public:
		ConvUInt(numt x, Convolution* father){X=x;master=father;}
		numt operator()(numt ksi){
			func1 A=master->A;
			func2 B=master->B;
			return A(ksi) * B(X-ksi);
		}
	};
public:
	Convolution(){}
	Convolution(Convolution &C){A=C.A;B=C.B;Ksi1=C.Ksi1;Ksi2=C.Ksi2;Step=C.Step;}
	Convolution(func1 a, func2 b){A=a;B=b;}
	void Init(numt ksi1, numt ksi2, numt step){Ksi1=ksi1;Ksi2=ksi2;Step=step;}
	virtual numt operator()(numt x){
		return Sympson(ConvUInt(x,this),Ksi1,Ksi2,Step);
	}
};
}
#endif // MATH_TEMPLATES_H
