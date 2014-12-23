// https://github.com/alexkernphysiker/math_h
#ifndef ___SYMPSON_H
#	define ___SYMPSON_H

// Sympson integral
template<class numX,class functype,class numY=numX>
numY Sympson(functype y,numX a, numX b, numX step){
	numX stp=step;
	if((stp*(b-a))<=0) stp=-stp;
	numX halfstep=stp/2;
	numY lastfunc=y(a);
	numY res=0;
	for(numX x=a+stp;(a-b)*(x-b)>0;x+=stp){
		numY midfunc=y(x-halfstep);
		numY nextfunc=y(x);
		res+=(lastfunc+4*midfunc+nextfunc)*stp/6;
		lastfunc=nextfunc;
	}
	return res;
}
// Sympson integral (output table)
template<class numX,class indexerX,class functype,class numY=numX>
numY* SympsonTable(int N,indexerX x,functype y){
	numY res=0;
	if(N<=0)return nullptr;
	numY *table=new numY[N];
	table[0]=0;
	numX lastx=x[0];numY lastfunc=y(lastx);
	for(int i=1;i<N;i++){
		numX thisarg=x[i];
		numY midfunc=y((thisarg+lastx)/2);
		numY nextfunc=y(thisarg);
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
			return master->A(ksi) * master->B(X-ksi);
		}
	};
public:
	Convolution(){}
	Convolution(Convolution &C){A=C.A;B=C.B;Ksi1=C.Ksi1;Ksi2=C.Ksi2;Step=C.Step;}
	Convolution(func1 a, func2 b){A=a;B=b;}
	void Init(numt ksi1, numt ksi2, numt step){Ksi1=ksi1;Ksi2=ksi2;Step=step;}
	virtual numt operator()(numt x){
		return Sympson<numt,ConvUInt>(ConvUInt(x,this),Ksi1,Ksi2,Step);
	}
};
#endif // ___SYMPSON_H
