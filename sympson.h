#ifndef ___SYMPSON_H
#	define ___SYMPSON_H

// Sympson integral
template<class numt,class functype>
numt Sympson(functype y,numt a, numt b, numt step){
	numt stp=step;
	if((stp*(b-a))<=0) stp=-stp;
	numt halfstep=stp/2;
	numt lastfunc=y(a);
	numt res=0;
	for(numt x=a+stp;(a-b)*(x-b)>0;x+=stp){
		numt midfunc=y(x-halfstep);
		numt nextfunc=y(x);
		res+=(lastfunc+4*midfunc+nextfunc)*stp/6;
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
			return master->A(ksi) * master->B(X-ksi);
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
#endif // ___SYMPSON_H
