#ifndef GRIVHOWXKUEHYGQF
#	define GRIVHOWXKUEHYGQF
#include <functional>
template<class numX,class numY=numX,class functype=std::function<numY(numX)>>
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
template<class numX,class indexerX,class numY=numX,class functype=std::function<numY(numX)>>
numY* SympsonTable(functype y,indexerX X,int N){
	numY res=0;
	if(N<=0)return nullptr;
	numY *table=new numY[N];
	table[0]=0;
	numX lastx=X[0];
	numY lastfunc=y(lastx);
	for(int i=1;i<N;i++){
		numX thisarg=X[i];
		numY midfunc=y((thisarg+lastx)/2);
		numY nextfunc=y(thisarg);
		res+=(lastfunc+4*midfunc+nextfunc)*(thisarg-lastx)/6;
		lastfunc=nextfunc;
		lastx=thisarg;
		table[i]=res;
	}
	return table;
}
template<class numX,class numY=numX,class func1=std::function<numY(numX)>,class func2=std::function<numY(numX)>>
class Convolution{
private:
	func1 A;
	func2 B;
	numX Ksi1;
	numX Ksi2;
	numX Step;
public:
	Convolution(){}
	Convolution(Convolution &C):A(C.A),B(C.B),Ksi1(C.Ksi1),Ksi2(C.Ksi2),Step(C.Step){}
	Convolution(func1 a, func2 b,numX ksi1, numX ksi2, numX step):A(a),B(b){Ksi1=ksi1;Ksi2=ksi2;Step=step;}
	numY operator()(numX x){
		return Sympson<numX,numY>([this,x](numX ksi){return A(ksi)*B(x-ksi);},Ksi1,Ksi2,Step);
	}
};
#endif
