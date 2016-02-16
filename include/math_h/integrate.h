// this file is distributed under 
// MIT license
#ifndef GRIVHOWXKUEHYGQF
#	define GRIVHOWXKUEHYGQF
#include <functional>
#include "interpolate.h"
#include "error.h"
namespace MathTemplates{
	using namespace std;
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
	template<class numX,class numY=numX>
	LinearInterpolation<numX,numY> Int_Trapez_Table(const LinearInterpolation<numX,numY>&source){
		LinearInterpolation<numX,numY> res;
		res<<make_pair(source.left().first,0);
		for(size_t i=1;i<source.size();i+=1)
			res<<make_pair(source[i].first,res.right().second+(source[i].second+source[i-1].second)*(source[i].first-source[i-1].first)/numY(2));
		return res;
	}
	template<class numX,class numY=numX>//Accepts only positive function
	LinearInterpolation<numX,numY>  Int_Trapez_Table_PositiveStrict(const LinearInterpolation<numX,numY>&source){
		LinearInterpolation<numX,numY> res;
		res<<make_pair(source.left().first,0);
		for(size_t i=1;i<source.size();i+=1){
			if(source[i-0].second<0)
				throw Exception<LinearInterpolation<numX,numY>>("SympsonTablePositiveStrict: negative value detected");
			res<<make_pair(source[i].first,res.right().second+(source[i].second+source[i-1].second)*(source[i].first-source[i-1].first)/numY(2));
		}
		return res;
	}
	template<class numX,class numY=numX>
	LinearInterpolation<numX,numY> Int_Trapez_Table_RV(LinearInterpolation<numX,numY>&&source){
		return Int_Trapez_Table(source);
	}
	template<class numX,class numY=numX>
	LinearInterpolation<numX,numY> Int_Trapez_Table_PositiveStrict_RV(LinearInterpolation<numX,numY>&&source){
		return Int_Trapez_Table_PositiveStrict(source);
	}
	template<class numX,class numY=numX,class func1=function<numY(numX)>,class func2=function<numY(numX)>>
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
		numY operator()(numX x)const{
			return Sympson<numX,numY>([this,x](numX ksi){return A(ksi)*B(x-ksi);},Ksi1,Ksi2,Step);
		}
	};
};
#endif
