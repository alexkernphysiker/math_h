// this file is distributed under
// LGPLv3 license
#ifndef GRIVHOWXKUEHYGQF
#	define GRIVHOWXKUEHYGQF
#include "error.h"
#include "functions.h"
#include "tabledata.h"
namespace MathTemplates
{
template<class numX, class numY = numX, class functype = std::function<numY(numX)>>
numY Sympson(const functype y, const numX &a, const numX &b, const numX &step)
{
    numX stp = step;
    if ((stp * (b - a)) <= 0) stp = -stp;
    numX halfstep = stp / 2;
    numY lastfunc = y(a);
    numY res = 0;
    for (numX x = a + stp; (a - b) * (x - b) > 0; x += stp) {
        numY midfunc = y(x - halfstep);
        numY nextfunc = y(x);
        res += (lastfunc + 4 * midfunc + nextfunc) * stp / 6;
        lastfunc = nextfunc;
    }
    return res;
}
template<class numX, class numY = numX>
SortedPoints<numX, numY> Int_Trapez_Table(const SortedPoints<numX, numY> &source)
{
    SortedPoints<numX, numY> res;
    res << point<numX, numY>(source.left().X(), 0);
    for (size_t i = 1; i < source.size(); i += 1) {
        res << point<numX, numY>(
                source[i].X(),
                res.right().Y() + (source[i].Y() + source[i - 1].Y())*numY(source[i].X() - source[i - 1].X()) / numY(2)
            );
    }
    return res;
}
template<class numX, class numY = numX> //Accepts only positive function
SortedPoints<numX, numY>  Int_Trapez_Table_PositiveStrict(const SortedPoints<numX, numY> &source)
{
    SortedPoints<numX, numY> res;
    res << point<numX, numY>(source.left().X(), 0);
    for (size_t i = 1; i < source.size(); i += 1) {
        if (source[i - 0].Y() < 0)
            throw Exception<SortedPoints<numX, numY>>("SympsonTablePositiveStrict: negative value detected");
        res << point<numX, numY>(
                source[i].X(),
                res.right().Y() + (source[i].Y() + source[i - 1].Y())*numY(source[i].X() - source[i - 1].X()) / numY(2)
            );
    }
    return std::move(res);
}

namespace AdaptiveQuadrature_details{
template<size_t mindepth,class numX, class numY = numX, class functype = std::function<numY(numX)>>
    numY AdaptiveQuadrature(const functype&y, const numX &a, const numX &b, const numY &epsilon,const numY&A,const numY&B,const numY&S, size_t depth)
    {
	    const numX m = (a+b)/numX(2);
	    const numY M=y(m);
	    numY Sl=(A+M)*numY(m-a)/numY(2);
	    numY Sr=(M+B)*numY(b-m)/numY(2);
	    if((depth<=mindepth)||(abs(Sl+Sr-S)>epsilon)){
	        Sl=AdaptiveQuadrature<mindepth>(y,a,m,epsilon,A,M,Sl,depth+1);
    	    Sr=AdaptiveQuadrature<mindepth>(y,m,b,epsilon,M,B,Sr,depth+1);
	    }
	    return Sl+Sr;
    }
};
template<size_t mindepth=5,class numX, class numY = numX, class functype = std::function<numY(numX)>>
inline numY AdaptiveQuadrature(const functype y, const numX &a, const numX &b, const numY &epsilon)
{
    const numY A=y(a),B=y(b),S=(A+B)*numY(b-a)/numY(2);
    return AdaptiveQuadrature_details::AdaptiveQuadrature<mindepth>(y,a,b,epsilon,A,B,S,0);
}

template<size_t mindepth=5,class numX = double, class numY = numX,
         class func1 = std::function<numY(const numX &)>, class func2 = std::function<numY(const numX &)>
>
class Convolution: public IFunction<numY, const numX &>
{
private:
    func1 A;
    func2 B;
    numX Ksi1;
    numX Ksi2;
    numY Eps;
public:
    Convolution(const Convolution &C)
        : A(C.A), B(C.B), Ksi1(C.Ksi1), Ksi2(C.Ksi2), Eps(C.Eps) {}
    Convolution(const func1 a, const func2 b, const numX &ksi1, const numX &ksi2, const numY &eps)
        : A(a), B(b)
    {
        Ksi1 = ksi1;
        Ksi2 = ksi2;
        Eps = eps;
    }
    virtual numY operator()(const numX &x)const override
    {
        return AdaptiveQuadrature<mindepth,numX, numY>([this,&x](const numX & ksi) {
            return A(ksi) * B(x - ksi);
        }, Ksi1, Ksi2, Eps);
    }
};
template<size_t mindepth=5,class numX = double, class numY = numX,
         class func1 = std::function<numY(const numX &)>, class func2 = std::function<numY(const numX &)>
         >
Convolution<mindepth,numX, numY, func1, func2> make_convolution(
    const func1 a, const func2 b, const numX &ksi1, const numX &ksi2, const numY &eps)
{
    return Convolution<mindepth,numX, numY, func1, func2>(a, b, ksi1, ksi2, eps);
}
};
#endif
