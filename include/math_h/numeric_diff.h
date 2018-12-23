// this file is distributed under
// LGPLv3 license
#ifndef _______NUMERIC_DIFF___H_________
#	define _______NUMERIC_DIFF___H_________
#include "error.h"
#include "functions.h"
#include "vectors.h"
namespace MathTemplates
{
    template<class numtx=double,class numty=numtx>
    class der1{
    private:
	numtx m_delta;
    public:
	der1(const numtx&delta):m_delta(delta){}
	template<class FUNC>
	numty operator()(FUNC F,const numtx&x){
	    const auto x1=x-m_delta,x2=x+m_delta;
	    const auto y1=F(x1),y2=F(x2);
	    return (y2-y1)/numty(x2-x1);
	}
    };
    template<class numtx=double,class numty=numtx>
    class der2{
    private:
	numtx m_delta;
    public:
	der2(const numtx&delta):m_delta(delta){}
	template<class FUNC>
	numty operator()(FUNC F,const numtx&x){
	    const auto x1=x-m_delta,x2=x,x3=x+m_delta;
	    const auto y1=F(x1),y2=F(x2),y3=F(x3);
	    return (y1-y2*numty(2)+y3)/numty(m_delta*m_delta);
	}
    };
};
#endif
