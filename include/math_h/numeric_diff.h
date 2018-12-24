// this file is distributed under
// LGPLv3 license
#ifndef _______NUMERIC_DIFF___H_________
#	define _______NUMERIC_DIFF___H_________
#include "error.h"
#include "functions.h"
#include "vectors.h"
namespace MathTemplates
{
    template<class numt=double>
    inline Opr<numt> num_der1(const numt&x,const numt&delta){
	return Opr<numt>([x,delta](const IFunction<numt,const numt&>&F){
	    const numt x1=x-delta,x2=x+delta;
	    const numt y1=F(x1),y2=F(x2);
	    return (y2-y1)/(x2-x1);
	});
    }
    template<class numt=double>
    inline Opr<numt> num_der2(const numt&x,const numt&delta){
	return Opr<numt>([x,delta](const IFunction<numt,const numt&>&F){
	    const numt x1=x-delta,x2=x,x3=x+delta;
	    const numt y1=F(x1),y2=F(x2),y3=F(x3);
	    return (y1-y2*numt(2)+y3)/(delta*delta);
	});
    }
};
#endif
