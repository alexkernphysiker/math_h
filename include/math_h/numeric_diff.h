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
    template<size_t size,class numt>
    inline Opr<Vector<size,numt>,numt> num_Pder1(const Vector<size,numt>&x,const Vector<size,numt>&delta){
	return Opr<Vector<size,numt>,numt>([x,delta](const IFunction<numt,const Vector<size,numt>&>&F){
	    const auto x1=x-delta,x2=x+delta;
	    const numt y1=F(x1),y2=F(x2);
	    return (y2-y1)/(x2-x1).M();
	});
    }
    namespace numeric_diff_details{
	template<size_t size,class numt,size_t index>
	Vector<index,Opr<Vector<size,numt>,numt>> ___nabla_partial(const Vector<size,numt>&X,const numt&delta){
	    if constexpr(index==1) return vec(num_Pder1(X,Vector<size,numt>::template basis_vector<1>()*delta));
	    else return ___nabla_partial<size,numt,index-1>(X,delta).template InsertComponent<index>(
		num_Pder1(X,Vector<size,numt>::template basis_vector<index>()*delta)
	    );
	}
    }
    template<size_t size,class numt>
    inline Vector<size,Opr<Vector<size,numt>,numt>> nabla(const Vector<size,numt>&X,const numt&delta){
	static_assert(size>0,"invalid space dimmensions count");
	return numeric_diff_details::___nabla_partial<size,numt,size>(X,delta);
    }
};
#endif
