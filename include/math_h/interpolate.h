// this file is distributed under 
// MIT license
#ifndef PPRNSANGJGXVGERD
#define PPRNSANGJGXVGERD
#include <vector>
#include <utility>
#include <functional>
#include "functions.h"
#include "tabledata.h"
namespace MathTemplates{
	template<class numX, class numY=numX>
	class LinearInterpolation:public SortedPoints<numX,numY>,public IFunction<numY,const numX&>{
	public:
	    typedef typename SortedPoints<numX,numY>::Func Func;
	    LinearInterpolation(){}
	    LinearInterpolation(const std::initializer_list<point<numX,numY>>&points)
	    :SortedPoints<numX,numY>(points){}
	    LinearInterpolation(const SortedChain<point<numX,numY>>&points)
	    :SortedPoints<numX,numY>(points){}
	    LinearInterpolation(const SortedChain<point<numX,numY>>&&points)
	    :SortedPoints<numX,numY>(points){}
	    LinearInterpolation(const Func f,const SortedChain<numX>&chain)
	    :SortedPoints<numX,numY>(f,chain){}
	    LinearInterpolation(const Func f,const SortedChain<numX>&&chain)
	    :SortedPoints<numX,numY>(f,chain){}
	    LinearInterpolation(const Func f,const std::initializer_list<numX>&chain)
	    :SortedPoints<numX,numY>(f,chain){}
	    LinearInterpolation(const SortedChain<point_editable_y<numX,numY>>&source)
	    :SortedPoints<numX,numY>(source){}
	    LinearInterpolation(const SortedChain<point_editable_y<numX,numY>>&&source)
	    :SortedPoints<numX,numY>(source){}
	    LinearInterpolation(const SortedPoints<numX,numY>&source)
	    :SortedPoints<numX,numY>(source){}
	    LinearInterpolation(const SortedPoints<numX,numY>&&source)
	    :SortedPoints<numX,numY>(source){}
	    virtual ~LinearInterpolation(){}
	    virtual numY operator()(const numX& x)const override{
		auto tbl=[this](size_t i){return this->operator[](i);};
		size_t sz=this->size();
		if(x==tbl(0).X())
		    return tbl(0).Y();
		if(x==tbl(sz-1).X())
		    return tbl(sz-1).Y();
		auto p=point<numX,numY>(x,numY(0));
		size_t i=details::WhereToInsert<point<numX,numY>>(0,sz-1,*this,p);
		if((i==0)||(i>=sz))
		    throw Exception<LinearInterpolation>("Attempt to interpolate outside the given region.");
		if((x-tbl(i-1).X())<(tbl(i).X()-x)){
		    numX k=(x-tbl(i-1).X())/(tbl(i).X()-tbl(i-1).X());
		    return tbl(i-1).Y()+(tbl(i).Y()-tbl(i-1).Y())*numY(k);
		}else{
		    numX k=(tbl(i).X()-x)/(tbl(i).X()-tbl(i-1).X());
		    return tbl(i).Y()+(tbl(i-1).Y()-tbl(i).Y())*numY(k);
		}
	    }
	};
	template<class numX, class numY=numX>
	class SplineInterpolation:public SortedPoints<numX,numY>,public IFunction<numY,const numX&>{
	private:
	    std::vector<std::vector<numY>> coefficients;
	public:
	    typedef typename SortedPoints<numX,numY>::Func Func;
	    SplineInterpolation(){}
	    SplineInterpolation(const std::initializer_list<point<numX,numY>>&points)
	    :SortedPoints<numX,numY>(points){}
	    SplineInterpolation(const SortedChain<point<numX,numY>>&points)
	    :SortedPoints<numX,numY>(points){}
	    SplineInterpolation(const SortedChain<point<numX,numY>>&&points)
	    :SortedPoints<numX,numY>(points){}
	    SplineInterpolation(const Func f,const SortedChain<numX>&chain)
	    :SortedPoints<numX,numY>(f,chain){}
	    SplineInterpolation(const Func f,const SortedChain<numX>&&chain)
	    :SortedPoints<numX,numY>(f,chain){}
	    SplineInterpolation(const Func f,const std::initializer_list<numX>&chain)
	    :SortedPoints<numX,numY>(f,chain){}
	    SplineInterpolation(const SortedChain<point_editable_y<numX,numY>>&source)
	    :SortedPoints<numX,numY>(source){}
	    SplineInterpolation(const SortedChain<point_editable_y<numX,numY>>&&source)
	    :SortedPoints<numX,numY>(source){}
	    SplineInterpolation(const SortedPoints<numX,numY>&source)
	    :SortedPoints<numX,numY>(source){}
	    SplineInterpolation(const SortedPoints<numX,numY>&&source)
	    :SortedPoints<numX,numY>(source){}
	    virtual ~SplineInterpolation(){}
	    virtual numY operator()(const numX& x)const override{
		auto tbl=[this](size_t i){return this->operator[](i);};
		size_t sz=this->size();
		if(x==tbl(0).X())
		    return tbl(0).Y();
		if(x==tbl(sz-1).X())
		    return tbl(sz-1).Y();
		auto p=point<numX,numY>(x,numY(0));
		size_t i=details::WhereToInsert<point<numX,numY>>(0,sz-1,*this,p);
		if((i==0)||(i>=sz))
		    throw Exception<SplineInterpolation>("Attempt to interpolate outside the given region.");
		if(coefficients.size()!=(sz-1)){
		    coefficients.clean();
		}
		const auto&pol=coefficients[i-1];
		numY pw=1,res=0;
		for(size_t i=0,k=pol.size();i<k;i++){
		    res+=pw*pol[i];
		    pw*=x;
		}
		return res;
	    }
	};

	template<class numtX,class numtY=numtX,class numtZ=numtY>
	class BiLinearInterpolation:public BiSortedPoints<numtX,numtY,numtZ>,public IFunction<numtZ,const numtX&,const numtY&>{
	public:
		BiLinearInterpolation(const std::initializer_list<numtX>&X,const std::initializer_list<numtY>&Y)
		:BiSortedPoints<numtX,numtY,numtZ>(X,Y){}
		BiLinearInterpolation(const std::initializer_list<numtX>&&X,const std::initializer_list<numtY>&&Y)
		:BiSortedPoints<numtX,numtY,numtZ>(X,Y){}
		BiLinearInterpolation(const SortedChain<numtX>&X,const SortedChain<numtY>&Y)
		:BiSortedPoints<numtX,numtY,numtZ>(X,Y){}
		BiLinearInterpolation(const SortedChain<numtX>&&X,const SortedChain<numtY>&&Y)
		:BiSortedPoints<numtX,numtY,numtZ>(X,Y){}
		BiLinearInterpolation()
		:BiLinearInterpolation({},{}){}
		BiLinearInterpolation(const BiSortedPoints<numtX,numtY,numtZ>&source)
		:BiSortedPoints<numtX,numtY,numtZ>(source){}
		virtual ~BiLinearInterpolation(){}
		virtual numtZ operator()(const numtX& x,const numtY& y)const override{
			int szx=this->X().size();
			int i=details::WhereToInsert<numtX>(0,szx-1,this->X(),x);
			if((i<=0)||(i>=szx))
				throw Exception<BiLinearInterpolation>("Attempt to interpolate outside the given region (x).");
			int szy=this->Y().size();
			int j=details::WhereToInsert<numtY>(0,szy-1,this->Y(),y);
			if((j<=0)||(j>=szy))
				throw Exception<BiLinearInterpolation>("Attempt to interpolate outside the given region (y).");
			numtZ U=numtZ(this->X()[i]-this->X()[i-1]);
			numtZ T=numtZ(this->Y()[j]-this->Y()[j-1]);
			numtZ u=numtZ(x-this->X()[i-1])/U;
			numtZ t=numtZ(y-this->Y()[j-1])/T;
			numtZ _u=numtZ(this->X()[i]-x)/U;
			numtZ _t=numtZ(this->Y()[j]-y)/T;
			return this->operator[](i-1)[j-1]*_u*_t
				+this->operator[](i)[j-1]*u*_t
				+this->operator[](i-1)[j]*_u*t
				+this->operator[](i)[j]*u*t;
		}
	};
};
#endif
