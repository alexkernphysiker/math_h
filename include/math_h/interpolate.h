// this file is distributed under 
// MIT license
#ifndef PPRNSANGJGXVGERD
#define PPRNSANGJGXVGERD
#include <vector>
#include <utility>
#include <functional>
#include "functions.h"
#include "structures.h"
namespace MathTemplates{
	using namespace std;
	template<class numX, class numY=numX>
	class LinearInterpolation:public SortedPoints<numX,numY>,public IFunction<numY,const numX&>{
	public:
		LinearInterpolation(){}
		LinearInterpolation(const SortedChain<point_editable_y<numX,numY>>&source):SortedPoints<numX,numY>(source){}
		LinearInterpolation(const SortedChain<point_editable_y<numX,numY>>&&source):SortedPoints<numX,numY>(source){}
		LinearInterpolation(const SortedPoints<numX,numY>&source):SortedPoints<numX,numY>(source){}
		LinearInterpolation(const SortedPoints<numX,numY>&&source):SortedPoints<numX,numY>(source){}
		virtual ~LinearInterpolation(){}
		virtual numY operator()(const numX& x)const override{
			auto tbl=[this](size_t i){return this->operator[](i);};
			size_t sz=this->size();
			if(x==tbl(0).X())
				return tbl(0).Y();
			if(x==tbl(sz-1).X())
				return tbl(sz-1).Y();
			auto p=point<numX,numY>(x,numY(0));
			int i=details::WhereToInsert<point<numX,numY>>(0,sz-1,*this,p);
			if((i<=0)||(i>=sz))
				throw Exception<LinearInterpolation>("Attempt to interpolate outside the given region.");
			numX k=(x-tbl(i-1).X())/(tbl(i).X()-tbl(i-1).X());
			return tbl(i-1).Y()+(tbl(i).Y()-tbl(i-1).Y())*numY(k);
		}
	};
	template<class numtX,class numtY=numtX,class numtZ=numtY>
	class BiLinearInterpolation:public BiSortedPoints<numtX,numtY,numtZ>,public IFunction<numtZ,const numtX&,const numtY&>{
	public:
		BiLinearInterpolation(const initializer_list<numtX>&X,const initializer_list<numtY>&Y):BiSortedPoints<numtX,numtY,numtZ>(X,Y){}
		BiLinearInterpolation(const initializer_list<numtX>&&X,const initializer_list<numtY>&&Y):BiSortedPoints<numtX,numtY,numtZ>(X,Y){}
		BiLinearInterpolation(const vector<numtX>&X,const vector<numtY>&Y):BiSortedPoints<numtX,numtY,numtZ>(X,Y){}
		BiLinearInterpolation(const vector<numtX>&&X,const vector<numtY>&&Y):BiSortedPoints<numtX,numtY,numtZ>(X,Y){}
		BiLinearInterpolation():BiLinearInterpolation({},{}){}
		BiLinearInterpolation(const BiSortedPoints<numtX,numtY,numtZ>&source):BiSortedPoints<numtX,numtY,numtZ>(source){}
		virtual ~BiLinearInterpolation(){}
		virtual numtZ operator()(const numtX& x,const numtY& y)const override{
			size_t szx=this->X().size();
			int i=details::WhereToInsert<numtX>(0,szx-1,this->X(),x);
			if((i<=0)||(i>=szx))
				throw Exception<BiLinearInterpolation>("Attempt to interpolate outside the given region (x).");
			size_t szy=this->Y().size();
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
