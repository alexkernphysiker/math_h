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
	class LinearInterpolation:public SortedPoints<numX,numY>,public IFunction<numX,numY>{
	public:
		typedef pair<numX,numY> Point;
	private:
		vector<Point> data;
	public:
		LinearInterpolation(){}
		LinearInterpolation(const SortedPoints<numX,numY>&source):SortedPoints<numX,numY>(source){}
		LinearInterpolation(const SortedPoints<numX,numY>&&source):SortedPoints<numX,numY>(source){}
		virtual ~LinearInterpolation(){}
		virtual numY operator()(const numX x)const override{
			auto tbl=[this](size_t i){return SortedPoints<numX,numY>::operator[](i);};
			size_t sz=SortedPoints<numX,numY>::size();
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
};
#endif
