// this file is distributed under 
// MIT license
#ifndef PPRNSANGJGXVGERD
#define PPRNSANGJGXVGERD
#include <vector>
#include <utility>
#include <functional>
#include "functions.h"
#include "structures.h"
#include "error.h"
namespace MathTemplates{
	using namespace std;
	template<class numX, class numY=numX, class PointIndexer=std::vector<point<numX,numY>>>
	numY  Interpolate_Linear(const int from,const int to,const PointIndexer tbl,const numX x){
		using namespace std;
		if(x==tbl[from].X())
			return tbl[from].Y();
		if(x==tbl[to].X())
			return tbl[to].Y();
		auto p=point<numX,numY>(x,numY(0));
		int i=details::WhereToInsert<point<numX,numY>>(from,to,tbl,p);
		if((i<=from)||(i>to))
			throw Exception<point<numX,numY>>("Attempt to interpolate outside the given region.");
		numX k=(x-tbl[i-1].X())/(tbl[i].X()-tbl[i-1].X());
		return tbl[i-1].Y()+(tbl[i].Y()-tbl[i-1].Y())*numY(k);
	}
	template<class numX, class numY=numX, class PointIndexer=std::vector<point<numX,numY>>,class Size=std::function<int()>>
	numY InterpolateLinear(const numX x,const PointIndexer tbl,const Size size){
		return Interpolate_Linear(0,size()-1,tbl,x);
	}
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
			return InterpolateLinear<numX,numY>(x,*this,field_size(*this));
		}
	};
};
#endif
