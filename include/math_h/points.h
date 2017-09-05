// this file is distributed under
// MIT license
#ifndef ______POINTS_H_______
#	define ______POINTS_H_______
#include <vector>
#include <functional>
#include "error.h"
namespace MathTemplates
{
template<class numtX = double, class numtY = numtX>
class point
{
private:
    numtX x;
    numtY y;
public:
    point(const numtX &pos): x(pos), y(numtY(0)) {}
    point(const numtX &pos, const numtY &val): x(pos), y(val) {}
    point(const point &source): x(source.x), y(source.y) {}
    template<class numt = numtY>
    point(const std::initializer_list<numt> &source)
    {
        if (source.size() == 0)
            throw Exception<point>("wrong initialization of point from emply list");
        if (source.size() > 2)
            throw Exception<point>("wrong initialization of value from list with more than two numbers");
        std::vector<numt> v;
        for (const numt &x : source)v.push_back(x);
        x = numtX(v[0]);
        y = numtY(v[1]);
    }
    const numtX &X()const
    {
        return x;
    }
    const numtY &Y()const
    {
        return y;
    }
    bool operator<(const point &b)const
    {
        return x < b.x;
    }
    bool operator>(const point &b)const
    {
        return x > b.x;
    }
    const point operator+(const numtY&val)const
    {
	return point(X(),Y()+val);
    }
    const point operator-(const numtY&val)const
    {
	return point(X(),Y()-val);
    }
    const point operator*(const numtY&val)const
    {
	return point(X(),Y()*val);
    }
    const point operator/(const numtY&val)const
    {
	return point(X(),Y()/val);
    }
    const point operator+(const point&other)const
    {
	if(other.X()==X())
	    return operator+(other.Y());
	else
	    throw Exception<point>("Cannot perform arithmetic operation with two points that have different X-coordinate");
    }
    const point operator-(const point&other)const
    {
	if(other.X()==X())
	    return operator-(other.Y());
	else
	    throw Exception<point>("Cannot perform arithmetic operation with two points that have different X-coordinate");
    }
    const point operator*(const point&other)const
    {
	if(other.X()==X())
	    return operator*(other.Y());
	else
	    throw Exception<point>("Cannot perform arithmetic operation with two points that have different X-coordinate");
    }
    const point operator/(const point&other)const
    {
	if(other.X()==X())
	    return operator/(other.Y());
	else
	    throw Exception<point>("Cannot perform arithmetic operation with two points that have different X-coordinate");
    }
protected:
    numtX &__X()
    {
        return x;
    }
    numtY &__Y()
    {
        return y;
    }
};
template<class numtX = double, class numtY = numtX>
class point_editable_x: public point<numtX, numtY>
{
public:
    point_editable_x(const point<numtX, numtY> &source): point<numtX, numtY>(source) {}
    numtX &varX()
    {
        return point<numtX, numtY>::__X();
    }
};
template<class numtX = double, class numtY = numtX, class numtZ = numtY>
class point3d
{
private:
    numtX x;
    numtY y;
    numtZ z;
public:
    point3d(const numtX &_x, const numtY &_y, const numtZ &_z)
        : x(_x), y(_y), z(_z) {}
    virtual ~point3d() {}
    const numtX &X()const
    {
        return x;
    }
    const numtY &Y()const
    {
        return y;
    }
    const numtZ &Z()const
    {
        return z;
    }
};

};
#endif
