// this file is distributed under
// MIT license
#ifndef ______POINTS_H_______
#	define ______POINTS_H_______
#include "chains.h"
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
    const numtX &X()const
    {
        return x;
    }
    const numtY &Y()const
    {
        return y;
    }
    point(const numtX &pos): x(pos), y(numtY(0)) {}
    point(const numtX &pos, const numtY &val): x(pos), y(val) {}
    template<class numtX2 = numtX, class numtY2 = numtY>
    point(const point<numtX2,numtY2> &source): x(source.X()), y(source.Y()) {}
    template<class numt = numtY>
    point(const std::initializer_list<numt> &source)
    {
        if (source.size() == 0)
            throw Exception<point>("wrong initialization of point from emply list");
        if (source.size() > 2)
            throw Exception<point>("wrong initialization of value from list with more than two numbers");
        Chain<numt> v;
        for (const numt &x : source)v.push_back(x);
        x = numtX(v[0]);
        y = numtY(v[1]);
    }
    template<class numtX2 = numtX, class numtY2 = numtY>
    point &operator=(const point<numtX2,numtY2> &source)
    {
        x = source.X();
        y = source.Y();
	return *this;
    }
    point &operator=(const numtY &s)
    {
        y = s;
	return *this;
    }
    bool operator<(const point &b)const
    {
        return x < b.x;
    }
    bool operator>(const point &b)const
    {
        return x > b.x;
    }
    const point operator+(const numtY &val)const
    {
        return {X(), Y() + val};
    }
    const point operator-(const numtY &val)const
    {
        return {X(), Y() - val};
    }
    const point operator*(const numtY &val)const
    {
        return {X(), Y() * val};
    }
    const point operator/(const numtY &val)const
    {
        return {X(), Y() / val};
    }
    const point operator+(const point &other)const
    {
        if (other.X() == X())
            return operator+(other.Y());
        else
            throw Exception<point>("Cannot perform arithmetic operation with two points that have different X-coordinate");
    }
    const point operator-(const point &other)const
    {
        if (other.X() == X())
            return operator-(other.Y());
        else
            throw Exception<point>("Cannot perform arithmetic operation with two points that have different X-coordinate");
    }
    const point operator*(const point &other)const
    {
        if (other.X() == X())
            return operator*(other.Y());
        else
            throw Exception<point>("Cannot perform arithmetic operation with two points that have different X-coordinate");
    }
    const point operator/(const point &other)const
    {
        if (other.X() == X())
            return operator/(other.Y());
        else
            throw Exception<point>("Cannot perform arithmetic operation with two points that have different X-coordinate");
    }
};
template<class numtX = double, class numtY = numtX>
inline const point<numtX,numtY>make_point(const numtX&x,const numtY&y){return point<numtX,numtY>(x,y);}
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
    template<class numtX2 = numtX, class numtY2 = numtY,class numtZ2=numtZ>
    point3d(const point3d<numtX2,numtY2,numtZ2>&source)
	: x(source.X()), y(source.Y()), z(source.Z()){}
};
template<class numtX = double, class numtY = numtX,class numtZ=numtY>
inline const point3d<numtX,numtY,numtZ>make_point(const numtX&x,const numtY&y,const numtZ&z){
    return point3d<numtX,numtY,numtZ>(x,y,z);
}

};
#endif
