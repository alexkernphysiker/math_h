// this file is distributed under
// LGPLv3 license
#ifndef ____table_data_H_____
#	define ____table_data_H_____
#if __cplusplus<201100L
#error c++>=11 is needed for using math_h headers
#endif
#include <vector>
#include <functional>
#include <map>
#include "error.h"
#include "sigma.h"
#include "functions.h"
namespace MathTemplates
{
    namespace table_data_details
    {
        template<class comparable, class indexer = Chain<comparable>>
        inline size_t  WhereToInsert(const size_t from, const size_t to, const indexer& X, const comparable& x)
        {
            if (from > to) return from;
            if (x < X[from]) return from;
            if (x > X[to]) return to + 1;
            size_t beg = from, end = to;
            if (x > X[end]) return end + 1;
            if (x < X[beg]) return beg;
            while (1 < (end - beg)) {
                size_t mid = (beg + end) / 2;
                if (x < X[mid]) end = mid;
                else if (x > X[mid]) beg = mid;
                else return mid;
            }
            return end;
        }
        template<class comparable, class indexer, class Size, class Insert>
        inline void InsertSorted(const comparable& x, indexer& X, const Size size, const Insert insert)
        {
            if (size() == 0) insert(0, x);
            else insert(WhereToInsert(0, size() - 1, X, x), x);
        }
    }
#define std_size(vector) [&vector](){return (vector).size();}
#define std_insert(vector,type) [&vector](int pos,type x){(vector).insert((vector).begin()+pos,x);}
#define field_size(vector)  [this](){return (vector).size();}
#define field_insert(vector,type)  [this](int pos,type x){(vector).insert((vector).begin()+pos,x);}
    template<class comparable = double>
    class SortedChain
    {
    private:
        Chain<comparable> data;
    public:
        const Chain<comparable>& operator()()const { return data; }
        SortedChain() {}
        SortedChain(const SortedChain& points) = delete;
        SortedChain& operator=(const SortedChain& points) = delete;
        SortedChain(SortedChain&& points) = default;
        SortedChain& operator=(SortedChain&& points) = default;
        SortedChain clone()const
        {
            SortedChain res;
            res.data = data;
            return res;
        }
        template<class comparable2>
        SortedChain(const SortedChain<comparable2>& points)
        {
            for (const auto& p : points.data)
                data.push_back(comparable(p));
        }
        SortedChain IndexRange(const size_t from, const size_t count)const
        {
            SortedChain result;
            size_t final = from + count;
            for (size_t i = from; (i < final) && (i < data.size()); i++)
                result.data.push_back(data[i]);
            return result;
        }
        SortedChain SelectByCondition(const std::function<bool(const comparable&)> condition)const
        {
            SortedChain result;
            for (size_t i = 0; i < data.size(); i++)
                if (condition(data[i]))
                    result.data.push_back(data[i]);
            return result;
        }
        SortedChain& operator<<(const comparable& p)
        {
            table_data_details::InsertSorted(p, data, field_size(data), field_insert(data, comparable));
            return *this;
        }
        template<class comparable2>
        SortedChain(const Chain<comparable2>& points)
        {
            for (const auto& p : points)
                operator<<(comparable(p));
        }
        virtual ~SortedChain() {}
        void clear()
        {
            data.clear();
        }

        size_t size()const
        {
            return data.size();
        }
        const comparable& operator[](const size_t i)const
        {
            if (size() <= i)
                throw Exception<SortedChain>("Range check error");
            return data[i];
        }
        const comparable& left()const
        {
            if (size() < 1)
                throw Exception<SortedChain>("Attempt to obtain empty properties.");
            return data[0];
        }
        const comparable& right()const
        {
            if (size() < 1)
                throw Exception<SortedChain>("Attempt to obtain empty properties.");
            return data[size() - 1];
        }
        typedef typename std::vector<comparable>::const_iterator const_iterator;
        const_iterator begin()const
        {
            return data.begin();
        }
        const_iterator cbegin()const
        {
            return data.cbegin();
        }
        const_iterator end() const
        {
            return data.end();
        }
        const_iterator cend() const
        {
            return data.cend();
        }
    protected:
        comparable& accessBin(const size_t i)
        {
            if (data.size() <= i)
                throw Exception<SortedChain>("range check error");
            return data[i];
        }
    };
    template<class numX>
    Chain<numX> ChainWithStep(const numX& from, const numX& step, const numX& to)
    {
        if (from >= to)throw Exception<Chain<numX>>("wrong binning ranges");
        if (step <= 0)throw Exception<Chain<numX>>("wrong binning step");
        Chain<numX> res;
        for (numX x = from; x <= to; x += step)res.push_back(x);
        return res;
    }
    template<class numX>
    Chain<numX> ChainWithCount(const size_t cont, const numX& from, const numX& to)
    {
        if (from >= to)throw Exception<Chain<numX>>("wrong binning ranges");
        if (0 == cont)throw Exception<Chain<numX>>("wrong bins count");
        numX step = (to - from) / numX(cont);
        Chain<numX> res;
        for (numX x = from; x <= to; x += step)res.push_back(x);
        return res;
    }
    template<class numt>
    Chain<value<numt>> BinsByStep(const numt from, const numt step, const numt to)
    {
        if (0 >= step)throw Exception<Chain<value<numt>>>("wrong bin width");
        if (to <= from)throw Exception<Chain<value<numt>>>("wrong range");
        numt delta = step / numt(2);
        Chain<value<numt>> res;
        for (numt x = from + delta; x < to; x += step)
            res.push_back(value<numt>(x, delta));
        return res;
    }
    template<class numt>
    Chain<value<numt>> BinsByCount(const size_t count, const numt from, const numt to)
    {
        if (0 == count)throw Exception<Chain<value<numt>>>("wrong bins count");
        return BinsByStep(from, (to - from) / numt(count), to);
    }


    template<class numtX = double, class numtY = numtX>
    class point
    {
    private:
        numtX x;
        numtY y;
    public:
        const numtX& X()const
        {
            return x;
        }
        const numtY& Y()const
        {
            return y;
        }
        numtY& Y_var()
        {
            return y;
        }
        point(const numtX& pos) : x(pos), y(numtY(0)) {}
        point(const numtX& pos, const numtY& val) : x(pos), y(val) {}
        template<class numtX2, class numtY2>
        point(const point<numtX2, numtY2>& source) : x(source.X()), y(source.Y()) {}
        point& operator=(const point& source) = default;
        inline point& operator=(const numtY& source)
        {
            y = source;
            return *this;
        }
        inline bool operator<(const point& b)const
        {
            return x < b.x;
        }
        inline bool operator>(const point& b)const
        {
            return x > b.x;
        }
        template<class numtY2>
        auto operator+(const numtY2& val)const->point<numtX, decltype(Y() + val)>
        {
            return { X(), Y() + val };
        }
        template<class numtY2>
        auto operator+(const point<numtX, numtY2>& val)const->point<numtX, decltype(Y() + val.Y())>
        {
            if (val.X() == X())
                return operator+(val.Y());
            else
                throw Exception<point>("Cannot perform arithmetic operation with two points that have different X-coordinate");
        }
        template<class numtY2>
        inline auto operator-(const numtY2& val)const->point<numtX, decltype(Y() - val)>
        {
            return { X(), Y() - val };
        }
        template<class numtY2>
        auto operator-(const point<numtX, numtY2>& val)const->point<numtX, decltype(Y() - val.Y())>
        {
            if (val.X() == X())
                return operator-(val.Y());
            else
                throw Exception<point>("Cannot perform arithmetic operation with two points that have different X-coordinate");
        }
        template<class numtY2>
        inline auto operator*(const numtY2& val)const->point<numtX, decltype(Y()* val)>
        {
            return { X(), Y() * val };
        }
        template<class numtY2>
        auto operator*(const point<numtX, numtY2>& val)const->point<numtX, decltype(Y()* val.Y())>
        {
            if (val.X() == X())
                return operator*(val.Y());
            else
                throw Exception<point>("Cannot perform arithmetic operation with two points that have different X-coordinate");
        }
        template<class numtY2>
        inline auto operator/(const numtY2& val)const->point<numtX, decltype(Y() / val)>
        {
            return { X(), Y() / val };
        }
        template<class numtY2>
        auto operator/(const point<numtX, numtY2>& val)const->point<numtX, decltype(Y() / val.Y())>
        {
            if (val.X() == X())
                return operator/(val.Y());
            else
                throw Exception<point>("Cannot perform arithmetic operation with two points that have different X-coordinate");
        }

        template<class numtY2>
        point& operator<<(const point<numtX, numtY2>& other)
        {
            if (other.X() == X())
                y << other.Y();
            else
                throw Exception<point>("Cannot perform arithmetic operation with two points that have different X-coordinate");
            return *this;
        }
    };
    template<class numtX = double, class numtY = numtX>
    inline point<numtX, numtY>make_point(const numtX& x, const numtY& y) { return point<numtX, numtY>(x, y); }
    template<class numtX = double, class numtY = numtX, class numtZ = numtY>
    class point3d
    {
    private:
        numtX x;
        numtY y;
        numtZ z;
    public:
        point3d(const numtX& _x, const numtY& _y, const numtZ& _z)
            : x(_x), y(_y), z(_z) {
        }
        virtual ~point3d() {}
        const numtX& X()const
        {
            return x;
        }
        const numtY& Y()const
        {
            return y;
        }
        const numtZ& Z()const
        {
            return z;
        }
        template<class numtX2, class numtY2, class numtZ2>
        point3d(const point3d<numtX2, numtY2, numtZ2>& source)
            : x(source.X()), y(source.Y()), z(source.Z()) {
        }
    };
    template<class numtX = double, class numtY = numtX, class numtZ = numtY>
    inline point3d<numtX, numtY, numtZ>make_point(const numtX& x, const numtY& y, const numtZ& z) {
        return point3d<numtX, numtY, numtZ>(x, y, z);
    }


    template<class numX = double, class numY = numX>
    using Points = Chain<point<numX, numY>>;
    template<class numX = double, class numY = numX>
    class SortedPoints : public SortedChain<point<numX, numY>>
    {
    public:
        typedef IFunction<numY, const numX&> Func;
        typedef point<numX, numY> Point;
        SortedPoints() {}
        SortedPoints(const SortedPoints&) = delete;
        SortedPoints& operator=(const SortedPoints&) = delete;
        SortedPoints(SortedChain<point<numX, numY>>&& points) :SortedChain<point<numX, numY>>(std::move(points)) {}
        SortedPoints(SortedPoints&& points) = default;
        SortedPoints& operator=(SortedChain<point<numX, numY>>&& points)
        {
            SortedChain<point<numX, numY>>::operator=(std::move(points));
            return *this;
        }
        template<class numY2>
        SortedPoints(const Points<numX, numY2>& chain)
        {
            for (const auto& x : chain)
                SortedChain<point<numX, numY>>::operator<<(point<numX, numY>(x));
        }
        template<class numX2, class numY2>
        SortedPoints(const SortedChain<point<numX2, numY2>>& chain)
        {
            for (const auto& p : chain) {
                numX a;numY b;
                a = p.X();b = p.Y();
                SortedChain<point<numX, numY>>::operator<<(make_point(a, b));
            }
        }
        SortedPoints(const Func& f, const SortedChain<numX>& chain)
        {
            for (const numX& x : chain)
                SortedChain<point<numX, numY>>::operator<<(point<numX, numY>(x, f(x)));
        }
        SortedPoints(const Func& f, const Chain<numX>& chain)
        {
            for (const numX& x : chain)
                SortedChain<point<numX, numY>>::operator<<(point<numX, numY>(x, f(x)));
        }
        template<class numY2>
        SortedPoints(const Func& f, const SortedChain<point<numX, numY2>>& chain)
        {
            for (const auto& p : chain)
                SortedChain<point<numX, numY>>::operator<<(point<numX, numY>(p.X(), f(p.X())));
        }
        template<class FUNC, class CHAIN>
        SortedPoints(FUNC F, const CHAIN& c) :SortedPoints(static_cast<const Func&>(FunctionWrap<numY, const numX&>(F)), c) {}
        virtual ~SortedPoints() {}
        inline SortedPoints Clone()const
        {
            return SortedPoints(SortedChain<point<numX, numY>>::clone());
        }
    protected:
        point<numX, numY>& Bin(const size_t i)
        {
            return SortedChain<point<numX, numY>>::accessBin(i);
        }
    public:
        Points<numY, numX> Transponate()const
        {
            Points<numY, numX> res;
            for (const auto& p : *this)
                res.push_back(point<numX, numY>(p.Y(), p.X()));
            return res;
        }
        SortedChain<point<numY, numX>> TransponateAndSort()const
        {
            SortedChain<point<numY, numX>> res;
            for (const auto& p : *this)
                res << point<numX, numY>(p.Y(), p.X());
            return res;
        }
        SortedPoints XRange(const numX& from, const numX& to)const
        {
            return SortedChain<point<numX, numY>>::SelectByCondition([from, to](const point<numX, numY>& P) {
                return (P.X() >= from) && (P.X() <= to);
                });
        }
        SortedPoints YRange(const numY& from, const numY& to)const
        {
            return SortedChain<point<numX, numY>>::SelectByCondition([from, to](const point<numX, numY>& P) {
                return (P.Y() >= from) && (P.Y() <= to);
                });
        }
        SortedPoints XExclude(const numX& from, const numX& to)const
        {
            return SortedChain<point<numX, numY>>::SelectByCondition([from, to](const point<numX, numY>& P) {
                return (P.X() < from) || (P.X() > to);
                });
        }
        SortedPoints YExclude(const numY& from, const numY& to)const
        {
            return SortedChain<point<numX, numY>>::SelectByCondition([from, to](const point<numX, numY>& P) {
                return (P.Y() < from) || (P.Y() > to);
                });
        }
        SortedPoints& FillWithValues(const numY& v)
        {
            for (size_t i = 0, n = this->size(); i < n; i++)
                Bin(i) = v;
            return *this;
        }
        SortedPoints& Transform(const std::function<numY(const numX&, const numY&)>& F)
        {
            for (size_t i = 0, n = this->size(); i < n; i++)
                Bin(i) = F(Bin(i).X(), Bin(i).Y());
            return *this;
        }
        SortedPoints& operator+=(const SortedPoints& second)
        {
            if (second.size() != this->size())
                throw Exception<SortedPoints>("Cannot perform arithmetic operation for two sets of points with different size");
            for (size_t i = 0, n = this->size(); i < n; i++) {
                Bin(i) = Bin(i) + second[i];
            }
            return *this;
        }
        SortedPoints& operator+=(const Func& f)
        {
            for (size_t i = 0, n = this->size(); i < n; i++)
                Bin(i) = Bin(i) + f(Bin(i).X());
            return *this;
        }
        SortedPoints& operator+=(const numY& c)
        {
            for (size_t i = 0, n = this->size(); i < n; i++)
                Bin(i) = Bin(i) + c;
            return *this;
        }

        SortedPoints& operator-=(const SortedPoints& second)
        {
            if (second.size() != this->size())
                throw Exception<SortedPoints>("Cannot perform arithmetic operation for two sets of points with different size");
            for (size_t i = 0, n = this->size(); i < n; i++) {
                Bin(i) = Bin(i) - second[i];
            }
            return *this;
        }
        SortedPoints& operator-=(const Func& f)
        {
            for (size_t i = 0, n = this->size(); i < n; i++)
                Bin(i) = Bin(i) - f(Bin(i).X());
            return *this;
        }
        SortedPoints& operator-=(const numY& c)
        {
            for (size_t i = 0, n = this->size(); i < n; i++)
                Bin(i) = Bin(i) - c;
            return *this;
        }

        SortedPoints& operator*=(const SortedPoints& second)
        {
            if (second.size() != this->size())
                throw Exception<SortedPoints>("Cannot perform arithmetic operation for two sets of points with different size");
            for (size_t i = 0, n = this->size(); i < n; i++) {
                Bin(i) = Bin(i) * second[i];
            }
            return *this;
        }
        SortedPoints& operator*=(const Func& f)
        {
            for (size_t i = 0, n = this->size(); i < n; i++)
                Bin(i) = Bin(i) * f(Bin(i).X());
            return *this;
        }
        SortedPoints& operator*=(const numY& c)
        {
            for (size_t i = 0, n = this->size(); i < n; i++)
                Bin(i) = Bin(i) * c;
            return *this;
        }

        SortedPoints& operator/=(const SortedPoints& second)
        {
            if (second.size() != this->size())
                throw Exception<SortedPoints>("Cannot perform arithmetic operation for two sets of points with different size");
            for (size_t i = 0, n = this->size(); i < n; i++) {
                Bin(i) = Bin(i) / second[i];
            }
            return *this;
        }
        SortedPoints& operator/=(const Func& f)
        {
            for (size_t i = 0, n = this->size(); i < n; i++)
                Bin(i) = Bin(i) / f(Bin(i).X());
            return *this;
        }
        SortedPoints& operator/=(const numY& c)
        {
            for (size_t i = 0, n = this->size(); i < n; i++)
                Bin(i) = Bin(i) / c;
            return *this;
        }
        template<class numtY2>
        SortedPoints& leftArrow(const SortedPoints<numX, numtY2>& other)
        {
            for (size_t i = 0, n = this->size(); i < n; i++) {
                Bin(i) << other[i];
            }
            return *this;
        }
        template<class numtY2>
        inline SortedPoints& leftArrow(const Points<numX, numtY2>& other)
        {
            return leftArrow(SortedPoints<numX, numtY2>(other));
        }

        template<class numtY2>
        auto operator+(const SortedPoints<numX, numtY2>& other)const->SortedPoints<numX, decltype(SortedChain<point<numX, numY>>::operator[](0).Y() + other[0].Y())>
        {
            SortedPoints<numX, decltype(SortedChain<point<numX, numY>>::operator[](0).Y() + other[0].Y())> res;
            for (size_t i = 0, n = this->size(); i < n; i++) {
                res << (SortedChain<point<numX, numY>>::operator[](i) + other[i]);
            }
            return res;
        }
        template<class numtY2>
        auto operator-(const SortedPoints<numX, numtY2>& other)const->SortedPoints<numX, decltype(SortedChain<point<numX, numY>>::operator[](0).Y() - other[0].Y())>
        {
            SortedPoints<numX, decltype(SortedChain<point<numX, numY>>::operator[](0).Y() - other[0].Y())> res;
            for (size_t i = 0, n = this->size(); i < n; i++) {
                res << (SortedChain<point<numX, numY>>::operator[](i) - other[i]);
            }
            return res;
        }
        template<class numtY2>
        auto operator*(const SortedPoints<numX, numtY2>& other)const->SortedPoints<numX, decltype(SortedChain<point<numX, numY>>::operator[](0).Y()* other[0].Y())>
        {
            SortedPoints<numX, decltype(SortedChain<point<numX, numY>>::operator[](0).Y()* other[0].Y())> res;
            for (size_t i = 0, n = this->size(); i < n; i++) {
                res << (SortedChain<point<numX, numY>>::operator[](i) * other[i]);
            }
            return res;
        }
        template<class numtY2>
        auto operator/(const SortedPoints<numX, numtY2>& other)const->SortedPoints<numX, decltype(SortedChain<point<numX, numY>>::operator[](0).Y() / other[0].Y())>
        {
            SortedPoints<numX, decltype(SortedChain<point<numX, numY>>::operator[](0).Y() / other[0].Y())> res;
            for (size_t i = 0, n = this->size(); i < n; i++) {
                res << (SortedChain<point<numX, numY>>::operator[](i) / other[i]);
            }
            return res;
        }


        template<class numtY2>
        auto operator+(const numtY2& other)const->SortedPoints<numX, decltype(SortedChain<point<numX, numY>>::operator[](0).Y() + other)>
        {
            SortedPoints<numX, decltype(SortedChain<point<numX, numY>>::operator[](0).Y() + other)> res;
            for (size_t i = 0, n = this->size(); i < n; i++) {
                res << (SortedChain<point<numX, numY>>::operator[](i) + other);
            }
            return res;
        }
        template<class numtY2>
        auto operator-(const numtY2& other)const->SortedPoints<numX, decltype(SortedChain<point<numX, numY>>::operator[](0).Y() - other)>
        {
            SortedPoints<numX, decltype(SortedChain<point<numX, numY>>::operator[](0).Y() - other)> res;
            for (size_t i = 0, n = this->size(); i < n; i++) {
                res << (SortedChain<point<numX, numY>>::operator[](i) - other);
            }
            return res;
        }
        template<class numtY2>
        auto operator*(const numtY2& other)const->SortedPoints<numX, decltype(SortedChain<point<numX, numY>>::operator[](0).Y()* other)>
        {
            SortedPoints<numX, decltype(SortedChain<point<numX, numY>>::operator[](0).Y()* other)> res;
            for (size_t i = 0, n = this->size(); i < n; i++) {
                res << (SortedChain<point<numX, numY>>::operator[](i) * other);
            }
            return res;
        }
        template<class numtY2>
        auto operator/(const numtY2& other)const->SortedPoints<numX, decltype(SortedChain<point<numX, numY>>::operator[](0).Y() / other)>
        {
            SortedPoints<numX, decltype(SortedChain<point<numX, numY>>::operator[](0).Y() / other)> res;
            for (size_t i = 0, n = this->size(); i < n; i++) {
                res << (SortedChain<point<numX, numY>>::operator[](i) / other);
            }
            return res;
        }
        SortedPoints operator+(const Func& other)const { return SortedPoints(*this) += other; }
        SortedPoints operator-(const Func& other)const { return SortedPoints(*this) -= other; }
        SortedPoints operator*(const Func& other)const { return SortedPoints(*this) *= other; }
        SortedPoints operator/(const Func& other)const { return SortedPoints(*this) /= other; }

        //In case of types with uncertainty
        SortedPoints CloneEmptyBins()const
        {
            SortedPoints res;
            for (const auto& P : *this)res << make_point(P.X(), std_error(numY(0).val()));
            return res;
        }
        numY TotalSum(std::function<void(numY&)> prepare = [](numY&) {})const
        {
            numY res = 0;
            prepare(res);
            for (const auto& P : *this)res += P.Y();
            return res;
        }
        SortedPoints Scale(const size_t sc_x)const
        {
            SortedChain<numX> new_x, sorted_x;
            for (const auto& item : *this)
                sorted_x << item.X();
            for (size_t i = sc_x - 1, n = sorted_x.size(); i < n; i += sc_x) {
                const auto min = sorted_x[i + 1 - sc_x].min(),
                    max = sorted_x[i].max();
                new_x << numX((max + min) / 2, (max - min) / 2);
            }
            SortedPoints res([](const numX&) {return numY(0);}, new_x);
            for (size_t i = 0; i < new_x.size(); i++) {
                auto v = numY(0).val();
                for (size_t ii = 0; ii < sc_x; ii++)
                    v += this->operator[](i* sc_x + ii).Y().val();
                res.Bin(i) = std_error(v);
            }
            return res;
        }
        SortedPoints& imbibe(const SortedPoints& second)
        {
            for (int i = 0, n = this->size(); i < n; i++) {
                if (this->operator[](i).X() == second[i].X()) {
                    this->Bin(i) = std_error(this->operator[](i).Y().val() + second[i].Y().val());
                }
                else
                    throw Exception<SortedPoints>("Cannot imbibe histogram. bins differ");
            }
            return *this;
        }
        auto toLine()const
        {
            SortedPoints<typename numX::NumberType, typename numY::NumberType> res;
            for (int i = 0, n = this->size(); i < n; i++) {
                const auto& P = this->operator[](i);
                res << make_point(P.X().val(), P.Y().val());
            }
            return res;
        }
        auto removeXerorbars()const
        {
            SortedPoints<typename numX::NumberType, numY> res;
            for (int i = 0, n = this->size(); i < n; i++) {
                const auto& P = this->operator[](i);
                res << make_point(P.X().val(), P.Y());
            }
            return res;
        }
        auto removeYerorbars()const
        {
            SortedPoints<numX, typename numY::NumberType> res;
            for (int i = 0, n = this->size(); i < n; i++) {
                const auto& P = this->operator[](i);
                res << make_point(P.X(), P.Y().val());
            }
            return res;
        }
    };
    template<class numtX = double, class numtY = numtX, class numtZ = numtY>
    class BiSortedPoints
    {
    private:
        SortedChain<numtX> m_x_axis;
        SortedChain<numtY> m_y_axis;
        Chain<Chain<numtZ>> m_data;
        template<typename...Args>
        void init(Args...args)
        {
            m_data.clear();
            for (size_t i = 0, I = m_x_axis.size(); i < I; i++) {
                m_data.push_back(Chain<numtZ>());
                for (size_t j = 0, J = m_y_axis.size(); j < J; j++)
                    m_data[m_data.size() - 1].push_back(numtZ(args...));
            }
        }
    public:
        const SortedChain<numtX>& X()const
        {
            return m_x_axis;
        }
        const SortedChain<numtY>& Y()const
        {
            return m_y_axis;
        }
        template<typename...Args>
        BiSortedPoints(SortedChain<numtX>&& X, SortedChain<numtY>&& Y, Args...args)
            : m_x_axis(std::move(X)), m_y_axis(std::move(Y))
        {
            init(args...);
        }
        BiSortedPoints() : BiSortedPoints({}, {}) {}
        BiSortedPoints(const BiSortedPoints&) = delete;
        BiSortedPoints& operator=(const BiSortedPoints&) = delete;
        BiSortedPoints(BiSortedPoints&&) = default;
        BiSortedPoints& operator=(BiSortedPoints&&) = default;
        BiSortedPoints clone() const
        {
            BiSortedPoints res;
            res.m_x_axis = m_x_axis.clone();
            res.m_y_axis = m_y_axis.clone();
            for (size_t i = 0, I = m_data.size(); i < I; i++) {
                res.m_data.push_back(m_data[i]);
            }
            return res;
        }
        virtual ~BiSortedPoints() {}
        typedef typename Chain<Chain<numtZ>>::const_iterator const_iterator;
        const_iterator begin()const
        {
            return m_data.cbegin();
        }
        const_iterator cbegin()const
        {
            return m_data.cbegin();
        }
        const_iterator end() const
        {
            return m_data.cend();
        }
        const_iterator cend() const
        {
            return m_data.cend();
        }
        size_t size()const
        {
            return m_data.size();
        }
        const Chain<numtZ>& operator[](const size_t i)const
        {
            if (size() <= i)throw Exception<BiSortedPoints>("X-range check error");
            return m_data[i];
        }
        numtZ& Bin(const size_t i, const size_t j)
        {
            if (size() <= i)throw Exception<BiSortedPoints>("X-range check error " + std::to_string(i) + ":" + std::to_string(j));
            if (m_data[i].size() <= j)throw Exception<BiSortedPoints>("Y-range check error " + std::to_string(i) + ":" + std::to_string(j));
            return m_data[i][j];
        }
        SortedPoints<numtX, numtZ>CutX(const size_t j)const
        {
            if (m_y_axis.size() <= j)throw Exception<BiSortedPoints>("range check error");
            SortedPoints<numtX, numtZ> res;
            for (size_t i = 0, I = m_x_axis.size(); i < I; i++) {
                res << point<numtX, numtZ>(m_x_axis[i], m_data[i][j]);
            }
            return res;
        }
        SortedPoints<numtY, numtZ>CutY(const size_t i)const
        {
            if (m_x_axis.size() <= i)throw Exception<BiSortedPoints>("range check error");
            SortedPoints<numtY, numtZ> res;
            for (size_t j = 0, J = m_y_axis.size(); j < J; j++) {
                res << point<numtY, numtZ>(m_y_axis[j], m_data[i][j]);
            }
            return res;
        }
        point3d<numtX, numtY, numtZ> operator()(const size_t i, const size_t j)const
        {
            if (size() <= i)throw Exception<BiSortedPoints>("range check error");
            if (m_y_axis.size() <= j)throw Exception<BiSortedPoints>("range check error");
            return point3d<numtX, numtY, numtZ>(m_x_axis[i], m_y_axis[j], m_data[i][j]);
        }
        void FullCycle(const std::function<void(const point3d<numtX, numtY, numtZ>&)>f)const
        {
            for (size_t i = 0, I = m_x_axis.size(); i < I; i++)
                for (size_t j = 0, J = m_y_axis.size(); j < J; j++) {
                    point3d<numtX, numtY, numtZ> P(m_x_axis[i], m_y_axis[j], m_data[i][j]);
                    f(P);
                }
        }
        void FullCycle(const std::function<void(const numtX&, const numtY&, const numtZ&)>f)const
        {
            for (size_t i = 0, I = m_x_axis.size(); i < I; i++)
                for (size_t j = 0, J = m_y_axis.size(); j < J; j++)
                    f(m_x_axis[i], m_y_axis[j], m_data[i][j]);
        }
        BiSortedPoints& FullCycleVar(const std::function<void(const numtX&, const numtY&, numtZ&)>f)
        {
            for (size_t i = 0, I = m_x_axis.size(); i < I; i++) {
                for (size_t j = 0, J = m_y_axis.size(); j < J; j++)
                    f(m_x_axis[i], m_y_axis[j], m_data[i][j]);
            }
            return *this;
        }
        BiSortedPoints& FullCycleVar(const std::function<void(numtZ&)>f)
        {
            for (size_t i = 0, I = m_x_axis.size(); i < I; i++) {
                for (size_t j = 0, J = m_y_axis.size(); j < J; j++)
                    f(m_data[i][j]);
            }
            return *this;
        }
        //In case of types with uncertainty
        BiSortedPoints Scale(const size_t sc_x, const size_t sc_y)const
        {
            SortedChain<numtX> new_x;
            for (size_t i = sc_x - 1, n = this->X().size(); i < n; i += sc_x) {
                const auto min = this->X()[i + 1 - sc_x].min(),
                    max = this->X()[i].max();
                new_x << (numtX((max + min) / 2, (max - min) / 2));
            }
            SortedChain<numtY> new_y;
            for (size_t i = sc_y - 1, n = this->Y().size(); i < n; i += sc_y) {
                const auto min = this->Y()[i + 1 - sc_y].min(),
                    max = this->Y()[i].max();
                new_y << (numtY((max + min) / 2, (max - min) / 2));
            }
            BiSortedPoints res(new_x.clone(), new_y.clone());
            for (size_t i = 0; i < new_x.size(); i++)for (size_t j = 0; j < new_y.size(); j++) {
                auto v = numtZ(0).val();
                for (size_t ii = 0; ii < sc_x; ii++) {
                    for (size_t jj = 0; jj < sc_y; jj++)
                        v += this->operator[](i* sc_x + ii)[j * sc_y + jj].val();
                }
                res.Bin(i, j) = std_error(v);
            }
            return res;
        }
        BiSortedPoints& imbibe(const BiSortedPoints& second)
        {
            if ((this->X().size() != second.X().size()) || (this->Y().size() != second.Y().size()))
                throw Exception<BiSortedPoints>("cannot imbibe second histogram: bins differ");
            for (int i = 0, n = this->size(); i < n; i++)for (int j = 0, m = this->operator[](i).size(); j < m; j++)
                this->Bin(i, j) = std_error(this->operator[](i)[j].val() + second[i][j].val());
            return *this;
        }
    };


    template<class numX, class numY>
    SortedPoints<typename numX::NumberType, typename numY::NumberType> toLine(const SortedPoints<numX, numY>& source)
    {
        SortedPoints<typename numX::NumberType, typename numY::NumberType> res;
        for (int i = 0, n = source.size(); i < n; i++) {
            const auto& P = source[i];
            res << make_point(P.X().val(), P.Y().val());
        }
        return res;
    }
    template<class numX, class numY>
    SortedPoints<typename numX::NumberType, numY> removeXerorbars(const SortedPoints<numX, numY>& source)
    {
        SortedPoints<typename numX::NumberType, numY> res;
        for (int i = 0, n = source.size(); i < n; i++) {
            const auto& P = source[i];
            res << make_point(P.X().val(), P.Y());
        }
        return res;
    }
    template<class numX, class numY>
    SortedPoints<numX, typename numY::NumberType> removeYerorbars(const SortedPoints<numX, numY>& source)
    {
        SortedPoints<numX, typename numY::NumberType> res;
        for (int i = 0, n = source.size(); i < n; i++) {
            const auto& P = source[i];
            res << make_point(P.X(), P.Y().val());
        }
        return res;
    }


    template<class numtX, class numtY>
    class hist_avr_calculator :public SortedPoints<numtX, WeightedAverage<numtY>>
    {
    public:
        inline hist_avr_calculator(const SortedPoints<numtX, value<numtY>>& source)
            :SortedPoints<numtX, WeightedAverage<numtY>>(source) {}
        template<typename...Args>
        inline hist_avr_calculator(const SortedPoints<numtX, value<numtY>>& source, const Args&...args)
            : hist_avr_calculator(args...) {
            SortedPoints<numtX, WeightedAverage<numtY>>::leftArrow(source);
        }
    };
    template<class numtX, class numtY, typename...Args>
    inline hist_avr_calculator<numtX, numtY> hist_avr(const SortedPoints<numtX, value<numtY>>& source, const Args&...args)
    {
        return hist_avr_calculator<numtX, numtY>(source, args...);
    }
    template<class numtX, class numtY, typename...Args>
    inline hist_avr_calculator<numtX, numtY> hist_avr(const Points<numtX, value<numtY>>& source, const Args&...args)
    {
        return hist_avr_calculator<numtX, numtY>(SortedPoints<numtX, value<numtY>>(source), args...);
    }
    template<class numtX, class numtY>
    class hist_stdev_calculator :public SortedPoints<numtX, StandardDeviation<numtY>>
    {
    public:
        hist_stdev_calculator(const SortedPoints<numtX, numtY>& source)
            :SortedPoints<numtX, StandardDeviation<numtY>>(source) {}
        template<class Arg, typename...Args>
        hist_stdev_calculator(const Arg& source, Args...args)
            : SortedPoints<numtX, StandardDeviation<numtY>>(args...) {
            SortedPoints<numtX, StandardDeviation<numtY>>::leftArrow(source);
        }
    };
    template<class numtX, class numtY, typename...Args>
    inline hist_stdev_calculator<numtX, numtY> hist_stdev(const SortedPoints<numtX, numtY>& source, Args...args)
    {
        return hist_stdev_calculator<numtX, numtY>(source, args...);
    }
    template<class numtX, class numtY, typename...Args>
    inline hist_stdev_calculator<numtX, numtY> hist_stdev(const Points<numtX, numtY>& source, Args...args)
    {
        return hist_stdev_calculator<numtX, numtY>(SortedPoints<numtX, numtY>(source), args...);
    }
    template<class numtX = double, class numtY = numtX>
    using hist = SortedPoints<value<numtX>, value<numtY>>;

    template<class numtX = double, class numtY = numtX, class numtZ = numtY>
    using hist2d = BiSortedPoints<value<numtX>, value<numtY>, value<numtZ>>;

    template<class numtX = double, class numtY = numtX>
    class Distribution1D : public SortedPoints<value<numtX>, value<numtY>>
    {
    private:
        unsigned long long counter;
    public:
        Distribution1D(const Chain<value<numtX>>& data)
            :SortedPoints<value<numtX>, value<numtY>>([](const value<numtX>&) {return std_error(numtY(0));}, data)
        {
            counter = 0;
        }
        Distribution1D(const SortedChain<value<numtX>>& data)
            :SortedPoints<value<numtX>, value<numtY>>([](const value<numtX>&) {return std_error(numtY(0));}, data)
        {
            counter = 0;
        }
        Distribution1D& Fill(const numtX& v)
        {
            counter++;
            for (size_t i = 0, n = this->size(); i < n; i++) {
                if (this->Bin(i).X().Contains(v))
                    this->Bin(i) = std_error(this->Bin(i).Y().val() + numtY(1));
            }
            return *this;
        }
        unsigned long long Entries()const
        {
            return counter;
        }
    };
    template<class numtX = double, class numtY = numtX, class numtZ = numtY>
    class Distribution2D : public hist2d<numtX, numtY, numtZ>
    {
    private:
        unsigned long long counter;
        void init()
        {
            counter = 0;
            this->FullCycleVar([](const value<numtX>&, const value<numtY>&, value<numtZ>& z) {
                z = std_error(numtZ(0));
                });
        }
    public:
        Distribution2D(const Chain<value<numtX>>& X, const Chain<value<numtY>>& Y)
            : hist2d<numtX, numtY, numtZ>(X, Y)
        {
            init();
        }
        Distribution2D(SortedChain<value<numtX>>&& X, SortedChain<value<numtY>>&& Y)
            : hist2d<numtX, numtY, numtZ>(std::move(X), std::move(Y))
        {
            init();
        }
        Distribution2D& Fill(const numtX& x, const numtY& y)
        {
            counter++;
            for (size_t i = 0, I = this->X().size(); i < I; i++)if (this->X()[i].Contains(x)) {
                for (size_t j = 0, J = this->Y().size(); j < J; j++)if (this->Y()[j].Contains(y))
                    this->Bin(i, j) = std_error(this->operator[](i)[j].val() + numtZ(1));
            }
            return *this;
        }
        unsigned long long Entries()const
        {
            return counter;
        }
    };

    template<class KEY, class DATA>
    class Cache {
    private:
        std::map<KEY, DATA> m_map;
    public:
        Cache() {}
        virtual ~Cache() {}
        template<class FUNC>
        inline const DATA& operator()(const KEY& key, FUNC F) {
            if (m_map.find(key) == m_map.end())m_map.insert(std::make_pair(key, F()));
            return m_map[key];
        }
    };

};
#endif
