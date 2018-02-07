// this file is distributed under
// LGPLv3 license
#ifndef ____HISTS_H_____
#	define ____HISTS_H_____
#include "error.h"
#include "sigma.h"
#include "tabledata.h"
namespace MathTemplates
{
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
template<class numtX, class numtY>
class hist_avr_calculator:public SortedPoints<numtX, WeightedAverage<numtY>>
{
public:
    inline hist_avr_calculator(const SortedPoints<numtX, value<numtY>>&source)
    :SortedPoints<numtX, WeightedAverage<numtY>>(source){}
    template<typename...Args>
    inline hist_avr_calculator(const SortedPoints<numtX, value<numtY>>&source,Args...args)
    :SortedPoints<numtX, WeightedAverage<numtY>>(args...){
	SortedPoints<numtX, WeightedAverage<numtY>>::leftArrow(source);
    }
};
template<class numtX, class numtY,typename...Args>
inline hist_avr_calculator<numtX,numtY> hist_avr(const SortedPoints<numtX, value<numtY>>&source,Args...args)
{
    return hist_avr_calculator<numtX,numtY>(source,args...);
}
template<class numtX, class numtY,typename...Args>
inline hist_avr_calculator<numtX,numtY> hist_avr(const Points<numtX, value<numtY>>&source,Args...args)
{
    return hist_avr_calculator<numtX,numtY>(SortedPoints<numtX, value<numtY>>(source),args...);
}
template<class numtX, class numtY>
class hist_stdev_calculator:public SortedPoints<numtX, StandardDeviation<numtY>>
{
public:
    hist_stdev_calculator(const SortedPoints<numtX, numtY>&source)
    :SortedPoints<numtX, StandardDeviation<numtY>>(source){}
    template<class Arg,typename...Args>
    hist_stdev_calculator(const Arg&source,Args...args)
    :SortedPoints<numtX, StandardDeviation<numtY>>(args...){
	SortedPoints<numtX, StandardDeviation<numtY>>::leftArrow(source);
    }
};
template<class numtX, class numtY,typename...Args>
inline hist_stdev_calculator<numtX,numtY> hist_stdev(const SortedPoints<numtX,numtY>&source,Args...args)
{
    return hist_stdev_calculator<numtX,numtY>(source,args...);
}
template<class numtX, class numtY,typename...Args>
inline hist_stdev_calculator<numtX,numtY> hist_stdev(const Points<numtX,numtY>&source,Args...args)
{
    return hist_stdev_calculator<numtX,numtY>(SortedPoints<numtX,numtY>(source),args...);
}
template<class numtX = double, class numtY = numtX>
using hist=SortedPoints<value<numtX>, value<numtY>>;

template<class numtX = double, class numtY = numtX, class numtZ = numtY>
using hist2d=BiSortedPoints<value<numtX>, value<numtY>, value<numtZ>>;

template<class numtX = double, class numtY = numtX>
class Distribution1D: public SortedPoints<value<numtX>, value<numtY>>
{
private:
    unsigned long long counter;
public:
    Distribution1D(const Chain<value<numtX>> &data)
	:SortedPoints<value<numtX>, value<numtY>>([](const value<numtX>&){return std_error(numtY(0));},data)
    {
        counter = 0;
    }
    Distribution1D(const SortedChain<value<numtX>> &data)
	:SortedPoints<value<numtX>, value<numtY>>([](const value<numtX>&){return std_error(numtY(0));},data)
    {
        counter = 0;
    }
    Distribution1D &Fill(const numtX &v)
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
class Distribution2D: public hist2d<numtX, numtY, numtZ>
{
private:
    unsigned long long counter;
    void init()
    {
        counter = 0;
        this->FullCycleVar([](const value<numtX> &, const value<numtY> &, value<numtZ> &z) {
            z = std_error(numtZ(0));
        });
    }
public:
    Distribution2D(const Chain<value<numtX>> &X, const Chain<value<numtY>> &Y)
        : hist2d<numtX, numtY, numtZ>(X, Y)
    {
        init();
    }
    Distribution2D(const SortedChain<value<numtX>> &X, const SortedChain<value<numtY>> &Y)
        : hist2d<numtX, numtY, numtZ>(X, Y)
    {
        init();
    }
    Distribution2D &Fill(const numtX &x, const numtY &y)
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
};
#endif
