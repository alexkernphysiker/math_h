// this file is distributed under
// LGPLv3 license
#ifndef ______STATISTICS_H_______
#define ______STATISTICS_H_______
#include <iostream>
#include <memory>
#include <utility>
#include <vector>
#include <functional>
#include <math.h>
#include "error.h"
#include "vectors.h"
#include "matrices.h"
namespace MathTemplates
{
template<class comparable = double>
using Chain=std::vector<comparable>;

template<class numt = double>
class AverageObtainer
{
private:
    Chain<numt> m_list;
    std::pair<size_t,numt> m_cache_avr;
public:
    typedef typename Chain<numt>::const_iterator const_iterator;
    inline const_iterator begin()const{return m_list.begin();}
    inline const_iterator end() const{return m_list.end();}
    inline const numt&operator[](const size_t index)const{
	if(index>=m_list.size())throw Exception<AverageObtainer>("range check error");
	return m_list[index];
    }
    inline size_t size()const{return m_list.size();}
    inline size_t count()const{return m_list.size();}
    template<typename... Args>
    AverageObtainer(Args...args):m_cache_avr(0,numt(args...)){}
    virtual ~AverageObtainer() {}
    inline AverageObtainer &operator<<(const numt &x)
    {
        m_list.push_back(x);
        return *this;
    }
    const numt&Average()const{
	const size_t sz = count();
	if(sz==0)throw Exception<AverageObtainer>("Cannot obtain the average of empty sample");
        if (m_cache_avr.first!=sz) {
            numt res = m_list[0];
	    for(size_t i=1;i<sz;i++)res+=m_list[i];
            const_cast<AverageObtainer&>(*this).m_cache_avr =
                std::make_pair(sz,res/sz);
        }
        return m_cache_avr.second;
    }
protected:
    void ForEach(const std::function<void(const numt&)> F)const{
	for(const numt&x:m_list)F(x);
    }
};

template<size_t N,class numt = double>
class Sampling:public AverageObtainer<Vector<N,numt>>{
private:
    std::pair<size_t,Matrix<N,Vector<N,numt>>> m_cache_cov;
public:
    enum {Dimensions = N};
    typedef numt NumberType;
    Sampling():AverageObtainer<Vector<N,numt>>(Vector<N,numt>::zero()),m_cache_cov(0,Matrix<N,Vector<N,numt>>::zero()){}
    virtual ~Sampling(){}
    const Matrix<N,Vector<N,numt>>&Cov()const{
	const size_t sz = AverageObtainer<Vector<N,numt>>::count();
	if(sz<2)throw Exception<Sampling>("Cannot obtain Cov matrix for sample less than 2 elements");
        if (m_cache_cov.first!=sz) {
	    auto res=Matrix<N,Vector<N,numt>>::zero();
	    const auto&avr=AverageObtainer<Vector<N,numt>>::Average();
	    AverageObtainer<Vector<N,numt>>::ForEach(
		[&res,&avr](const Vector<N,numt>&x){
		    const auto d=x-avr;
		    res=(res+(columns(d)*rows(d)));
		}
	    );
            const_cast<Sampling&>(*this).m_cache_cov =
                std::make_pair(sz,res/(sz-1));
        }
        return m_cache_cov.second;
    }
};
}
#endif
