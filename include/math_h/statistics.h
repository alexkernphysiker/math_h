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
    const numt&operator[](const size_t index)const{
	if(index>=m_list.size())throw Exception<AverageObtainer>("range check error");
	return m_list[index];
    }
    inline size_t size()const{return m_list.size();}
    inline size_t count()const{return m_list.size();}
    template<typename... Args>
    AverageObtainer(Args...args):m_cache_avr(0,numt(args...)){}
    virtual ~AverageObtainer() {}
    inline void Fill(const numt &x){m_list.push_back(x);}
    const numt&Average()const{
	const size_t sz = count();
	if(sz==0)throw Exception<AverageObtainer>("Cannot obtain the average of empty sample");
        if (m_cache_avr.first!=sz) {
            numt res = m_list[0];
	    for(size_t i=1;i<sz;i++)res=res+m_list[i];
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
    template<typename...Args>
    inline void Fill(Args...args){AverageObtainer<Vector<N,numt>>::Fill(vec(args...));}
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

template<size_t Nx,size_t Ny,class numt=double>
using DiVector=std::pair<Vector<Nx,numt>,Vector<Ny,numt>>;
template<size_t Nx,size_t Ny,class numt>
inline DiVector<Nx,Ny,numt> operator+(const DiVector<Nx,Ny,numt>&A,const DiVector<Nx,Ny,numt>&B){
    return DiVector<Nx,Ny,numt>(A.first+B.first,A.second+B.second);
}
template<size_t Nx,size_t Ny,class numt>
inline DiVector<Nx,Ny,numt> operator/(const DiVector<Nx,Ny,numt>&A,const size_t b){
    return DiVector<Nx,Ny,numt>(A.first/b,A.second/b);
}

template<size_t Nx,size_t Ny,class numt=double>
class SamplingXY:public AverageObtainer<DiVector<Nx,Ny,numt>>{
private:
    std::pair<size_t,Matrix<Nx,Vector<Nx,numt>>> m_cache_cov_x;
    std::pair<size_t,Matrix<Ny,Vector<Ny,numt>>> m_cache_cov_y;
    std::pair<size_t,Matrix<Nx,Vector<Ny,numt>>> m_cache_cov_xy;
public:
    enum {DimensionsX = Nx};
    enum {DimensionsY = Ny};
    typedef numt NumberType;
    SamplingXY():AverageObtainer<DiVector<Nx,Ny,numt>>(Vector<Nx,numt>::zero(),Vector<Ny,numt>::zero())
	,m_cache_cov_x(0,Matrix<Nx,Vector<Nx,numt>>::zero())
	,m_cache_cov_y(0,Matrix<Ny,Vector<Ny,numt>>::zero())
	,m_cache_cov_xy(0,Matrix<Nx,Vector<Ny,numt>>::zero())
	{}
    virtual ~SamplingXY(){}
    template<typename...Args>
    inline void Fill(Args...args){AverageObtainer<DiVector<Nx,Ny,numt>>::Fill(DiVector<Nx,Ny,numt>(args...));}
    const Matrix<Nx,Vector<Nx,numt>>&CovX()const{
	const size_t sz = AverageObtainer<DiVector<Nx,Ny,numt>>::count();
	if(sz<2)throw Exception<SamplingXY>("Cannot obtain Cov matrix for sample less than 2 elements");
        if (m_cache_cov_x.first!=sz) {
	    auto res=Matrix<Nx,Vector<Nx,numt>>::zero();
	    const auto&avr=AverageObtainer<DiVector<Nx,Ny,numt>>::Average();
	    AverageObtainer<DiVector<Nx,Ny,numt>>::ForEach(
		[&res,&avr](const DiVector<Nx,Ny,numt>&item){
		    const auto d=item.first-avr.first;
		    res=(res+(columns(d)*rows(d)));
		}
	    );
            const_cast<SamplingXY&>(*this).m_cache_cov_x =
                std::make_pair(sz,res/(sz-1));
        }
        return m_cache_cov_x.second;
    }
    const Matrix<Ny,Vector<Ny,numt>>&CovY()const{
	const size_t sz = AverageObtainer<DiVector<Nx,Ny,numt>>::count();
	if(sz<2)throw Exception<SamplingXY>("Cannot obtain Cov matrix for sample less than 2 elements");
        if (m_cache_cov_y.first!=sz) {
	    auto res=Matrix<Ny,Vector<Ny,numt>>::zero();
	    const auto&avr=AverageObtainer<DiVector<Nx,Ny,numt>>::Average();
	    AverageObtainer<DiVector<Nx,Ny,numt>>::ForEach(
		[&res,&avr](const DiVector<Nx,Ny,numt>&item){
		    const auto d=item.second-avr.second;
		    res=(res+(columns(d)*rows(d)));
		}
	    );
            const_cast<SamplingXY&>(*this).m_cache_cov_y =
                std::make_pair(sz,res/(sz-1));
        }
        return m_cache_cov_y.second;
    }
    const Matrix<Nx,Vector<Ny,numt>>&CovXY()const{
	const size_t sz = AverageObtainer<DiVector<Nx,Ny,numt>>::count();
	if(sz<2)throw Exception<SamplingXY>("Cannot obtain Cov matrix for sample less than 2 elements");
        if (m_cache_cov_xy.first!=sz) {
	    auto res=Matrix<Nx,Vector<Ny,numt>>::zero();
	    const auto&avr=AverageObtainer<DiVector<Nx,Ny,numt>>::Average();
	    AverageObtainer<DiVector<Nx,Ny,numt>>::ForEach(
		[&res,&avr](const DiVector<Nx,Ny,numt>&item){
		    const auto dx=item.first-avr.first;
		    const auto dy=item.second-avr.second;
		    res=(res+(columns(dx)*rows(dy)));
		}
	    );
            const_cast<SamplingXY&>(*this).m_cache_cov_xy =
                std::make_pair(sz,res/(sz-1));
        }
        return m_cache_cov_xy.second;
    }
};
}
#endif
