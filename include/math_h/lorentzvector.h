// this file is distributed under
// LGPLv3 license
#ifndef ___________LORENTZVECTOR_H_____
#	define ___________LORENTZVECTOR_H_____
#include <tuple>
#include <math.h>
#include <type_traits>
#include "error.h"
#include "vectors.h"
#include "vectortransformations.h"
namespace MathTemplates
{

template<class Time = double, class Space = Vector<3, double>>
class LorentzVector
{
public:
    typedef Space SpaceVectorType;
    typedef Time TimeCoordinateType;

private:
    TimeCoordinateType m_time;
    SpaceVectorType m_space;

public:
    virtual ~LorentzVector() {}

    inline LorentzVector(const Time &t, const Space &S)
        : m_time(t), m_space(S) {}

    inline const auto &E()const
    {
        return m_time;
    }

    inline const auto &P()const
    {
        return m_space;
    }

    template<class Time2, class Space2>
    inline LorentzVector(const LorentzVector<Time2, Space2> &source)
        : m_time(source.E()), m_space(source.P()) {}

    LorentzVector &operator=(const LorentzVector &source)
    {
        m_space = source.m_space;
        m_time = source.m_time;
        return *this;
    }

    bool operator==(const LorentzVector &second)const
    {
        return (m_time == second.m_time) && (m_space == second.m_space);
    }
    bool operator!=(const LorentzVector &second)const
    {
        return !operator==(second);
    }

    LorentzVector &operator+=(const LorentzVector &second)
    {
        m_time += second.m_time;
        m_space += second.m_space;
        return *this;
    }

    LorentzVector operator+(const LorentzVector &second)const
    {
        return LorentzVector(m_time + second.m_time, m_space + second.m_space);
    }

    LorentzVector &operator-=(const LorentzVector &second)
    {
        m_time -= second.m_time;
        m_space -= second.m_space;
        return *this;
    }

    LorentzVector operator-(const LorentzVector &second)const
    {
        return LorentzVector(m_time - second.m_time, m_space - second.m_space);
    }

    auto operator*(const LorentzVector &second)const
    {
        return (m_time * second.m_time) - (m_space * second.m_space);
    }

    inline auto M_sqr()const
    {
        return operator*(*this);
    }

    inline auto M()const
    {
        return sqrt(M_sqr());
    }

    inline auto Ekin()const
    {
	    return E()-M();
    }

    inline static auto zero()
    {
        return LorentzVector(Time(0), Space::zero());
    }

    template<class...Args>
    inline auto Rotate(Args...args)const
    {
        return LorentzVector(E(), Rotation(args...) * P());
    }

    auto Transform(const Space &Beta)const
    {
        const Time beta = Beta.length();
        if (beta == 0.0)return *this;
        if (beta >= Time(1))throw Exception<LorentzVector>("Bad Lorentz transformation");
        const auto bn = Beta / beta;
        const Time gamma = Time(1) / sqrt(Time(1) - beta * beta);
        const auto ST = ONE<Space::Dimensions,typename Space::NumberType>() + TensorProduct(bn, bn) * (gamma - Time(1));
        const auto TT = -Beta * gamma;
        return LorentzVector((E() * gamma) + (P() * TT), (ST * P()) + (TT * E()));
    }

    auto Beta()const
    {
        return P() / E();
    }

};

template<class Time = double, class Space = Vector<3, Time>>
inline auto lorentzVector(const Time &t, const Space &s)
{
    return LorentzVector<Time, Space>(t, s);
}

template<size_t dim,class numt = double>
inline auto lorentz_rest(const numt &l4)
{
    return LorentzVector<numt,Vector<dim,numt>>(l4 * l4,Vector<dim,numt>::zero());
}

template<class numt = double>
inline auto lorentz_Rest(const numt &l4)
{
    return lorentz_rest<3,numt>(l4);
}

template<class Time = double, class Space = Vector<3, double>>
inline auto lorentz_byPM(const Space &s, const Time &l4)
{
    return LorentzVector<Time, Space>(sqrt(s.length_sqr() + l4 * l4), s);
}

template<class Time = double, class Space = Vector<3, double>>
inline auto lorentz_byEM(const Time &t, const Time &l4, const typename Space::DType &dir)
{
    Time Sp = sqrt(t * t - l4 * l4);
    return LorentzVector<Time, Space>(t, dir * Sp);
}

template<class Time = double, class Space = Vector<3, double>>
inline auto lorentz_byEM(const Time &t, const Time &l4, const Space &Dir)
{
    return lorentz_byEM(t, l4, direction(Dir));
}

template<class Time = double, class Space = Vector<3, double>>
inline auto lorentz_byEkM(const Time &e, const Time &l4, const typename Space::DType &dir)
{
    return lorentz_byEM(e+l4, l4, dir);
}

template<class Time = double, class Space = Vector<3, double>>
inline auto lorentz_byEkM(const Time &e, const Time &l4, const Space &Dir)
{
    return lorentz_byEM(e+l4, l4, direction(Dir));
}

template<class numt>
auto binaryDecay(const numt &IM, const numt &m1, const numt &m2)
{
    if (m1 < 0)throw Exception<std::pair<LorentzVector<numt, Vector<1, numt>>, LorentzVector<numt, Vector<1, numt>>>>("Negative mass1 error");
    if (m2 < 0)throw Exception<std::pair<LorentzVector<numt, Vector<1, numt>>, LorentzVector<numt, Vector<1, numt>>>>("Negative mass2 error");
    if (IM < (m1 + m2))throw Exception<std::pair<LorentzVector<numt, Vector<1, numt>>, LorentzVector<numt, Vector<1, numt>>>>("Invariant mass of decaying system is less then masses of products");
    const auto E1 = (IM * IM - m2 * m2 + m1 * m1) / (IM * numt(2));
    const auto p = sqrt(E1 * E1 - m1 * m1);
    return std::make_pair(lorentz_byPM(vec(p), m1), lorentz_byPM(-vec(p), m2));
}

template<size_t size, class numt>
inline auto binaryDecay(const numt &IM, const numt &m1, const numt &m2, const Direction<size, numt> &dir)
{
    const auto D = binaryDecay(IM,m1,m2);
    return std::make_pair(lorentzVector(D.first.E(),dir*D.first.P().x()),lorentzVector(D.second.E(),dir*D.second.P().x()));
}

};
#endif
