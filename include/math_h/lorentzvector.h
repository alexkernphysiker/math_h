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
template<class numt = double, class Space = Vector<3, numt>>
class LorentzVector
{
public:
    typedef Space SpaceVectorType;
    typedef numt TimeCoordinateType;
private:
    TimeCoordinateType m_time;
    SpaceVectorType m_space;
public:
    virtual ~LorentzVector() {}
    inline LorentzVector(const numt &t, const Space &S): m_time(t), m_space(S) {}
    inline const TimeCoordinateType &E()const{return m_time;}
    inline const SpaceVectorType &P()const{return m_space;}
    template<class numt2, class Space2>
    inline LorentzVector(const LorentzVector<numt2, Space2> &source): m_time(source.E()), m_space(source.P()) {}
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
    numt operator*(const LorentzVector &second)const
    {
        return (m_time * second.m_time) - (m_space * second.m_space);
    }
    inline numt M_sqr()const
    {
        return operator*(*this);
    }
    inline numt M()const
    {
        return sqrt(M_sqr());
    }
    inline numt Ekin()const
    {
	return E()-M();
    }

    inline static LorentzVector zero()
    {
        return LorentzVector(numt(0), Space::zero());
    }
    template<class...Args>
    inline LorentzVector Rotate(Args...args)const
    {
        return LorentzVector(E(), Rotation(args...) * P());
    }
    LorentzVector Transform(const Space &Beta)const
    {
        const numt beta = Beta.length();
        if (beta == 0.0)return *this;
        if (beta >= numt(1))throw Exception<LorentzVector>("Bad Lorentz transformation");
        const auto bn = Beta / beta;
        const numt gamma = numt(1) / sqrt(numt(1) - beta * beta);
        const auto ST = ONE<Space::Dimensions,typename Space::NumberType>() + TensorProduct(bn, bn) * (gamma - numt(1));
        const auto TT = -Beta * gamma;
        return LorentzVector((E() * gamma) + (P() * TT), (ST * P()) + (TT * E()));
    }
    Space Beta()const
    {
        return P() / E();
    }
};
template<class numt = double, class Space = Vector<3, numt>>
inline LorentzVector<numt, Space> lorentzVector(const numt &t, const Space &s)
{
    return LorentzVector<numt, Space>(t, s);
}
template<size_t dim,class numt = double>
inline LorentzVector<numt, Vector<dim,numt>> lorentz_rest(const numt &l4)
{
    return LorentzVector<numt,Vector<dim,numt>>(l4 * l4,Vector<dim,numt>::zero());
}
template<class numt = double>
inline LorentzVector<numt, Vector<3,numt>> lorentz_Rest(const numt &l4)
{
    return lorentz_rest<3,numt>(l4);
}
template<class numt = double, class Space = Vector<3, numt>>
inline LorentzVector<numt, Space> lorentz_byPM(const Space &s, const numt &l4)
{
    return LorentzVector<numt, Space>(sqrt(s.length_sqr() + l4 * l4), s);
}
template<class numt = double, class Space = Vector<3, numt>>
inline LorentzVector<numt, Space> lorentz_byEM(const numt &t, const numt &l4, const typename Space::DType &dir)
{
    numt Sp = sqrt(t * t - l4 * l4);
    return LorentzVector<numt, Space>(t, dir * Sp);
}
template<class numt = double, class Space = Vector<3, numt>>
inline LorentzVector<numt, Space> lorentz_byEM(const numt &t, const numt &l4, const Space &Dir)
{
    return lorentz_byEM(t, l4, direction(Dir));
}
template<class numt = double, class Space = Vector<3, numt>>
inline LorentzVector<numt, Space> lorentz_byEkM(const numt &e, const numt &l4, const typename Space::DType &dir)
{
    return lorentz_byEM(e+l4, l4, dir);
}
template<class numt = double, class Space = Vector<3, numt>>
inline LorentzVector<numt, Space> lorentz_byEkM(const numt &e, const numt &l4, const Space &Dir)
{
    return lorentz_byEM(e+l4, l4, direction(Dir));
}
template<class numt>
std::pair<LorentzVector<numt, Vector<1, numt>>, LorentzVector<numt, Vector<1, numt>>>
        binaryDecay(const numt &IM, const numt &m1, const numt &m2)
{
    if (m1 < 0)throw Exception<std::pair<LorentzVector<numt, Vector<1, numt>>, LorentzVector<numt, Vector<1, numt>>>>("Negative mass1 error");
    if (m2 < 0)throw Exception<std::pair<LorentzVector<numt, Vector<1, numt>>, LorentzVector<numt, Vector<1, numt>>>>("Negative mass2 error");
    if (IM < (m1 + m2))throw Exception<std::pair<LorentzVector<numt, Vector<1, numt>>, LorentzVector<numt, Vector<1, numt>>>>("Invariant mass of decaying system is less then masses of products");
    const auto E1 = (IM * IM - m2 * m2 + m1 * m1) / (IM * numt(2));
    const auto p = sqrt(E1 * E1 - m1 * m1);
    return std::make_pair(lorentz_byPM(vec(p), m1), lorentz_byPM(-vec(p), m2));
}
template<size_t size, class numt>
inline std::pair<LorentzVector<numt, Vector<size, numt>>, LorentzVector<numt, Vector<size, numt>>>
        binaryDecay(const numt &IM, const numt &m1, const numt &m2, const Direction<size, numt> &dir)
{
    const auto D = binaryDecay(IM,m1,m2);
    return std::make_pair(lorentzVector(D.first.E(),dir*D.first.P().x()),lorentzVector(D.second.E(),dir*D.second.P().x()));
}

};
#endif
