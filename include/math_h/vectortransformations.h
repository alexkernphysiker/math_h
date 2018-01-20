// this file is distributed under
// LGPLv3 license
#ifndef ___________VECTORTRANSFORMATIONS_H_____
#	define ___________VECTORTRANSFORMATIONS_H_____
#if __cplusplus<201100L
#error c++>=11 is needed for using math_h headers
#endif
#include <tuple>
#include <math.h>
#include <type_traits>
#include "error.h"
#include "vectors.h"
#include "matrices.h"
#if __cplusplus>201700L
#define ____optimized_version_of_vectortransformations_h_____
#else
#warning compiler does not support "if constexpr(...)". c++>=17 is needed. classes from vectortransformations.h will work slower
#endif
namespace MathTemplates
{

template<class numt>
class Direction<1, numt>
{
    template<size_t sizef, class n>friend class Vector;
    template<size_t sizef, class n>friend class Direction;
    template<size_t sizef, class n>friend class Matrix;
private:
    bool sign;
    static numt PHI(RANDOM &r)
    {
        static const RandomUniform<numt> res(-PI<numt>(), +PI<numt>());
        return res(r);
    }
protected:
    inline const bool &___sign()const
    {
        return sign;
    }
public:
    enum {Dimensions = 1};
    enum {Thetas = 0};
    typedef numt NumberType;
    typedef Vector<1, numt> VType;
    virtual ~Direction() {}
    inline Direction(RANDOM &RG): sign(PHI(RG) >= 0) {}
    template<class numt2>
    inline Direction(const Direction<Dimensions, numt2> &source): sign(source.___sign()) {}
    inline Direction(const VType &v): sign(v.x() >= 0) {}
    template<class... Args>
    inline Direction(const std::tuple<Args...> &args): sign(std::get<0>(args)) {}
    inline VType operator*(const numt &rho)const
    {
        if (sign)return desCartes(rho);
        else return desCartes(-rho);
    }
    inline bool operator==(const Direction &second)const
    {
        return (sign == second.sign);
    }
    inline numt dir()const
    {
        return sign ? numt(1) : numt(-1);
    }
    inline Matrix<Dimensions, VType> Rotations()const
    {
        return sign ? lines(desCartes(numt(1))) : lines(desCartes(numt(-1)));
    }
    inline Matrix<Dimensions, VType> AntiRotations()const
    {
        return sign ? lines(desCartes(numt(1))) : lines(desCartes(numt(-1)));
    }
};

template<class numt>
class Direction<2, numt>
{
    template<size_t sizef, class n>friend class Vector;
    template<size_t sizef, class n>friend class Direction;
    template<size_t sizef, class n>friend class Matrix;
public:
    enum {Dimensions = 2};
    enum {Thetas = 0};
    typedef numt NumberType;
    typedef Vector<2, numt> VType;
private:
    numt m_phi;
    void ___phi_norm()
    {
        while (m_phi > PI<numt>())m_phi -= PI<numt>() * numt(2);
        while (m_phi < -PI<numt>())m_phi += PI<numt>() * numt(2);
    }
    static numt PHI(RANDOM &r)
    {
        static const RandomUniform<numt> res(-PI<numt>(), +PI<numt>());
        return res(r);
    }
    inline const numt &___last_component()const
    {
        return m_phi;
    }
public:
    virtual ~Direction() {}
    template<class... Args>
    inline Direction(const std::tuple<Args...> &v): m_phi(std::get<0>(v))
    {
        ___phi_norm();
    }
    inline Direction(RANDOM &RG): m_phi(PHI(RG)) {}
    template<class numt2>
    inline Direction(const Direction<Dimensions, numt2> &source): m_phi(source.___last_component())
    {
        ___phi_norm();
    }
    inline const numt &phi()const
    {
        return m_phi;
    }
    Direction(const VType &V): m_phi(atan2(V.y(), V.x())) {}
    VType operator*(const numt &rho)const
    {
        return desCartes(cos(m_phi), sin(m_phi)) * rho;
    }
    bool operator==(const Direction &second)const
    {
        return (m_phi == second.m.phi);
    }
    Matrix<Dimensions, VType> Rotations()const
    {
        return Matrix<Dimensions, VType>::template RotationInPlane<1, 2>(m_phi);
    }
    Matrix<Dimensions, VType> AntiRotations()const
    {
        return Matrix<Dimensions, VType>::template RotationInPlane<1, 2>(-m_phi);
    }
};
template<class numt>
class Direction<3, numt>
{
    template<size_t sizef, class n>friend class Vector;
    template<size_t sizef, class n>friend class Direction;
    template<size_t sizef, class n>friend class Matrix;
public:
    enum {Dimensions = 3};
    enum {Thetas = 1};
    typedef numt NumberType;
    typedef Vector<3, numt> VType;
    typedef Direction<2, numt> DirectionN;
private:
    DirectionN m_ld;
    numt m_theta;
    void ___theta_check()const
    {
        if (m_theta < numt(0))throw Exception<Direction>("wrong theta value <0");
        if (m_theta > PI<numt>())throw Exception<Direction>("wrong theta value >pi");
    }
    static numt CTHETA(RANDOM &r)
    {
        static const RandomUniform<numt> res(numt(-1), numt(+1));
        return res(r);
    }
protected:
    inline const DirectionN &___recursive()const
    {
        return m_ld;
    }
    inline const numt &___last_component()const
    {
        return m_theta;
    }
public:
    virtual ~Direction() {}
    template<class numt2>
    inline Direction(const Direction<Dimensions, numt2> &source): m_ld(source.___recursive()), m_theta(source.___last_component())
    {
        ___theta_check();
    }
    template<class... Args>
    inline Direction(const std::tuple<Args...> &v): m_ld(v), m_theta(std::get<1>(v))
    {
        ___theta_check();
    }
    inline Direction(RANDOM &RG): m_ld(RG), m_theta(acos(CTHETA(RG))) {}
    inline const numt &phi()const
    {
        return m_ld.phi();
    }
#ifdef ____optimized_version_of_vectortransformations_h_____
    template<size_t index = 1>
    inline const numt &th()const
    {
	static_assert(index ==1,"dimension index is out of range");
        return m_theta;
    }
#else
    template<size_t index = 1>
    const numt &th()const
    {
	static_assert(index >0,"dimension index is out of range");
	if(index>1)throw Exception<Direction>("dimension index is out of range");
        return m_theta;
    }
#endif
    Direction(const VType &V): m_ld(V.___recursive()), m_theta(acos(V.template component<Dimensions>() / V.M())) {}
    VType operator*(const numt &rho)const
    {
        return VType(m_ld * (rho * sin(m_theta)), rho * cos(m_theta));
    }
    bool operator==(const Direction &second)const
    {
        return (m_theta == second.m_theta) && (m_ld == second.m_ld);
    }
    Matrix<Dimensions, VType> Rotations()const
    {
        return
            Matrix<Dimensions, VType>(
                m_ld.Rotations().___add_column(DirectionN::VType::zero()),
                VType::template basis_vector<Dimensions>())
            * Matrix<Dimensions, VType>::template RotationInPlane<3, 1>(m_theta);
    }
    Matrix<Dimensions, VType> AntiRotations()const
    {
        return
            Matrix<Dimensions, VType>::template RotationInPlane<3, 1>(-m_theta) *
        Matrix<Dimensions, VType>(
            m_ld.AntiRotations().___add_column(DirectionN::VType::zero()),
            VType::template basis_vector<Dimensions>());
    }
};
template<size_t size, class numt>
class Direction
{
    template<size_t sizef, class n>friend class Vector;
    template<size_t sizef, class n>friend class Direction;
    template<size_t sizef, class n>friend class Matrix;
public:
    enum {Dimensions = size};
    enum {Thetas = size - 2};
    typedef numt NumberType;
    typedef Vector<size, numt> VType;
    typedef Direction < size - 1, numt > DirectionN;
private:
    DirectionN m_ld;
    numt m_theta;
    void ___theta_check()const
    {
        if (m_theta < numt(0))throw Exception<Direction>("wrong theta value <0");
        if (m_theta > PI<numt>())throw Exception<Direction>("wrong theta value >pi");
    }
    static numt CTHETA(RANDOM &r)
    {
        static const RandomUniform<numt> res(numt(-1), numt(+1));
        return res(r);
    }
protected:
    inline const DirectionN &___recursive()const
    {
        return m_ld;
    }
    inline const numt &___last_component()const
    {
        return m_theta;
    }
public:
    virtual ~Direction() {}
    template<class numt2>
    inline Direction(const Direction<Dimensions, numt2> &source): m_ld(source.___recursive()), m_theta(source.___last_component())
    {
        ___theta_check();
    }
    template<class... Args>
    inline Direction(const std::tuple<Args...> &v): m_ld(v), m_theta(std::get<Thetas>(v))
    {
        ___theta_check();
    }
    inline Direction(RANDOM &RG): m_ld(RG), m_theta(acos(CTHETA(RG))) {}
    inline const numt &phi()const
    {
        return m_ld.phi();
    }
#ifdef ____optimized_version_of_vectortransformations_h_____
    template<size_t index>
    inline const numt &th()const
    {
	static_assert(index > 0,"dimension index is out of range");
	static_assert(index <= Thetas,"dimension index is out of range");
        if constexpr(index == Thetas) return m_theta;
	else return m_ld.template th<index>();
    }
#else
    template<size_t index>
    const numt &th()const
    {
	static_assert(index > 0,"dimension index is out of range");
	if(index > Thetas)throw Exception<Direction>("dimension index is out of range");
        if(index == Thetas) return m_theta;
	else return m_ld.template th<index>();
    }
#endif
    Direction(const VType &V): m_ld(V.___recursive()), m_theta(acos(V.template component<Dimensions>() / V.M())) {}
    VType operator*(const numt &rho)const
    {
        return VType(m_ld * (rho * sin(m_theta)), rho * cos(m_theta));
    }
    Matrix<Dimensions, VType> Rotations()const
    {
        return
            Matrix<Dimensions, VType>(
                m_ld.Rotations().___add_column(DirectionN::VType::zero()),
                VType::template basis_vector<Dimensions>())
            * Matrix<Dimensions, VType>::template RotationInPlane < Dimensions, Dimensions - 1 > (m_theta);
    }
    Matrix<Dimensions, VType> AntiRotations()const
    {
        return
            Matrix<Dimensions, VType>::template RotationInPlane < Dimensions, Dimensions - 1 > (-m_theta) *
        Matrix<Dimensions, VType>(
            m_ld.AntiRotations().___add_column(DirectionN::VType::zero()),
            VType::template basis_vector<Dimensions>());
    }
};
template<class numt = double>
inline Direction< 1, numt > direction()
{
    return Direction < 1, numt > (std::make_tuple(true));
}
template<class numt = double, class... Args>
inline Direction < 2 + sizeof...(Args), numt > direction(const numt &phi, Args... other)
{
    return Direction < 2 + sizeof...(Args), numt > (std::make_tuple(phi, other...));
}
template<size_t size, class numt = double>
inline Direction< size, numt > direction(const Vector<size, numt> &V)
{
    return Direction < size, numt > (V);
}
template<size_t size = 3, class numt = double, class RG = RANDOM>
inline Direction<size, numt> randomIsotropic(RG &generator)
{
    return Direction<size, numt>(generator);
}
template<class VType>
struct VectorDecomposition {
    VType tau;
    VType n;
};
template<class VType>
inline VectorDecomposition<VType> decompose_by_direction(const VType &source, const typename VType::DType &dir)
{
    typedef typename VType::NumberType numt;
    const auto t = dir * ((dir * numt(1)) * source);
    return {.tau = t, .n = source - t};
}
template<class VType>
inline VectorDecomposition<VType> decompose_by_plane_normale(const VType &source, const typename VType::DType &pn)
{
    typedef typename VType::NumberType numt;
    const auto t = pn * ((pn * numt(1)) * source);
    return {.tau = source - t, .n = t};
}
template<class numt = double>
inline Matrix<2, Vector<2, numt>> Rotation(const numt &theta)
{
    const numt cost = cos(theta), sint = sin(theta);
    return lines(
               desCartes(cost, -sint),
               desCartes(sint, cost)
           );
}
template<class numt = double>
inline Matrix<3, Vector<3, numt>> Rotation(const Direction<3, numt> &axis, const numt &theta)
{
    const auto n = axis * numt(1);
    const numt cost = cos(theta), sint = sin(theta), one = 1;
    return lines(
               desCartes(cost + (one - cost) * n.x() * n.x(),	(one - cost) * n.x() * n.y() - sint * n.z(),	(one - cost) * n.x() * n.z() + sint * n.y()),
               desCartes((one - cost) * n.y() * n.x() + sint * n.z(),	    cost + (one - cost) * n.y() * n.y(),	(one - cost) * n.y() * n.z() - sint * n.x()),
               desCartes((one - cost) * n.z() * n.x() - sint * n.y(),	(one - cost) * n.z() * n.y() + sint * n.x(),	    cost + (one - cost) * n.z() * n.z())
           );
}


#define AC(n) (A.template component<n>())
#define ZeRo numt(0)
template<class numt = double>
inline Matrix<1, Vector<2, numt>> SkewM(const Vector<2, numt> &A)
{
    return lines(desCartes(-AC(2), AC(1)));
}
template<class numt = double>
inline Matrix<3, Vector<3, numt>> SkewM(const Vector<3, numt> &A)
{
    return lines(
               desCartes(ZeRo, -AC(3), AC(2)),
               desCartes(AC(3),  ZeRo, -AC(1)),
               desCartes(-AC(2), AC(1),  ZeRo)
           );
}
template<class numt = double>
inline Matrix<7, Vector<7, numt>> SkewM(const Vector<7, numt> &A)
{
    return lines(
               desCartes(ZeRo, -AC(4), -AC(7), AC(2), -AC(6), AC(5), AC(3)),
               desCartes(AC(4),  ZeRo, -AC(5), -AC(1), AC(3), -AC(7), AC(6)),
               desCartes(AC(7), AC(5),  ZeRo, -AC(6), -AC(2), AC(4), -AC(1)),
               desCartes(-AC(2), AC(1), AC(6),  ZeRo, -AC(7), -AC(3), AC(5)),
               desCartes(AC(6), -AC(3), AC(2), AC(7),  ZeRo, -AC(1), -AC(4)),
               desCartes(-AC(5), AC(7), -AC(4), AC(3), AC(1),  ZeRo, -AC(2)),
               desCartes(-AC(3), -AC(6), AC(1), -AC(5), AC(4), AC(2), ZeRo)
           );
}
#undef ZeRo
#undef AC
template<class A, class B>
inline auto operator^(const A &a, const B &b)->decltype(SkewM(a)*b)
{
    return SkewM(a) * b;
}


};
#endif