// this file is distributed under
// MIT license
#ifndef ___________VECTORS_H_____
#	define ___________VECTORS_H_____
#include <math.h>
#include "error.h"
#include "chains.h"
#include "functions.h"
#include "randomfunc.h"
namespace MathTemplates
{
template<class numt, class... Args>
const numt &MakeException()
{
    throw Exception<numt>("Vector components indexing error");
    static numt x = 0;
    return x;
}
template<size_t size = 3, class numt = double>class Vector;
template<class numt>
class Vector<1, numt>
{
public:
    enum {Dimensions = 1};
    typedef numt NumberType;
private:
    numt m_x;
public:
    virtual ~Vector() {}
    Vector(const numt &x): m_x(x) {}
    Vector(const Chain<numt> &v, const size_t i): m_x(v[i]) {}
    Vector(const Chain<numt> &v): Vector(v,0) {}
    static inline const Vector zero()
    {
        return Vector(numt(0));
    }
    template<size_t index>
    static inline const Vector basis_vector()
    {
        return (index == 1) ? numt(1) : MakeException<numt>();
    }
    template<size_t index>
    const numt &component()const
    {
        return (index == 1) ? m_x : MakeException<numt>();
    }
    const numt &x()const
    {
        return component<1>();
    }
    const numt mag_sqr()const
    {
        return (m_x * m_x);
    }
    inline const numt mag()const
    {
        return sqrt(mag_sqr());
    }
    Vector(const Vector &source): m_x(source.m_x) {}
    Vector &operator=(const Vector &source)
    {
        m_x = source.m_x;
        return *this;
    }
    Vector &operator+=(const Vector &second)
    {
        m_x += second.m_x;
        return *this;
    }
    const Vector operator+(const Vector &second)const
    {
        return Vector(m_x + second.m_x);
    }
    Vector &operator-=(const Vector &second)
    {
        m_x -= second.m_x;
        return *this;
    }
    const Vector operator-(const Vector &second)const
    {
        return Vector(m_x - second.m_x);
    }

    Vector &operator*=(const numt &second)
    {
        m_x *= second;
        return *this;
    }
    const Vector operator*(const numt &second)const
    {
        return Vector(m_x * second);
    }
    Vector &operator/=(const numt &second)
    {
        m_x /= second;
        return *this;
    }
    const Vector operator/(const numt &second)const
    {
        return Vector(m_x / second);
    }
    const numt operator*(const Vector &second)const
    {
        return (m_x * second.m_x);
    }
    const bool operator==(const Vector &second)const
    {
        return (m_x == second.m_x);
    }
    const Chain<numt> to_vector()const
    {
        return {m_x};
    }
};
template<size_t size, class numt>
class Vector
{
public:
    enum {Dimensions = size};
    typedef numt NumberType;
    typedef Vector < size - 1, numt > VectorN;
private:
    numt m_x;
    VectorN m_other;
    Vector(const numt &x, const VectorN &other): m_x(x), m_other(other) {}
public:
    const VectorN &__last_coordinates()const
    {
        return m_other;
    }
    virtual ~Vector() {}
    template<class... Args>
    Vector(const numt &x, Args... other): m_x(x), m_other(other...) {}
    Vector(const Chain<numt> &v, const size_t i): m_x(v[i]), m_other(v, i + 1) {}
    Vector(const Chain<numt> &v): Vector(v,0) {}
    static inline const Vector zero()
    {
        return Vector(numt(0), Vector < size - 1 >::zero());
    }
    template<size_t index>
    static inline const Vector basis_vector()
    {
        return (index == 1) ? Vector(numt(1), VectorN::zero()) : Vector(numt(0), VectorN::template basis_vector < index - 1 > ());
    }
    template<size_t index>
    const numt &component()const
    {
        return (index == 1) ? m_x : m_other.template component < index - 1 > ();
    }
    const numt &x()const
    {
        return component<1>();
    }
    const numt &y()const
    {
        return component<2>();
    }
    const numt &z()const
    {
        return component<3>();
    }
    const numt mag_sqr()const
    {
        return (m_x * m_x) + m_other.mag_sqr();
    }
    inline const numt mag()const
    {
        return sqrt(mag_sqr());
    }
    Vector(const Vector &source): m_x(source.m_x), m_other(source.m_other) {}
    Vector &operator=(const Vector &source)
    {
        m_x = source.m_x;
        m_other = source.m_other;
        return *this;
    }
    Vector &operator+=(const Vector &second)
    {
        m_x += second.m_x;
        m_other += second.m_other;
        return *this;
    }
    const Vector operator+(const Vector &second)const
    {
        return Vector(m_x + second.m_x, m_other + second.m_other);
    }
    Vector &operator-=(const Vector &second)
    {
        m_x -= second.m_x;
        m_other -= second.m_other;
        return *this;
    }
    const Vector operator-(const Vector &second)const
    {
        return Vector(m_x - second.m_x, m_other - second.m_other);
    }

    Vector &operator*=(const numt &second)
    {
        m_x *= second;
        m_other *= second.m_other;
        return *this;
    }
    const Vector operator*(const numt &second)const
    {
        return Vector(m_x * second, m_other * second);
    }
    Vector &operator/=(const numt &second)
    {
        m_x /= second;
        m_other /= second.m_other;
        return *this;
    }
    const Vector operator/(const numt &second)const
    {
        return Vector(m_x / second, m_other / second);
    }
    const numt operator*(const Vector &second)const
    {
        return (m_x * second.m_x) + (m_other * second.m_other);
    }
    const bool operator==(const Vector &second)const
    {
        return (m_x == second.m_x) && (m_other == second.m_other);
    }
    const Chain<numt> to_vector()const
    {
        Chain<numt> res = m_other.to_vector();
        res.insert(res.begin(), m_x);
        return res;
    }
};
template<class numt = double, class... Args>
const Vector < 1 + sizeof...(Args), numt > DesCartes(const numt &x, Args... other)
{
    return Vector < 1 + sizeof...(Args), numt > (x, other...);
}
template<class numt = double>
const Vector<2, numt> I()
{
    return DesCartes(numt(1), numt(0));
}
template<class numt = double>
const Vector<2, numt> J()
{
    return DesCartes(numt(0), numt(1));
}
template<class numt = double>
const Vector<2, numt> zero()
{
    return DesCartes(numt(0), numt(0));
}

template<class numt = double>
const Vector<3, numt> X()
{
    return DesCartes(numt(1), numt(0), numt(0));
}
template<class numt = double>
const Vector<3, numt> Y()
{
    return DesCartes(numt(0), numt(1), numt(0));
}
template<class numt = double>
const Vector<3, numt> Z()
{
    return DesCartes(numt(0), numt(0), numt(1));
}
template<class numt = double>
const Vector<3, numt> Zero()
{
    return DesCartes(numt(0), numt(0), numt(0));
}

template<size_t i, class numt = double>
inline const Vector<i, numt> operator-(const Vector<i, numt> &V)
{
    return V * numt(-1);
}
template<size_t i, class numt = double>
std::ostream &operator<<(std::ostream &str, const Vector<i, numt> &V)
{
    return str << V.x() << " " << V.__last_coordinates();
}
template<class numt = double>
std::ostream &operator<<(std::ostream &str, const Vector<1, numt> &V)
{
    return str << V.x() << " ";
}

template<size_t i,class numt = double>
static const Vector<i,numt> Direction(const Vector<i,numt>&v)
{
    return v/v.mag();
}
template<class numt = double>
static const Vector<1, numt> Direction()
{
    return Vector<1, numt>(1);
}
template<class numt = double>
static const Vector<2, numt> Direction(const numt &phi)
{
    return Vector<2, numt>(cos(phi), sin(phi));
}
template<class numt = double>
static const Vector<3, numt> Direction(const numt &theta, const numt &phi)
{
    return Vector<3, numt>(cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta));
}
template<class numt = double, class... Args>
static inline const Vector < 1 + sizeof...(Args), numt > PolarCoordinates(const numt &x, Args... other)
{
    return Direction(other...) * x;
}

template<class numt = double>
const Vector<3, numt> operator^(const Vector<3, numt> &first, const Vector<3, numt> &second)
{
    return Vector<3, numt>(
               (first.y() * second.z()) - (second.y() * first.z()),
               (first.z() * second.x()) - (second.z() * first.x()),
               (first.x() * second.y()) - (second.x() * first.y())
           );
}
template<class numt = double>
const Vector<3, numt> Rotate(const Vector<3, numt> &src, const Vector<3, numt> &axis, const numt &theta)
{
    const auto n = axis / axis.mag();
    const Chain<numt> source = {src.x(), src.y(), src.z()};
    const numt cost = cos(theta), sint = sin(theta), one = 1;
    Chain<Chain<numt>> M = {
        {    cost + (one - cost) *n.x() *n.x(),	(one - cost) *n.x() *n.y() - sint * n.z(),	(one - cost) *n.x() *n.z() + sint * n.y()},
        {(one - cost) *n.y() *n.x() + sint * n.z(),	    cost + (one - cost) *n.y() *n.y(),	(one - cost) *n.y() *n.z() - sint * n.x()},
        {(one - cost) *n.z() *n.x() - sint * n.y(),	(one - cost) *n.z() *n.y() + sint * n.x(),	    cost + (one - cost) *n.z() *n.z()}
    };
    Chain<numt> dest;
    for (size_t i = 0; i < 3; i++) {
        numt v = 0;
        for (size_t j = 0; j < 3; j++)v += M[i][j] * source[j];
        dest.push_back(v);
    }
    return DesCartes(dest[0], dest[1], dest[2]);
}
template<class numt = double, class RG = RANDOM>
static const Vector<3, numt> RandomIsotropicDirection3(RG &generator)
{
    static RandomUniform<numt> Z(-1, 1);
    static RandomUniform<numt> Phi(0, PI<numt>() * 2.0);
    const numt z = Z(generator);
    const numt x_y = sqrt(1.0 - pow(z, 2));
    const numt phi = Phi(generator);
    return DesCartes(cos(phi) * x_y, sin(phi) * x_y, z);
}

template<class numt = double>
const numt operator^(const Vector<2, numt> &first, const Vector<2, numt> &second)
{
    return (first.x() * second.y()) - (second.x() * first.y());
}
template<class numt = double>
const Vector<2, numt> Rotate(const Vector<2, numt> &v, const numt &theta)
{
    const Chain<numt> source = {v.x(), v.y()};
    const numt cost = cos(theta), sint = sin(theta);
    Chain<Chain<numt>> M = {
        {cost, -sint},
        {sint, cost}
    };
    Chain<numt> dest;
    for (size_t i = 0; i < 2; i++) {
        numt v = 0;
        for (size_t j = 0; j < 2; j++)v += M[i][j] * source[j];
        dest.push_back(v);
    }
    return Vector<2, numt>(dest[0], dest[1]);
}
template<class numt = double, class RG = RANDOM>
const Vector<2, numt> RandomIsotropicDirection2(RG &generator)
{
    static RandomUniform<numt> Phi(0, PI<numt>() * 2.0);
    const numt phi = Phi(generator);
    return Vector<2, numt>(cos(phi), sin(phi));
}



template<class numt = double, class Space = Vector<3, numt>>
class LorentzVector
{
private:
    numt m_time;
    Space m_space;
public:
    virtual ~LorentzVector() {}
    LorentzVector(const numt &t, const Space &S): m_time(t), m_space(S) {}
    typedef Space SpaceVectorType;
    typedef numt TimeCoordinateType;
    LorentzVector(const LorentzVector &source): LorentzVector(source.m_time, source.m_space) {}
    LorentzVector &operator=(const LorentzVector &source)
    {
        m_space = source.m_space;
        m_time = source.m_time;
        return *this;
    }
    const numt &time_component()const
    {
        return m_time;
    }
    const Space &space_component()const
    {
        return m_space;
    }
    const numt length4_sqr()const
    {
        return (m_time * m_time) - m_space.mag_sqr();
    }
    const numt length4()const
    {
        return sqrt(length4_sqr());
    }
    const bool operator==(const LorentzVector &second)const
    {
        return (m_time == second.m_time) && (m_space == second.m_space);
    }
    LorentzVector &operator+=(const LorentzVector &second)
    {
        m_time += second.m_time;
        m_space += second.m_space;
        return *this;
    }
    const LorentzVector operator+(const LorentzVector &second)const
    {
        return LorentzVector(m_time + second.m_time, m_space + second.m_space);
    }
    LorentzVector &operator-=(const LorentzVector &second)
    {
        m_time -= second.m_time;
        m_space -= second.m_space;
        return *this;
    }
    const LorentzVector operator-(const LorentzVector &second)const
    {
        return LorentzVector(m_time - second.m_time, m_space - second.m_space);
    }
    const numt operator*(const LorentzVector &second)const
    {
        return (m_time * second.m_time) - (m_space * second.m_space);
    }
    static const LorentzVector zero()
    {
        return LorentzVector(numt(0), Space::zero());
    }
    LorentzVector(const numt &t): LorentzVector(t, Space::zero()) {}
    template<class...Args>
    const LorentzVector Rotate(Args...args)const
    {
        return LorentzVector(time_component(), MathTemplates::Rotate(space_component(), args...));
    }
    const LorentzVector Lorentz(const Space &Beta)const
    {
        const numt beta = Beta.mag();
        if (beta == 0.0)return *this;
        if (beta >= numt(1))throw Exception<LorentzVector>("Bad Lorentz transformation");
        const auto n = (Beta / beta).to_vector();
        const numt gamma = numt(1) / sqrt(numt(1) - beta * beta);
        auto source = space_component().to_vector();
        source.insert(source.begin(), time_component());
        Chain<Chain<numt>> M;
        for (size_t i = 0; i < source.size(); i++) {
            Chain<numt> V;
            for (size_t j = 0; j < source.size(); j++) {
                if ((0 == i) && (0 == j)) {
                    V.push_back(gamma);
                } else {
                    if (0 == i) {
                        V.push_back(-gamma * beta * n[j - 1]);
                    } else {
                        if (0 == j) {
                            V.push_back(-gamma * beta * n[i - 1]);
                        } else {
                            if (i == j) {
                                V.push_back(numt(1) + (gamma - numt(1))*n[i - 1]*n[j - 1]);
                            } else {
                                V.push_back((gamma - numt(1))*n[i - 1]*n[j - 1]);
                            }
                        }
                    }
                }
            }
            M.push_back(V);
        }
        Chain<numt> dest;
        for (size_t i = 0; i < source.size(); i++) {
            numt v = 0;
            for (size_t j = 0; j < source.size(); j++)v += M[i][j] * source[j];
            dest.push_back(v);
        }
        const auto T = dest[0];
        dest.erase(dest.begin());
        return LorentzVector(T, dest);
    }
    const Space Beta()const //Makes physical sense only if it's a lorentz-momentum
    {
        return space_component() / time_component();
    }
};
template<class numt = double, class Space = Vector<3, numt>>
const LorentzVector<numt, Space> lorentzVector(const numt &t, const Space &s)
{
    return LorentzVector<numt, Space>(t, s);
}
template<class numt = double, class Space = Vector<3, numt>>
const LorentzVector<numt, Space> lorentz_byPM(const Space &s, const numt &l4)
{
    return LorentzVector<numt, Space>(sqrt(s.mag_sqr() + l4 * l4), s);
}
template<class numt = double, class Space = Vector<3, numt>, class...Args>
const LorentzVector<numt, Space> lorentz_byEM(const numt &t, const numt &l4, Args...args)
{
    numt Sp = sqrt(t * t - l4 * l4);
    return LorentzVector<numt, Space>(t, Direction(args...) * Sp);
}
template<class numt = double, class...Args>
const auto binaryDecay(const numt &IM, const numt &m1, const numt &m2, Args...args)->decltype(std::make_pair(lorentz_byPM(Direction(args...),IM),lorentz_byPM(Direction(args...),IM)))
{
    if (m1 < 0)throw Exception<std::pair<LorentzVector<numt, Vector<1, numt>>, LorentzVector<numt, Vector<1, numt>>>>("Negative mass error");
    if (m2 < 0)throw Exception<std::pair<LorentzVector<numt, Vector<1, numt>>, LorentzVector<numt, Vector<1, numt>>>>("Negative mass error");
    if (IM < (m1 + m2))throw Exception<std::pair<LorentzVector<numt, Vector<1, numt>>, LorentzVector<numt, Vector<1, numt>>>>("Invariant mass of decaying system is less then masses of products");
    const auto E1=(IM*IM-m2*m2+m1*m1)/(IM*numt(2));
    const auto p=sqrt(E1*E1-m1*m1);
    const auto P=Direction(args...)*p;
    return std::make_pair(lorentz_byPM(P,m1),lorentz_byPM(-P,m2));
}

template<class numt = double>
class Plane3D
{
public:
    typedef Vector<3,numt> Space;
    typedef Vector<2,numt> Plane;
private:
    Space basis_x, basis_y;
public:
    const Space operator()(const Plane &v)const
    {
        return basis_x * v.x() + basis_y * v.y();
    }
    virtual ~Plane3D() {}
    Plane3D(const Space &x, const Space &y): basis_x(x), basis_y(y)
    {
        if ((x ^ y).mag_sqr() == numt(0))throw Exception<Plane3D>("Invalid basis vectors for converting 2d vectors to 3d");
    }
    template<size_t i, size_t j>
    static const Plane3D basis()
    {
        return Plane3D(Space::template basis_vector<i>(), Space::template basis_vector<j>());
    }
    static const Plane3D ByNormalVectorAndTheta(const Space &N, const numt &theta)
    {
        if (N.mag() == 0)throw Exception<Plane3D>("Invalid normale vector");
        const auto n = N / N.mag();
        if ((n ^ Space::template basis_vector<3>()).mag_sqr() == 0) {
            const auto X = Rotate(Space::template basis_vector<1>(), n, theta);
            return Plane3D(X, n ^ X);
        } else {
            const auto X = (n ^ Space::template basis_vector<3>());
            return Plane3D(X / X.mag(), (n ^ (X / X.mag())));
        }
        throw Exception<Plane3D>("This line should not be reached");
    }
};

};
#endif
