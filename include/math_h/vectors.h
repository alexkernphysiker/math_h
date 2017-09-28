// this file is distributed under
// MIT license
#ifndef ___________VECTORS_H_____
#	define ___________VECTORS_H_____
#include <tuple>
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
template<size_t size = 3, class numt = double>class Direction;
template<class numt>
class Vector<1, numt>
{
public:
    enum {Dimensions = 1};
    typedef numt NumberType;
    typedef Direction<1,numt> DType;
private:
    numt m_x;
public:
    virtual ~Vector() {}
    Vector(const Vector &source): m_x(source.m_x) {}
    Vector(const Chain<numt> &v): m_x(v[0]) {}
    template<class... Args>
    Vector(const std::tuple<Args...>&v): m_x(std::get<0>(v)) {}
    const Chain<numt> chain()const
    {
        return {m_x};
    }
    static const Vector zero()
    {
        return Vector(std::make_tuple(numt(0)));
    }
    template<size_t index>
    static const Vector basis_vector()
    {
        return (index == Dimensions) ? Vector(std::make_tuple(numt(1))) : Vector(std::make_tuple(MakeException<numt>()));
    }
    template<size_t index>
    const numt &component()const
    {
        return (index == Dimensions) ? m_x : MakeException<numt>();
    }
    const numt &x()const
    {
        return m_x;
    }
    const numt M_sqr()const
    {
        return (m_x * m_x);
    }
    inline const numt M()const
    {
        return sqrt(M_sqr());
    }
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
        return Vector(std::make_tuple(m_x + second.m_x));
    }
    Vector &operator-=(const Vector &second)
    {
        m_x -= second.m_x;
        return *this;
    }
    const Vector operator-(const Vector &second)const
    {
        return Vector(std::make_tuple(m_x - second.m_x));
    }

    Vector &operator*=(const numt &second)
    {
        m_x *= second;
        return *this;
    }
    const Vector operator*(const numt &second)const
    {
        return Vector(std::make_tuple(m_x * second));
    }
    Vector &operator/=(const numt &second)
    {
        m_x /= second;
        return *this;
    }
    const Vector operator/(const numt &second)const
    {
        return Vector(std::make_tuple(m_x / second));
    }
    const numt operator*(const Vector &second)const
    {
        return m_x * second.m_x;
    }
    const bool operator==(const Vector &second)const
    {
        return (m_x == second.m_x);
    }
};
template<size_t size, class numt>
class Vector
{
    friend class Direction<size,numt>;
public:
    enum {Dimensions = size};
    typedef numt NumberType;
    typedef Vector < size - 1, numt > VectorN;
    typedef Direction<size,numt> DType;
private:
    VectorN m_other;
    numt m_x;
    Vector(const VectorN&other,const numt&x): m_other(other),m_x(x) {}
protected:
    const VectorN&___recursive()const{return m_other;}
public:
    virtual ~Vector() {}
    Vector(const Vector &source): m_other(source.m_other),m_x(source.m_x) {}
    Vector(const Chain<numt> &v): m_other(v),m_x(v[Dimensions-1]) {}
    template<class... Args>
    Vector(const std::tuple<Args...>&v): m_other(v),m_x(std::get<Dimensions-1>(v)) {}
    const Chain<numt> chain()const
    {
        Chain<numt> res = m_other.chain();
        res.push_back(m_x);
        return res;
    }
    static const Vector zero()
    {
        return Vector(VectorN::zero(),numt(0));
    }
    template<size_t index>
    static const Vector basis_vector()
    {
        return (index == Dimensions) ? Vector(VectorN::zero(),numt(1)) : Vector(VectorN::template basis_vector<index>(),numt(0));
    }
    template<size_t index>
    const numt &component()const
    {
        return (index == Dimensions) ? m_x : m_other.template component<index>();
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
    const numt M_sqr()const
    {
        return m_other.M_sqr()+(m_x * m_x);
    }
    inline const numt M()const
    {
        return sqrt(M_sqr());
    }
    Vector &operator+=(const Vector &second)
    {
        m_x += second.m_x;
        m_other += second.m_other;
        return *this;
    }
    const Vector operator+(const Vector &second)const
    {
        return Vector(m_other + second.m_other,m_x + second.m_x);
    }
    Vector &operator-=(const Vector &second)
    {
        m_x -= second.m_x;
        m_other -= second.m_other;
        return *this;
    }
    const Vector operator-(const Vector &second)const
    {
        return Vector(m_other - second.m_other,m_x - second.m_x);
    }

    Vector &operator*=(const numt &second)
    {
        m_x *= second;
        m_other *= second.m_other;
        return *this;
    }
    const Vector operator*(const numt &second)const
    {
        return Vector(m_other * second,m_x * second);
    }
    Vector &operator/=(const numt &second)
    {
        m_x /= second;
        m_other /= second.m_other;
        return *this;
    }
    const Vector operator/(const numt &second)const
    {
        return Vector(m_other / second,m_x / second);
    }
    const numt operator*(const Vector &second)const
    {
        return (m_other * second.m_other)+(m_x * second.m_x);
    }
    const bool operator==(const Vector &second)const
    {
        return (m_x == second.m_x) && (m_other == second.m_other);
    }
};
template<class numt,class... Args>
const Vector < sizeof...(Args)+1,numt> desCartes(const numt&x,Args... args)
{
    return Vector<sizeof...(Args)+1,numt>(std::make_tuple(x,args...));
}
template<class numt = double>
const Vector<2, numt> x()
{
    return desCartes(numt(1),numt(0));
}
template<class numt = double>
const Vector<2, numt> y()
{
    return desCartes(numt(0), numt(1));
}
template<class numt=double>
const Vector<2,numt> zero()
{
    return desCartes(numt(0),numt(0));
}
template<class numt = double>
const Vector<3, numt> X()
{
    return desCartes(numt(1), numt(0),numt(0));
}
template<class numt = double>
const Vector<3, numt> Y()
{
    return desCartes(numt(0), numt(1),numt(0));
}
template<class numt = double>
const Vector<3, numt> Z()
{
    return desCartes(numt(0), numt(0), numt(1));
}
template<class numt = double>
const Vector<3, numt> Zero()
{
    return desCartes(numt(0),numt(0),numt(0));
}
template<size_t i, class numt = double>
inline const Vector<i, numt> operator-(const Vector<i, numt> &V)
{
    return V * numt(-1);
}
template<class numt = double>
const Vector<3, numt> operator^(const Vector<3, numt> &first, const Vector<3, numt> &second)
{
    return desCartes(
               (first.y() * second.z()) - (second.y() * first.z()),
               (first.z() * second.x()) - (second.z() * first.x()),
               (first.x() * second.y()) - (second.x() * first.y())
           );
}
template<class numt = double>
const numt operator^(const Vector<2, numt> &first, const Vector<2, numt> &second)
{
    return (first.x() * second.y()) - (second.x() * first.y());
}

template<class numt>
class Direction<1,numt>
{
public:
    enum {Dimensions = 1};
    enum {Thetas = 0};
    typedef numt NumberType;
    typedef Vector<1,numt> VType;
public:
    virtual ~Direction() {}
    Direction(){}
    Direction(RANDOM&){}
    Direction(const Direction&){}
    Direction(const VType&){}
    template<class... Args>
    Direction(const std::tuple<Args...>&){}
    const VType operator*(const numt&rho)const
    {
	return desCartes(rho);
    }
};

template<class numt>
class Direction<2,numt>
{
public:
    enum {Dimensions = 2};
    enum {Thetas = 0};
    typedef numt NumberType;
    typedef Vector<2,numt> VType;
private:
    numt m_phi;
    static const numt PHI(RANDOM&r){
	static const RandomUniform<numt> res(-PI<numt>(),+PI<numt>());
	return res(r);
    }
public:
    virtual ~Direction() {}
    template<class... Args>
    Direction(const std::tuple<Args...>&v): m_phi(std::get<0>(v)) {
	while(m_phi>PI<numt>())m_phi-=PI<numt>()*numt(2);
	while(m_phi<-PI<numt>())m_phi+=PI<numt>()*numt(2);
    }
    Direction(RANDOM&RG):m_phi(PHI(RG)){}
    Direction(const Direction&source):m_phi(source.m_phi){}
    const numt &phi()const
    {
        return m_phi;
    }
    Direction(const VType&V):m_phi(atan2(V.y(),V.x())){}
    const VType operator*(const numt&rho)const{
	return desCartes(cos(m_phi),sin(m_phi))*rho;
    }

};
template<class numt>
class Direction<3,numt>
{
public:
    enum {Dimensions = 3};
    enum {Thetas = 1};
    typedef numt NumberType;
    typedef Vector<3,numt> VType;
    typedef Direction<2,numt> DirectionN;
private:
    DirectionN m_ld;
    numt m_theta;
    static const numt CTHETA(RANDOM&r){
	static const RandomUniform<numt> res(numt(-1),numt(+1));
	return res(r);
    }
public:
    virtual ~Direction() {}
    template<class... Args>
    Direction(const std::tuple<Args...>&v):m_ld(v),m_theta(std::get<1>(v)) {
	if(m_theta<numt(0))throw Exception<Direction>("wrong theta value <0");
	if(m_theta>PI<numt>())throw Exception<Direction>("wrong theta value >pi");
    }
    Direction(RANDOM&RG):m_ld(RG),m_theta(acos(CTHETA(RG))){}
    const numt &phi()const
    {
        return m_ld.phi();
    }
    template<size_t index=1>
    const numt &th()const
    {
        return (index == 1) ? m_theta : MakeException<numt>();
    }
    Direction(const VType&V):m_ld(V.___recursive()),m_theta(acos(V.template component<Dimensions>()/V.M())){}
    const VType operator*(const numt&rho)const{
	return VType(m_ld*(rho*sin(m_theta)),rho*cos(m_theta));
    }
};
template<size_t size, class numt>
class Direction
{
public:
    enum {Dimensions = size};
    enum {Thetas = size-2};
    typedef numt NumberType;
    typedef Vector<size,numt> VType;
    typedef Direction<size-1,numt> DirectionN;
private:
    DirectionN m_ld;
    numt m_theta;
    static const numt CTHETA(RANDOM&r){
	static const RandomUniform<numt> res(numt(-1),numt(+1));
	return res(r);
    }
public:
    virtual ~Direction() {}
    template<class... Args>
    Direction(const std::tuple<Args...>&v):m_ld(v),m_theta(std::get<Thetas>(v)) {
	if(m_theta<numt(0))throw Exception<Direction>("wrong theta value <0");
	if(m_theta>PI<numt>())throw Exception<Direction>("wrong theta value >pi");
    }
    Direction(RANDOM&RG):m_ld(RG),m_theta(acos(CTHETA(RG))){}
    inline const numt &phi()const
    {
        return m_ld.phi();
    }
    template<size_t index>
    inline const numt &th()const
    {
        return (index == Thetas) ? m_theta : m_ld.template th<index>();
    }
    Direction(const VType&V):m_ld(V.___recursive()),m_theta(acos(V.template component<Dimensions>()/V.M())){}
    const VType operator*(const numt&rho)const{
	return VType(m_ld*(rho*sin(m_theta)),rho*cos(m_theta));
    }
};
template<class numt = double, class... Args>
const Direction< 2 + sizeof...(Args), numt > direction(const numt &phi, Args... other)
{
    return Direction < 2 + sizeof...(Args), numt > (std::make_tuple(phi,other...));
}
template<size_t size,class numt = double>
const Direction< size, numt > direction(const Vector<size,numt>&V)
{
    return Direction < size, numt > (V);
}
template<class numt = double, class... Args>
const auto fromPolar(const numt &rho, Args... other)->typename decltype(direction(other...))::VType
{
    return direction(other...)*rho;
}
template<size_t size=3,class numt = double, class RG = RANDOM>
const Direction<size, numt> randomIsotropic(RG &generator)
{
    return Direction<size,numt>(generator);
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
    return desCartes(dest[0], dest[1]);
}
template<class numt = double>
const Vector<3, numt> Rotate(const Vector<3, numt> &src, const Direction<3, numt> &axis, const numt &theta)
{
    const auto n = axis*numt(1);
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
    return desCartes(dest[0], dest[1], dest[2]);

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
    const numt &T()const
    {
        return m_time;
    }
    const Space &S()const
    {
        return m_space;
    }
    const numt M_sqr()const
    {
        return (m_time * m_time) - m_space.M_sqr();
    }
    const numt M()const
    {
        return sqrt(M_sqr());
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
    template<class...Args>
    const LorentzVector Rotate(Args...args)const
    {
        return LorentzVector(T(), MathTemplates::Rotate(S(), args...));
    }
    const LorentzVector Transform(const Space &Beta)const
    {
        const numt beta = Beta.M();
        if (beta == 0.0)return *this;
        if (beta >= numt(1))throw Exception<LorentzVector>("Bad Lorentz transformation");
        const auto n = (Beta / beta).chain();
        const numt gamma = numt(1) / sqrt(numt(1) - beta * beta);
        auto source = S().chain();
        source.insert(source.begin(),T());
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
        return S()/T();
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
    return LorentzVector<numt, Space>(sqrt(s.M_sqr() + l4 * l4), s);
}
template<class numt = double, class Space = Vector<3, numt>, class...Args>
const LorentzVector<numt, Space> lorentz_byEM(const numt &t, const numt &l4, Args...args)
{
    numt Sp = sqrt(t * t - l4 * l4);
    return LorentzVector<numt, Space>(t, direction(args...) * Sp);
}
template<size_t size,class numt>
const std::pair<LorentzVector<numt,Vector<size,numt>>,LorentzVector<numt,Vector<size,numt>>> 
binaryDecay(const numt &IM, const numt &m1, const numt &m2,const Direction<size,numt>&dir)
{
    if (m1 < 0)throw Exception<std::pair<LorentzVector<numt, Vector<1, numt>>, LorentzVector<numt, Vector<1, numt>>>>("Negative mass1 error");
    if (m2 < 0)throw Exception<std::pair<LorentzVector<numt, Vector<1, numt>>, LorentzVector<numt, Vector<1, numt>>>>("Negative mass2 error");
    if (IM < (m1 + m2))throw Exception<std::pair<LorentzVector<numt, Vector<1, numt>>, LorentzVector<numt, Vector<1, numt>>>>("Invariant mass of decaying system is less then masses of products");
    const auto E1=(IM*IM-m2*m2+m1*m1)/(IM*numt(2));
    const auto p=sqrt(E1*E1-m1*m1);
    const auto P=dir*p;
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
        if ((x ^ y).M_sqr() == numt(0))throw Exception<Plane3D>("Invalid basis vectors for converting 2d vectors to 3d");
    }
    Plane3D(const typename Space::DType &x, const typename Space::DType &y):Plane3D(x*numt(1),y*numt(1)){}
    static const Plane3D ByNormalVectorAndTheta(const typename Space::DType &N, const numt &theta)
    {
        if (((N*numt(1)) ^ Space::template basis_vector<3>()).M() == 0) {
            const auto X = Rotate(Space::template basis_vector<1>(), N, theta);
            return Plane3D(X, (N*numt(1)) ^ X);
        } else {
            const auto X = ((N*numt(1)) ^ Space::template basis_vector<3>());
            return Plane3D(X / X.M(), ((N*numt(1)) ^ (X / X.M())));
        }
        throw Exception<Plane3D>("This line should not be reached");
    }
};




};
#endif
