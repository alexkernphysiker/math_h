// this file is distributed under 
// MIT license
#ifndef ___________VECTORS_H_____
#	define ___________VECTORS_H_____
#include <math.h>
#include <vector>
#include "error.h"
#include "functions.h"
#include "randomfunc.h"
namespace MathTemplates{
    template<class numt=double>
    class Vector2{
    private:
	numt m_x,m_y;
	Vector2(const numt&x,const numt&y):m_x(x),m_y(y){}
    public:
	static inline const Vector2 zero(){return Vector2(numt(0),numt(0));}
	static inline const Vector2 basis_x(){return Vector2(numt(1),numt(0));}
	static inline const Vector2 basis_y(){return Vector2(numt(0),numt(1));}
	static inline const Vector2 DesCartes(const numt&x,const numt&y){return Vector3(x,y);}
	static const Vector2 Polar(const numt&mag,const numt&phi){
	    return Vector2(mag*cos(phi),mag*sin(phi));
	}
	static inline const Vector2 Direction(const numt&phi){return Polar(numt(1),phi);}
	const numt&x()const{return m_x;}
	const numt&y()const{return m_y;}
	const numt mag_sqr()const{
	    return (m_x*m_x)+(m_y*m_y);
	}
	const numt mag()const{
	    return sqrt(mag_sqr());
	}
	Vector2(const Vector2&source):m_x(source.m_x),m_y(source.m_y){}
	Vector2&operator=(const Vector2&source){
	    m_x=source.m_x;m_y=source.m_y;
	    return *this;
	}
	Vector2&operator+=(const Vector2&second){
	    m_x+=second.m_x;m_y+=second.m_y;
	    return *this;
	}
	const Vector2 operator+(const Vector2&second)const{
	    return Vector3(m_x+second.m_x,m_y+second.m_y);
	}
	Vector2&operator-=(const Vector2&second){
	    m_x-=second.m_x;m_y-=second.m_y;
	    return *this;
	}
	const Vector2 operator-(const Vector2&second)const{
	    return Vector3(m_x-second.m_x,m_y-second.m_y);
	}
		
	Vector2&operator*=(const numt&second){
	    m_x*=second;m_y*=second;
	    return *this;
	}
	const Vector2 operator*(const numt&second)const{
	    return Vector3(m_x*second,m_y*second);
	}
	Vector2&operator/=(const numt&second){
	    m_x/=second;m_y/=second;
	    return *this;
	}
	const Vector2 operator/(const numt&second)const{
	    return Vector3(m_x/second,m_y/second);
	}
	const numt operator*(const Vector2&second)const{
	    return (m_x*second.m_x)+(m_y*second.m_y);
	}
	const bool operator==(const Vector2&second)const{
	    return (m_x==second.m_x)&&(m_y==second.m_y);
	}
	const Vector2 Rotate(const numt&theta)const{
	    const std::vector<numt> source={x(),y()};
	    const numt cost=cos(theta),sint=sin(theta);
	    std::vector<std::vector<numt>> M={
		{cost,sint},
		{-sint,cost}
	    };
	    std::vector<numt> dest;
	    for(size_t i=0;i<2;i++){
		numt v=0;
		for(size_t j=0;j<2;j++)v+=M[i][j]*source[j];
		dest.push_back(v);
	    }
	    return Vector2(dest[0],dest[1]);

	}

	template<class RG=RANDOM>
	static inline const Vector2 RandomIsotropicDirection(RG&generator){
	    static RandomUniform<numt> Phi(0,PI<numt>()*2.0);
	    const numt phi=Phi(generator);
	    return Vector2(cos(phi),sin(phi));
	}
    };
    template<class numt>
    inline const Vector2<numt> operator-(const Vector2<numt>&V){return V*numt(-1);}
    template<typename numt>
    std::ostream&operator<<(std::ostream&str,const Vector2<numt>&V){
	return str<<V.x()<<" "<<V.y()<<" ";
    }
    template<class numt=double>
    class Vector3{
    private:
	numt m_x,m_y,m_z;
	Vector3(const numt&x,const numt&y,const numt&z):m_x(x),m_y(y),m_z(z){}
    public:
	static inline const Vector3 zero(){return Vector3(numt(0),numt(0),numt(0));}
	static inline const Vector3 basis_x(){return Vector3(numt(1),numt(0),numt(0));}
	static inline const Vector3 basis_y(){return Vector3(numt(0),numt(1),numt(0));}
	static inline const Vector3 basis_z(){return Vector3(numt(0),numt(0),numt(1));}
	static inline const Vector3 DesCartes(const numt&x,const numt&y,const numt&z){return Vector3(x,y,z);}
	static const Vector3 Polar(const numt&mag,const numt&theta,const numt&phi){
	    return Vector3(mag*cos(phi)*sin(theta),mag*sin(phi)*sin(theta),mag*cos(theta));
	}
	static inline const Vector3 Direction(const numt&theta,const numt&phi){
	    return Polar(numt(1),theta,phi);
	}
	const numt&x()const{return m_x;}
	const numt&y()const{return m_y;}
	const numt&z()const{return m_z;}
	const numt mag_sqr()const{
	    return (m_x*m_x)+(m_y*m_y)+(m_z*m_z);
	}
	const numt mag()const{
	    return sqrt(mag_sqr());
	}
	const numt cos_theta()const{
	    return m_z/mag();
	}
	const numt sin_theta()const{
	    return sqrt((m_x*m_x)+(m_y*m_y))/mag();
	}
	const numt cos_phi()const{
	    return m_x/sqrt((m_x*m_x)+(m_y*m_y));
	}
	const numt sin_phi()const{
	    return m_x/sqrt((m_x*m_x)+(m_y*m_y));
	}
		
	Vector3(const Vector3&source):m_x(source.m_x),m_y(source.m_y),m_z(source.m_z){}
	Vector3&operator=(const Vector3&source){
	    m_x=source.m_x;m_y=source.m_y;m_z=source.m_z;
	    return *this;
	}
	Vector3&operator+=(const Vector3&second){
	    m_x+=second.m_x;m_y+=second.m_y;m_z+=second.m_z;
	    return *this;
	}
	const Vector3 operator+(const Vector3&second)const{
	    return Vector3(m_x+second.m_x,m_y+second.m_y,m_z+second.m_z);
	}
	Vector3&operator-=(const Vector3&second){
	    m_x-=second.m_x;m_y-=second.m_y;m_z-=second.m_z;
	    return *this;
	}
	const Vector3 operator-(const Vector3&second)const{
	    return Vector3(m_x-second.m_x,m_y-second.m_y,m_z-second.m_z);
	}
		
	Vector3&operator*=(const numt&second){
	    m_x*=second;m_y*=second;m_z*=second;
	    return *this;
	}
	const Vector3 operator*(const numt&second)const{
	    return Vector3(m_x*second,m_y*second,m_z*second);
	}
	Vector3&operator/=(const numt&second){
	    m_x/=second;m_y/=second;m_z/=second;
	    return *this;
	}
	const Vector3 operator/(const numt&second)const{
	    return Vector3(m_x/second,m_y/second,m_z/second);
	}
	const numt operator*(const Vector3&second)const{
	    return (m_x*second.m_x)+(m_y*second.m_y)+(m_z*second.m_z);
	}
		
	const Vector3 VecP(const Vector3&second)const{
	    return Vector3<numt>::DesCartes(
		(y()*second.z())-(second.y()*z()),
		(z()*second.x())-(second.z()*x()),
		(x()*second.y())-(second.x()*y())
	    );
	}
	const bool operator==(const Vector3&second)const{
	    return (m_x==second.m_x)&&(m_y==second.m_y)&&(m_z==second.m_z);
	}
	const Vector3 Rotate(const Vector3&axis,const numt&theta)const{
	    const auto n=axis/axis.mag();
	    const std::vector<numt> source={x(),y(),z()};
	    const numt cost=cos(theta),sint=sin(theta),one=1;
	    std::vector<std::vector<numt>> M={
		{    cost+(one-cost)*n.x()*n.x(),	(one-cost)*n.x()*n.y()-sint*n.z(),	(one-cost)*n.x()*n.z()+sint*n.y()},
		{(one-cost)*n.y()*n.x()+sint*n.z(),	    cost+(one-cost)*n.y()*n.y(),	(one-cost)*n.y()*n.z()-sint*n.x()},
		{(one-cost)*n.z()*n.x()-sint*n.y(),	(one-cost)*n.z()*n.y()+sint*n.x(),	    cost+(one-cost)*n.z()*n.z()}
	    };
	    std::vector<numt> dest;
	    for(size_t i=0;i<3;i++){
		numt v=0;
		for(size_t j=0;j<3;j++)v+=M[i][j]*source[j];
		dest.push_back(v);
	    }
	    return Vector3(dest[0],dest[1],dest[2]);

	}

	template<class RG=RANDOM>
	static inline const Vector3 RandomIsotropicDirection(RG&generator){
	    static RandomUniform<numt> Z(-1,1);
	    static RandomUniform<numt> Phi(0,PI<numt>()*2.0);
	    const numt z=Z(generator);
	    const numt x_y=sqrt(1.0-pow(z,2));
	    const numt phi=Phi(generator);
	    return Vector3(cos(phi)*x_y,sin(phi)*x_y,z);
	}
	template<class RG=RANDOM>
	static inline const Vector3 RandomIsotropicXYDirection(RG&generator){
	    static RandomUniform<numt> Phi(0,PI<numt>()*2.0);
	    const numt phi=Phi(generator);
	    return Vector3(cos(phi),sin(phi),0);
	}
	template<class RG=RANDOM>
	static inline const Vector3 RandomIsotropicYZDirection(RG&generator){
	    static RandomUniform<numt> Phi(0,PI<numt>()*2.0);
	    const numt phi=Phi(generator);
	    return Vector3(0,cos(phi),sin(phi));
	}
	template<class RG=RANDOM>
	static inline const Vector3 RandomIsotropicXZDirection(RG&generator){
	    static RandomUniform<numt> Phi(0,PI<numt>()*2.0);
	    const numt phi=Phi(generator);
	    return Vector3(cos(phi),0,sin(phi));
	}
    };
    template<class numt>
    inline const Vector3<numt> operator-(const Vector3<numt>&V){return V*numt(-1);}
    template<typename numt>
    std::ostream&operator<<(std::ostream&str,const Vector3<numt>&V){
	return str<<V.x()<<" "<<V.y()<<" "<<V.z()<<" ";
    }
    template<class numt=double>
    class Vector4{
    public:
	typedef Vector3<numt> Space;
    private:
	numt m_time;
	Space m_space;
	Vector4(const numt&t,const Space&S):m_time(t),m_space(S){}
    public:
	Vector4(const Vector4&source):Vector4(source.m_time,source.m_space){}
	Vector4&operator=(const Vector4&source){
	    m_space=source.m_space;
	    m_time=source.m_time;
	    return *this;
	}
	const numt&time_component()const{return m_time;}
	const Space&space_component()const{return m_space;}
	const numt Sqr4()const{
	    return (m_time*m_time)-m_space.mag_sqr();
	}
	const numt length4()const{
	    return sqrt(Sqr4());
	}
	const bool operator==(const Vector4&second)const{
	    return (m_time==second.m_time)&&(m_space==second.m_space);
	}
	Vector4&operator+=(const Vector4&second){
	    m_time+=second.m_time;
	    m_space+=second.m_space;
	    return *this;
	}
	const Vector4 operator+(const Vector4&second)const{
	    return Vector4(m_time+second.m_time,m_space+second.m_space);
	}
	Vector4&operator-=(const Vector4&second){
	    m_time-=second.m_time;
	    m_space-=second.m_space;
	    return *this;
	}
	const Vector4 operator-(const Vector4&second)const{
	    return Vector4(m_time-second.m_time,m_space-second.m_space);
	}
	const numt operator*(const Vector4&second)const{
	    return (m_time*second.m_time)-(m_space*second.m_space);
	}
	static const Vector4 zero(){return Vector4(numt(0),Space::zero());}
	static const Vector4 byComponents(const numt&t,const numt&x,const numt&y,const numt&z){
	    return Vector4(t,Space::DesCartes(x,y,z));
	}
	static const Vector4 byComponents(const numt&t,const Space&s){
	    return Vector4(t,s);
	}
	Vector4(const numt&t):Vector4(t,Space::zero()){}
	static const Vector4 byOnlyTimeC(const numt&t){
	    return Vector4(t);
	}
	static const Vector4 bySpaceC_and_Length4(const Space&s,const numt&l4){
	    return Vector4(sqrt(s.mag_sqr()+l4*l4),s);
	}
	static const Vector4 byTime_Dir_and_Length4(const numt&t,const numt&theta,const numt&phi,const numt&l4){
	    numt Sp=sqrt(t*t-l4*l4);
	    return Vector4(t,Space::Direction(theta,phi)*Sp);
	}
	const Vector4 Rotate(const Space&axis,const numt&theta)const{
	    return Vector4(time_component(),space_component().Rotate(axis,theta));
	}
	const Vector4 Lorentz(const Space&Beta)const{
	    const numt beta=Beta.mag();
	    if(beta==0.0)return *this;
	    if(beta>=numt(1))throw Exception<Vector4>("Bad Lorentz transformation");
	    const Space n=Beta/beta;
	    const numt gamma=numt(1)/sqrt(numt(1)-beta*beta);
	    const std::vector<numt> source={
		time_component(),
		space_component().x(),
		space_component().y(),
		space_component().z()
	    };
	    std::vector<std::vector<numt>> M={
		{ gamma           ,-gamma*Beta.x(),-gamma*Beta.y(),-gamma*Beta.z()},
		{-gamma*Beta.x()  ,numt(1)+(gamma-numt(1))*n.x()*n.x(),numt(0)+(gamma-numt(1))*n.x()*n.y(),numt(0)+(gamma-numt(1))*n.x()*n.z()},
		{-gamma*Beta.y()  ,numt(0)+(gamma-numt(1))*n.y()*n.x(),numt(1)+(gamma-numt(1))*n.y()*n.y(),numt(0)+(gamma-numt(1))*n.y()*n.z()},
		{-gamma*Beta.z()  ,numt(0)+(gamma-numt(1))*n.z()*n.x(),numt(0)+(gamma-numt(1))*n.z()*n.y(),numt(1)+(gamma-numt(1))*n.z()*n.z()}
	    };
	    std::vector<numt> dest;
	    for(size_t i=0;i<4;i++){
		numt v=0;
		for(size_t j=0;j<4;j++)v+=M[i][j]*source[j];
		dest.push_back(v);
	    }
	    return Vector4(dest[0],Space::DesCartes(dest[1],dest[2],dest[3]));
	}
	const Space Beta()const{//Makes physical sense only if it's a 4-momentum
	    return space_component()/time_component();
	}
    };
};
#endif
