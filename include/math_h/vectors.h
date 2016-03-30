// this file is distributed under 
// MIT license
#ifndef ___________VECTORS_H_____
#	define ___________VECTORS_H_____
#include <math.h>
#include "error.h"
namespace MathTemplates{
	template<class numt>
	class Vector3{
	private:
		numt m_x,m_y,m_z;
		Vector3(const numt&x,const numt&y,const numt&z):m_x(x),m_y(y),m_z(z){}
	public:
		static inline const Vector3 zero(){return Vector3(numt(0),numt(0),numt(0));}
		static inline const Vector3 basis_x(){return Vector3(numt(1),numt(0),numt(0));}
		static inline const Vector3 basis_y(){return Vector3(numt(0),numt(1),numt(0));}
		static inline const Vector3 basis_z(){return Vector3(numt(0),numt(0),numt(1));}
		static const Vector3 DesCartes(const numt&x,const numt&y,const numt&z){return Vector3(x,y,z);}
		static const Vector3 DesCartes(const numt&&x,const numt&&y,const numt&&z){return Vector3(x,y,z);}
		static const Vector3 Polar(const numt&mag,const numt&theta,const numt&phi){
			return Vector3(mag*cos(phi)*sin(theta),mag*sin(phi)*sin(theta),mag*cos(theta));
		}
		static inline const Vector3 Polar(const numt&mag,const numt&&theta,const numt&&phi){return Polar(mag,theta,phi);}
		static inline const Vector3 Polar(const numt&&mag,const numt&theta,const numt&phi){return Polar(mag,theta,phi);}
		static inline const Vector3 Polar(const numt&&mag,const numt&&theta,const numt&&phi){return Polar(mag,theta,phi);}
		static inline const Vector3 Direction(const numt&theta,const numt&phi){return Polar(numt(1),theta,phi);}
		static inline const Vector3 Direction(const numt&&theta,const numt&phi){return Polar(numt(1),theta,phi);}
		static inline const Vector3 Direction(const numt&theta,const numt&&phi){return Polar(numt(1),theta,phi);}
		static inline const Vector3 Direction(const numt&&theta,const numt&&phi){return Polar(numt(1),theta,phi);}
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
		Vector3&operator+=(const Vector3&&second){
			return operator+=(second);
		}
		const Vector3 operator+(const Vector3&second)const{
			return Vector3(m_x+second.m_x,m_y+second.m_y,m_z+second.m_z);
		}
		const Vector3 operator+(const Vector3&&second)const{
			return operator+(second);
		}
		Vector3&operator-=(const Vector3&second){
			m_x-=second.m_x;m_y-=second.m_y;m_z-=second.m_z;
			return *this;
		}
		Vector3&operator-=(const Vector3&&second){
			return operator-=(second);
		}
		const Vector3 operator-(const Vector3&second)const{
			return Vector3(m_x-second.m_x,m_y-second.m_y,m_z-second.m_z);
		}
		const Vector3 operator-(const Vector3&&second)const{
			return operator-(second);
		}
		
		Vector3&operator*=(const numt&second){
			m_x*=second;m_y*=second;m_z*=second;
			return *this;
		}
		Vector3&operator*=(const numt&&second){
			return operator*=(second);
		}
		const Vector3 operator*(const numt&second)const{
			return Vector3(m_x*second,m_y*second,m_z*second);
		}
		const Vector3 operator*(const numt&&second)const{
			return operator*(second);
		}
		
		Vector3&operator/=(const numt&second){
			m_x/=second;m_y/=second;m_z/=second;
			return *this;
		}
		Vector3&operator/=(const numt&&second){
			return operator/=(second);
		}
		const Vector3 operator/(const numt&second)const{
			return Vector3(m_x/second,m_y/second,m_z/second);
		}
		const Vector3 operator/(const numt&&second)const{
			return operator/(second);
		}
		
		const numt operator*(const Vector3&second)const{
			return (m_x*second.m_x)+(m_y*second.m_y)+(m_z*second.m_z);
		}
		const numt operator*(const Vector3&&second)const{
			return operator*(second);
		}
		
		const Vector3 VecP(const Vector3&second)const{
			return Vector3<numt>::DesCartes(
				(y()*second.z())-(second.y()*z()),
							(z()*second.x())-(second.z()*x()),
							(x()*second.y())-(second.x()*y())
			);
		}
		const Vector3 VecP(const Vector3&&second)const{
			return VecP(second);
		}
		
		bool operator==(const Vector3&second)const{
			return (m_x==second.m_x)&&(m_y==second.m_y)&&(m_z==second.m_z);
		}
		bool operator==(const Vector3&&second)const{return operator==(second);}
	};
	template<class numt>
	inline const Vector3<numt> operator-(const Vector3<numt>&V){return V*numt(-1);}
	template<class numt>
	inline const Vector3<numt> operator-(const Vector3<numt>&&V){return V*numt(-1);}
	template<typename numt>
	std::ostream&operator<<(std::ostream&str,const Vector3<numt>&V){
		return str<<V.x()<<" "<<V.y()<<" "<<V.z()<<" ";
	}
	template<typename numt>
	std::ostream&operator<<(std::ostream&str,const Vector3<numt>&&V){return str<<V;}
	
	
	template<class numt>
	class Vector4{
	public:
		typedef Vector3<numt> Space;
	private:
		numt m_time;
		Space m_space;
	public:
		Vector4(const numt&t,const Space&S):m_time(t),m_space(S){}
		Vector4(const numt&&t,const Space&S):Vector4(t,S){}
		Vector4(const numt&t,const Space&&S):Vector4(t,S){}
		Vector4(const numt&&t,const Space&&S):Vector4(t,S){}
		Vector4(const Vector4&source):Vector4(source.m_time,source.m_space){}
		Vector4&operator=(const Vector4&source){
			m_space=source.m_space;
			m_time=source.m_time;
			return *this;
		}
		const numt&time_component()const{return m_time;}
		const Space&space_component()const{return m_space;}
		const numt length4()const{
			return sqrt((m_time*m_time)-m_space.mag_sqr());
		}
		
		Vector4&operator+=(const Vector4&second){
			m_time+=second.m_time;
			m_space+=second.m_space;
			return *this;
		}
		Vector4&operator+=(const Vector4&&second){
			return operator+=(second);
		}
		const Vector4 operator+(const Vector4&second)const{
			return Vector4(m_time+second.m_time,m_space+second.m_space);
		}
		const Vector4 operator+(const Vector4&&second)const{
			return operator+(second);
		}
		Vector4&operator-=(const Vector4&second){
			m_time-=second.m_time;
			m_space-=second.m_space;
			return *this;
		}
		Vector4&operator-=(const Vector4&&second){
			return operator-=(second);
		}
		const Vector4 operator-(const Vector4&second)const{
			return Vector4(m_time-second.m_time,m_space-second.m_space);
		}
		const Vector4 operator-(const Vector4&&second)const{
			return operator-(second);
		}
		
		
		const numt operator*(const Vector4&second)const{
			return (m_time*second.m_time)-m_space*second.m_space;
		}
		const numt operator*(const Vector4&&second)const{
			return operator*(second);
		}
		
		static const Vector4 zero(){return Vector4(numt(0),Space::zero());}
		static const Vector4 SpaceLength4(const Space&s,const numt&l4){
			return Vector4(sqrt(s.mag_sqr()+l4*l4),s);
		}
		static inline const Vector4 SpaceLength4(const Space&s,const numt&&l4){return SpaceLength4(s,l4);}
		static inline const Vector4 SpaceLength4(const Space&&s,const numt&l4){return SpaceLength4(s,l4);}
		static inline const Vector4 SpaceLength4(const Space&&s,const numt&&l4){return SpaceLength4(s,l4);}
		static const Vector4 TimeDirLength4(const numt&t,const numt&theta,const numt&phi,const numt&l4){
			numt Sp=sqrt(t*t-l4*l4);
			return Vector4(t,Space::Direction(theta,phi)*Sp);
		}
		static const Vector4 TimeDirLength4(const numt&&t,const numt&theta,const numt&phi,const numt&l4){return TimeDirLength4(t,theta,phi,l4);}
		static const Vector4 TimeDirLength4(const numt&t,const numt&&theta,const numt&&phi,const numt&l4){return TimeDirLength4(t,theta,phi,l4);}
		static const Vector4 TimeDirLength4(const numt&t,const numt&theta,const numt&phi,const numt&&l4){return TimeDirLength4(t,theta,phi,l4);}
		static const Vector4 TimeDirLength4(const numt&&t,const numt&&theta,const numt&&phi,const numt&l4){return TimeDirLength4(t,theta,phi,l4);}
		static const Vector4 TimeDirLength4(const numt&t,const numt&&theta,const numt&&phi,const numt&&l4){return TimeDirLength4(t,theta,phi,l4);}
		static const Vector4 TimeDirLength4(const numt&&t,const numt&theta,const numt&phi,const numt&&l4){return TimeDirLength4(t,theta,phi,l4);}
		static const Vector4 TimeDirLength4(const numt&&t,const numt&&theta,const numt&&phi,const numt&&l4){return TimeDirLength4(t,theta,phi,l4);}
	};
};
#endif