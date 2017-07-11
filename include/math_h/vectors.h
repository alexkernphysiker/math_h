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
    public:
	typedef numt NumberType;
	virtual ~Vector2(){}
	Vector2(const numt&x,const numt&y):m_x(x),m_y(y){}
	const std::vector<numt> to_vector()const{return {m_x,m_y};}
	Vector2(const std::vector<numt>&v){
	    if(v.size()!=2)throw Exception<Vector2>("Bad Vector2 initializing by vector");
	    m_x=v[0];m_y=v[1];
	}
	static inline const Vector2 zero(){return Vector2(numt(0),numt(0));}
	static inline const Vector2 basis_x(){return Vector2(numt(1),numt(0));}
	static inline const Vector2 basis_y(){return Vector2(numt(0),numt(1));}
	static inline const Vector2 DesCartes(const numt&x,const numt&y){return Vector2(x,y);}
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
	    return Vector2(m_x+second.m_x,m_y+second.m_y);
	}
	Vector2&operator-=(const Vector2&second){
	    m_x-=second.m_x;m_y-=second.m_y;
	    return *this;
	}
	const Vector2 operator-(const Vector2&second)const{
	    return Vector2(m_x-second.m_x,m_y-second.m_y);
	}
		
	Vector2&operator*=(const numt&second){
	    m_x*=second;m_y*=second;
	    return *this;
	}
	const Vector2 operator*(const numt&second)const{
	    return Vector2(m_x*second,m_y*second);
	}
	Vector2&operator/=(const numt&second){
	    m_x/=second;m_y/=second;
	    return *this;
	}
	const Vector2 operator/(const numt&second)const{
	    return Vector2(m_x/second,m_y/second);
	}
	const numt operator*(const Vector2&second)const{
	    return (m_x*second.m_x)+(m_y*second.m_y);
	}
	const bool operator==(const Vector2&second)const{
	    return (m_x==second.m_x)&&(m_y==second.m_y);
	}
	const numt VecP(const Vector2&second)const{
	    return (x()*second.y())-(second.x()*y());
	}
	const Vector2 Rotate(const numt&theta)const{
	    const std::vector<numt> source={x(),y()};
	    const numt cost=cos(theta),sint=sin(theta);
	    std::vector<std::vector<numt>> M={
		{cost,-sint},
		{sint,cost}
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
    public:
	virtual ~Vector3(){}
	Vector3(const numt&x,const numt&y,const numt&z):m_x(x),m_y(y),m_z(z){}
	typedef numt NumberType;
	const std::vector<numt> to_vector()const{return {m_x,m_y,m_z};}
	Vector3(const std::vector<numt>&v){
	    if(v.size()!=3)throw Exception<Vector3>("Bad Vector3 initializing by vector");
	    m_x=v[0];m_y=v[1];m_z=v[2];
	}
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
    template<class Space=Vector3<double>,class numt=typename Space::NumberType>
    class LorentzVector{
    private:
	numt m_time;
	Space m_space;
    public:
	virtual ~LorentzVector(){}
	LorentzVector(const numt&t,const Space&S):m_time(t),m_space(S){}
	typedef Space SpaceVectorType;
	typedef numt TimeCoordinateType;
	LorentzVector(const LorentzVector&source):LorentzVector(source.m_time,source.m_space){}
	LorentzVector&operator=(const LorentzVector&source){
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
	const bool operator==(const LorentzVector&second)const{
	    return (m_time==second.m_time)&&(m_space==second.m_space);
	}
	LorentzVector&operator+=(const LorentzVector&second){
	    m_time+=second.m_time;
	    m_space+=second.m_space;
	    return *this;
	}
	const LorentzVector operator+(const LorentzVector&second)const{
	    return LorentzVector(m_time+second.m_time,m_space+second.m_space);
	}
	LorentzVector&operator-=(const LorentzVector&second){
	    m_time-=second.m_time;
	    m_space-=second.m_space;
	    return *this;
	}
	const LorentzVector operator-(const LorentzVector&second)const{
	    return LorentzVector(m_time-second.m_time,m_space-second.m_space);
	}
	const numt operator*(const LorentzVector&second)const{
	    return (m_time*second.m_time)-(m_space*second.m_space);
	}
	static const LorentzVector zero(){return LorentzVector(numt(0),Space::zero());}
	static const LorentzVector byComponents(const numt&t,const numt&x,const numt&y,const numt&z){
	    return LorentzVector(t,Space::DesCartes(x,y,z));
	}
	static const LorentzVector byComponents(const numt&t,const Space&s){
	    return LorentzVector(t,s);
	}
	LorentzVector(const numt&t):LorentzVector(t,Space::zero()){}
	static const LorentzVector byOnlyTimeC(const numt&t){
	    return LorentzVector(t);
	}
	static const LorentzVector bySpaceC_and_Length4(const Space&s,const numt&l4){
	    return LorentzVector(sqrt(s.mag_sqr()+l4*l4),s);
	}
	static const LorentzVector byTime_Dir_and_Length4(const numt&t,const numt&theta,const numt&phi,const numt&l4){
	    numt Sp=sqrt(t*t-l4*l4);
	    return LorentzVector(t,Space::Direction(theta,phi)*Sp);
	}
	const LorentzVector Rotate(const Space&axis,const numt&theta)const{
	    return LorentzVector(time_component(),space_component().Rotate(axis,theta));
	}
	const LorentzVector Lorentz(const Space&Beta)const{
	    const numt beta=Beta.mag();
	    if(beta==0.0)return *this;
	    if(beta>=numt(1))throw Exception<LorentzVector>("Bad Lorentz transformation");
	    const auto n=(Beta/beta).to_vector();
	    const numt gamma=numt(1)/sqrt(numt(1)-beta*beta);
	    auto source=space_component().to_vector();
	    source.insert(source.begin(),time_component());
	    std::vector<std::vector<numt>> M;
	    for(size_t i=0;i<source.size();i++){
		std::vector<numt> V;
		for(size_t j=0;j<source.size();j++){
		    if((0==i)&&(0==j)){
			V.push_back(gamma);
		    }else{
			if(0==i){
			    V.push_back(-gamma*beta*n[j-1]);
			}else{
			    if(0==j){
				V.push_back(-gamma*beta*n[i-1]);
			    }else{
				if(i==j){
				    V.push_back(numt(1)+(gamma-numt(1))*n[i-1]*n[j-1]);
				}else{
				    V.push_back((gamma-numt(1))*n[i-1]*n[j-1]);
				}
			    }
			}
		    }
		}
		M.push_back(V);
	    }
	    std::vector<numt> dest;
	    for(size_t i=0;i<source.size();i++){
		numt v=0;
		for(size_t j=0;j<source.size();j++)v+=M[i][j]*source[j];
		dest.push_back(v);
	    }
	    const auto T=dest[0];
	    dest.erase(dest.begin());
	    return LorentzVector(T,dest);
	}
	const Space Beta()const{//Makes physical sense only if it's a 4-momentum
	    return space_component()/time_component();
	}
    };
    template<class numt=double>
    using Vector4=LorentzVector<Vector3<numt>,numt>;

    template<class numt=double>
    class Plane3D{
    public:
	typedef Vector3<numt> Space;
	typedef Vector2<numt> Plane;
    private:
	Space basis_x,basis_y;
    public:
	const Space operator()(const Plane&v)const{
	    return basis_x*v.x()+basis_y*v.y();
	}
	const LorentzVector<Space> operator()(const LorentzVector<Plane>&v)const{
	    return LorentzVector<Space>(v.time_component(),basis_x*v.space_component().x()+basis_y*v.space_component().y());
	}
	virtual ~Plane3D(){}
	Plane3D(const Space&x,const Space&y):basis_x(x),basis_y(y){
	    if(x.VecP(y).mag_sqr()==numt(0))throw Exception<Plane3D>("Invalid basis vectors for converting 2d vectors to 3d");
	}
	static const Plane3D XY(){return Plane3D(Space::basis_x(),Space::basis_y());}
	static const Plane3D YX(){return Plane3D(Space::basis_y(),Space::basis_x());}
	static const Plane3D XZ(){return Plane3D(Space::basis_x(),Space::basis_z());}
	static const Plane3D ZX(){return Plane3D(Space::basis_z(),Space::basis_x());}
	static const Plane3D YZ(){return Plane3D(Space::basis_y(),Space::basis_z());}
	static const Plane3D ZY(){return Plane3D(Space::basis_z(),Space::basis_y());}
	static const Plane3D ByNormalVectorAndTheta(const Space&N,const numt&theta){
	    if(N.mag()==0)throw Exception<Plane3D>("Invalid normale vector");
	    const auto n=N/N.mag();
	    if(n.VecP(Space::basis_z()).mag_sqr()==0){
		const auto X=Space::basis_x().Rotate(n,theta);
		return Plane3D(X,n.VecP(X));
	    }else{
		const auto X=n.VecP(Space::basis_z());
		return Plane3D(X/X.mag(),n.VecP(X/X.mag()));
	    }
	    throw Exception<Plane3D>("This line should not be reached");
	}
    };

};
#endif
