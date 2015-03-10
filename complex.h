#ifndef LECTGMFGFWVHIFQT
#	define LECTGMFGFWVHIFQT
#include <math.h>
template<typename numt>
class Complex{
private:
	numt m_re;
	numt m_im;
public:
	Complex():m_re(0),m_im(0){}
	Complex(const numt re):m_re(re),m_im(0){}
	Complex(const numt re, const numt im):m_re(re),m_im(im){}
	Complex(const Complex& parent):m_re(parent.m_re),m_im(parent.m_im){}
	virtual ~Complex(){}
	numt re()const{
		return m_re;
	}
	numt im()const{
		return m_im;
	}
	numt sqrabs()const{
		return m_re*m_re+m_im*m_im;
	}
	numt abs()const{
		return sqrt(sqrabs());
	}
	numt arg()const{numt r=sqrabs();
		if(r>0){
			return atan2(m_im,m_re);
		}else{
			return 0;
		}
	}
	Complex conj()const{
		return Complex(m_re,-m_im);
	}
	Complex& operator+=(const Complex z){
		m_re+=z.m_re;
		m_im+=z.m_im;
		return *this;
	}
	Complex& operator+=(const numt c){
		m_re+=c;
		return *this;
	}
	Complex& operator-=(const Complex z){
		m_re-=z.m_re;
		m_im-=z.m_im;
		return *this;
	}
	Complex& operator-=(const numt c){
		m_re-=c;
		return *this;
	}
	Complex& operator*=(const Complex z){
		numt r=m_re*z.m_re-m_im*z.m_im;
		numt i=m_im*z.m_re+m_re*z.m_im;
		m_re=r;
		m_im=i;
		return *this;
	}
	Complex& operator*=(const numt c){
		m_re*=c;
		m_im*=c;
		return *this;
	}
	Complex& operator/=(const numt c){
		if(c==0)throw;
		m_re/=c;
		m_im/=c;
		return *this;
	}
	Complex& operator/=(const Complex z){
		numt below=z.absSqr();
		if(below==0) throw;
		numt r=m_re*z.m_re+m_im*z.m_im;
		numt i=m_im*z.m_re-m_re*z.m_im;
		m_re=r;
		m_im=i;
		operator/=(below);
		return *this;
	}
};

template<typename numt>
inline Complex<numt> I(){
	return Complex<numt>(0,1);
}
template<typename numt>
inline numt re(Complex<numt> c){
	return c.re();
}
template<typename numt>
inline numt im(Complex<numt> c){
	return c.im();
}
template<typename numt>
inline numt abs(const Complex<numt> z){
	return z.abs();
}
template<typename numt>
inline numt arg(const Complex<numt> z){
	return z.arg();
}
template<typename numt>
inline numt conj(const Complex<numt> z){
	return z.conj();
}
template<typename numt, typename B>
inline Complex<numt> operator+(Complex<numt> a, B b){
	Complex<numt> res=a;
	res+=b;
	return res;
}
template<typename numt, typename B>
inline Complex<numt> operator-(Complex<numt> a, B b){
	Complex<numt> res=a;
	res-=b;
	return res;
}
template<typename numt, typename B>
inline Complex<numt> operator*(Complex<numt> a, B b){
	Complex<numt> res=a;
	res*=b;
	return res;
}
template<typename numt, typename B>
inline Complex<numt> operator/(Complex<numt> a, B b){
	Complex<numt> res=a;
	res/=b;
	return res;
}
template<typename numt>
inline Complex<numt> operator-(const Complex<numt> a){
	Complex<numt> res=0;
	return res-=a;
}
template<typename numt>
Complex<numt> exp(const Complex<numt> z){
	return Complex<numt>(cos(z.im()),sin(z.im()))*exp(z.re());
}
template<typename numt>
Complex<numt> log(const Complex<numt> z){
	return Complex<numt>(log(z.abs()),z.arg());
}
template <typename numt, typename B>
Complex<numt> pow(Complex<numt> a, B b){
	Complex<numt> r=log(a);
	r*=b;
	return exp(r);
}
#endif
