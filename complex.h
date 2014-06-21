#ifndef COMPLEX____
#	define COMPLEX____
template<class numt=double>
class Complex{
private:
	numt m_re;
	numt m_im;
public:
	Complex(){m_re=0;m_im=0;}
	Complex(const numt& re){m_re=re;m_im=0;}
	Complex(const numt& re, const numt& im){m_re=re;m_im=im;}
	Complex(const Complex& parent){m_re=parent.re();m_im=parent.im();}
	virtual ~Complex(){}
	//properties
	numt re()const{return m_re;}	numt im()const{return m_im;}
	numt absSqr()const{return m_re*m_re+m_im*m_im;}	numt abs()const{return sqrt(absSqr());}
	numt arg()const{numt r=absSqr();
		if(r>0){return atan2(m_im,m_re);}else{return 0;}
	}
	Complex conj()const{return Complex(m_re,-m_im);}
	//changing this
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
		double r=m_re*z.m_re-m_im*z.m_im;
		double i=m_im*z.m_re+m_re*z.m_im;
		m_re=r;m_im=i;
		return *this;
	}
	Complex& operator*=(const numt c){
		m_re*=c;
		m_im*=c;
		return *this;
	}
	Complex& operator/=(const Complex z){
		numt below=z.absSqr();
		if(below==0) throw;//division by zero
		numt r=m_re*z.m_re+m_im*z.m_im;
		numt i=m_im*z.m_re-m_re*z.m_im;
		m_re=r;m_im=i;
		operator*=(1/below);
		return *this;
	}
	Complex& operator/=(const numt c){
		m_re/=c;
		m_im/=c;
		return *this;
	}
};
template<class numt=double>
Complex<numt> I(){return Complex<numt>(0,1);}
// operators
template<class A, class B>
A operator+(A a, B b){A res=a; res+=A(b);return res;}
template<class A, class B>
A operator-(A a, B b){A res=a; res-=A(b);return res;}
template<class A, class B>
A operator*(A a, B b){A res=a; res*=A(b);return res;}
template<class A, class B>
A operator/(A a, B b){A res=a; res/=A(b);return res;}
template<class numt=double>
Complex<numt> operator-(const Complex<numt> a){Complex<numt> res=0;return res-=a;}

//functions
template<class numt=double>
numt abs(const Complex<numt> z){return z.abs();}
template<class numt=double>
Complex<numt> exp(const Complex<numt> z){
	return Complex<numt>(cos(z.im()),sin(z.im()))*exp(z.re());
}
template<class numt=double>
Complex<numt> log(const Complex<numt> z){
	return Complex<numt>(log(abs(z)),z.arg());
}
template <class A, class B>  A pow(A a, B b){
	A r=log(a);r*=b;return exp(r);
}
#endif // COMPLEX____
