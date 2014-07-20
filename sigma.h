//https://github.com/alexkernphysiker/math_h
#ifndef ___SIGMA_H___
#	define ___SIGMA_H___
#include <list>
template<typename numt>
class Sigma{
private:
	std::list<numt> m_list;
public:
	Sigma(){}
	virtual ~Sigma(){}
	void AddValue(numt value){m_list.push_back(value);}
	numt getAverage(){
		numt res=0;
		foreach(double value,m_list)
			res+=value;
		return res/m_list.size();
	}
	numt getSigma(){
		numt avr=getAverage();
		numt res=0;
		foreach(double value,m_list)
			res+=pow(value-avr,2);
		res/=m_list.size()-1;
		return sqrt(res);
	}
};
#endif
