//https://github.com/alexkernphysiker/math_h
#ifndef ___SIGMA_H___
#	define ___SIGMA_H___
#include <list>
template<typename numt>
class Sigma{
private:
	std::list<numt> m_list;
	numt m_avr;bool cached=false;
public:
	Sigma(){cached=false;}
	virtual ~Sigma(){}
	void AddValue(numt value){m_list.push_back(value);cached=false;}
	numt getAverage(){
		if(!cached){
			numt res=0;
			foreach(double value,m_list)
				res+=value;
			m_avr=res/m_list.size();
			cached=true;
		}
		return m_avr;
	}
	numt operator[](int i){return m_list[i];}
	int count(){return m_list.size();}
	numt getSigma(){
		getAverage();
		numt res=0;
		foreach(double value,m_list)
			res+=pow(value-m_avr,2);
		res/=m_list.size()-1;
		return sqrt(res);
	}
};
#endif
