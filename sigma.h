//https://github.com/alexkernphysiker/math_h
#ifndef ___SIGMA_H___
#	define ___SIGMA_H___
#include <list>
template<typename numt>
class Sigma{
private:
	std::list<numt> m_list;numt m_sum;
	bool cache; numt m_sigsqr;
public:
	Sigma(){m_sum=0;cache=false;m_sigsqr=0;}
	virtual ~Sigma(){}
	void AddValue(numt value){m_list.push_back(value);m_sum+=value;cache=false;}
	numt getAverage(){return m_sum/m_list.size();}
	int count(){return m_list.size();}
	numt getSigmaSqr(){
		if(!cache){
			double average=getAverage();
			m_sigsqr=0;
			foreach(double value,m_list)
				m_sigsqr+=pow(value-average,2);
			m_sigsqr/=m_list.size()-1;
			cache=true;
		}
		return m_sigsqr;
	}
	numt getSigma(){return sqrt(getSigmaSqr());}
};
#endif
