#ifndef YJIOGPKSIIYJVKND
#define YJIOGPKSIIYJVKND
#include <list>
#include <math.h>
#include <exception>
template<typename numt>
class Sigma{
private:
	std::list<numt> m_list;numt m_sum;
	bool cache; numt m_sigsqr;
public:
	Sigma(){m_sum=0;cache=false;m_sigsqr=0;}
	virtual ~Sigma(){}
	void AddValue(numt value){
		m_list.push_back(value);
		m_sum+=value;
		cache=false;
	}
	numt getAverage(){
		int sz=m_list.size();
		if(sz<=0)throw std::exception();
		return m_sum/sz;
	}
	int count(){return m_list.size();}
	numt getSigmaSqr(){
		int sz=m_list.size();
		if(sz<=1)throw std::exception();
		if(!cache){
			numt average=getAverage();
			m_sigsqr=0;
			for(auto value:m_list)
				m_sigsqr+=pow(value-average,2);
			m_sigsqr/=sz-1;
			cache=true;
		}
		return m_sigsqr;
	}
	numt getSigma(){return sqrt(getSigmaSqr());}
};
#endif
