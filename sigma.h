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
	void AddValue(numt value){
		m_list.push_back(value);
		m_sum+=value;
		cache=false;
	}
	numt getAverage(){
		int sz=m_list.size();
		if(sz<=0)throw;
		return m_sum/sz;
	}
	int count(){return m_list.size();}
	numt getSigmaSqr(){
		if(!cache){
			numt average=getAverage();
			m_sigsqr=0;
			auto it=m_list.begin();
			while(it!=m_list.end()){
				numt value=*it;
				numt d=abs(value-average);
				m_sigsqr+=d*d;
				it++;
			}
			m_sigsqr/=m_list.size()-1;
			cache=true;
		}
		return m_sigsqr;
	}
	numt getSigma(){return sqrt(getSigmaSqr());}
};
#endif
