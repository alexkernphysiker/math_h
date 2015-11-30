// this file is distributed under 
// GPL v 3.0 license
#ifndef zVoOnNfd
#define zVoOnNfd
#include <exception>
#include <string>
template<class SOURCE>
class math_h_error:public std::exception{
private:
	std::string m_msg;
public:
	math_h_error(std::string msg){m_msg=msg;}
	virtual ~math_h_error() throw(){}
	virtual const char* what() const throw(){return m_msg.c_str();}
};

#endif  
