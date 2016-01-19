// this file is distributed under 
// MIT license
#ifndef zVoOnNfd
#define zVoOnNfd
#include <exception>
#include <string>
template<class SOURCE, int code=0>
class Error:public std::exception{
private:
	std::string m_msg;
public:
	Error(std::string msg){m_msg=msg;}
	virtual ~Error() throw(){}
	virtual const char* what() const throw(){return m_msg.c_str();}
};

#endif  
