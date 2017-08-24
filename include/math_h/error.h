// this file is distributed under
// MIT license
#ifndef zVoOnNfd
#define zVoOnNfd
#include <exception>
#include <string>
namespace MathTemplates
{
template<typename SOURCE, int code = 0>
class Exception: public std::exception
{
private:
    std::string m_msg;
public:
    Exception(std::string &&msg)
    {
        m_msg = msg;
    }
    virtual ~Exception() throw() {}
    virtual const char *what() const throw()
    {
        return m_msg.c_str();
    }
};
};
#endif
