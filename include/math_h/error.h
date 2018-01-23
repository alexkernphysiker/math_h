// this file is distributed under
// LGPLv3 license
#ifndef zVoOnNfd
#define zVoOnNfd
#if __cplusplus<201100L
#error c++>=11 is needed for using math_h headers
#endif
#if __cplusplus>201700L
#define ____full_version_of_math_h_____
#else
#warning c++17 is needed to compile full featured version of math_h
#endif
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
