#ifndef ___SINGLE_PARAM_H
#	define ___SINGLE_PARAM_H

#include <utility>
namespace detail{// implementation details
	template <
		typename ret,typename F, typename Tuple,
		bool Done,
		 int Total, int... N
	>
	struct call_impl{
		static ret call(F f, Tuple && t){
			return call_impl<
					ret, F, Tuple,
					Total == 1 + sizeof...(N),
					Total, N..., sizeof...(N)
			>::call(f, std::forward<Tuple>(t));
		}
	};
	template <
			typename ret,typename F, typename Tuple,
			int Total, int... N
	>
	struct call_impl<ret, F, Tuple, true, Total, N...>{
		static ret call(F f, Tuple && t){
			return f(std::get<N>(std::forward<Tuple>(t))...);
		}
	};
}
//allows to call multi-parameter function as a single-parameter one
template<typename restype,int parn, typename... Args>
class SingleParam{
private:
	typedef restype (*F)(Args...);
	typedef std::tuple<Args...> Tuple;
	typedef typename std::decay<Tuple>::type ttype;
	typedef typename std::tuple_element<parn,Tuple>::type partype;
	F m_func;
	Tuple m_args;
public:
	SingleParam(){}
	SingleParam(SingleParam &S){m_args=S.m_args;m_func=S.m_func;}
	SingleParam(F func,Args... args):m_args(args...){m_func=func;}
	restype operator()(partype x){
		std::get<parn>(m_args)=x;
		return detail::call_impl<
				restype,F, Tuple,
				0 == std::tuple_size<ttype>::value,
				std::tuple_size<ttype>::value
		>::call(m_func, std::forward<Tuple>(m_args));
	}
};
#endif //___SINGLE_PARAM_H