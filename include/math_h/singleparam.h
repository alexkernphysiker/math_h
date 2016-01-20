// this file is distributed under 
// MIT license
#ifndef NEZOWVHXCAGGDHQC
#	define NEZOWVHXCAGGDHQC
#include <utility>
#include <tuple>
namespace MathTemplates{
	using namespace std;
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
				>::call(f, forward<Tuple>(t));
			}
		};
		template <
		typename ret,typename F, typename Tuple,
		int Total, int... N
		>
		struct call_impl<ret, F, Tuple, true, Total, N...>{
			static ret call(F f, Tuple && t){
				return f(get<N>(forward<Tuple>(t))...);
			}
		};
	}
	template<typename restype,int parn, typename... Args>
	class SingleParam{
	private:
		typedef restype (*F)(Args...);
		typedef std::tuple<Args...> Tuple;
		typedef typename decay<Tuple>::type ttype;
		typedef typename tuple_element<parn,Tuple>::type partype;
		F m_func;
		Tuple m_args;
	public:
		SingleParam(){}
		SingleParam(SingleParam &S){m_args=S.m_args;m_func=S.m_func;}
		SingleParam(F func,Args... args):m_args(args...){m_func=func;}
		restype operator()(partype x){
			get<parn>(m_args)=x;
			return detail::call_impl<
				restype,F, Tuple,
				0 == tuple_size<ttype>::value,
				tuple_size<ttype>::value
			>::call(m_func, forward<Tuple>(m_args));
		}
	};
};
#endif
