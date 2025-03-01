// this file is distributed under
// LGPLv3 license
#ifndef _______NUMERIC_DIFF___H_________
#	define _______NUMERIC_DIFF___H_________
#include "error.h"
#include "functions.h"
#include "vectors.h"
namespace MathTemplates
{
	template<class numt>
	inline auto num_der1(const numt& delta) {
		return Opr<FunctionWrap<numt, const numt&>, FunctionWrap<numt, const numt&>>(
			[delta](const FunctionWrap<numt, const numt&>& F) {
				return FunctionWrap<numt, const numt&>([delta, &F](const numt& x) {
					const numt x1 = x - delta, x2 = x + delta;
					const numt y1 = F(x1), y2 = F(x2);
					return (y2 - y1) / (x2 - x1);
					});
			}
		);
	}
	template<class numt>
	inline auto num_der2(const numt& delta) {
		return Opr<FunctionWrap<numt, const numt&>, FunctionWrap<numt, const numt&>>(
			[delta](const FunctionWrap<numt, const numt&>& F) {
				return FunctionWrap<numt, const numt&>([delta, &F](const numt& x) {
					const numt x1 = x - delta, x2 = x, x3 = x + delta;
					const numt y1 = F(x1), y2 = F(x2), y3 = F(x3);
					return (y1 - y2 * numt(2) + y3) / (delta * delta);
					});
			}
		);
	}
	template<size_t index, size_t size, class numt>
	inline auto num_Pder1(const numt& delta) {
		return Opr<FunctionWrap<numt, const Vector<size, numt>&>, FunctionWrap<numt, const Vector<size, numt>&>>(
			[delta](const FunctionWrap<numt, const Vector<size, numt>&>& F) {
				return FunctionWrap<numt, const Vector<size, numt>&>([delta, F](const Vector<size, numt>& x) {
					const auto D = Vector<size, numt>::template basis_vector<index>()* delta;
					const auto x1 = x - D, x2 = x + D;
					const numt y1 = F(x1), y2 = F(x2);
					return (y2 - y1) / (delta * numt(2));
					});
			}
		);
	}
	namespace numeric_diff_details {
		template<size_t size, class numt, size_t index>
		Vector<index, Opr<FunctionWrap<numt, const Vector<size, numt>&>, FunctionWrap<numt, const Vector<size, numt>&>>>
			___nabla_partial(const numt& delta) {
			if constexpr (index == 1) return vec(num_Pder1<index, size, numt>(delta));
			else return ___nabla_partial<size, numt, index - 1>(delta).template InsertComponent<index>(num_Pder1<index, size, numt>(delta));
		}
	}
	template<size_t size, class numt>
	inline auto nabla(const numt& delta) {
		static_assert(size > 0, "invalid space dimmensions count");
		return numeric_diff_details::___nabla_partial<size, numt, size>(delta);
	}
	template<size_t size, class numt>
	inline auto laplace(const numt& delta) {
		static_assert(size > 0, "invalid space dimmensions count");
		return nabla<size, numt>(delta)* nabla<size, numt>(delta);
	}
};
#endif
