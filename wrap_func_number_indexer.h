// https://github.com/alexkernphysiker/math_h
// Be carefull with using this header file!!!!!!!!!!
// It requires two defined expressions and can be included more than once!
#ifdef use_num_type
#ifdef use_indexer_type
#include "functions.h"
namespace wrap_func_number_indexer{
	inline use_num_type arg(use_num_type x,use_indexer_type){return x;}
	template<int p_ind>inline use_num_type par(use_num_type,use_indexer_type P){return P[p_ind];}
#define ___p_decl___ use_num_type,use_indexer_type
#define ___p_decl2___ use_num_type x,use_indexer_type P
#define ___p_ x,P
#include "wrap.cc"
#undef ___p_
#undef ___p_decl2___
#undef ___p_decl___
	template<int par_ind, int p>
	use_num_type polynom(use_num_type x,use_indexer_type P){return Polynom(x,P,p,par_ind);}
}
#endif
#endif
