// https://github.com/alexkernphysiker/math_h
// Be carefull with using this header file!!!!!!!!!!
// It requires two defined expressions and can be included more than once!
#ifdef use_num_type
#ifdef use_indexer_type
#include "functions.h"
namespace wrap_func_indexer_indexer{
	template<int x_ind>inline use_num_type arg(use_indexer_type X,use_indexer_type){return X[x_ind];}
	template<int p_ind>inline use_num_type par(use_indexer_type,use_indexer_type P){return P[p_ind];}
#define ___p_decl___ use_indexer_type,use_indexer_type
#define ___p_decl2___ use_indexer_type X,use_indexer_type P
#define ___p_ X,P
#include "wrap.cc"
#undef ___p_
#undef ___p_decl2___
#undef ___p_decl___
	template<int x_ind, int par_ind, int p>
	use_num_type polynom(use_indexer_type X,use_indexer_type P){return Polynom(X[x_ind],P,p,par_ind);}
}
#endif
#endif
