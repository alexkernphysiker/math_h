// https://github.com/alexkernphysiker/math_h
// Be carefull with using this header file!!!!!!!!!!
// It requires two defined expressions and can be included more than once!
#ifdef use_num_type
#ifdef use_indexer_type
namespace wrap_func_indexer{
	template<int i>inline use_num_type par(use_indexer_type P){return P[i];}
#define ___p_decl___ use_indexer_type
#define ___p_decl2___ use_indexer_type P
#define ___p_ P
#include "wrap.cc"
#undef ___p_
#undef ___p_decl2___
#undef ___p_decl___
}
#endif
#endif
