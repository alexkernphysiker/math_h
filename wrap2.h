// Be carefull with using this header file!!!!!!!!!!
#ifdef use_num_type
#ifdef use_indexer_type
#include "functions.h"
namespace FuncWrappers2{
	template<int x_ind>
	use_num_type arg(use_indexer_type X,use_indexer_type){return X[x_ind];}
	template<int p_ind>
	use_num_type par(use_indexer_type,use_indexer_type P){return P[p_ind];}
	template<use_num_type (f)(use_num_type),use_num_type(F)(use_indexer_type,use_indexer_type)>
	use_num_type func(use_indexer_type X,use_indexer_type P){return f(F(X,P));}
	template<use_num_type (f)(use_num_type,use_num_type),use_num_type(F1)(use_indexer_type,use_indexer_type),use_num_type(F2)(use_indexer_type,use_indexer_type)>
	use_num_type func2(use_indexer_type X,use_indexer_type P){return f(F1(X,P),F2(X,P));}
	template<use_num_type (f)(use_num_type,use_num_type,use_num_type),use_num_type(F1)(use_indexer_type,use_indexer_type),use_num_type(F2)(use_indexer_type,use_indexer_type),use_num_type(F3)(use_indexer_type,use_indexer_type)>
	use_num_type func3(use_indexer_type X,use_indexer_type P){return f(F1(X,P),F2(X,P),F3(X,P));}
	template<use_num_type(F1)(use_indexer_type,use_indexer_type),use_num_type(F2)(use_indexer_type,use_indexer_type)>
	use_num_type add(use_indexer_type X,use_indexer_type P){return F1(X,P)+F2(X,P);}
	template<use_num_type(F1)(use_indexer_type,use_indexer_type),use_num_type(F2)(use_indexer_type,use_indexer_type)>
	use_num_type sub(use_indexer_type X,use_indexer_type P){return F1(X,P)-F2(X,P);}
	template<use_num_type(F1)(use_indexer_type,use_indexer_type),use_num_type(F2)(use_indexer_type,use_indexer_type)>
	use_num_type mul(use_indexer_type X,use_indexer_type P){return F1(X,P)*F2(X,P);}
	template<use_num_type(F1)(use_indexer_type,use_indexer_type),use_num_type(F2)(use_indexer_type,use_indexer_type)>
	use_num_type div(use_indexer_type X,use_indexer_type P){return F1(X,P)/F2(X,P);}
	template<use_num_type(F1)(use_indexer_type,use_indexer_type),use_num_type(F2)(use_indexer_type,use_indexer_type)>
	use_num_type power(use_indexer_type X,use_indexer_type P){return pow(F1(X,P),F2(X,P));}
	template<int x_ind, int par_ind, int p>
	use_num_type polynom(use_indexer_type X,use_indexer_type P){return Polynom(X[x_ind],P,p,par_ind);}
}
#endif
#endif
