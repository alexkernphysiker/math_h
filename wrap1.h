// Be carefull with using this header file!!!!!!!!!!
#ifdef use_num_type
#ifdef use_indexer_type
namespace FuncWrappers{
	template<int i>
	use_num_type par(use_indexer_type P){return P[i];}
	template<use_num_type (f)(use_num_type),use_num_type(F)(use_indexer_type)>
	use_num_type func(use_indexer_type P){return f(F(P));}
	template<use_num_type (f)(use_num_type,use_num_type),use_num_type(F1)(use_indexer_type),use_num_type(F2)(use_indexer_type)>
	use_num_type func2(use_indexer_type P){return f(F1(P),F2(P));}
	template<use_num_type (f)(use_num_type,use_num_type,use_num_type),use_num_type(F1)(use_indexer_type),use_num_type(F2)(use_indexer_type),use_num_type(F3)(use_indexer_type)>
	use_num_type func3(use_indexer_type P){return f(F1(P),F2(P),F3(P));}
	template<use_num_type(F1)(use_indexer_type),use_num_type(F2)(use_indexer_type)>
	use_num_type add(use_indexer_type P){return F1(P)+F2(P);}
	template<use_num_type(F1)(use_indexer_type),use_num_type(F2)(use_indexer_type)>
	use_num_type sub(use_indexer_type P){return F1(P)-F2(P);}
	template<use_num_type(F1)(use_indexer_type),use_num_type(F2)(use_indexer_type)>
	use_num_type mul(use_indexer_type P){return F1(P)*F2(P);}
	template<use_num_type(F1)(use_indexer_type),use_num_type(F2)(use_indexer_type)>
	use_num_type div(use_indexer_type P){return F1(P)/F2(P);}
	template<use_num_type(F1)(use_indexer_type),use_num_type(F2)(use_indexer_type)>
	use_num_type power(use_indexer_type P){return pow(F1(P),F2(P));}
}
#endif
#endif
