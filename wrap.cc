#ifdef use_num_type
#ifdef use_indexer_type
#ifdef ___p_decl___
#ifdef ___p_decl2___
#ifdef ___p_
//wrap global variable
template<use_num_type& value>
inline use_num_type var(___p_decl2___){return value;}
//wrap functions requiring use_num_type as parameters
template<use_num_type (f)(use_num_type),use_num_type(F)(___p_decl___)>
inline use_num_type func(___p_decl2___){return f(F(___p_));}
template<use_num_type (f)(use_num_type,use_num_type),use_num_type(F1)(___p_decl___),use_num_type(F2)(___p_decl___)>
inline use_num_type func2(___p_decl2___){return f(F1(___p_),F2(___p_));}
template<use_num_type (f)(use_num_type,use_num_type,use_num_type),use_num_type(F1)(___p_decl___),use_num_type(F2)(___p_decl___),use_num_type(F3)(___p_decl___)>
inline use_num_type func3(___p_decl2___){return f(F1(___p_),F2(___p_),F3(___p_));}
template<use_num_type (f)(use_num_type,use_num_type,use_num_type,use_num_type),use_num_type(F1)(___p_decl___),use_num_type(F2)(___p_decl___),use_num_type(F3)(___p_decl___),use_num_type(F4)(___p_decl___)>
inline use_num_type func4(___p_decl2___){return f(F1(___p_),F2(___p_),F3(___p_),F4(___p_));}
template<use_num_type (f)(use_num_type,use_num_type,use_num_type,use_num_type,use_num_type),use_num_type(F1)(___p_decl___),use_num_type(F2)(___p_decl___),use_num_type(F3)(___p_decl___),use_num_type(F4)(___p_decl___),use_num_type(F5)(___p_decl___)>
inline use_num_type func5(___p_decl2___){return f(F1(___p_),F2(___p_),F3(___p_),F4(___p_),F5(___p_));}
// combine functions
template<use_num_type(F)(___p_decl___)>
inline use_num_type minus(___p_decl2___){return -F(___p_);}
template<use_num_type(F1)(___p_decl___),use_num_type(F2)(___p_decl___)>
inline use_num_type add(___p_decl2___){return F1(___p_)+F2(___p_);}
template<use_num_type(F1)(___p_decl___),use_num_type(F2)(___p_decl___)>
inline use_num_type sub(___p_decl2___){return F1(___p_)-F2(___p_);}
template<use_num_type(F1)(___p_decl___),use_num_type(F2)(___p_decl___)>
inline use_num_type mul(___p_decl2___){return F1(___p_)*F2(___p_);}
template<use_num_type(F1)(___p_decl___),use_num_type(F2)(___p_decl___)>
inline use_num_type div(___p_decl2___){return F1(___p_)/F2(___p_);}
template<use_num_type(F1)(___p_decl___),use_num_type(F2)(___p_decl___)>
inline use_num_type power(___p_decl2___){return pow(F1(___p_),F2(___p_));}
//we also need to provide constants
namespace details{
template<unsigned char digits=0> struct movecoef{enum{value=movecoef<digits-1>::value * 10};};
template<> struct movecoef<0>{enum{value=1};};
};
template<int basic,unsigned char move_point=0>inline use_num_type num(___p_decl2___){return use_num_type(basic)/use_num_type(details::movecoef<move_point>::value);}
template<use_num_type(F)(___p_decl___)>
use_num_type wrap(___p_decl2___){return F(___p_);}
#endif
#endif
#endif
#endif
#endif
