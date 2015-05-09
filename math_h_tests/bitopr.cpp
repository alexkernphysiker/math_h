#include <gtest/gtest.h>
#include <functional>
#include <bit_opr.h>
using namespace std;
using namespace BitOperations;
template<int bn>
inline numtype bit_test(){
	numtype v=bit_test<bn-1>();
	EXPECT_EQ(v,bit<bn>::set);
	EXPECT_EQ(!v,bit<bn>::unset);
	return v*2;
}
template<>inline numtype bit_test<-1>(){return 1;}
template<int bn>
inline numtype bits_test(){
	numtype v=bits_test<bn-1>();
	numtype set=bits_in<bn,0>::set;
	numtype unset=bits_in<bn,0>::unset;
	EXPECT_EQ(v,set);
	EXPECT_EQ(!v,unset);
	return v*2+1;
}
template<>inline numtype bits_test<-1>(){return 1;}
template<int bn>
inline numtype o_bits_test(){
	numtype v=o_bits_test<bn-1>();
	numtype V=occupy_bits<bn,bn>(1);
	EXPECT_EQ(v,V);
	return v*2;
}
template<>inline numtype o_bits_test<-1>(){return 1;}
TEST(bit,BasicTest){bit_test<sizeof(numtype)*8-1>();}
TEST(bits_in,BasicTest){bits_test<sizeof(numtype)*8-1>();}
TEST(occupy_bits,BasicTest){o_bits_test<sizeof(numtype)*8-1>();}
TEST(occupy_bits,Throws){
	auto must_throw=[](){return occupy_bits<7,4>(255);};
	EXPECT_THROW(must_throw(),exception);
	auto must_not_throw=[](){return occupy_bits<7,4>(15);};
	EXPECT_NO_THROW(must_not_throw());
}