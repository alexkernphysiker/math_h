// this file is distributed under
// LGPLv3 license
#include <gtest/gtest.h>
#include <math_h/vectors.h>
#include <math_h/hists.h>
#include <math_h/sigma.h>
using namespace std;
using namespace MathTemplates;
TEST(Vector, base1d)
{
    
    RandomUniform<> X(-50., 50);
    for (size_t i = 0; i < 100; i++) {
        const auto x = X();
        const auto V = vec(x);
        EXPECT_EQ(x, V.x());
        EXPECT_EQ(x, V.component<1>());
        const auto x2 = X();
        const auto V2 = vec(x2);
        EXPECT_EQ(x2, V2.x());
        EXPECT_EQ(x2, V2.component<1>());
	EXPECT_EQ(x+x2,(V+V2).x());
	EXPECT_EQ(x-x2,(V-V2).x());
	EXPECT_EQ(x*x2,V*V2);
	const auto a=X();
	EXPECT_EQ(x*a,(V*a).x());
	if(a!=0)
	EXPECT_EQ(x/a,(V/a).x());
    }
}
TEST(Vector, base2d)
{
    
    RandomUniform<> X(-50., 50);
    for (size_t i = 0; i < 100; i++) {
        const auto x = X(), y = X();
        const auto V = vec(x, y);
        EXPECT_EQ(x, V.x());
        EXPECT_EQ(x, V.component<1>());
        EXPECT_EQ(y, V.y());
        EXPECT_EQ(y, V.component<2>());
        const auto x2 = X(), y2 = X();
        const auto V2 = vec(x2, y2);
        EXPECT_EQ(x2, V2.x());
        EXPECT_EQ(x2, V2.component<1>());
        EXPECT_EQ(y2, V2.y());
        EXPECT_EQ(y2, V2.component<2>());
	EXPECT_EQ(x+x2,(V+V2).x());
	EXPECT_EQ(x-x2,(V-V2).x());
	EXPECT_EQ(y+y2,(V+V2).y());
	EXPECT_EQ(y-y2,(V-V2).y());
	EXPECT_EQ(x*x2+y*y2,V*V2);
	const auto a=X();
	EXPECT_EQ(x*a,(V*a).x());
	EXPECT_EQ(y*a,(V*a).y());
	if(a!=0){
	    EXPECT_EQ(x/a,(V/a).x());
	    EXPECT_EQ(y/a,(V/a).y());
	}
    }
}
TEST(Vector, base3d)
{
    
    RandomUniform<> X(-50., 50);
    for (size_t i = 0; i < 100; i++) {
        const auto x = X(), y = X(), z = X();
        const auto V = vec(x, y, z);
        EXPECT_EQ(x, V.x());
        EXPECT_EQ(x, V.component<1>());
        EXPECT_EQ(y, V.y());
        EXPECT_EQ(y, V.component<2>());
        EXPECT_EQ(z, V.z());
        EXPECT_EQ(z, V.component<3>());
        const auto x2 = X(), y2 = X(), z2 = X();
        const auto V2 = vec(x2, y2, z2);
        EXPECT_EQ(x2, V2.x());
        EXPECT_EQ(x2, V2.component<1>());
        EXPECT_EQ(y2, V2.y());
        EXPECT_EQ(y2, V2.component<2>());
        EXPECT_EQ(z2, V2.z());
        EXPECT_EQ(z2, V2.component<3>());
	EXPECT_EQ(x+x2,(V+V2).x());
	EXPECT_EQ(x-x2,(V-V2).x());
	EXPECT_EQ(y+y2,(V+V2).y());
	EXPECT_EQ(y-y2,(V-V2).y());
	EXPECT_EQ(z+z2,(V+V2).z());
	EXPECT_EQ(z-z2,(V-V2).z());
	EXPECT_EQ(x*x2+y*y2+z*z2,V*V2);
	const auto a=X();
	EXPECT_EQ(x*a,(V*a).x());
	EXPECT_EQ(y*a,(V*a).y());
	EXPECT_EQ(z*a,(V*a).z());
	if(a!=0){
	    EXPECT_EQ(x/a,(V/a).x());
	    EXPECT_EQ(y/a,(V/a).y());
	    EXPECT_EQ(z/a,(V/a).z());
	}
    }
}
TEST(Vector, base4d)
{
    
    RandomUniform<> X(-50., 50);
    for (size_t i = 0; i < 100; i++) {
        const auto x = X(), y = X(), z = X(), zz = X();
        const auto V = vec(x, y, z, zz);
        EXPECT_EQ(x, V.x());
        EXPECT_EQ(x, V.component<1>());
        EXPECT_EQ(y, V.y());
        EXPECT_EQ(y, V.component<2>());
        EXPECT_EQ(z, V.z());
        EXPECT_EQ(z, V.component<3>());
        EXPECT_EQ(zz, V.component<4>());
        const auto x2 = X(), y2 = X(), z2 = X(), zz2 = X();
        const auto V2 = vec(x2, y2, z2, zz2);
        EXPECT_EQ(x2, V2.x());
        EXPECT_EQ(x2, V2.component<1>());
        EXPECT_EQ(y2, V2.y());
        EXPECT_EQ(y2, V2.component<2>());
        EXPECT_EQ(z2, V2.z());
        EXPECT_EQ(z2, V2.component<3>());
        EXPECT_EQ(zz2, V2.component<4>());
	EXPECT_EQ(x+x2,(V+V2).x());
	EXPECT_EQ(x-x2,(V-V2).x());
	EXPECT_EQ(y+y2,(V+V2).y());
	EXPECT_EQ(y-y2,(V-V2).y());
	EXPECT_EQ(z+z2,(V+V2).z());
	EXPECT_EQ(z-z2,(V-V2).z());
	EXPECT_EQ(zz+zz2,(V+V2).component<4>());
	EXPECT_EQ(zz-zz2,(V-V2).component<4>());
	EXPECT_EQ(x*x2+y*y2+z*z2+zz*zz2,V*V2);
	const auto a=X();
	EXPECT_EQ(x*a,(V*a).x());
	EXPECT_EQ(y*a,(V*a).y());
	EXPECT_EQ(z*a,(V*a).z());
	EXPECT_EQ(zz*a,(V*a).component<4>());
	if(a!=0){
	    EXPECT_EQ(x/a,(V/a).x());
	    EXPECT_EQ(y/a,(V/a).y());
	    EXPECT_EQ(z/a,(V/a).z());
	    EXPECT_EQ(zz/a,(V/a).component<4>());
	}
    }
}
TEST(Vector, scalar_prod2)
{
    EXPECT_EQ(0, x()*y());
    EXPECT_EQ(0, y()*x());
    EXPECT_EQ(1, x()*x());
    EXPECT_EQ(1, y()*y());

    EXPECT_EQ(0, x() * (-y()));
    EXPECT_EQ(0, y() * (-x()));
    EXPECT_EQ(0, (-x())*y());
    EXPECT_EQ(0, (-y())*x());

    EXPECT_EQ(-1, x() * (-x()));
    EXPECT_EQ(-1, y() * (-y()));
    EXPECT_EQ(-1, (-x())*x());
    EXPECT_EQ(-1, (-y())*y());
}
TEST(Vector, scalar_prod3)
{
    EXPECT_EQ(0, X()*Y());
    EXPECT_EQ(0, Y()*Z());
    EXPECT_EQ(0, X()*Z());
    EXPECT_EQ(1, X()*X());
    EXPECT_EQ(1, Y()*Y());
    EXPECT_EQ(1, Z()*Z());

    EXPECT_EQ(0, X() * (-Y()));
    EXPECT_EQ(0, (-X())*Z());

    EXPECT_EQ(-1, (-X())*X());
    EXPECT_EQ(-1, Y() * (-Y()));
    EXPECT_EQ(1, (-Z()) * (-Z()));
}

#ifdef ____optimized_version_of_vectors_h_____
TEST(Vector, remove_component)
{
    EXPECT_EQ(vec(1,2).RemoveComponent<1>(),vec(2));
    EXPECT_EQ(vec(1,2).RemoveComponent<2>(),vec(1));

    EXPECT_EQ(vec(1,2,3).RemoveComponent<1>(),vec(2,3));
    EXPECT_EQ(vec(1,2,3).RemoveComponent<2>(),vec(1,3));
    EXPECT_EQ(vec(1,2,3).RemoveComponent<3>(),vec(1,2));

    EXPECT_EQ(vec(1,2,3,4).RemoveComponent<1>(),vec(2,3,4));
    EXPECT_EQ(vec(1,2,3,4).RemoveComponent<2>(),vec(1,3,4));
    EXPECT_EQ(vec(1,2,3,4).RemoveComponent<3>(),vec(1,2,4));
    EXPECT_EQ(vec(1,2,3,4).RemoveComponent<4>(),vec(1,2,3));

}
TEST(Vector, insert_component)
{
    EXPECT_EQ(vec(1).InsertComponent<1>(0),vec(0,1));
    EXPECT_EQ(vec(1).InsertComponent<2>(0),vec(1,0));

    EXPECT_EQ(vec(1,2).InsertComponent<1>(0),vec(0,1,2));
    EXPECT_EQ(vec(1,2).InsertComponent<2>(0),vec(1,0,2));
    EXPECT_EQ(vec(1,2).InsertComponent<3>(0),vec(1,2,0));

    EXPECT_EQ(vec(1,2,3).InsertComponent<1>(0),vec(0,1,2,3));
    EXPECT_EQ(vec(1,2,3).InsertComponent<2>(0),vec(1,0,2,3));
    EXPECT_EQ(vec(1,2,3).InsertComponent<3>(0),vec(1,2,0,3));
    EXPECT_EQ(vec(1,2,3).InsertComponent<4>(0),vec(1,2,3,0));
}
#endif
