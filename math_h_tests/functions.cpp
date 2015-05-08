#include <gtest/gtest.h>
#include <functional>
#include <functions.h>
using namespace std;
template<class numt>
void Test_Peak_Shape(function<numt(numt)> f,numt from,numt peak,numt to){
	for(numt x=from;x<=to;x+=0.05){
		ASSERT_TRUE(f(x)>=0);
		if(x<=peak)
			ASSERT_TRUE(f(x)>=f(x-0.01));
		else
			ASSERT_TRUE(f(x)>=f(x+0.01));
	}
}
TEST(Gaussian,Shape){
	for(double sigma=0.5;sigma<10;sigma+=0.5)
		for(double X=-5;X<=5;X+=1)
			Test_Peak_Shape<double>([sigma,X](double x){return Gaussian(x,X,sigma);},X-sigma*5,X,X+sigma*5);
}
TEST(BreitWigner,Shape){
	for(double sigma=0.5;sigma<10;sigma+=0.5)
		for(double X=-5;X<=5;X+=1)
			Test_Peak_Shape<double>([sigma,X](double x){return BreitWigner(x,X,sigma);},X-sigma*2,X,X+sigma*2);
}
TEST(Novosibirsk,Shape){
	for(double sigma=0.5;sigma<10;sigma+=0.5)
		for(double X=-5;X<=5;X+=1)for(double asym=-1;asym<=1;asym+=0.1)
			Test_Peak_Shape<double>([sigma,X,asym](double x){return Novosibirsk<>(x,X,sigma,asym);},X-sigma*2,X,X+sigma*2);
}
template<class numt>
void Test_stair_Shape(function<numt(numt)> f,numt from,numt middle,numt to){
	numt a=f(from);
	numt b=f(to);
	numt m=f(middle);
	ASSERT_TRUE(a!=b);
	ASSERT_TRUE(pow((m-a)*(b-a)-0.5,2)<0.4);
	for(numt x=from;x<=to;x+=0.05){
		ASSERT_TRUE((f(x)-f(x-0.01))*(a-b)<=0);
	}
}
TEST(FermiFunc,Shape){
	for(double diffuse=0.5;diffuse<10;diffuse+=0.5)
		for(double X=-5;X<=5;X+=1)
			Test_stair_Shape<double>([diffuse,X](double x){return FermiFunc(x,X,diffuse);},X-diffuse*5,X,X+diffuse*5);
}

TEST(Polynom,Base){
	for(int p=0;p<10;p++){
		double C[p+1];
		for(int i=0;i<p;i++)C[i]=0;
		C[p]=1;
		for(double x=-2;x<=2;x+=0.01)for(int P=p;P>=0;P--){
			EXPECT_EQ(true,pow(pow(x,P)-Polynom(x,C,P,p-P),2)<0.01);
		}
	}
}
TEST(Polynom,Extended){
	double C[9+1];
	for (int i=0;i<=9;i++)C[i]=1;
	for(double x=-2;x<=2;x+=0.01){
		double V=0;
		for(int p=0;p<=9;p++){
			V+=pow(x,p);
			EXPECT_EQ(true,pow(V-Polynom(x,C,p),2)<0.01);
		}
	}
}
