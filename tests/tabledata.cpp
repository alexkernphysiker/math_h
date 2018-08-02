// this file is distributed under
// LGPLv3 license
#include <gtest/gtest.h>
#include <math_h/tabledata.h>
#include <math_h/sigma.h>
#include <math_h/randomfunc.h>
//this file contains unit tests for tabledata.h
using namespace std;
using namespace MathTemplates;
double TestArray[] = { -2, -1, 0, 1, 2, 3, 4, 5, 6, 7};
TEST(WhereToInsert, BorderConditions)
{
    for (int beg = 0; beg < 10; beg++)
        for (double V = TestArray[beg] - 2; V <= TestArray[beg] + 2; V += 0.5) {
            int index = 0;
            if (beg > 0) {
                index = table_data_details::WhereToInsert(beg, beg - 1, TestArray, V);
                EXPECT_EQ(beg, index);
            }
            index = table_data_details::WhereToInsert(beg, beg, TestArray, V);
            if (V < TestArray[beg]) {
                EXPECT_EQ(beg, index);
            }
            if (V > TestArray[beg]) {
                EXPECT_EQ(beg + 1, index);
            }
            if (V == TestArray[beg]) {
                EXPECT_TRUE((index == beg) || (index == (beg + 1)));
            }
        }
}
TEST(WhereToInsert, NormalConditions)
{
    for (double x = TestArray[0] - 0.5; x <= TestArray[9] + 0.5; x += 0.5)
        for (int beg = 0; beg < 10; beg++)for (int end = beg + 1; end < 10; end++) {
                int index = table_data_details::WhereToInsert(beg, end, TestArray, x);
                if (x < TestArray[beg])EXPECT_EQ(beg, index);
                else if (x > TestArray[end])EXPECT_EQ(end + 1, index);
                else EXPECT_TRUE((x >= TestArray[index - 1]) && (x <= TestArray[index]));
            }
}
TEST(InsertSorted, BasicTest)
{
    Chain<int> X;
    for (int i = 0; i < 50; i++) {
        table_data_details::InsertSorted(rand() % 10, X,std_size(X),std_insert(X,int));
        for (int j = 0; j < i; j++)
            EXPECT_TRUE(X[j] <= X[j + 1]);
    }
}
TEST(SortedChain, basetest)
{
    SortedChain<> chain;
    RandomUniform<> get(-10, 10);
    for (size_t i = 0; i < 100; i++)
        chain << get();
    for (size_t i = 1; i < 100; i++)
        EXPECT_TRUE(chain[i] > chain[i - 1]);
}
TEST(point, test2d){
{
    double x(25), y(78);
    point<double> p(x, y);
    EXPECT_EQ(x, p.X());
    EXPECT_EQ(y, p.Y());
}
{
    double x(25), y(78);
    point<> p = {x, y};
    EXPECT_EQ(x, p.X());
    EXPECT_EQ(y, p.Y());
}
{
    int x(25);
    double y(78);
    point<int, double> p = make_point(25, 78);
    EXPECT_EQ(x, p.X());
    EXPECT_EQ(y, p.Y());
    point<int, double> p2 = make_point(25., 78.);
    EXPECT_EQ(x, p2.X());
    EXPECT_EQ(y, p2.Y());
}
}
TEST(point, test3d)
{
    double x(25), y(78), z(55);
    point3d<double> p(x, y, z);
    EXPECT_EQ(x, p.X());
    EXPECT_EQ(y, p.Y());
    EXPECT_EQ(z, p.Z());
}
TEST(point,test_value)
{
    value<double> x(25), y(78);
    point<value<double>> p(x, y);
    EXPECT_EQ(x.val(), p.X().val());
    EXPECT_EQ(x.uncertainty(), p.X().uncertainty());
    EXPECT_EQ(y.val(), p.Y().val());
    EXPECT_EQ(y.uncertainty(), p.Y().uncertainty());
}
TEST(point, test_value_3d)
{
    value<double> x(25), y(78), z(55);
    point3d<value<double>> p(x, y, z);
    EXPECT_EQ(x.val(), p.X().val());
    EXPECT_EQ(x.uncertainty(), p.X().uncertainty());
    EXPECT_EQ(y.val(), p.Y().val());
    EXPECT_EQ(y.uncertainty(), p.Y().uncertainty());
    EXPECT_EQ(z.val(), p.Z().val());
    EXPECT_EQ(z.uncertainty(), p.Z().uncertainty());
}
TEST(SortedPoints, size){
{
    SortedPoints<double> chain;
    EXPECT_EQ(0, chain.size());
    chain << point<double>(0, 1);
    EXPECT_EQ(1, chain.size());
    chain << point<double>(2, 1);
    EXPECT_EQ(2, chain.size());
    chain << point<double>(1, 1);
    EXPECT_EQ(3, chain.size());
    chain << point<double>(3, 1);
    EXPECT_EQ(4, chain.size());
    chain << point<double>(5, 1);
    EXPECT_EQ(5, chain.size());
    chain << point<double>(4, 1);
    EXPECT_EQ(6, chain.size());
    EXPECT_EQ(0, chain.left().X());
    EXPECT_EQ(5, chain.right().X());
    for (size_t i = 0; i < chain.size(); i++) {
        EXPECT_EQ(1, chain[i].Y());
        EXPECT_EQ(i, chain[i].X());
    }
}
{
    SortedPoints<value<>> chain;
    EXPECT_EQ(0, chain.size());
    chain << make_point(0, 1);
    EXPECT_EQ(1, chain.size());
    chain << make_point(2, 1);
    EXPECT_EQ(2, chain.size());
    chain << make_point(1, 1);
    EXPECT_EQ(3, chain.size());
    chain << make_point(3, 1);
    EXPECT_EQ(4, chain.size());
    chain << make_point(5, 1);
    EXPECT_EQ(5, chain.size());
    chain << make_point(4, 1);
    EXPECT_EQ(6, chain.size());
    EXPECT_EQ(0, chain.left().X().val());
    EXPECT_EQ(5, chain.right().X().val());
    for (size_t i = 0; i < chain.size(); i++) {
        EXPECT_EQ(1, chain[i].Y().val());
        EXPECT_EQ(i, chain[i].X().val());
    }
}
}
TEST(SortedPoints, range)
{
    SortedPoints<> chain=Points<>{{0, 1},{2, 1.1},{3, 1.2},{4, 1.3},{5, 1.4}};
    auto xcut = chain.XRange(0.5, 4.5);
    EXPECT_EQ(3, xcut.size());
    EXPECT_EQ(2, xcut.left().X());
    EXPECT_EQ(4, xcut.right().X());
    auto ycut = chain.YRange(1.15, 1.35).TransponateAndSort();
    EXPECT_EQ(2, ycut.size());
    EXPECT_EQ(1.2, ycut.left().X());
    EXPECT_EQ(1.3, ycut.right().X());
}
TEST(SortedPoints, WeightedAvr){
{
    SortedPoints<double,WeightedAverage<>> chain;
    chain
    << make_point(0, WeightedAverage<>())
    << make_point(1, WeightedAverage<>())
    << make_point(2, WeightedAverage<>());
    for (size_t i = 0; i < chain.size(); i++) {
        EXPECT_EQ(i, chain[i].X());
        EXPECT_THROW(chain[i].Y().val(),Exception<WeightedAverage<>>);
    }
    SortedPoints<double,value<>> opr1=Points<double,value<>>{{0,{0,1}},{1,{1,1}},{2,{2,1}}};
    chain.leftArrow(opr1);
    for (size_t i = 0; i < chain.size(); i++) {
        EXPECT_EQ(i, chain[i].X());
        EXPECT_EQ(i, chain[i].Y().val());
    }
    chain.leftArrow(Points<double,value<>>{{0,{0,1}},{1,{0,1}},{2,{0,1}}});
    for (size_t i = 0; i < chain.size(); i++) {
        EXPECT_EQ(i, chain[i].X());
        EXPECT_EQ(double(i)/2, chain[i].Y().val());
    }
#ifdef ____full_version_of_tabledata_h_____
    const auto chain2=chain+2;
    for (size_t i = 0; i < chain2.size(); i++) {
        EXPECT_EQ(i, chain2[i].X());
        EXPECT_EQ(double(i)/2 + 2, chain2[i].Y().val());
    }
    const auto chain3=chain-2;
    for (size_t i = 0; i < chain3.size(); i++) {
        EXPECT_EQ(i, chain3[i].X());
        EXPECT_EQ(double(i)/2 - 2, chain3[i].Y().val());
    }
    const auto chain4=chain*2;
    for (size_t i = 0; i < chain4.size(); i++) {
        EXPECT_EQ(i, chain4[i].X());
        EXPECT_EQ(double(i), chain4[i].Y().val());
    }
    const auto chain5=chain/2;
    for (size_t i = 0; i < chain5.size(); i++) {
        EXPECT_EQ(i, chain5[i].X());
        EXPECT_EQ(double(i)/4, chain5[i].Y().val());
    }
#endif
}
{
    SortedPoints<double,WeightedAverage<>> chain=Points<double,value<>>{{0,{0,1}},{1,{1,1}},{2,{2,1}}};
    for (size_t i = 0; i < chain.size(); i++) {
        EXPECT_EQ(i, chain[i].X());
        EXPECT_EQ(i, chain[i].Y().val());
    }
    chain.leftArrow(Points<double,value<>>{{0,{0,1}},{1,{0,1}},{2,{0,1}}});
    for (size_t i = 0; i < chain.size(); i++) {
        EXPECT_EQ(i, chain[i].X());
        EXPECT_EQ(double(i)/2, chain[i].Y().val());
    }
}
}
TEST(hist, scale_norm)
{
    hist<> H([](const value<>&x)->const value<>{return 10.0 + 4.0 * sin(x.val());},BinsByCount(10, 0.0, 1.0));
    double s1 = 0;
    for (const auto &p : H) s1+= p.Y().val();
    double s2 = 0;
    auto H2 = H.Scale(2);
    for (const auto &p : H2)s2 += p.Y().val();
    EXPECT_TRUE(pow(s1 - s2, 2) < 0.000001);
    for (size_t x = 0; x < H2.size(); x++) {
        EXPECT_EQ(sqrt(H2[x].Y().val()), H2[x].Y().uncertainty());
        EXPECT_EQ(H2[x].Y().val(), H[x * 2].Y().val() + H[x * 2 + 1].Y().val());
    }
}
TEST(hist, scale_norm_2d)
{
    hist2d<> H(BinsByCount(10, 0.0, 1.0), BinsByCount(10, 0.0, 1.0));
    H.FullCycleVar([](const value<> &x, const value<> &y, value<> &z) {
        z = value<>(10.0 + 4.0 * sin(x.val() + y.val()));
    });
    double s1 = 0;
    H.FullCycle([&s1](const value<> &, const value<> &, const value<> &z) {
        s1 += z.val();
    });
    double s2 = 0;
    auto H2 = H.Scale(2, 2);
    H2.FullCycle([&s2](const value<> &, const value<> &, const value<> &z) {
        s2 += z.val();
    });
    EXPECT_TRUE(pow(s1 - s2, 2) < 0.000001);
    for (size_t x = 0; x < H2.X().size(); x++)
        for (size_t y = 0; y < H2.Y().size(); y++) {
            EXPECT_EQ(sqrt(H2[x][y].val()), H2[x][y].uncertainty());
            EXPECT_EQ(H2[x][y].val(), H[x * 2][y * 2].val() + H[x * 2 + 1][y * 2].val() + H[x * 2][y * 2 + 1].val() + H[x * 2 + 1][y * 2 + 1].val());
        }
}
TEST(hist, cut2d)
{
    hist2d<> H(BinsByCount(10, 0.0, 1.0), BinsByCount(20, 0.0, 1.0));
    H.FullCycleVar([](const value<> &x, const value<> &y, value<> &z) {
        z = value<>(10.0 + 4.0 * sin(x.val() + y.val()));
    });
    for (size_t x = 0; x < H.X().size(); x++) {
        const auto Y = H.CutY(x);
        EXPECT_EQ(Y.size(), H.Y().size());
        for (size_t y = 0; y < H.Y().size(); y++)
            EXPECT_EQ(H.Bin(x, y), Y[y].Y());
    }
    for (size_t y = 0; y < H.Y().size(); y++) {
        const auto X = H.CutX(y);
        EXPECT_EQ(X.size(), H.X().size());
        for (size_t x = 0; x < H.X().size(); x++)
            EXPECT_EQ(H.Bin(x, y), X[x].Y());
    }
}
TEST(hist,Distribution1D)
{
    Distribution1D<> D=Chain<value<>>{{-1.0, 0.5},{-0.0, 0.5},{1.0, 0.5}};
    ASSERT_EQ(3, D.size());
    ASSERT_EQ(-1, D[0].X().val());
    ASSERT_EQ(0, D[1].X().val());
    ASSERT_EQ(1, D[2].X().val());
    ASSERT_EQ(0.5, D[0].X().uncertainty());
    ASSERT_EQ(0.5, D[1].X().uncertainty());
    ASSERT_EQ(0.5, D[2].X().uncertainty());
    EXPECT_EQ(0, D[0].Y().val());
    EXPECT_EQ(1, D[0].Y().uncertainty());
    EXPECT_EQ(0, D[1].Y().val());
    EXPECT_EQ(1, D[1].Y().uncertainty());
    EXPECT_EQ(0, D[2].Y().val());
    EXPECT_EQ(1, D[2].Y().uncertainty());
    D.Fill(0.0);
    EXPECT_EQ(0, D[0].Y().val());
    EXPECT_EQ(1, D[0].Y().uncertainty());
    EXPECT_EQ(1, D[1].Y().val());
    EXPECT_EQ(1, D[1].Y().uncertainty());
    EXPECT_EQ(0, D[2].Y().val());
    EXPECT_EQ(1, D[2].Y().uncertainty());
    D.Fill(1.0);
    EXPECT_EQ(0, D[0].Y().val());
    EXPECT_EQ(1, D[0].Y().uncertainty());
    EXPECT_EQ(1, D[1].Y().val());
    EXPECT_EQ(1, D[1].Y().uncertainty());
    EXPECT_EQ(1, D[2].Y().val());
    EXPECT_EQ(1, D[2].Y().uncertainty());
    D.Fill(-0.7);
    EXPECT_EQ(1, D[0].Y().val());
    EXPECT_EQ(1, D[0].Y().uncertainty());
    EXPECT_EQ(1, D[1].Y().val());
    EXPECT_EQ(1, D[1].Y().uncertainty());
    EXPECT_EQ(1, D[2].Y().val());
    EXPECT_EQ(1, D[2].Y().uncertainty());
    D.Fill(-0.2);
    EXPECT_EQ(1, D[0].Y().val());
    EXPECT_EQ(1, D[0].Y().uncertainty());
    EXPECT_EQ(2, D[1].Y().val());
    EXPECT_EQ(sqrt(2.0), D[1].Y().uncertainty());
    EXPECT_EQ(1, D[2].Y().val());
    EXPECT_EQ(1, D[2].Y().uncertainty());
}
TEST(hist,Distribution2D)
{
    Distribution2D<> D(
    {value<>(-1.0, 0.5), value<>(-0.0, 0.5), value<>(1.0, 0.5)},
    {value<>(-0.0, 0.5), value<>(1.0, 0.5)}
    );
    ASSERT_EQ(3, D.X().size());
    ASSERT_EQ(2, D.Y().size());
    ASSERT_EQ(3, D.size());
    ASSERT_EQ(2, D[0].size());
    ASSERT_EQ(2, D[1].size());
    ASSERT_EQ(2, D[2].size());
    EXPECT_EQ(0, D[0][0].val());
    EXPECT_EQ(0, D[0][1].val());
    EXPECT_EQ(0, D[1][0].val());
    EXPECT_EQ(0, D[1][1].val());
    EXPECT_EQ(0, D[2][0].val());
    EXPECT_EQ(0, D[2][1].val());
    D.Fill(0.0, 0.0).Fill(0.1, 0.2).Fill(0.9, 1.1).Fill(-0.9, 1.1).Fill(0.6, 0.6).Fill(0.6, 0.4);
    EXPECT_EQ(0, D[0][0].val());
    EXPECT_EQ(1, D[0][1].val());
    EXPECT_EQ(2, D[1][0].val());
    EXPECT_EQ(0, D[1][1].val());
    EXPECT_EQ(1, D[2][0].val());
    EXPECT_EQ(2, D[2][1].val());
    Chain<point3d<value<>>> dbg;
    D.FullCycle([&dbg](const point3d<value<>> &p) {
        dbg.push_back(p);
    });
    ASSERT_EQ(6, dbg.size());
    EXPECT_EQ(-1, dbg[0].X().val());
    EXPECT_EQ(0, dbg[0].Y().val());
    EXPECT_EQ(0, dbg[0].Z().val());
    EXPECT_EQ(-1, dbg[1].X().val());
    EXPECT_EQ(1, dbg[1].Y().val());
    EXPECT_EQ(1, dbg[1].Z().val());
    EXPECT_EQ(0, dbg[2].X().val());
    EXPECT_EQ(0, dbg[2].Y().val());
    EXPECT_EQ(2, dbg[2].Z().val());
    EXPECT_EQ(0, dbg[3].X().val());
    EXPECT_EQ(1, dbg[3].Y().val());
    EXPECT_EQ(0, dbg[3].Z().val());
    EXPECT_EQ(1, dbg[4].X().val());
    EXPECT_EQ(0, dbg[4].Y().val());
    EXPECT_EQ(1, dbg[4].Z().val());
    EXPECT_EQ(1, dbg[5].X().val());
    EXPECT_EQ(1, dbg[5].Y().val());
    EXPECT_EQ(2, dbg[5].Z().val());
}
TEST(hist,hist_avr)
{
    auto chain=hist_avr(
	Points<value<>>{{{0,0.5},{0,1}},{{1,0.5},{1,1}},{{2,0.5},{2,1}}},
	Points<value<>>{{{0,0.5},{0,1}},{{1,0.5},{0,1}},{{2,0.5},{0,1}}}
    );
    for (size_t i = 0; i < chain.size(); i++) {
        EXPECT_EQ(i, chain[i].X().val());
        EXPECT_EQ(double(i)/2, chain[i].Y().val());
    }
    const hist<>
    a=Points<value<>>{{{0,0.5},{0,1}},{{1,0.5},{1,1}},{{2,0.5},{2,1}}},
    b=Points<value<>>{{{0,0.5},{0,1}},{{1,0.5},{0,1}},{{2,0.5},{0,1}}};
    chain=hist_avr(a,b);
    for (size_t i = 0; i < chain.size(); i++) {
        EXPECT_EQ(i, chain[i].X().val());
        EXPECT_EQ(double(i)/2, chain[i].Y().val());
    }
}
TEST(hist,hist_stdev)
{
    const auto chain=hist_stdev(
	Points<>{{0,-1},{1,0},{2,1}},
	Points<>{{0,1},{1,2},{2,3}}
    );
    for (size_t i = 0; i < chain.size(); i++) {
        EXPECT_EQ(i, chain[i].X());
        EXPECT_EQ(int(chain[i].Y().val()*10),i*10);
    }
}
TEST(Cache,test_stored_value1){
    Cache<size_t,string> C;
    for(size_t k=0;k<10;k++){
	for(size_t c=k;c<20;c++){
	    const auto&val=C(k,[&c](){return to_string(c);});
	    EXPECT_EQ(val,to_string(k));
	}
    }
}
class CacheTest{
private:
    size_t m_val1;
    string m_val2;
public:
    //We need a class with two parameters constructor
    CacheTest(size_t v1,const string&v2):m_val1(v1),m_val2(v2){}
    CacheTest(const CacheTest&other):m_val1(other.m_val1),m_val2(other.m_val2){}
    CacheTest():m_val1(0),m_val2(""){}
    ~CacheTest(){}
    bool operator==(const CacheTest&other)const{return (m_val1==other.m_val1)&&(m_val2==other.m_val2);}
};
TEST(Cache,test_stored_value2){
    Cache<size_t,CacheTest> C;
    for(size_t k=0;k<10;k++){
	for(size_t c=k;c<20;c++){
	    EXPECT_EQ(C(k,c,"test"),CacheTest(k,"test"));
	}
    }
}
