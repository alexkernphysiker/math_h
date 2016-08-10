// this file is distributed under 
// MIT license
#include <gtest/gtest.h>
#include <vector>
#include <math.h>
#include <math_h/matrices.h>
using namespace std;
using namespace MathTemplates;
typedef Exception<vector<vector<double>>,0> invalid_matr;
typedef Exception<vector<vector<double>>,1> size_mismatch;
TEST(Matrix,MatrixValid){
	EXPECT_FALSE(MatrixValid<double>({}));
	EXPECT_FALSE(MatrixValid<double>({{}}));
	EXPECT_FALSE(MatrixValid<double>({{},{}}));
	EXPECT_TRUE(MatrixValid<double>({{1}}));
	EXPECT_TRUE(MatrixValid<double>({{1},{1}}));
	EXPECT_TRUE(MatrixValid<double>({{1},{1},{1}}));
	EXPECT_FALSE(MatrixValid<double>({{1},{1},{}}));
	EXPECT_FALSE(MatrixValid<double>({{1},{},{1}}));
	EXPECT_FALSE(MatrixValid<double>({{},{1},{1}}));
	EXPECT_TRUE(MatrixValid<double>({{1,1}}));
	EXPECT_TRUE(MatrixValid<double>({{1,1},{1,1}}));
	EXPECT_TRUE(MatrixValid<double>({{1,1},{1,1},{1,1}}));
	EXPECT_FALSE(MatrixValid<double>({{1,1},{1,1},{1}}));
	EXPECT_FALSE(MatrixValid<double>({{1,1},{1},{1,1}}));
	EXPECT_FALSE(MatrixValid<double>({{1},{1,1},{1,1}}));
	EXPECT_FALSE(MatrixValid<double>({{1,1},{1,1},{1,1,1}}));
	EXPECT_FALSE(MatrixValid<double>({{1,1},{1,1,1},{1,1}}));
	EXPECT_FALSE(MatrixValid<double>({{1,1,1},{1,1},{1,1}}));
	EXPECT_FALSE(MatrixValid<double>({{1,1,1},{1,1},{1}}));
	EXPECT_FALSE(MatrixValid<double>({{1},{1,1},{1,1,1}}));
	EXPECT_TRUE(MatrixValid<double>({{1,1,1}}));
	EXPECT_TRUE(MatrixValid<double>({{1,1,1},{1,1,1}}));
	EXPECT_TRUE(MatrixValid<double>({{1,1,1},{1,1,1},{1,1,1}}));
	EXPECT_FALSE(MatrixValid<double>({{1,1,1},{1,1,1},{1,1}}));
	EXPECT_FALSE(MatrixValid<double>({{1,1,1},{1,1},{1,1,1}}));
	EXPECT_FALSE(MatrixValid<double>({{1,1},{1,1,1},{1,1,1}}));
	EXPECT_FALSE(MatrixValid<double>({{1,1,1},{1,1,1},{1}}));
	EXPECT_FALSE(MatrixValid<double>({{1,1,1},{1},{1,1,1}}));
	EXPECT_FALSE(MatrixValid<double>({{1},{1,1,1},{1,1,1}}));	
}
TEST(Matrix,MatrixSquare){
	EXPECT_FALSE(MatrixSquare<double>({}));
	EXPECT_FALSE(MatrixSquare<double>({{}}));
	EXPECT_FALSE(MatrixSquare<double>({{},{}}));
	EXPECT_TRUE(MatrixSquare<double>({{1}}));
	EXPECT_FALSE(MatrixSquare<double>({{1},{1}}));
	EXPECT_FALSE(MatrixSquare<double>({{1},{1},{1}}));
	EXPECT_FALSE(MatrixSquare<double>({{1},{1},{}}));
	EXPECT_FALSE(MatrixSquare<double>({{1},{},{1}}));
	EXPECT_FALSE(MatrixSquare<double>({{},{1},{1}}));
	EXPECT_FALSE(MatrixSquare<double>({{1,1}}));
	EXPECT_TRUE(MatrixSquare<double>({{1,1},{1,1}}));
	EXPECT_FALSE(MatrixSquare<double>({{1,1},{1,1},{1,1}}));
	EXPECT_FALSE(MatrixSquare<double>({{1,1},{1,1},{1}}));
	EXPECT_FALSE(MatrixSquare<double>({{1,1},{1},{1,1}}));
	EXPECT_FALSE(MatrixSquare<double>({{1},{1,1},{1,1}}));
	EXPECT_FALSE(MatrixSquare<double>({{1,1},{1,1},{1,1,1}}));
	EXPECT_FALSE(MatrixSquare<double>({{1,1},{1,1,1},{1,1}}));
	EXPECT_FALSE(MatrixSquare<double>({{1,1,1},{1,1},{1,1}}));
	EXPECT_FALSE(MatrixSquare<double>({{1,1,1},{1,1},{1}}));
	EXPECT_FALSE(MatrixSquare<double>({{1},{1,1},{1,1,1}}));
	EXPECT_FALSE(MatrixSquare<double>({{1,1,1}}));
	EXPECT_FALSE(MatrixSquare<double>({{1,1,1},{1,1,1}}));
	EXPECT_TRUE(MatrixSquare<double>({{1,1,1},{1,1,1},{1,1,1}}));
	EXPECT_FALSE(MatrixSquare<double>({{1,1,1},{1,1,1},{1,1}}));
	EXPECT_FALSE(MatrixSquare<double>({{1,1,1},{1,1},{1,1,1}}));
	EXPECT_FALSE(MatrixSquare<double>({{1,1},{1,1,1},{1,1,1}}));
	EXPECT_FALSE(MatrixSquare<double>({{1,1,1},{1,1,1},{1}}));
	EXPECT_FALSE(MatrixSquare<double>({{1,1,1},{1},{1,1,1}}));
	EXPECT_FALSE(MatrixSquare<double>({{1},{1,1,1},{1,1,1}}));	
}
TEST(Matrix,Compare){
	EXPECT_THROW(MatrixSizeEqual<double>({},{}),invalid_matr);
	EXPECT_THROW(MatricesEqual<double>({},{}),invalid_matr);

	EXPECT_TRUE(MatrixSizeEqual<double>({{1}},{{1}}));
	EXPECT_TRUE(MatricesEqual<double>({{1}},{{1}}));
	EXPECT_TRUE(MatrixSizeEqual<double>({{1}},{{0}}));
	EXPECT_FALSE(MatricesEqual<double>({{1}},{{0}}));

	EXPECT_TRUE(MatrixSizeEqual<double>({{1,1}},{{1,1}}));
	EXPECT_TRUE(MatricesEqual<double>({{1,1}},{{1,1}}));
	EXPECT_TRUE(MatrixSizeEqual<double>({{1,0}},{{1,1}}));
	EXPECT_FALSE(MatricesEqual<double>({{1,0}},{{1,1}}));
	
	EXPECT_FALSE(MatrixSizeEqual<double>({{1,1},{1,1}},{{1,1}}));
	EXPECT_FALSE(MatricesEqual<double>({{1,1},{1,1}},{{1,1}}));
	EXPECT_FALSE(MatrixSizeEqual<double>({{1,0},{1,1}},{{1,1},{1,1},{1,1}}));
	EXPECT_FALSE(MatricesEqual<double>({{1,0},{1,1}},{{1,1},{1,1},{1,1}}));
}
TEST(Matrix,MulMatrices){
	EXPECT_THROW(MulMatrices<double>({{}},{{1}}),invalid_matr);
	EXPECT_THROW(MulMatrices<double>({{1}},{{}}),invalid_matr);
	EXPECT_THROW(MulMatrices<double>({{1,1}},{{1,1}}),size_mismatch);
	EXPECT_TRUE(MatricesEqual(MulMatrices<double>({{1}},{{1}}),{{1}}));
	EXPECT_TRUE(MatricesEqual(MulMatrices<double>({{2}},{{3}}),{{6}}));
	EXPECT_TRUE(MatricesEqual(MulMatrices<double>({{1,0},{0,1}},{{1,0},{0,1}}),{{1,0},{0,1}}));
	EXPECT_TRUE(MatricesEqual(MulMatrices<double>({{2,0},{0,4}},{{3,0},{0,5}}),{{6,0},{0,20}}));
	EXPECT_TRUE(MatricesEqual(MulMatrices<double>({{1,0,0},{0,1,0},{0,0,1}},{{1,0,0},{0,1,0},{0,0,1}}),{{1,0,0},{0,1,0},{0,0,1}}));
	EXPECT_TRUE(MatricesEqual(MulMatrices<double>({{1,0,0},{0,2,0},{0,0,3}},{{4,0,0},{0,5,0},{0,0,6}}),{{4,0,0},{0,10,0},{0,0,18}}));

	EXPECT_TRUE(MatricesEqual(MulMatrices<double>({{1,2},{3,4}},{{5,6},{7,8}}),{{19,22},{43,50}}));
	EXPECT_TRUE(MatricesEqual(MulMatrices<double>({{1,2}},{{3},{4}}),{{11}}));
}