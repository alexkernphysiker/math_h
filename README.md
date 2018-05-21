Header files with mathematical template functions and classes (C++17)
=====================================================================

This library is distributed under LGPL v.3 license.
Don't be afraid, you can use headers with template classes/functions in your proprietary software if you don't modify them.


What features does this library provide
=======================================

- Using physical magnitudes values given with uncertainties. Arithmetical operations with them (uncertainties are calculated as well)

- Standard deviation, weighted average and covariance matrices calculations

- Numeric integrating of function using Sympson folmula. Convolution integral calculation. Integrating function given as table data using trapeze method

- Linear interpolation

- Generating random values with different distributions. Even ones given by table of density function.

- Vectors and matrices with arbitrary dimensions count. Vectors can be rotated and represented in polar coordinates.

- Lorentz vectors and relativistic kinematics facilities.

- Determinant calculation and solving systems of linear equations using Cramer's method

- Histograms, arithmetic operations with them.

- Plotting data with gnuplot


How to add into your project
============================

If you have your own git repository with cmake project you can add this repository as a submodule:

	git submodule add https://github.com/alexkernphysiker/math_h.git
	git submodule update --init --recursive
	
Then add to your project's CMakeLists.txt

	add_definitions(--std=c++17) #the most recommended compiler mode
	#set(GTEST ON) #uncomment for compiling unit-tests. Requires gtest
	add_subdirectory(math_h)
	include_directories(${MATH_H_INC})
	
Then commit your changes.

This library still can be compiled with c++11 or c++14 compiler but some features will be absent.

Header files
============

	include/math_h/*.h
	include/*.h


How to use in your program
==========================

Typical patterns of using the template classes and functions from this library can be learnt from:

Examples

	Examples/*.cpp

Unit tests (are instead of documentation)

	tests/*.cpp


Includes and classes provided
=============================

<math_h/error.h>

	Exception - template class for exceptions that are thrown by objects in this library

<math_h/functions.h>

	IFunction - template class providing a function-like interface with operator()

	PI,E - template functions returning pi and e numbers

	Gaussian, Lorentzian, BreitWigner, Novosibirsk - peak shaped functions often used in physical calculations

	FermiFunc - slope function often used in physical calculations

	Polynom - two template functions for calculating polinomials

<math_h/vectors.h>

	Vector<x> - template class that represents fixed but arbitrary sized vector. Supports adding subtracting, multiplying by a number.  Also supports scalar(*) multiplication. Axial-vector or pseudoscalar multiplication (^) are declared in math_h/vectortransformations.h

	vec(...) - create vector by coordinates

	zero() - 2d zero vector.

	x(),y() - 2d basis vectors

	Zero() - 3d zero vector

	X(),Y(),Z() - 3d basis vector

<math_h/matrices.h>

	Matrix<y,x> - template class that represents fixed but arbitrary sized matrix. Supports adding subtracting matrices with the same size, multiplying by a number. Supports matrix multiplication, determinant calculation and solving system of linear equation using Cramer method.

	rows(...) - create matrix from several vectors concidering them as rows

	columns(...) - create matrix from several vectors concidering them as columns

<math_h/vectortransformations.h>

	Direction - template class representing vectors direction (phi and theta angles). Supports converting to vector (Direction*length = vector)

<math_h/statistics.h>

	Sampling - template class that calculates covariances for a sample of vectors

	SamplingXY - the same but conciders a sample of vector pairs

<math_h/sigma.h>

	value<type> - template class representing a value with uncertainty. Supports arithmetic operations.

	std_error(...) - creates value N with uncertainty of sqrt(N)

	func_with_uncertainty(...) - applies function to value with uncertainty

	StandardDeviation - class that is used for calculation of average and standard deviation. use operator<< to input single measurements.

	WeightedAverage - weighted average calculator. Use operator<< to input single measurements.

<math_h/sigma2.h>

	abstract_value_with_uncertainty_numeric - abstract class representing value with uncertainty that has some distribution that can be different from normal. Also can be used for adding magnitudes that have dependent uncertainties.

	value_numeric_const - inherits abstract_value_with_uncertainty_numeric and implements constant value

	value_numeric_distr - inherits abstract_value_with_uncertainty_numeric and implements uncertainty with some distribution (template parameter)

	value_numeric(...) - converts value into value_numeric_distr with normal distribution

	SQRT, SQR, EXP, LOG, SIN, COS, TAN, ATAN, ATAN2, POW, operator+, operator-, operator*, operator/ - arithmetic actions over such values with uncertainties

	FUNC(...) - applying a function to values with uncertainties.

<math_h/lorentzvector.h>

	LorentzVector - template class that represents energy-momenta lorentz 4-vectors. Arithmetical actions of adding and subtracting are supported. Lorentz transformation is implemented as well.

	lorentz_Rest - creates lorentz vector of particle in rest by it's mass

	lorentz_byPM - creates lorents vector my momentum and mass

	lorentz_byEM - by full energy and mass

	lorentz_byEkM - by kinetic energy and mass

	binaryDecay - simulates decay of compound system into two parts

<math_h/tabledata.h>

	SortedChain - template class providing vector-like interface to a sorted set of values of any type that can be compared. Operator << provides element adding and it's inserted in proper place immediately.

	ChainWithStep, ChainWithCount - functions that generate chains of values uniformly filling the given range

	BinsByStep, BinsByCount - generate chains of values with uncertainties filling the given range

	point - template class representing data point with x and y-values

	point3d - template class representing data point with x,y and z-values

	SortedPoints - template class inheriting SortedChain that is used to represent data points and histograms. Supports arithmetic operations

	BiSortedPoints - template class that is used to represent 3-d data points (2-d historgams)

	hist<x,y> - shortcut for SortedPoints<value<x>,value<y>>

	hist_avr(...) - calculates wieghted average for several histograms

	Distribution1D, Distribution2D - template classes that are used for remembering distribution of a magnitude. Fill(...) method is used for adding a single value to the distribution. Then is can be operated as a histogram.

<math_h/sigma3.h>

	Uncertainties - template class representing value with several uncertainties (for example statistical and systematic uncertainties) that are independent. Supports arithmetic actions.

	uncertainties(...) - creates value with several uncertainties

	extend_value<i,n>(...) - converts value<> itno Uncertainties<n> using value's uncertainty as i-th one. All others are set to zero.

	ext_hist<n> - is like hist<> but y-values are Uncertainties<n>

	extend_hist<i,n>(...) - convert hist<> to ext_hist<n>

	wrap_hist(...) - convert ext_hist<n> to hist<> (total uncertainties are calculated)

<math_h/interpolate.h>

	LinearInterpolation - inherits SortedPoints and provides function-like interface to linear interpolation algorithm

	BiLinearInterpolation - inherits BiSortedPoints and provides bilinear interpolation

<math_h/integrate.h>

	Sympson - function that performs numeric integrating with Sympson formula. Argument function is given like an std::function or IFunction

	Int_Trapez_Table - integrates using trapez formula. Argument function is given as a table data (for example SortedPoints)

	Convolution - template class that calculates convolution integral and provides function-like interface

<math_h/randomfunc.h>

	RandomUniform - random number generator implementing uniform distribution

	RandomGauss - random number generator implementing Gaussian distribution

	RandomValueTableDistr - generator implementing distribution given by a table of density values (SortedPoints)

	Poisson() - creates generator implementing Poisson distribution

