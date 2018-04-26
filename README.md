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

