Header files with mathematical template functions and classes (C++11)
====================================================================
This repository also contains a small static library which is an interface between these classes and gnuplot




Compiling
=========
If you have your own git repository with cmake project you can add this repo as a submodule:

	git submodule add https://github.com/alexkernphysiker/math_h.git
	git submodule update --init

Then add to your CMakeLists.txt

	set(gnuplot_wrap ON CACHE BOOL "") #if you want to use gnuplot interface library
	add_subdirectory(math_h)
	include_directories(math_h/include)


CMake Options
=============
	debug - if ON the project is compiled in debug mode
	test - if ON the tests are compiled
	gnuplot_wrap - if ON the library which is an interface for gnuplot is compiled

	
Header files
============
	math_h/bit_opr.h - useful templates for constructing bit masks
	math_h/error.h - template exception class that requires "what"-message as constructor argument and the template parameter is the name of class that throws it.
	math_h/interpolate.h - templates for searching in sorted arrays and for linear interpolated functions based on table data.
	math_h/integrate.h - templates for calculating numeric integrals
	math_h/sigma.h - templates for calculations of dispersion and wieghted average.
	math_h/hist.h - template classes for creating histograms
	math_h/randomfunc.h - templates for generating random values distributed by given function.
	math_h/functions.h - templates for functions that I frequently need in my calculations.

	gnuplot_wrap.h - header for gnuplot interface library (libgpwrap)



Tests
=====
The directory FitGen_tests contains cmake project with unit tests.
The tests require GoogleTest framework to be installed
