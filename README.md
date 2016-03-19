Header files with mathematical template functions and classes (C++11)
====================================================================
This repository also contains a small static library which is an interface between these classes and gnuplot




Compiling
=========
If you have your own git repository with cmake project you can add this repo as a submodule:

	git submodule add https://github.com/alexkernphysiker/math_h.git
	git submodule update --init

Then add to your CMakeLists.txt

	add_subdirectory(math_h)
	include_directories(math_h/include)

Then commit your changes :)

CMake Options
=============
	debug - if ON the project is compiled in debug mode
	test - if ON the tests are compiled
	example - if ON the examples are compiled
	
Header files
============
	math_h/bit_opr.h - useful templates for constructing bit masks
	math_h/error.h - template exception class that requires "what"-message as constructor argument and the template parameter is the name of class that throws it.
	math_h/functions.h - interface IFunction and some useful template functions.
	math_h/sigma.h - class for representing value with uncertainty and templates for calculations of dispersion and wieghted average.
	math_h/structures.h - classes for representing data structures usually useful in calculations
	math_h/interpolate.h - interpolation algorithm(s).
	math_h/integrate.h - calculating numeric integrals
	math_h/randomfunc.h - templates for generating random values distributed by given function.

	gnuplot_wrap.h - header for gnuplot interface library (libgpwrap)


Examples
========
	Examples/*.cpp - simple programs that show how to use this template classes
	

Tests
=====
The directory tests contains cmake project with unit tests.
The tests require GoogleTest framework to be installed
