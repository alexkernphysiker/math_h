Header files with mathematical template functions and classes (C++17)
====================================================================

The library is distributed under LGPL v.3 license.

    Question of how LGPL v.3 deals with c++ templates is explained here:
    http://eigen.tuxfamily.org/index.php?title=Licensing_FAQ&oldid=1117
    That's a page of C++ library that also contains templates



How to use it
=============
If you have your own git repository with cmake project you can add this repo as a submodule:

	git submodule add https://github.com/alexkernphysiker/math_h.git
	git submodule update --init --recursive
	
Then add to your CMakeLists.txt

	add_definitions(-std=c++17)
	add_subdirectory(math_h)
	include_directories(${MATH_H_INC})
	
Then commit your changes

CMake Options
=============

	tests
if ON the tests are compiled

Header files
============

	include/math_h/*.h
	include/*.h

Contain template classes and can be included in your project

Examples
========
	Examples/*.cpp
Simple programs that show how to use this template classes

Tests
=====
	tests/*.cpp
The unit tests (require GoogleTest framework to be installed)
