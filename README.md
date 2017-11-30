Header files with mathematical template functions and classes (C++17)
=====================================================================

The library is distributed under LGPL v.3 license.
Don't be afraid, you still can use these headers in your proprietary software if you don't modify them.

Question of how LGPL v.3 deals with c++ templates is explained here:

    http://eigen.tuxfamily.org/index.php?title=Licensing_FAQ&oldid=1117

That's a page of C++ library that also consists of templates



How to use it
=============

If you have your own git repository with cmake project you can add this repository as a submodule:

	git submodule add https://github.com/alexkernphysiker/math_h.git
	git submodule update --init --recursive
	
Then add to your project's CMakeLists.txt

	add_definitions(-std=c++17)
	add_subdirectory(math_h)
	include_directories(${MATH_H_INC})
	
Then commit your changes.

This library still can be compiled with c++11 or c++14 compiler but some features will be less optimized and work slower.


CMake Options
=============

Option that turns on unit tests compilation

	tests

Header files (that are used)
============================

	include/math_h/*.h
	include/*.h


Examples
========

The examples that show the main facilities of the library and good using patterns are here:

	Examples/*.cpp

Please look through them before you start using this code

Unit tests (require GoogleTest framework)
=========================================

The unit tests for this library are here

	tests/*.cpp

For compilation use "tests" option for cmake project.
This requires GoogleTest framework installed into your system.
You can look into these cpp files and use them as documentation for classes declared in the library
