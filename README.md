Header files with mathematical template functions and classes (C++11)
====================================================================

this repository also contains a small static library which is an interface between these classes and gnuplot




Compiling
=========

If you are going to use only headers with mathematical templates then you don't need to compile anything else outside your project.

But this repository also contains one library (gnuplot interface) and unit tests for headers with templates.

For compiling this stuff you need to run:

git clone https://github.com/alexkernphysiker/math_h.git

cd math_h

cmake .

make

For using this project as a submodule in your git repository please run:

git submodule add https://github.com/alexkernphysiker/math_h.git

git submodule update --init

If you have your own cmake project, you can use this library by adding add_subdirectory instruction.
But you will have to add path math_h/include to your include directories.



CMake Options
=============

debug - if ON the project is compiled in debug mode

test - if ON the tests are compiled

gnuplot_wrap - if ON the library which is an interface for gnuplot is compiled


Header files
============

math_h/bit_opr.h - useful templates for constructing bit masks

math_h/error.h - template exception class that requires "what"-message as constructor argument 
and the template parameter is the name of class that throws it.

math_h/functions.h - templates for functions that I frequently need in my calculations.

math_h/interpolate.h - templates for searching in sorted arrays and for linear interpolated functions based on table data.

math_h/randomfunc.h - templates for generating random values distributed by given function.

math_h/sigma.h - templates for calculations of dispersion and wieghted average.

math_h/sympson.h - templates for calculating numeric integrals by Sympson's formula

gnuplot_wrap.h - header for gnuplot interface library (libgpwrap)



Tests
=====

The directory FitGen_tests contains cmake project with unit tests.
The tests require GoogleTest framework to be installed
