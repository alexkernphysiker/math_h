cmake_minimum_required(VERSION 2.6)
project(math_h)

add_definitions(-std=c++11 -Wall -Wunused-parameter)
set(CMAKE_CXX_FLAGS "-pthread -fPIC")
get_directory_property(has_parent PARENT_DIRECTORY)
if(has_parent)
	set(MATH_H_INC 
		${CMAKE_CURRENT_SOURCE_DIR}/include
		PARENT_SCOPE)
endif()

file(GLOB SOURCES gnuplot_wrap/*.cpp)
file(GLOB TESTS tests/*.cpp)
file(GLOB EXAMPLES Example/*.cpp)

include_directories(include)
add_library(gpwrap STATIC ${SOURCES})
foreach(example ${EXAMPLES})
	get_filename_component(ex_app ${example} NAME_WE)
	add_executable(${ex_app} ${example})
	target_link_libraries(${ex_app} gpwrap)
endforeach()

option(tests "Compile tests." OFF) 
if(tests)
	find_package(GTest REQUIRED)
	include_directories(${GTEST_INCLUDE_DIRS})
	add_executable(math_h_test.exe ${TESTS})
	target_link_libraries(math_h_test.exe ${GTEST_BOTH_LIBRARIES})
endif()
