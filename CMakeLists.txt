cmake_minimum_required(VERSION 3.0)
project(math_h)
if(NOT DEFINED MATH_H_INC)
    option(INSTALL OFF)
    if(INSTALL)
	add_compile_options(-std=c++17)
    endif()
    get_directory_property(has_parent PARENT_DIRECTORY)
    if(has_parent)
	set(MATH_H_INC 
	    ${CMAKE_CURRENT_SOURCE_DIR}/include
	    PARENT_SCOPE)
    endif()

    include_directories(include)
    install(DIRECTORY include DESTINATION .)
    set(inst_targets "")

    file(GLOB EXAMPLES Example/*.cpp)
    foreach(example ${EXAMPLES})
	get_filename_component(ex_app ${example} NAME_WE)
	add_executable(${ex_app} ${example})
    endforeach()
    if(DEFINED GTEST)
	file(GLOB TESTS tests/*.cpp)
	add_executable(math_h_test.exe ${TESTS})
	target_link_libraries(math_h_test.exe gtest gtest_main)
    endif()
endif()
