cmake_minimum_required(VERSION 3.5)
project(Viability_IK)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(BUILD_CXX_EXAMPLES "Build example_cxx" ON)

add_subdirectory(src)

if(BUILD_CXX_EXAMPLES)
    add_subdirectory(examples) 
endif()

