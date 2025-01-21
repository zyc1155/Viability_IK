cmake_minimum_required(VERSION 3.10.0)
project(viability_ik VERSION 1.1.0 LANGUAGES C CXX)
set(CMAKE_CXX_STANDARD 11)

option(BUILD_IN_CONDA "Build program in conda enviroment" ON)
option(BUILD_PYTHON "Build python wrapper" ON)
option(BUILD_CXX_EXAMPLE "Build cxx example" ON)

if(BUILD_IN_CONDA)
    # Get CONDA_PREFIX from the environment
    if(DEFINED ENV{CONDA_PREFIX})
        set(CONDA_PREFIX "$ENV{CONDA_PREFIX}")
    else()
        message(FATAL_ERROR "CONDA_PREFIX is not set. Activate Conda environment first.")
    endif()
    if(WIN32)
        # Windows: Install to ${CONDA_PREFIX_UNIX}/Library
        file(TO_CMAKE_PATH "${CONDA_PREFIX}" CONDA_PREFIX_UNIX)
        set(CMAKE_INSTALL_PREFIX "${CONDA_PREFIX_UNIX}/Library")
    else()
        # macOS and Linux: Install to ${CONDA_PREFIX}
        set(CMAKE_INSTALL_PREFIX "${CONDA_PREFIX}")
    endif()
endif()

include(CMakePackageConfigHelpers)

add_subdirectory(cddlib)

add_subdirectory(src)

if(BUILD_CXX_EXAMPLE)
    add_subdirectory(examples)
endif()

