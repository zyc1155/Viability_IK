find_package(Eigen3 REQUIRED NO_MODULE)
find_package(CasADi REQUIRED)
find_package(osqp REQUIRED)
find_package(SCIP REQUIRED)

set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

set(HEADERS_VIK MILPSolver.h scip_exception.h viability_ik.h)
set(HEADERS_EIGENCDD eigen-cddlib/Polyhedron.h eigen-cddlib/typedefs.h)

add_library(${PROJECT_NAME} SHARED viability_ik.cpp MILPSolver.cpp eigen-cddlib/Polyhedron.cpp)

target_include_directories(${PROJECT_NAME} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<INSTALL_INTERFACE:include>)
target_include_directories(${PROJECT_NAME} PRIVATE ${EIGEN3_INCLUDE_DIR})

target_link_libraries(${PROJECT_NAME} PUBLIC
        Eigen3::Eigen
        casadi
        libscip
        cdd
)

set_target_properties(${PROJECT_NAME} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR}/Debug
        RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR}/Release
        LIBRARY_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR}/Debug
        LIBRARY_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR}/Release)

install(FILES ${HEADERS_VIK}
        DESTINATION include/${PROJECT_NAME})
install(FILES ${HEADERS_EIGENCDD}
        DESTINATION include/${PROJECT_NAME}/eigen-cddlib)

install(TARGETS ${PROJECT_NAME}
        EXPORT ${PROJECT_NAME}Targets
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        RUNTIME DESTINATION bin)

configure_package_config_file(
        ${CMAKE_SOURCE_DIR}/cmake/Config.cmake.in
        ${CMAKE_BINARY_DIR}/${PROJECT_NAME}Config.cmake
        INSTALL_DESTINATION lib/cmake/${PROJECT_NAME}
)

install(FILES ${CMAKE_BINARY_DIR}/${PROJECT_NAME}Config.cmake
        DESTINATION lib/cmake/${PROJECT_NAME}
)

install(EXPORT ${PROJECT_NAME}Targets
        FILE ${PROJECT_NAME}Targets.cmake
        DESTINATION lib/cmake/${PROJECT_NAME})


set(BUILD_PYTHON_PATH "lib" CACHE PATH "Installation path for python wrapper")
if(BUILD_PYTHON)
        if(BUILD_IN_CONDA)
                set(Python3_FIND_VIRTUALENV ONLY)
                find_package(Python3 MODULE COMPONENTS Interpreter Development)
                if(WIN32)
                        file(TO_CMAKE_PATH "${Python3_SITEARCH}" Python3_SITEARCH_UNIX)
                        set(CMAKE_INSTALL_PREFIX "${Python3_SITEARCH_UNIX}/Library")
                else()
                        set(BUILD_PYTHON_PATH "${Python3_SITEARCH}")
                endif()
        else()
                set(PYBIND11_FINDPYTHON ON)
        endif()
        find_package(pybind11 REQUIRED)

        pybind11_add_module(py_${PROJECT_NAME} MODULE bindings.cpp)
        target_include_directories(py_${PROJECT_NAME} PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<INSTALL_INTERFACE:include>)
        target_include_directories(py_${PROJECT_NAME} PRIVATE ${EIGEN3_INCLUDE_DIR})
        target_link_libraries(py_${PROJECT_NAME} PRIVATE ${PROJECT_NAME})
        install(TARGETS py_${PROJECT_NAME}
                LIBRARY DESTINATION ${BUILD_PYTHON_PATH})
endif()

