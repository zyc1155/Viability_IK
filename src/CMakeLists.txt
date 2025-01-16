
option(BUILD_PYTHON_LIB "Build python library" ON)

if(BUILD_PYTHON_LIB)
        set(PYBIND11_FINDPYTHON ON)
        find_package(pybind11 REQUIRED)  
endif()

find_package (Eigen3 REQUIRED NO_MODULE)
find_package(eigen-cdd REQUIRED)
find_package(CasADi REQUIRED)


if(TARGET SCIP::SCIP)
  find_package(SCIP CONFIG PATHS ${SCIP_BINARY_DIR} REQUIRED)
else()
  find_package(SCIP REQUIRED)
endif()


add_library(viability_ik_shared SHARED viability_ik.cpp MILPSolver.cpp)
add_library(viability_ik_static STATIC viability_ik.cpp MILPSolver.cpp)

target_include_directories(viability_ik_shared PRIVATE 
        ${EIGEN3_INCLUDE_DIR}
        ${SCIP_INCLUDE_DIRS}
        ${CMAKE_SOURCE_DIR}/include
        )
target_include_directories(viability_ik_static PRIVATE 
        ${EIGEN3_INCLUDE_DIR}
        ${SCIP_INCLUDE_DIRS}
        ${CMAKE_SOURCE_DIR}/include
        )

target_link_libraries(viability_ik_shared PUBLIC
  Eigen3::Eigen
  eigen-cdd::eigen-cdd
  casadi
  ${SCIP_LIBRARIES}
)
target_link_libraries(viability_ik_static PUBLIC
  Eigen3::Eigen
  eigen-cdd::eigen-cdd
  casadi
  ${SCIP_LIBRARIES}
)

set_target_properties(viability_ik_shared PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR}/Debug
    RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR}/Release
    LIBRARY_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR}/Debug
    LIBRARY_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR}/Release)

set(VIABILITY_IK_LIB viability_ik_shared PARENT_SCOPE)
set(VIABILITY_IK_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include PARENT_SCOPE)

install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/
        DESTINATION include)
install(TARGETS viability_ik_shared
        LIBRARY DESTINATION lib)

if(BUILD_PYTHON_LIB)
pybind11_add_module(pyviability_ik MODULE viability_ik.cpp MILPSolver.cpp bindings.cpp)
target_include_directories(pyviability_ik PRIVATE ${EIGEN3_INCLUDE_DIR} ${SCIP_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(pyviability_ik PRIVATE Eigen3::Eigen eigen-cdd::eigen-cdd casadi ${SCIP_LIBRARIES})
install(TARGETS pyviability_ik
        LIBRARY DESTINATION lib)
endif()

