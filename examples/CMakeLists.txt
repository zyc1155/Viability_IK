
find_package (Eigen3 REQUIRED NO_MODULE)

add_executable(example_cxx example_cxx.cpp)
target_include_directories(example_cxx PRIVATE 
        ${EIGEN3_INCLUDE_DIR}
        ${VIABILITY_IK_INCLUDE_DIR}
        )
target_link_libraries(example_cxx PRIVATE ${VIABILITY_IK_LIB} Eigen3::Eigen)

set_target_properties(example_cxx PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR}/Debug
    RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR}/Release
    LIBRARY_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR}/Debug
    LIBRARY_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR}/Release)

install(TARGETS example_cxx
    RUNTIME DESTINATION ${CMAKE_SOURCE_DIR}/examples)