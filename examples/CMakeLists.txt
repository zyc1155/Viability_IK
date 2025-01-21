
find_package(Eigen3 REQUIRED NO_MODULE)

add_executable(${PROJECT_NAME}_example_cxx example_cxx.cpp)
target_include_directories(${PROJECT_NAME}_example_cxx PRIVATE ${EIGEN3_INCLUDE_DIR})
target_include_directories(${PROJECT_NAME}_example_cxx PRIVATE $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src> $<INSTALL_INTERFACE:include>)
target_link_libraries(${PROJECT_NAME}_example_cxx PRIVATE ${PROJECT_NAME} Eigen3::Eigen)

set_target_properties(${PROJECT_NAME}_example_cxx PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR}/Debug
    RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR}/Release
    LIBRARY_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR}/Debug
    LIBRARY_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR}/Release)

install(TARGETS ${PROJECT_NAME}_example_cxx
    RUNTIME DESTINATION bin)
