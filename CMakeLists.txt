cmake_minimum_required(VERSION 3.10)
project(WaveSimC VERSION 0.8)



set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fconcepts")
set (source_dir "${PROJECT_SOURCE_DIR}/src")
find_package(Boost 1.40 COMPONENTS program_options REQUIRED)

add_executable(Main "${source_dir}/main.cpp")
target_include_directories(Main PRIVATE "${source_dir}")

add_library(cxxopts STATIC "${source_dir}/ExternalLibraries/cxxopts.hpp")
set_target_properties(cxxopts PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(Main cxxopts Boost::program_options)

add_executable(VariadicTest "${source_dir}/tests/variadic.cpp")
target_link_libraries(VariadicTest Boost::program_options)
target_include_directories(VariadicTest PRIVATE "${source_dir}")

#configure_file(Config.h.in Config.h)