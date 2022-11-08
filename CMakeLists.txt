cmake_minimum_required(VERSION 3.10)
project(WaveSimC VERSION 0.8)



set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
set (source_dir "${PROJECT_SOURCE_DIR}/src")


add_executable(Main "${source_dir}/main.cpp")
target_include_directories(Main PRIVATE "${source_dir}")

add_library(cxxopts STATIC "${source_dir}/ExternalLibraries/cxxopts.hpp")
set_target_properties(cxxopts PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(Main PRIVATE cxxopts)

#configure_file(Config.h.in Config.h)