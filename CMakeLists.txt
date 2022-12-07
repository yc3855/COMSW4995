cmake_minimum_required(VERSION 3.10)
project(WaveSimC VERSION 0.8)



set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fconcepts")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fmodules-ts")
set (source_dir "${PROJECT_SOURCE_DIR}/src")

#
#include(FetchContent)
#FetchContent_Declare(matplotplusplus
#        GIT_REPOSITORY https://github.com/alandefreitas/matplotplusplus
#        GIT_TAG origin/master) # or whatever tag you want
#FetchContent_GetProperties(matplotplusplus)
#        if(NOT matplotplusplus_POPULATED)
#            FetchContent_Populate(matplotplusplus)
#            add_subdirectory(${matplotplusplus_SOURCE_DIR} ${matplotplusplus_BINARY_DIR} EXCLUDE_FROM_ALL)
#        endif()

find_package(Matplot++ REQUIRED)
find_package(Boost 1.40 COMPONENTS program_options REQUIRED)
#add_subdirectory("${source_dir}/ExternalLibraries/matplotplusplus")

add_executable(Main "${source_dir}/main.cpp")
target_include_directories(Main PRIVATE "${source_dir}")

add_library(cxxopts STATIC "${source_dir}/ExternalLibraries/cxxopts.hpp")
set_target_properties(cxxopts PROPERTIES LINKER_LANGUAGE CXX)
#target_link_libraries(Main cxxopts matplot Boost::program_options)
target_link_libraries(Main cxxopts Matplot++::matplot Boost::program_options)
#target_link_libraries(Main cxxopts Boost::program_options)

add_executable(VariadicTest "${source_dir}/tests/variadic.cpp")
target_link_libraries(VariadicTest Boost::program_options)
target_include_directories(VariadicTest PRIVATE "${source_dir}")

add_executable(CoreTests "${source_dir}/tests/CoreTests.cpp")
target_link_libraries(CoreTests Matplot++::matplot Boost::program_options)
target_include_directories(CoreTests PRIVATE "${source_dir}")

add_executable(MatPlotTest "${source_dir}/tests/MatPlotTest.cpp")
target_link_libraries(MatPlotTest Matplot++::matplot)
target_include_directories(MatPlotTest PRIVATE "${source_dir}")

#configure_file(Config.h.in Config.h)