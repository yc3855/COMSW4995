cmake_minimum_required(VERSION 3.10)
project(WaveSimC VERSION 1.0)



set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fconcepts")
set (source_dir "${PROJECT_SOURCE_DIR}/src")


find_package(Matplot++ REQUIRED)
find_package(Boost 1.40 COMPONENTS program_options REQUIRED)

add_executable(WaveSimCExec "${source_dir}/main.cpp")
target_include_directories(WaveSimCExec PRIVATE "${source_dir}")

add_library(cxxopts STATIC "${source_dir}/ExternalLibraries/cxxopts.hpp")
set_target_properties(cxxopts PROPERTIES LINKER_LANGUAGE CXX)

target_link_libraries(WaveSimCExec cxxopts Matplot++::matplot Boost::program_options)


add_executable(VariadicTest "${source_dir}/tests/variadic.cpp")
target_link_libraries(VariadicTest Boost::program_options)
target_include_directories(VariadicTest PRIVATE "${source_dir}")

add_executable(CoreTests "${source_dir}/tests/CoreTests.cpp")
target_link_libraries(CoreTests Matplot++::matplot Boost::program_options)
target_include_directories(CoreTests PRIVATE "${source_dir}")

add_executable(MatPlotTest "${source_dir}/tests/MatPlotTest.cpp")
target_link_libraries(MatPlotTest Matplot++::matplot)
target_include_directories(MatPlotTest PRIVATE "${source_dir}")

add_executable(ExampleWithAnimation "${source_dir}/examples/wave_solver_with_animation.cpp")
target_link_libraries(ExampleWithAnimation Matplot++::matplot Boost::program_options)
target_include_directories(ExampleWithAnimation PRIVATE "${source_dir}")

#configure_file(Config.h.in Config.h)