#include <iostream>
#include <string>
#include "ExternalLibraries/cxxopts.hpp"
#include "CustomLibraries/np.hpp"
#include <matplot/matplot.h>

// Command line arguments
cxxopts::Options options("WaveSimC", "A wave propagation simulator written in C++ for seismic data processing.");
int main(int argc, char *argv[])
{
    // Parse command line arguments
    options.add_options()("d,debug", "Enable debugging")("i,input_file", "Input file path", cxxopts::value<std::string>())("o,output_file", "Output file path", cxxopts::value<std::string>())("v,verbose", "Verbose output", cxxopts::value<bool>()->default_value("false"));
    auto result = options.parse(argc, argv);

    std::cout << "Hello World"
              << "\n";
}