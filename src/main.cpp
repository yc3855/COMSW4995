#include <iostream>
#include <string>
#include "ExternalLibraries/cxxopts.hpp"

#include <boost/multi_array.hpp>
#include <boost/array.hpp>
#include "CustomLibraries/np.hpp"
#include "CustomLibraries/np_to_matplot.hpp"
#include "CustomLibraries/wavePlotter.hpp"
#include <matplot/matplot.h>
#include <cassert>
#include <sstream>

#include "CoreAlgorithm/helper_func.hpp"
#include "CoreAlgorithm/coeff.hpp"
#include "CoreAlgorithm/source.hpp"
#include "CoreAlgorithm/computational.hpp"
#include "CoreAlgorithm/solver.hpp"

// Command line arguments
cxxopts::Options options("WaveSimPP", "A wave propagation simulator written in C++.");
int main(int argc, char *argv[])
{
    // Parse command line arguments
    options.add_options()("d,debug", "Enable debugging", cxxopts::value<bool>()->default_value("false"));
    options.add_options()("animate", "Render an animation at the end", cxxopts::value<bool>()->default_value("true"));
    options.add_options()("render", "Render each of the frames at the end", cxxopts::value<bool>()->default_value("false"));
    options.add_options()("export", "Export the data to a series of csv files", cxxopts::value<bool>()->default_value("false"));
    options.add_options()("o,output_dir", "Output directory path", cxxopts::value<std::string>()->default_value("output"));
    options.add_options()("output_filename", "Output filename", cxxopts::value<std::string>()->default_value("output.mp4"));
    options.add_options()("framerate", "Framerate of output video", cxxopts::value<int>()->default_value("30"));
    options.add_options()("v,verbose", "Verbose output", cxxopts::value<bool>()->default_value("false"));
    options.add_options()("source_i", "Source i position", cxxopts::value<int>()->default_value("50"));
    options.add_options()("source_j", "Source j position", cxxopts::value<int>()->default_value("50"));
    options.add_options()("nt", "Number of time steps", cxxopts::value<int>()->default_value("1000"));
    options.add_options()("nx", "Number of x steps", cxxopts::value<int>()->default_value("100"));
    options.add_options()("nz", "Number of z steps", cxxopts::value<int>()->default_value("100"));
    options.add_options()("dt", "Time step size", cxxopts::value<double>()->default_value("0.001"));
    options.add_options()("dx", "x step size", cxxopts::value<double>()->default_value("0.01"));
    options.add_options()("dz", "z step size", cxxopts::value<double>()->default_value("0.01"));
    options.add_options()("f,frequency", "Frequency of source", cxxopts::value<double>()->default_value("10.0"));
    options.add_options()("a,amplitude", "Amplitude of source", cxxopts::value<double>()->default_value("1.0"));
    options.add_options()("s,shift", "Shift of source", cxxopts::value<double>()->default_value("0.1"));
    options.add_options()("r,radius", "Radius of velocity profile", cxxopts::value<double>()->default_value("150.0"));
    options.add_options()("l,num_levels", "Number of levels in the filled contour plot", cxxopts::value<int>()->default_value("100"));
    options.add_options()("h,help", "Print help");
    cxxopts::ParseResult result;
    try
    {
        result = options.parse(argc, argv);
    }
    catch (const cxxopts::exceptions::exception &e)
    {
        std::cerr << "WaveSimPP: " << e.what() << '\n';
        std::cerr << "usage: WaveSimPP [options] ...\n";
        return EXIT_FAILURE;
    }
    if (result.count("help"))
    {
        std::cout << options.help({"", "Group"}) << std::endl;
        return true;
    }

    // Define the constants for the simulation

    // Number of x and z grid points

    int nx = result["nx"].as<int>();

    int nz = result["nz"].as<int>();
    // Number of time steps
    int nt = result["nt"].as<int>();

    // Differentiation values
    double dx = result["dx"].as<double>();
    double dz = result["dz"].as<double>();
    double dt = result["dt"].as<double>();

    // Define the domain
    double xmin = 0.0;
    double xmax = nx * dx;
    double zmin = 0.0;
    double zmax = nz * dz;
    double tmin = 0.0;
    double tmax = nt * dt;

    // Define the source parameters
    double f_M = result["frequency"].as<double>();
    double amp = result["amplitude"].as<double>();
    double shift = result["shift"].as<double>();

    // Source location
    int source_is = result["source_i"].as<int>();
    int source_js = result["source_j"].as<int>();
    // Create the source
    boost::multi_array<double, 3> f = waveSimCore::ricker(source_is, source_js, f_M, amp, shift, tmin, tmax, nt, nx, nz);

    // Create the velocity profile
    double r = result["radius"].as<double>();

    boost::multi_array<double, 2> vel = waveSimCore::get_profile(xmin, xmax, zmin, zmax, nx, nz, r);

    // Solve the wave equation
    boost::multi_array<double, 3> u = waveSimCore::wave_solver(vel, dt, dx, dz, nt, nx, nz, f);

    // Define the number of different levels for the contour plot
    int num_levels = result["num_levels"].as<int>();
    // Create the levels for the contour plot based on the min and max values of u
    double min_u = np::min(u);
    double max_u = np::max(u);
    std::vector<double> levels = matplot::linspace(min_u, max_u, num_levels);

    // Create the x and z axis for the contour plot and convert them to matplot format
    boost::multi_array<double, 1> x = np::linspace(xmin, xmax, nx);
    boost::multi_array<double, 1> z = np::linspace(zmin, zmax, nz);
    const boost::multi_array<double, 1> axis[2] = {x, z};
    std::vector<boost::multi_array<double, 2>> XcZ = np::meshgrid(axis, false, np::xy);

    matplot::vector_2d Xp = np::convert_to_matplot(XcZ[0]);
    matplot::vector_2d Zp = np::convert_to_matplot(XcZ[1]);

    // Create the plotter object and animate the results
    wavePlotter::Plotter my_plotter(u, Xp, Zp, num_levels, nt);
    my_plotter.setSaveDirectory(result["output_filename"].as<std::string>());

    // If you want to render a specific frame, use this:
    // my_plotter.renderFrame(int frame_index);

    // Renders the entire animation from start_frame to end_frame
    int start_frame = 0;
    int end_frame = nt;
    int fps = result["framerate"].as<int>();
    if (result["export"].as<bool>())
    {
        my_plotter.renderAllFrames(start_frame, end_frame);
    }
    else if (result["animate"].as<bool>())
    {
        my_plotter.animate(result["output_filename"].as<std::string>(), start_frame, end_frame, fps);
    }
    else
    {
        my_plotter.renderAllFrames(start_frame, end_frame);
    }
}