#include <boost/multi_array.hpp>
#include <boost/array.hpp>
#include "CustomLibraries/np.hpp"
#include "CustomLibraries/np_to_matplot.hpp"
#include "CustomLibraries/wavePlotter.hpp"
#include <matplot/matplot.h>
#include <cassert>
#include <iostream>
#include <sstream>

#include "CoreAlgorithm/helper_func.hpp"
#include "CoreAlgorithm/coeff.hpp"
#include "CoreAlgorithm/source.hpp"
#include "CoreAlgorithm/computational.hpp"
#include "CoreAlgorithm/solver.hpp"

int main()
{
    // Define the constants for the simulation

    // Number of x and z grid points
    int nx = 100;
    int nz = 100;
    // Number of time steps
    int nt = 1000;

    // Differentiation values
    double dx = 0.01;
    double dz = 0.01;
    double dt = 0.001;

    // Define the domain
    double xmin = 0.0;
    double xmax = nx * dx;
    double zmin = 0.0;
    double zmax = nz * dz;
    double tmin = 0.0;
    double tmax = nt * dt;

    // Define the source parameters
    double f_M = 10.0;
    double amp = 1e0;
    double shift = 0.1;

    // Source location
    int source_is = 50;
    int source_js = 50;

    // Create the source
    boost::multi_array<double, 3> f = waveSimCore::ricker(source_is, source_js, f_M, amp, shift, tmin, tmax, nt, nx, nz);

    // Create the velocity profile
    double r = 150.0;
    boost::multi_array<double, 2> vel = waveSimCore::get_profile(xmin, xmax, zmin, zmax, nx, nz, r);

    // Solve the wave equation
    boost::multi_array<double, 3> u = waveSimCore::wave_solver(vel, dt, dx, dz, nt, nx, nz, f);

    // Define the number of different levels for the contour plot
    int num_levels = 100;
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

    // If you want to render a specific frame, use this:
    // my_plotter.renderFrame(int frame_index);

    // Renders the entire animation from start_frame to end_frame
    int start_frame = 20;
    int end_frame = nt - 1;
    int fps = 30;
    my_plotter.animate("example-wave.mp4", start_frame, end_frame, fps);

    // The animation will be saved in .
    // Frames will be saved to ./output
}