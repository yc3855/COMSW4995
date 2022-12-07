//
// Created by Yan Cheng on 12/2/22.
//

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

std::string format_num(int num, int length = 8)
{
    std::string str_num = std::to_string(num);

    int str_length = str_num.length();
    for (int i = 0; i < length - str_length; i++)
        str_num = "0" + str_num;
    return str_num;
}

void test_()
{
    int num_levels = 100;
    int nx = 100;
    int nz = 100;
    int nt = 1000;

    double dx = 0.01;
    double dz = 0.01;
    double dt = 0.001;

    double xmin = 0.0;
    double xmax = nx * dx;
    double zmin = 0.0;
    double zmax = nz * dz;
    double tmin = 0.0;
    double tmax = nt * dt;

    double f_M = 10.0;
    double amp = 1e0;
    double shift = 0.1;

    boost::multi_array<double, 3> f = waveSimCore::ricker(50, 50, f_M, amp, shift, tmin, tmax, nt, nx, nz);

    double r = 150.0;
    boost::multi_array<double, 2> vel = waveSimCore::get_profile(xmin, xmax, zmin, zmax, nx, nz, r);

    boost::multi_array<double, 3> u = waveSimCore::wave_solver(vel, dt, dx, dz, nt, nx, nz, f);
    double min_u = np::min(u);
    double max_u = np::max(u);
    std::cout << "min_u = " << min_u << " max_u = " << max_u << "\n";
    std::vector<double> levels = matplot::linspace(min_u, max_u, num_levels);
    std::cout << u[20][10][10] << "\n";

    // boost::multi_array<double, 2> u20(boost::extents[nx][nz]);
    // for (int i = 0; i < nx; i++)
    // {
    //     for (int j = 0; j < nz; j++)
    //     {
    //         u20[i][j] = u[10][i][j];
    //         std::cout << u20[i][j] << " ";
    //     }
    //     std::cout << "\n";
    // }
    // std::cout << "\n";

    boost::multi_array<double, 1> x = np::linspace(xmin, xmax, nx);
    boost::multi_array<double, 1> z = np::linspace(zmin, zmax, nz);
    const boost::multi_array<double, 1> axis[2] = {x, z};
    std::vector<boost::multi_array<double, 2>> XcZ = np::meshgrid(axis, false, np::xy);

    matplot::vector_2d Xp = np::convert_to_matplot(XcZ[0]);
    matplot::vector_2d Zp = np::convert_to_matplot(XcZ[1]);
    // matplot::contourf(Xp, Zp, vel, levels, "", num_levels);
    // matplot::show();

    wavePlotter::Plotter my_plotter(u, Xp, Zp, num_levels, nt);
    // my_plotter.exportAllFrames(0, nt - 1);
    my_plotter.animate("output-test.mp4", 20, nt - 1, 30);

    // for (int i = 10; i < nt - 1; i++)
    // {
    //     matplot::vector_2d Up = np::convert_to_matplot(u[i]);
    //     matplot::contourf(Xp, Zp, Up, levels);
    //     matplot::save("output/contourf_" + format_num(i) + ".png");
    // }
}

int main()
{
    test_();
}