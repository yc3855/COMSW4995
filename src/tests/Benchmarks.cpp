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
#include <chrono>

#include "CoreAlgorithm/helper_func.hpp"
#include "CoreAlgorithm/coeff.hpp"
#include "CoreAlgorithm/source.hpp"
#include "CoreAlgorithm/computational.hpp"
#include "CoreAlgorithm/solver.hpp"

void benchmark_solver()
{
    int nx = 101;
    int nz = 101;
    int nt = 750;

    double dx = 10;
    double dz = 10;
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
    auto start = std::chrono::high_resolution_clock::now();
    boost::multi_array<double, 3> u = waveSimCore::wave_solver(vel, dt, dx, dz, nt, nx, nz, f);
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "Time taken by the solving function: " << duration.count() << " microseconds" << std::endl;
}

void benchmark_gradient()
{
    int sizex = 100;
    int sizey = 100;
    int sizez = 100;
    int sizet = 10;
    double dx = 0.1;
    double dy = 0.1;
    double dz = 0.1;
    double dt = 0.1;
    int axis[4] = {sizex, sizey, sizez, sizet};
    boost::multi_array<double, 4> nd_array = np::zeros<double>(axis);
    auto start = std::chrono::high_resolution_clock::now();
    auto grad = np::gradient(nd_array, {dx, dy, dz, dt});
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "Time taken by the gradient function: " << duration.count() << " microseconds" << std::endl;
}

int main()
{
    benchmark_solver();
    benchmark_gradient();
}