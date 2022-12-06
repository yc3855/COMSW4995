//
// Created by Yan Cheng on 12/2/22.
//

#include <boost/multi_array.hpp>
#include <boost/array.hpp>
#include "CustomLibraries/np.hpp"
#include <cassert>
#include <iostream>

#include "CoreAlgorithm/helper_func.hpp"
#include "CoreAlgorithm/coeff.hpp"
#include "CoreAlgorithm/source.hpp"
#include "CoreAlgorithm/computational.hpp"
#include "CoreAlgorithm/solver.hpp"

void test_(){

    int nx = 101;
    int nz = 101;
    int nt = 750;

    double dx = 10.0;
    double dz = 10.0;
    double dt = 1e-3;

    double xmin = 0.0;
    double xmax = nx * dx;
    double zmin = 0.0;
    double zmax = nz * dz;
    double tmin = 0.0;
    double tmax = nt * dt;

    boost::multi_array<double, 2> sigma_1 = get_sigma_1(np::linspace(xmin, xmax, nx), dx, nx, nz, 3000.0);
    boost::multi_array<double, 2> sigma_2 = get_sigma_2(np::linspace(zmin, zmax, nz), dz, nx, nz, 3000.0);

//    for (int i = 0; i < nx; i++)
//    {
//        for (int j = 0; j < nz; j++)
//            std::cout << sigma_2[i][j] << " ";
//        std::cout << "\n";
//    }

    double f_M = 10.0;
    double amp = 1.0;
    double shift = 0.1;

    boost::multi_array<double, 3> f = ricker(10, 10, f_M, amp, shift, tmin, tmax, nt, nx, nz);
//    for (int i = 0; i < nt; i++)
//    {
//        std::cout << f[i][10][10] << " ";
//    }

    double r = 150;
    boost::multi_array<double, 2> vel = get_profile(xmin, xmax, zmin, zmax, nx, nz, r);

    boost::multi_array<double, 3> u = wave_solver(vel, dt, dx, dz, nt, nx, nz, f, sigma_1, sigma_2);
    std::cout << u[20][10][10] << "\n";

}


int main(){
    test_();
}