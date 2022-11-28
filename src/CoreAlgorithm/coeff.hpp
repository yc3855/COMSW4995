//
// Created by Yan Cheng on 11/28/22.
//

#ifndef WAVESIMC_COEFF_HPP
#define WAVESIMC_COEFF_HPP

#include "CustomLibraries/np.hpp"
#include <math.h>

boost::multi_array<double, 2> get_sigma_1(boost::multi_array<double, 1> x, double dx, int nx, int nz,
                                          double c_max, int n=10, double R=1e-3, int m=2)
{
    boost::multi_array<double, 2> sigma_1;
    const double PML_width = n * dx;
    const double sigma_max = - c_max * log(R) * (m+1) / (PML_width**(m+1));

    // TODO: max: find the maximum element in 1D array
    const double x_0 = max(x) - PML_width;
    sigma_1 = np::zeros(nx, nz);


    return sigma_1;
}

#endif //WAVESIMC_COEFF_HPP
