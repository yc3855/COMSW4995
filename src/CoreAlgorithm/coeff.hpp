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
    boost::multi_array<double, 2> sigma_1 = np::zeros(nx, nz);
    const double PML_width = n * dx;
    const double sigma_max = - c_max * log(R) * (m+1) / (PML_width**(m+1));

    // TODO: max: find the maximum element in 1D array
    const double x_0 = max(x) - PML_width;

    // each column of sigma_1 is a 1D array named "polynomial"
    boost::multi_array<double, 1> polynomial = np::zeros(nx);
    for (int i=0; i<nx; i++)
    {
        if (x[i] > x_0)
        {
            // TODO: Does math.h have an absolute value function?
            polynomial[i] = sigma_max * abs(x[i] - x_0)**m;
            polynomial[nx-i] = polynomial[i];
        }
        else
        {
            polynomial[i] = 0;
        }
    }

    // Copy 1D array into each column of 2D array
    for (int i=0; i<nx; i++)
        for (int j=0; i<nz; j++)
            sigma_1[i][j] = polynomial[i];

    return sigma_1;
}



boost::multi_array<double, 2> get_sigma_2(boost::multi_array<double, 1> z, double dz, int nx, int nz,
                                          double c_max, int n=10, double R=1e-3, int m=2)
{
    boost::multi_array<double, 2> sigma_2 = np::zeros(nx, nz);
    const double PML_width = n * dz;
    const double sigma_max = - c_max * log(R) * (m+1) / (PML_width**(m+1));

    // TODO: max: find the maximum element in 1D array
    const double z_0 = max(z) - PML_width;

    // each column of sigma_1 is a 1D array named "polynomial"
    boost::multi_array<double, 1> polynomial = np::zeros(nz);
    for (int j=0; j<nz; j++)
    {
        if (z[j] > z_0)
        {
            // TODO: Does math.h have an absolute value function?
            polynomial[j] = sigma_max * abs(z[j] - z_0)**m;
            polynomial[nz-j] = polynomial[j];
        }
        else
        {
            polynomial[j] = 0;
        }
    }

    // Copy 1D array into each column of 2D array
    for (int i=0; i<nx; i++)
        for (int j=0; i<nz; j++)
            sigma_1[i][j] = polynomial[j];

    return sigma_2;
}

#endif //WAVESIMC_COEFF_HPP
