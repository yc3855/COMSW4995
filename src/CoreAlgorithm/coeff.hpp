//
// Created by Yan Cheng on 11/28/22.
//

#ifndef WAVESIMC_COEFF_HPP
#define WAVESIMC_COEFF_HPP

#include "CustomLibraries/np.hpp"
#include <math.h>

/*!
 *  \addtogroup waveSimCore
 *  @{
 */

namespace waveSimCore
{
    boost::multi_array<double, 2> get_sigma_1(boost::multi_array<double, 1> x, double dx, int nx, int nz,
                                              double c_max, int n = 15, double R = 1e-3, double m = 2.0)
    {
        boost::multi_array<double, 2> sigma_1(boost::extents[nx][nz]);
        const double PML_width = n * dx;

        const double sigma_max = -c_max * log(R) * (m + 1.0) / np::pow(PML_width, m + 1.0);

        const double x_0 = np::max(x) - PML_width;

        boost::multi_array<double, 1> polynomial(boost::extents[nx]);

        for (int i = 0; i < nx; i++)
        {
            if (x[i] > x_0)
            {
                polynomial[i] = sigma_max * np::pow(np::abs(x[i] - x_0), m);
                polynomial[nx - 1 - i] = polynomial[i];
            }
            else
            {
                polynomial[i] = 0;
            }
        }
        // Copy 1D array into each column of 2D array
        for (int i = 0; i < nx; i++)
            for (int j = 0; j < nz; j++)
                sigma_1[i][j] = polynomial[i];

        return sigma_1;
    }

    boost::multi_array<double, 2> get_sigma_2(boost::multi_array<double, 1> z, double dz, int nx, int nz,
                                              double c_max, int n = 10, double R = 1e-3, double m = 2.0)
    {
        boost::multi_array<double, 2> sigma_2(boost::extents[nx][nz]);
        const double PML_width = n * dz;
        const double sigma_max = -c_max * log(R) * (m + 1.0) / np::pow(PML_width, m + 1.0);

        const double z_0 = np::max(z) - PML_width;

        boost::multi_array<double, 1> polynomial(boost::extents[nz]);
        for (int j = 0; j < nz; j++)
        {
            if (z[j] > z_0)
            {
                // TODO: Does math.h have an absolute value function?
                polynomial[j] = sigma_max * np::pow(np::abs(z[j] - z_0), m);
                polynomial[nz - 1 - j] = polynomial[j];
            }
            else
            {
                polynomial[j] = 0;
            }
        }

        // Copy 1D array into each column of 2D array
        for (int i = 0; i < nx; i++)
            for (int j = 0; j < nz; j++)
                sigma_2[i][j] = polynomial[j];

        return sigma_2;
    }
}
#endif // WAVESIMC_COEFF_HPP
