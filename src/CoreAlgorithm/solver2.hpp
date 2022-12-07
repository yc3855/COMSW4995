//
// Created by Yan Cheng on 12/6/22.
//

#ifndef CORETESTS_CPP_SOLVER2_HPP
#define CORETESTS_CPP_SOLVER2_HPP

#include "CustomLibraries/np.hpp"

boost::multi_array<double, 3> wave_solver(boost::multi_array<double, 2> c,
                                          double dt, double dx, double dz, int nt, int nx, int nz,
                                          boost::multi_array<double, 3> f) {

    const boost::multi_array<double, 2> CX = np::pow(c * dt / dx, 2.0);
    const boost::multi_array<double, 2> CZ = np::pow(c * dt / dz, 2.0);
    const double Cf = np::pow(dt, 2.0);
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < nz; j++)
        {
            std::cout << CX[i][j] << " ";
        }
        std::cout << "\n";
    }

    int dimensions_1[] = {nt, nx, nz};
    boost::multi_array<double, 3> u = np::zeros<double>(dimensions_1);
    int dimensions_4[] = {nx, nz};
    boost::multi_array<double, 2> u_xx = np::zeros<double>(dimensions_4);
    int dimensions_5[] = {nx, nz};
    boost::multi_array<double, 2> u_zz = np::zeros<double>(dimensions_5);

    for (int n = 1; n < nt-1; n++)
    {
        for (int i = 1; i < nx-1; i++)
        {
            for (int j = 1; j < nz-1; j++)
            {
                u_xx[i][j] = u[n][i+1][j] + u[n][i-1][j] - 2.0 * u[n][i][j];
                u_zz[i][j] = u[n][i][j+1] + u[n][i][j-1] - 2.0 * u[n][i][j];
            }
        }

        // ! Update u
        for (int i = 0; i < nx; i++)
        {
            for (int j = 0; j < nz; j++)
            {
                u[n+1][i][j] = 2.0 * u[n][i][j] + CX[i][j] * u_xx[i][j] + CZ[i][j] * u_zz[i][j] + Cf * f[n][i][j] -u[n-1][i][j];
            }
        }

        // Dirichlet boundary condition
        for (int i = 0; i < nx; i++)
        {
            u[n+1][i][0] = 0.0;
            u[n+1][i][nz-1] = 0.0;
        }
        for (int j = 0; j < nz; j++)
        {
            u[n+1][0][j] = 0.0;
            u[n+1][nx-1][j] = 0.0;
        }
    }
    return u;
}

#endif //CORETESTS_CPP_SOLVER2_HPP
