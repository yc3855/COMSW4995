//
// Created by Yan Cheng on 11/28/22.
//

#ifndef WAVESIMC_SOLVER_HPP
#define WAVESIMC_SOLVER_HPP

#include "CustomLibraries/np.hpp"
#include "helper_func.hpp"

#include <cmath>

boost::multi_array<double, 3> wave_solver(boost::multi_array<double, 2> c,
                                          double dt, double dx, double dz, int nt, int nx, int nz,
                                          boost::multi_array<double, 3> f,
                                          boost::multi_array<double, 2> sigma_1, boost::multi_array<double, 2> sigma_2) {

    const boost::multi_array<double, 2> C1 = 1.0 + dt * (sigma_1 + sigma_2) * 0.5;
    const boost::multi_array<double, 2> C2 = sigma_1 * sigma_2 * np::pow(dt, 2.0) - 2.0;
    const boost::multi_array<double, 2> C3 = 1.0 - dt * (sigma_1 + sigma_2) * 0.5;
    const boost::multi_array<double, 2> C4 = np::pow(dt * c, 2.0);
    const boost::multi_array<double, 2> C5 = 1.0 + dt * sigma_1 * 0.5;
    const boost::multi_array<double, 2> C6 = 1.0 + dt * sigma_2 * 0.5;
    const boost::multi_array<double, 2> C7 = 1.0 - dt * sigma_1 * 0.5;
    const boost::multi_array<double, 2> C8 = 1.0 - dt * sigma_2 * 0.5;



    int dimensions_1[] = {nt, nx, nz};
    boost::multi_array<double, 3> u = np::zeros<double>(dimensions_1);

    int dimensions_2[] = {nx, nz};
    boost::multi_array<double, 2> q_1 = np::zeros<double>(dimensions_2);
    int dimensions_3[] = {nx, nz};
    boost::multi_array<double, 2> q_2 = np::zeros<double>(dimensions_3);

    int dimensions_4[] = {nx, nz};
    boost::multi_array<double, 2> u_xx = np::zeros<double>(dimensions_4);
    int dimensions_5[] = {nx, nz};
    boost::multi_array<double, 2> u_zz = np::zeros<double>(dimensions_5);

    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < nz; j++)
        {
            q_1[i][j] = 0.0;
            q_2[i][j] = 0.0;
        }
    }

    boost::multi_array<double, 2> f_n(boost::extents[nx][nz]);
    boost::multi_array<double, 2> u_n(boost::extents[nx][nz]);
    boost::multi_array<double, 2> u_n_1(boost::extents[nx][nz]);

    for (int n = 1; n < nt-1; n++)
    {

        for (int i = 0; i < nx; i++)
        {
            for (int j = 0; j < nz; j++)
            {
                f_n[i][j] = f[n][i][j];
                u_n[i][j] = u[n][i][j];
                u_n_1[i][j] = u[n-1][i][j];
            }
        }

        boost::multi_array<double, 2> div = divergence(q_1 * sigma_1, q_2 * sigma_2, dx, dz);
        boost::multi_array<double, 2> dq_1dx = dfdx(q_1, dx);
        boost::multi_array<double, 2> dq_2dz = dfdz(q_2, dz);
        u_xx = d2fdx2(u_n, dx); // (nx, nz)
        u_zz = d2fdz2(u_n, dz); // (nx, nz)

        for (int i = 0; i < nx; i++)
        {
            for (int j = 0; j < nz; j++)
            {
                u[n+1][i][j] = (C4[i][j] * ((u_xx[i][j] / np::pow(dx, 2.0)) + (u_zz[i][j] / np::pow(dz, 2.0))
                                      - div[i][j] + sigma_2[i][j] * dq_1dx[i][j] + sigma_1[i][j] * dq_2dz[i][j]
                                      + f_n[i][j]) - (C2[i][j] * u_n[i][j]) - (C3[i][j] * u_n_1[i][j])) / C1[i][j];

            }
        }

        boost::multi_array<double, 2> dudx = dfdx(u_n, dx);
        boost::multi_array<double, 2> dudz = dfdx(u_n, dz);

        for (int i = 0; i < nx; i++)
        {
            for (int j = 0; j < nz; j++)
            {
                q_1[i][j] = (dt * dudx[i][j] + C7[i][j] * q_1[i][j]) / C5[i][j];
                q_2[i][j] = (dt * dudz[i][j] + C8[i][j] * q_2[i][j]) / C6[i][j];
            }
        }
//
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

#endif //WAVESIMC_SOLVER_HPP
