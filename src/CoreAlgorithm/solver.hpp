//
// Created by Yan Cheng on 11/28/22.
//

#ifndef WAVESIMC_SOLVER_HPP
#define WAVESIMC_SOLVER_HPP

#include "CustomLibraries/np.hpp"
#include "helper_func.hpp"

boost::multi_array<double, 3> wave_solver(boost::multi_array<double, 2> c,
                                          double dt, double dx, double dz, int nt, int nx, int nz,
                                          boost::multi_array<double, 3> f,
                                          boost::multi_array<double, 2> sigma_1, boost::multi_array<double, 2> sigma_2)
{
    // TODO: "same shape" functionality of np::zeros
    boost::multi_array<double, 3> u(boost::extents[nt][nx][nz]);
    boost::multi_array<double, 2> u_xx(boost::extents[nx][nz]);
    boost::multi_array<double, 2> u_zz(boost::extents[nx][nz]);
    boost::multi_array<double, 2> q_1(boost::extents[nx][nz]);
    boost::multi_array<double, 2> q_2(boost::extents[nx][nz]);

    const boost::multi_array<double, 2> C1 = 1.0 + dt * (sigma_1 + sigma_2) / 2.0;
    const boost::multi_array<double, 2> C2 = sigma_1 * sigma_2 * np::pow(dt, 2.0) - 2.0;
    const boost::multi_array<double, 2> C3 = 1.0 - dt * (sigma_1 + sigma_2) / 2.0;
    const boost::multi_array<double, 2> C4 = np::pow(dt * c, 2.0);
    const boost::multi_array<double, 2> C5 = 1.0 + dt * sigma_1 / 2.0;
    const boost::multi_array<double, 2> C6 = 1.0 + dt * sigma_2 / 2.0;
    const boost::multi_array<double, 2> C7 = 1.0 - dt * sigma_1 / 2.0;
    const boost::multi_array<double, 2> C8 = 1.0 - dt * sigma_2 / 2.0;

    for (int n = 0; n < nt; n++)
    {
        u_xx = d2fdx2(u[n], dx);
        u_zz = d2fdz2(u[n], dz);

        u[n+1] = (C4 * ((u_xx / np::pow(dx, 2.0)) + (u_zz / np::pow(dz, 2.0))
                    - divergence(q_1 * sigma_1, q_2 * sigma_2, dx, dz)
                    + (sigma_2 * dfdx(q_1, dx)) + (sigma_1 * dfdz(q_2, dz) + f[n]))
                    - (C2 * u[n]) - (C3 * u[n-1])) / C1;

        q_1 = (dt * dfdx(u[n], dx) + C7 * q_1) / C5;
        q_2 = (dt * dfdz(u[n], dx) + C8 * q_2) / C6;

        // Dirichlet boundary condition
        for (int i = 0; i < nx; i++)
        {
            u[n+1][i][0] = 0;
            u[n+1][i][nx-1] = 0;
        }
        for (int j = 0; j < nz; j++)
        {
            u[n+1][0][j] = 0;
            u[n+1][nz-1][j] = 0;
        }
    }
    return u;
}

#endif //WAVESIMC_SOLVER_HPP
