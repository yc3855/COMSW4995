//
// Created by Yan Cheng on 11/28/22.
//

#ifndef WAVESIMC_COMPUTATIONAL_HPP
#define WAVESIMC_COMPUTATIONAL_HPP

boost::multi_array<double, 2> get_profile(double xmin, double xmax, double zmin, double zmax, int nx, int nz, double r)
{
    boost::multi_array<double, 2> c(boost::extents[nx][nz]);

    boost::multi_array<double, 1> x = np::linspace(xmin, xmax, nx);
    boost::multi_array<double, 1> z = np::linspace(zmin, zmax, nz);

    const boost::multi_array<double, 1> axis[2] = {x, z};
    std::vector<boost::multi_array<double, 2>> XZ = np::meshgrid(axis, false, np::ij);

    double x_0 = xmax / 2.0;
    double z_0 = zmax / 2.0;

    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < nz; j++)
        {
            if (np::pow(XZ[0][i][j]-x_0, 2.0) + np::pow(XZ[1][i][j]-z_0, 2.0) <= np::pow(r, 2.0))
                c[i][j] = 3.0;
            else
                c[i][j] = 3.0;
        }
    }

    return c;
}

#endif //WAVESIMC_COMPUTATIONAL_HPP
