//
// Created by Yan Cheng on 11/28/22.
//

#ifndef WAVESIMC_SOURCE_HPP
#define WAVESIMC_SOURCE_HPP


boost::multi_array<double, 3> ricker(int i_s, int j_s, double f=10, double amp=1e0, double shift=0.1)
{
    const double pi = 3.141592654;

    boost::multi_array<double, 1> t = np::linspace(tmin, tmax, nt);

    // TODO: element-wise operators
    boost::multi_array<double, 1> pft2 = (pi * f * (t - shift))**2;
    boost::multi_array<double, 1> r = amp * (1 - 2 * pft2) * exp(-pft2);

    boost::multi_array<double, 1> x = np.zeros(nx);
    boost::multi_array<double, 1> z = np.zeros(nz);
    x[i_s] = 1.0;
    z[j_s] = 1.0;
    boost::multi_array<double, 3> TXZ = np::meshgrid(r, x, z, sparse=True, indexing='ij');

    return TXZ;
}

#endif //WAVESIMC_SOURCE_HPP
