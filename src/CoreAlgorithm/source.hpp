//
// Created by Yan Cheng on 11/28/22.
//

#ifndef WAVESIMC_SOURCE_HPP
#define WAVESIMC_SOURCE_HPP


boost::multi_array<double, 3> ricker(int i_s, int j_s, double f, double amp, double shift,
                                     double tmin, double tmax, int nt, int nx, int nz)
{
    const double pi = 3.141592654;

    boost::multi_array<double, 1> t = np::linspace(tmin, tmax, nt);
    boost::multi_array<double, 1> pft2 = np::pow(pi * f * (t - shift), 2.0);
    boost::multi_array<double, 1> r = amp * (1.0 - 2.0 * pft2) * np::exp(-1.0 * pft2);


    int dimensions_x[] = {nx};
    boost::multi_array<double, 1> x = np::zeros<double>(dimensions_x);

    int dimensions_z[] = {nz};
    boost::multi_array<double, 1> z = np::zeros<double>(dimensions_z);

    x[i_s] = 1.0;
    z[j_s] = 1.0;

    const boost::multi_array<double, 1> axis[3] = {r, x, z};
    std::vector<boost::multi_array<double, 3>> RXZ = np::meshgrid(axis, false, np::ij);
    boost::multi_array<double, 3> source = RXZ[0] * RXZ[1] * RXZ[2];

//    boost::multi_array<double, 3> source(boost::extents[nt][nx][nz]);
//    for (int i=0; i<nt; i++)
//    {
//        for (int j=0; i<nx; j++)
//        {
//            for (int k=0; k<nz; k++)
//            {
//                if (j==i_s && k==j_s)
//                    source[i][j][k] = r[i];
//            }
//        }
//    }

    return source;
}

#endif //WAVESIMC_SOURCE_HPP
