//
// Created by Yan Cheng on 11/28/22.
//

#ifndef WAVESIMC_HELPER_FUNC_HPP
#define WAVESIMC_HELPER_FUNC_HPP

#include "CustomLibraries/np.hpp"

boost::multi_array<double, 2> dfdx(boost::multi_array<double, 2> f, double dx)
{
    std::vector<boost::multi_array<double, 2>> grad_f = np::gradient(f, {dx, dx});
    return grad_f[0];
}

boost::multi_array<double, 2> dfdz(boost::multi_array<double, 2> f, double dz)
{
    std::vector<boost::multi_array<double, 2>> grad_f = np::gradient(f, {dz, dz});
    return grad_f[1];
}

boost::multi_array<double, 2> d2fdx2(boost::multi_array<double, 2> f, double dx)
{
    boost::multi_array<double, 2> df = dfdx(f, dx);
    boost::multi_array<double, 2> df2 = dfdx(df, dx);
    return df2;
}

boost::multi_array<double, 2> d2fdz2(boost::multi_array<double, 2> f, double dz)
{
    boost::multi_array<double, 2> df = dfdz(f, dz);
    boost::multi_array<double, 2> df2 = dfdz(df, dz);
    return df2;
}

boost::multi_array<double, 2> divergence(boost::multi_array<double, 2> f1, boost::multi_array<double, 2> f2,
                                         double dx, double dz)
{
    boost::multi_array<double, 2> f_x = dfdx(f1, dx);
    boost::multi_array<double, 2> f_z = dfdz(f2, dz);
    // TODO: use element-wize add
    boost::multi_array<double, 2> div = f_x + f_z;
    return div;
}


#endif //WAVESIMC_HELPER_FUNC_HPP
