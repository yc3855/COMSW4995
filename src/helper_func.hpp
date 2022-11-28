//
// Created by Yan Cheng on 11/28/22.
//

#ifndef WAVESIMC_HELPER_FUNC_HPP
#define WAVESIMC_HELPER_FUNC_HPP

// Helper functions
boost::multi_array<double, 2> d2fdx2(boost::multi_array<double, 2> f, double dx)
{
    return f;
}

boost::multi_array<double, 2> d2fdz2(boost::multi_array<double, 2> f, double dz)
{
    return f;
}

boost::multi_array<double, 2> dfdx(boost::multi_array<double, 2> f, double dx)
{
    return f;
}

boost::multi_array<double, 2> dfdz(boost::multi_array<double, 2> f, double dz)
{
    return f;
}

boost::multi_array<double, 2> divergence(boost::multi_array<double, 2> f1, boost::multi_array<double, 2> f2,
                                         double dx, double dz)
{
    return f1;
}


#endif //WAVESIMC_HELPER_FUNC_HPP
