//
// Created by Yan Cheng on 12/2/22.
//

#include <boost/multi_array.hpp>
#include <boost/array.hpp>
#include "CustomLibraries/np.hpp"
#include <cassert>
#include <iostream>

#include "CoreAlgorithm/helper_func.hpp"
#include "CoreAlgorithm/coeff.hpp"
#include "CoreAlgorithm/source.hpp"
#include "CoreAlgorithm/computational.hpp"
//#include "CoreAlgorithm/solver.hpp"

void test_(){
    boost::multi_array<double, 2> sigma_1 = get_sigma_1(np::linspace(0.0, 1.0, 100), 1.0 / 100.0, 100, 100, 3000.0);

        int nx = 100;
        int nz  = 100;
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < nz; j++)
            std::cout << sigma_1[i][j] << " ";
        std::cout << "\n";
    }
}

int main(){
    test_();
}