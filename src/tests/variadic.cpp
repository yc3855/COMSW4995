#include "boost/multi_array.hpp"
#include "boost/array.hpp"
#include "CustomLibraries/np.hpp"
#include <cassert>
#include <iostream>

int main()
{
    // Create a 4D array that is 3 x 4 x 2 x 1
    typedef boost::multi_array<double, 4>::index index;
    boost::multi_array<double, 4> A(boost::extents[3][4][2][2]);

    // Assign values to the elements
    int values = 0;
    for (index i = 0; i != 3; ++i)
        for (index j = 0; j != 4; ++j)
            for (index k = 0; k != 2; ++k)
                for (index l = 0; l != 2; ++l)
                    A[i][j][k][l] = values++;

    // Verify values
    int verify = 0;
    for (index i = 0; i != 3; ++i)
        for (index j = 0; j != 4; ++j)
            for (index k = 0; k != 2; ++k)
                for (index l = 0; l != 2; ++l)
                    assert(A[i][j][k][l] == verify++);

    double dx, dy, dz, dt;
    dx = 1.0;
    dy = 1.0;
    dz = 1.0;
    dt = 1.0;
    std::vector<boost::multi_array<double, 4>> my_arrays = np::gradient(A, {dx, dy, dz, dt});
    // np::print(std::cout, my_arrays[0]);
    return 0;
}