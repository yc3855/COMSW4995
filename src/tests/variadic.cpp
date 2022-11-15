#include "boost/multi_array.hpp"
#include "boost/array.hpp"
#include "CustomLibraries/np.hpp"
#include <cassert>
#include <iostream>

void test_gradient()
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

    boost::multi_array<double, 1> x = np::linspace(0, 1, 5);
    std::vector<boost::multi_array<double, 1>> gradf = np::gradient(x, {1.0});
    for (int i = 0; i < 5; i++)
    {
        std::cout << gradf[0][i] << ",";
    }
    std::cout << "\n";
    // np::print(std::cout, my_arrays[0]);
}

void test_meshgrid()
{
    boost::multi_array<double, 1> x = np::linspace(0, 1, 5);
    boost::multi_array<double, 1> y = np::linspace(0, 1, 5);
    boost::multi_array<double, 1> z = np::linspace(0, 1, 5);
    boost::multi_array<double, 1> t = np::linspace(0, 1, 5);
    const boost::multi_array<double, 1> axis[4] = {x, y, z, t};
    std::vector<boost::multi_array<double, 4>> my_arrays = np::meshgrid(axis, false, np::xy);
    // np::print(std::cout, my_arrays[0]);
    int nx = 3;
    int ny = 2;
    boost::multi_array<double, 1> x2 = np::linspace(0, 1, nx);
    boost::multi_array<double, 1> y2 = np::linspace(0, 1, ny);
    const boost::multi_array<double, 1> axis2[2] = {x2, y2};
    std::vector<boost::multi_array<double, 2>> my_arrays2 = np::meshgrid(axis2, false, np::xy);
    std::cout << "xv\n";
    for (int i = 0; i < ny; i++)
    {
        for (int j = 0; j < nx; j++)
        {
            std::cout << my_arrays2[0][i][j] << " ";
        }
        std::cout << "\n";
    }
    std::cout << "yv\n";
    for (int i = 0; i < ny; i++)
    {
        for (int j = 0; j < nx; j++)
        {
            std::cout << my_arrays2[1][i][j] << " ";
        }
        std::cout << "\n";
    }
}

int main()
{
    test_gradient();
    test_meshgrid();
}