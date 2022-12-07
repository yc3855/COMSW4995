// Standard IO libraries
#include <iostream>
#include <fstream>

#include "CustomLibraries/np.hpp"

#include <math.h>

#include "solver.hpp"
#include "computational.hpp"
#include "coeff.hpp"
#include "source.hpp"
#include "helper_func.hpp"

int main()
{
    double dx, dy, dz, dt;
    dx = 1.0;
    dy = 1.0;
    dz = 1.0;
    dt = 1.0;
    std::vector<boost::multi_array<double, 4>> my_arrays = np::gradient(A, {dx, dy, dz, dt});
    return 0;
}