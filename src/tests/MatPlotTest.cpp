#include <matplot/matplot.h>
#include <thread>
#include "boost/multi_array.hpp"
#include "boost/array.hpp"
#include "CustomLibraries/np.hpp"

using namespace matplot;
void test_simple_plot()
{
    std::vector<double> x = linspace(-2 * pi, 2 * pi);
    std::vector<double> y = linspace(0, 4 * pi);
    auto [X, Y] = meshgrid(x, y);
    vector_2d Z =
        transform(X, Y, [](double x, double y)
                  { return sin(x) + cos(y); });
    contourf(X, Y, Z, 10);

    show();
}

matplot::vector_2d convert_to_matplot(const boost::multi_array<double, 2> &arr)
{
    std::vector<double> x = matplot::linspace(0, 0, arr.shape()[0]);
    std::vector<double> y = matplot::linspace(00, 0, arr.shape()[1]);
    matplot::vector_2d result = std::get<0>(matplot::meshgrid(x, y));
    for (size_t i = 0; i < arr.shape()[0]; i++)
    {
        for (size_t j = 0; j < arr.shape()[1]; j++)
        {
            result[i][j] = arr[i][j];
        }
    }
    return result;
}

void test_conversion()
{
    boost::multi_array<double, 1> x = np::linspace(0, 1, 100);
    boost::multi_array<double, 1> y = np::linspace(0, 1, 100);
    // x = np::pow(x, 2.0);
    // y = np::pow(y, 3.0);

    const boost::multi_array<double, 1> axis[2] = {x, y};
    std::vector<boost::multi_array<double, 2>> XcY = np::meshgrid(axis, false, np::xy);

    double dx, dy;
    dx = 1.0 / 100.0;
    dy = 1.0 / 100.0;

    boost::multi_array<double, 2> f = np::pow(XcY[0], 2.0) + XcY[0] * np::pow(XcY[1], 1.0);

    // g.push_back(np::gradient(XcY[0], {dx, dy}));
    // g.push_back(np::gradient(XcY[1], {dx, dy}));
    std::vector<boost::multi_array<double, 2>> gradf = np::gradient(f, {dx, dy});
    contourf(convert_to_matplot(XcY[0]), convert_to_matplot(XcY[1]), convert_to_matplot(gradf[0]), 10);
    show();
}

int main()
{
    // test_simple_plot();
    test_conversion();
    return 0;
}