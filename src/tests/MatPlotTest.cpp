#include <matplot/matplot.h>
#include <thread>
#include "boost/multi_array.hpp"
#include "boost/array.hpp"
#include "CustomLibraries/np.hpp"
#include "CustomLibraries/np_to_matplot.hpp"

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
    matplot::vector_2d Xp = np::convert_to_matplot(XcY[0]);
    matplot::vector_2d Yp = np::convert_to_matplot(XcY[1]);
    matplot::vector_2d X = matplot::meshgrid(linspace(0, 1, 100), linspace(0, 1, 100)).first;
    matplot::vector_2d Y = matplot::meshgrid(linspace(0, 1, 100), linspace(0, 1, 100)).second;

    matplot::vector_2d Z = np::convert_to_matplot(gradf[1]);
    std::cout << "X.size() = " << X.size() << std::endl;
    std::cout << "Y.size() = " << Y.size() << std::endl;
    std::cout << "Z.size() = " << Z.size() << std::endl;
    for (size_t i = 0; i < X.size(); i++)
    {
        for (size_t j = 0; j < X[i].size(); j++)
        {
            std::cout << "X[" << i << "][" << j << "] = " << X[i][j] << " | XP[" << i << "][" << j << "] = " << Xp[i][j] << " | YP[" << i << "][" << j << "] = " << Yp[i][j] << std::endl;
        }
    }
    contourf(Xp, Yp, Z, 10);
    show();
}

int main()
{
    // test_simple_plot();
    test_conversion();
    return 0;
}