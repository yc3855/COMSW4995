#ifndef NPTOMATPLOT_H_
#define NPTOMATPLOT_H_

#include <matplot/matplot.h>
#include <thread>
#include "boost/multi_array.hpp"
#include "boost/array.hpp"
/*!
 *  \addtogroup np
 *  @{
 */

namespace np
{
    //! Convert a 2D boost::multi_array to a matplot::vector_2d
    inline matplot::vector_2d convert_to_matplot(const boost::multi_array<double, 2> &arr)
    {
        std::vector<double> x = matplot::linspace(0, 0, arr.shape()[0]);
        std::vector<double> y = matplot::linspace(00, 0, arr.shape()[1]);
        matplot::vector_2d result = std::get<0>(matplot::meshgrid(x, y));
        // std::cout << "arr.shape()[0] = " << arr.shape()[0] << " arr.shape()[1] = " << arr.shape()[1] << std::endl;
        for (size_t i = 0; i < arr.shape()[0]; i++)
        {
            for (size_t j = 0; j < arr.shape()[1]; j++)
            {
                result[i][j] = arr[i][j];
            }
        }
        return result;
    }
}
#endif