#ifndef NP_H_
#define NP_H_

#include "boost/multi_array.hpp"
#include "boost/array.hpp"
#include "boost/cstdlib.hpp"
#include <type_traits>
#include <cassert>
#include <iostream>
#include <functional>

/*!
 *  \addtogroup np
 *  @{
 */

//! Custom implementation of numpy in C++
namespace np
{

    typedef double ndArrayValue;

    //! Gets the index of one element in a multi_array in one axis
    template <std::size_t ND>
    inline boost::multi_array<ndArrayValue, ND>::index
    getIndex(const boost::multi_array<ndArrayValue, ND> &m, const ndArrayValue *requestedElement, const unsigned short int direction)
    {
        int offset = requestedElement - m.origin();
        return (offset / m.strides()[direction] % m.shape()[direction] + m.index_bases()[direction]);
    }

    //! Gets the index of one element in a multi_array
    template <std::size_t ND>
    inline boost::array<typename boost::multi_array<ndArrayValue, ND>::index, ND>
    getIndexArray(const boost::multi_array<ndArrayValue, ND> &m, const ndArrayValue *requestedElement)
    {
        using indexType = boost::multi_array<ndArrayValue, ND>::index;
        boost::array<indexType, ND> _index;
        for (unsigned int dir = 0; dir < ND; dir++)
        {
            _index[dir] = getIndex(m, requestedElement, dir);
        }

        return _index;
    }

    //! Function to apply a function to all elements of a multi_array
    //! Simple overload
    template <typename Array, typename Element, typename Functor>
    inline void for_each(const boost::type<Element> &type_dispatch,
                         Array A, Functor &xform)
    {
        for_each(type_dispatch, A.begin(), A.end(), xform);
    }

    //! Function to apply a function to all elements of a multi_array
    template <typename Element, typename Functor>
    inline void for_each(const boost::type<Element> &, Element &Val, Functor &xform)
    {
        Val = xform(Val);
    }

    //! Function to apply a function to all elements of a multi_array
    template <typename Element, typename Iterator, typename Functor>
    inline void for_each(const boost::type<Element> &type_dispatch,
                         Iterator begin, Iterator end,
                         Functor &xform)
    {
        while (begin != end)
        {
            for_each(type_dispatch, *begin, xform);
            ++begin;
        }
    }

    //! Function to apply a function to all elements of a multi_array
    //! Simple overload
    template <typename Array, typename Functor>
    inline void for_each(Array &A, Functor xform)
    {
        // Dispatch to the proper function
        for_each(boost::type<typename Array::element>(), A.begin(), A.end(), xform);
    }

    //! Takes the gradient of a n-dimensional multi_array
    //! Todo: Actually implement the gradient calculation
    //! template <long unsigned int ND, typename... Args>
    template <long unsigned int ND>
    inline constexpr std::vector<boost::multi_array<double, ND>> gradient(boost::multi_array<double, ND> inArray, std::initializer_list<double> args)
    {
        // static_assert(args.size() == ND, "Number of arguments must match the number of dimensions of the array");
        using arrayIndex = boost::multi_array<double, ND>::index;

        using ndIndexArray = boost::array<arrayIndex, ND>;

        // constexpr std::size_t n = sizeof...(Args);
        std::size_t n = args.size();
        // std::tuple<Args...> store(args...);
        std::vector<double> arg_vector = args;
        boost::multi_array<double, ND> my_array;
        std::vector<boost::multi_array<double, ND>> output_arrays;
        for (std::size_t i = 0; i < n; i++)
        {
            boost::multi_array<double, ND> dfdh = inArray;
            output_arrays.push_back(dfdh);
        }

        ndArrayValue *p = inArray.data();
        ndIndexArray index;
        for (std::size_t i = 0; i < inArray.num_elements(); i++)
        {
            index = getIndexArray(inArray, p);
            /*
            std::cout << "Index: ";
            for (std::size_t j = 0; j < n; j++)
            {
                std::cout << index[j] << " ";
            }
            std::cout << "\n";
            */
            // Calculating the gradient now
            // j is the axis/dimension
            for (std::size_t j = 0; j < n; j++)
            {
                ndIndexArray index_high = index;
                double dh_high;
                if ((long unsigned int)index_high[j] < inArray.shape()[j] - 1)
                {
                    index_high[j] += 1;
                    dh_high = arg_vector[j];
                }
                else
                {
                    dh_high = 0;
                }
                ndIndexArray index_low = index;
                double dh_low;
                if (index_low[j] > 0)
                {
                    index_low[j] -= 1;
                    dh_low = arg_vector[j];
                }
                else
                {
                    dh_low = 0;
                }

                double dh = dh_high + dh_low;
                double gradient = (inArray(index_high) - inArray(index_low)) / dh;
                // std::cout << gradient << "\n";
                output_arrays[j](index) = gradient;
            }
            // std::cout << " value = " << inArray(index) << "  check = " << *p << std::endl;
            ++p;
        }
        return output_arrays;
    }

    //! Implements the numpy linspace function
    inline boost::multi_array<double, 1> linspace(double start, double stop, long unsigned int num)
    {
        double step = (stop - start) / (num - 1);
        boost::multi_array<double, 1> output(boost::extents[num]);
        for (std::size_t i = 0; i < num; i++)
        {
            output[i] = start + i * step;
        }
        return output;
    }

    //! Implements the numpy zeros function
    //! Todo: make it work for any number of dimensions
    inline boost::multi_array<double, 1> zeros(long unsigned int num)
    {
        boost::multi_array<double, 1> output(boost::extents[num]);
        for (std::size_t i = 0; i < num; i++)
        {
            output[i] = 0;
        }
        return output;
    }

    enum indexing
    {
        xy,
        ij
    };

    //! Implementation of meshgrid
    //! TODO: Implement sparsing=true
    //! If the indexing type is xx, then reverse the order of the first two elements of ci
    //! if the number of dimensions is 2 or 3
    //! In accordance with the numpy implementation
    template <long unsigned int ND>
    inline std::vector<boost::multi_array<double, ND>> meshgrid(const boost::multi_array<double, 1> (&cinput)[ND], bool sparsing = false, indexing indexing_type = xy)
    {
        using arrayIndex = boost::multi_array<double, ND>::index;
        using ndIndexArray = boost::array<arrayIndex, ND>;
        std::vector<boost::multi_array<double, ND>> output_arrays;
        boost::multi_array<double, 1> ci[ND];
        // Copy elements of cinput to ci, do the proper inversions
        for (std::size_t i = 0; i < ND; i++)
        {
            std::size_t source = i;
            if (indexing_type == xy && (ND == 3 || ND == 2))
            {
                switch (i)
                {
                case 0:
                    source = 1;
                    break;
                case 1:
                    source = 0;
                    break;
                default:
                    break;
                }
            }
            ci[i] = boost::multi_array<double, 1>();
            ci[i].resize(boost::extents[cinput[source].num_elements()]);
            ci[i] = cinput[source];
        }
        // Deducing the extents of the N-Dimensional output
        boost::detail::multi_array::extent_gen<ND> output_extents;
        std::vector<size_t> shape_list;
        for (std::size_t i = 0; i < ND; i++)
        {
            shape_list.push_back(ci[i].shape()[0]);
        }
        std::copy(shape_list.begin(), shape_list.end(), output_extents.ranges_.begin());

        // Creating the output arrays
        for (std::size_t i = 0; i < ND; i++)
        {
            boost::multi_array<double, ND> output_array(output_extents);
            ndArrayValue *p = output_array.data();
            ndIndexArray index;
            // Looping through the elements of the output array
            for (std::size_t j = 0; j < output_array.num_elements(); j++)
            {
                index = getIndexArray(output_array, p);
                boost::multi_array<double, 1>::index index_1d;
                index_1d = index[i];
                output_array(index) = ci[i][index_1d];
                ++p;
            }
            output_arrays.push_back(output_array);
        }
        return output_arrays;
    }

    //! Cretes a new array and fills it with the values of the result of the function called on the input array element-wise
    template <class T, long unsigned int ND>
    inline boost::multi_array<T, ND> element_wise_apply(const boost::multi_array<T, ND> &input_array, std::function<T(T)> func)
    {

        // Create output array copying extents
        using arrayIndex = boost::multi_array<double, ND>::index;
        using ndIndexArray = boost::array<arrayIndex, ND>;
        boost::detail::multi_array::extent_gen<ND> output_extents;
        std::vector<size_t> shape_list;
        for (std::size_t i = 0; i < ND; i++)
        {
            shape_list.push_back(input_array.shape()[i]);
        }
        std::copy(shape_list.begin(), shape_list.end(), output_extents.ranges_.begin());
        boost::multi_array<T, ND> output_array(output_extents);

        // Looping through the elements of the output array
        const T *p = input_array.data();
        ndIndexArray index;
        for (std::size_t i = 0; i < input_array.num_elements(); i++)
        {
            index = getIndexArray(input_array, p);
            output_array(index) = func(input_array(index));
            ++p;
        }
        return output_array;
    }

    // Complex operations

    template <class T, long unsigned int ND>
    inline boost::multi_array<T, ND> sqrt(const boost::multi_array<T, ND> &input_array)
    {
        std::function<T(T)> func = (T(*)(T))std::sqrt;
        return element_wise_apply(input_array, func);
    }
    template <class T>
    inline T sqrt(const T input)
    {
        return std::sqrt(input);
    }

    template <class T, long unsigned int ND>
    inline boost::multi_array<T, ND> exp(const boost::multi_array<T, ND> &input_array)
    {
        std::function<T(T)> func = (T(*)(T))std::exp;
        return element_wise_apply(input_array, func);
    }
    template <class T>
    inline T exp(const T input)
    {
        return std::exp(input);
    }

    template <class T, long unsigned int ND>
    inline boost::multi_array<T, ND> log(const boost::multi_array<T, ND> &input_array)
    {
        std::function<T(T)> func = std::log<T>();
        return element_wise_apply(input_array, func);
    }
    template <class T>
    inline T log(const T input)
    {
        return std::log(input);
    }
    template <class T, long unsigned int ND>
    inline boost::multi_array<T, ND> pow(const boost::multi_array<T, ND> &input_array, const T exponent)
    {
        std::function<T(T)> pow_func = [exponent](T input)
        { return std::pow(input, exponent); };
        return element_wise_apply(input_array, pow_func);
    }
    template <class T>
    inline T pow(const T input, const T exponent)
    {
        return std::pow(input, exponent);
    }

    //! Creates a new array in which the value at each index is the
    //! the result of the input function applied to an element of the left hand side array and one on the righ hand side array in the same index
    //! Outputs a copy of the result
    template <class T, long unsigned int ND>
    boost::multi_array<T, ND> element_wise_duo_apply(boost::multi_array<T, ND> const &lhs, boost::multi_array<T, ND> const &rhs, std::function<T(T, T)> func)
    {
        // Create output array copying extents
        using arrayIndex = boost::multi_array<double, ND>::index;
        using ndIndexArray = boost::array<arrayIndex, ND>;
        boost::detail::multi_array::extent_gen<ND> output_extents;
        std::vector<size_t> shape_list;
        for (std::size_t i = 0; i < ND; i++)
        {
            shape_list.push_back(lhs.shape()[i]);
        }
        std::copy(shape_list.begin(), shape_list.end(), output_extents.ranges_.begin());
        boost::multi_array<T, ND> output_array(output_extents);

        // Looping through the elements of the output array
        const T *p = lhs.data();
        ndIndexArray index;
        for (std::size_t i = 0; i < lhs.num_elements(); i++)
        {
            index = getIndexArray(lhs, p);
            output_array(index) = func(lhs(index), rhs(index));
            ++p;
        }
        return output_array;
    }
}

// Basic operators

template <class T, long unsigned int ND>
inline boost::multi_array<T, ND> operator*(boost::multi_array<T, ND> const &lhs, boost::multi_array<T, ND> const &rhs)
{
    std::function<T(T, T)> func = std::multiplies<T>();
    return np::element_wise_duo_apply(lhs, rhs, func);
}
template <class T, long unsigned int ND>
boost::multi_array<T, ND> operator+(boost::multi_array<T, ND> const &lhs, boost::multi_array<T, ND> const &rhs)
{
    std::function<T(T, T)> func = std::plus<T>();
    return np::element_wise_duo_apply(lhs, rhs, func);
}
template <class T, long unsigned int ND>
boost::multi_array<T, ND> operator-(boost::multi_array<T, ND> const &lhs, boost::multi_array<T, ND> const &rhs)
{
    std::function<T(T, T)> func = std::minus<T>();
    return np::element_wise_duo_apply(lhs, rhs, func);
}
template <class T, long unsigned int ND>
boost::multi_array<T, ND> operator/(boost::multi_array<T, ND> const &lhs, boost::multi_array<T, ND> const &rhs)
{
    std::function<T(T, T)> func = std::divides<T>();
    return np::element_wise_duo_apply(lhs, rhs, func);
}
/*! @} End of Doxygen Groups*/
#endif