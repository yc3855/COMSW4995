#ifndef NP_H_
#define NP_H_

#include "boost/multi_array.hpp"
#include "boost/array.hpp"
#include "boost/cstdlib.hpp"
#include <type_traits>
#include <cassert>
#include <iostream>
#include <functional>
#include <type_traits>

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
    template <typename T, long unsigned int ND>
    requires std::is_floating_point<T>::value inline constexpr std::vector<boost::multi_array<T, ND>> gradient(boost::multi_array<T, ND> inArray, std::initializer_list<T> args)
    {
        // static_assert(args.size() == ND, "Number of arguments must match the number of dimensions of the array");
        using arrayIndex = boost::multi_array<T, ND>::index;

        using ndIndexArray = boost::array<arrayIndex, ND>;

        // constexpr std::size_t n = sizeof...(Args);
        std::size_t n = args.size();
        // std::tuple<Args...> store(args...);
        std::vector<T> arg_vector = args;
        boost::multi_array<T, ND> my_array;
        std::vector<boost::multi_array<T, ND>> output_arrays;
        for (std::size_t i = 0; i < n; i++)
        {
            boost::multi_array<T, ND> dfdh = inArray;
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
                T dh_high;
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
                T dh_low;
                if (index_low[j] > 0)
                {
                    index_low[j] -= 1;
                    dh_low = arg_vector[j];
                }
                else
                {
                    dh_low = 0;
                }

                T dh = dh_high + dh_low;
                T gradient = (inArray(index_high) - inArray(index_low)) / dh;
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
    template <typename T, long unsigned int ND>
    requires std::is_arithmetic<T>::value inline constexpr std::vector<boost::multi_array<T, ND>> meshgrid(const boost::multi_array<T, 1> (&cinput)[ND], bool sparsing = false, indexing indexing_type = xy)
    {
        using arrayIndex = boost::multi_array<T, ND>::index;
        using oneDArrayIndex = boost::multi_array<T, 1>::index;
        using ndIndexArray = boost::array<arrayIndex, ND>;
        std::vector<boost::multi_array<T, ND>> output_arrays;
        boost::multi_array<T, 1> ci[ND];
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
            ci[i] = boost::multi_array<T, 1>();
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
            boost::multi_array<T, ND> output_array(output_extents);
            ndArrayValue *p = output_array.data();
            ndIndexArray index;
            // Looping through the elements of the output array
            for (std::size_t j = 0; j < output_array.num_elements(); j++)
            {
                index = getIndexArray(output_array, p);
                oneDArrayIndex index_1d;
                index_1d = index[i];
                output_array(index) = ci[i][index_1d];
                ++p;
            }
            output_arrays.push_back(output_array);
        }
        return output_arrays;
    }

    //! Creates a new array and fills it with the values of the result of the function called on the input array element-wise
    template <class T, long unsigned int ND>
    requires std::is_arithmetic<T>::value inline constexpr boost::multi_array<T, ND> element_wise_apply(const boost::multi_array<T, ND> &input_array, std::function<T(T)> func)
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

    //! Implements the numpy sqrt function on multi arrays
    template <class T, long unsigned int ND>
    requires std::is_arithmetic<T>::value inline constexpr boost::multi_array<T, ND> sqrt(const boost::multi_array<T, ND> &input_array)
    {
        std::function<T(T)> func = (T(*)(T))std::sqrt;
        return element_wise_apply(input_array, func);
    }

    //! Implements the numpy sqrt function on scalars
    template <class T>
    requires std::is_arithmetic<T>::value inline constexpr T sqrt(const T input)
    {
        return std::sqrt(input);
    }

    //! Implements the numpy exp function on multi arrays
    template <class T, long unsigned int ND>
    requires std::is_arithmetic<T>::value inline constexpr boost::multi_array<T, ND> exp(const boost::multi_array<T, ND> &input_array)
    {
        std::function<T(T)> func = (T(*)(T))std::exp;
        return element_wise_apply(input_array, func);
    }

    //! Implements the numpy exp function on scalars
    template <class T>
    requires std::is_arithmetic<T>::value inline constexpr T exp(const T input)
    {
        return std::exp(input);
    }

    //! Implements the numpy log function on multi arrays
    template <class T, long unsigned int ND>
    requires std::is_arithmetic<T>::value inline constexpr boost::multi_array<T, ND> log(const boost::multi_array<T, ND> &input_array)
    {
        std::function<T(T)> func = std::log<T>();
        return element_wise_apply(input_array, func);
    }

    //! Implements the numpy log function on scalars
    template <class T>
    requires std::is_arithmetic<T>::value inline constexpr T log(const T input)
    {
        return std::log(input);
    }

    //! Implements the numpy pow function on multi arrays
    template <class T, long unsigned int ND>
    requires std::is_arithmetic<T>::value inline constexpr boost::multi_array<T, ND> pow(const boost::multi_array<T, ND> &input_array, const T exponent)
    {
        std::function<T(T)> pow_func = [exponent](T input)
        { return std::pow(input, exponent); };
        return element_wise_apply(input_array, pow_func);
    }

    //! Implements the numpy pow function on scalars
    template <class T>
    requires std::is_arithmetic<T>::value inline constexpr T pow(const T input, const T exponent)
    {
        return std::pow(input, exponent);
    }

    //! Creates a new array in which the value at each index is the
    //! the result of the input function applied to an element of the left hand side array and one on the righ hand side array in the same index
    //! Outputs a copy of the result
    template <class T, long unsigned int ND>
    inline constexpr boost::multi_array<T, ND> element_wise_duo_apply(boost::multi_array<T, ND> const &lhs, boost::multi_array<T, ND> const &rhs, std::function<T(T, T)> func)
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

    //! Implements the numpy zeros function for an n-dimensionl multi array
    template <typename T, typename inT, long unsigned int ND>
    requires std::is_integral<inT>::value && std::is_arithmetic<T>::value inline constexpr boost::multi_array<T, ND> zeros(inT (&dimensions_input)[ND])
    {
        // Deducing the extents of the N-Dimensional output
        boost::detail::multi_array::extent_gen<ND> output_extents;
        std::vector<size_t> shape_list;
        for (std::size_t i = 0; i < ND; i++)
        {
            shape_list.push_back(dimensions_input[i]);
        }
        std::copy(shape_list.begin(), shape_list.end(), output_extents.ranges_.begin());
        // Applying a function to return zero always to all of its elements
        boost::multi_array<T, ND> output_array(output_extents);
        std::function<T(T)> zero_func = [](T input)
        { return 0; };
        return element_wise_apply(output_array, zero_func);
    }

    //! Implements the numpy max function for an n-dimensionl multi array
    template <typename T, long unsigned int ND>
    requires std::is_arithmetic<T>::value inline constexpr T max(boost::multi_array<T, ND> const &input_array)
    {
        T max = 0;
        bool max_not_set = true;
        const T *data_pointer = input_array.data();
        for (std::size_t i = 0; i < input_array.num_elements(); i++)
        {
            T element = *data_pointer;
            if (max_not_set || element > max)
            {
                max = element;
                max_not_set = false;
            }
            ++data_pointer;
        }
        return max;
    }

    //! Implements the numpy max function for an variadic number of arguments
    template <class T, class... Ts, class = std::enable_if_t<(std::is_same_v<T, Ts> && ...)>>
    requires std::is_arithmetic<T>::value inline constexpr T max(T input1, Ts... inputs)
    {
        T max = input1;
        for (T input : {inputs...})
        {
            if (input > max)
            {
                max = input;
            }
        }
        return max;
    }

    //! Implements the numpy min function for an n-dimensionl multi array
    template <typename T, long unsigned int ND>
    requires std::is_arithmetic<T>::value inline constexpr T min(boost::multi_array<T, ND> const &input_array)
    {
        T min = 0;
        bool min_not_set = true;
        const T *data_pointer = input_array.data();
        for (std::size_t i = 0; i < input_array.num_elements(); i++)
        {
            T element = *data_pointer;
            if (min_not_set || element < min)
            {
                min = element;
                min_not_set = false;
            }
            ++data_pointer;
        }
        return min;
    }

    //! Implements the numpy min function for an variadic number of arguments
    template <class T, class... Ts, class = std::enable_if_t<(std::is_same_v<T, Ts> && ...)>>
    inline constexpr T min(T input1, Ts... inputs) requires std::is_arithmetic<T>::value
    {
        T min = input1;
        for (T input : {inputs...})
        {
            if (input < min)
            {
                min = input;
            }
        }
        return min;
    }

    //! Implements the numpy abs function for an n-dimensionl multi array
    template <typename T, long unsigned int ND>
    requires std::is_arithmetic<T>::value inline constexpr boost::multi_array<T, ND> abs(boost::multi_array<T, ND> const &input_array)
    {
        std::function<T(T)> abs_func = [](T input)
        { return std::abs(input); };
        return element_wise_apply(input_array, abs_func);
    }

    //! Implements the numpy abs function for a scalar
    template <typename T>
    requires std::is_arithmetic<T>::value inline constexpr T abs(T input)
    {
        return std::abs(input);
    }

    //! Slices the array through one dimension and returns a ND - 1 dimensional array
    template <typename T, long unsigned int ND>
    requires std::is_arithmetic<T>::value inline constexpr boost::multi_array<T, ND - 1> slice(boost::multi_array<T, ND> const &input_array, std::size_t slice_index)
    {

        // Deducing the extents of the N-Dimensional output
        boost::detail::multi_array::extent_gen<ND - 1> output_extents;
        std::vector<size_t> shape_list;
        for (std::size_t i = 1; i < ND; i++)
        {
            shape_list.push_back(input_array.shape()[i]);
        }
        std::copy(shape_list.begin(), shape_list.end(), output_extents.ranges_.begin());

        boost::multi_array<T, ND - 1> output_array(output_extents);

        const T *p = input_array.data();
        boost::array<std::size_t, ND> index;
        for (std::size_t i = 0; i < input_array.num_elements(); i++)
        {
            index = getIndexArray(input_array, p);
            output_array(index) = input_array[slice_index](index);
            p++;
        }
        return output_array;
    }

}

// Override of operators in the boost::multi_array class to make them more np-like
// Basic operators
// All of the are element-wise

// Multiplication operator
//! Multiplication operator between two multi arrays, element-wise
template <class T, long unsigned int ND>
inline boost::multi_array<T, ND> operator*(boost::multi_array<T, ND> const &lhs, boost::multi_array<T, ND> const &rhs)
{
    std::function<T(T, T)> func = std::multiplies<T>();
    return np::element_wise_duo_apply(lhs, rhs, func);
}

//! Multiplication operator between a multi array and a scalar
template <class T, long unsigned int ND>
inline boost::multi_array<T, ND> operator*(T const &lhs, boost::multi_array<T, ND> const &rhs)
{
    std::function<T(T)> func = [lhs](T item)
    { return lhs * item; };
    return np::element_wise_apply(rhs, func);
}
//! Multiplication operator between a multi array and a scalar
template <class T, long unsigned int ND>
inline boost::multi_array<T, ND> operator*(boost::multi_array<T, ND> const &lhs, T const &rhs)
{
    return rhs * lhs;
}

// Plus operator
//! Addition operator between two multi arrays, element wise
template <class T, long unsigned int ND>
boost::multi_array<T, ND> operator+(boost::multi_array<T, ND> const &lhs, boost::multi_array<T, ND> const &rhs)
{
    std::function<T(T, T)> func = std::plus<T>();
    return np::element_wise_duo_apply(lhs, rhs, func);
}

//! Addition operator between a multi array and a scalar
template <class T, long unsigned int ND>
inline boost::multi_array<T, ND> operator+(T const &lhs, boost::multi_array<T, ND> const &rhs)
{
    std::function<T(T)> func = [lhs](T item)
    { return lhs + item; };
    return np::element_wise_apply(rhs, func);
}

//! Addition operator between a scalar and a multi array
template <class T, long unsigned int ND>
inline boost::multi_array<T, ND> operator+(boost::multi_array<T, ND> const &lhs, T const &rhs)
{
    return rhs + lhs;
}

// Subtraction operator
//! Minus operator between two multi arrays, element-wise
template <class T, long unsigned int ND>
boost::multi_array<T, ND> operator-(boost::multi_array<T, ND> const &lhs, boost::multi_array<T, ND> const &rhs)
{
    std::function<T(T, T)> func = std::minus<T>();
    return np::element_wise_duo_apply(lhs, rhs, func);
}

//! Minus operator between a scalar and a multi array, element-wise
template <class T, long unsigned int ND>
inline boost::multi_array<T, ND> operator-(T const &lhs, boost::multi_array<T, ND> const &rhs)
{
    std::function<T(T)> func = [lhs](T item)
    { return lhs - item; };
    return np::element_wise_apply(rhs, func);
}

//! Minus operator between a multi array and a scalar, element-wise
template <class T, long unsigned int ND>
inline boost::multi_array<T, ND> operator-(boost::multi_array<T, ND> const &lhs, T const &rhs)
{
    return rhs - lhs;
}

// Division operator
//! Division between two multi arrays, element wise
template <class T, long unsigned int ND>
boost::multi_array<T, ND> operator/(boost::multi_array<T, ND> const &lhs, boost::multi_array<T, ND> const &rhs)
{
    std::function<T(T, T)> func = std::divides<T>();
    return np::element_wise_duo_apply(lhs, rhs, func);
}

//! Division between a scalar and a multi array, element wise
template <class T, long unsigned int ND>
inline boost::multi_array<T, ND> operator/(T const &lhs, boost::multi_array<T, ND> const &rhs)
{
    std::function<T(T)> func = [lhs](T item)
    { return lhs / item; };
    return np::element_wise_apply(rhs, func);
}

//! Division between a multi array and a scalar, element wise
template <class T, long unsigned int ND>
inline boost::multi_array<T, ND> operator/(boost::multi_array<T, ND> const &lhs, T const &rhs)
{
    return rhs / lhs;
}

/*! @} End of Doxygen Groups*/
#endif