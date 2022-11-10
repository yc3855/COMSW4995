#ifndef NP_H_
#define NP_H_

#include "boost/multi_array.hpp"
#include "boost/array.hpp"
#include <cassert>
#include <iostream>

namespace np
{

    typedef double ndArrayValue;

    // Gets the index of one element in a multi_array in one axis
    template <std::size_t ND>
    inline boost::multi_array<ndArrayValue, ND>::index
    getIndex(const boost::multi_array<ndArrayValue, ND> &m, const ndArrayValue *requestedElement, const unsigned short int direction)
    {
        int offset = requestedElement - m.origin();
        return (offset / m.strides()[direction] % m.shape()[direction] + m.index_bases()[direction]);
    }

    // Gets the index of one element in a multi_array
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

    // Function to apply a function to all elements of a multi_array
    // Simple overload
    template <typename Array, typename Element, typename Functor>
    inline void for_each(const boost::type<Element> &type_dispatch,
                         Array A, Functor &xform)
    {
        for_each(type_dispatch, A.begin(), A.end(), xform);
    }

    // Function to apply a function to all elements of a multi_array
    template <typename Element, typename Functor>
    inline void for_each(const boost::type<Element> &, Element &Val, Functor &xform)
    {
        Val = xform(Val);
    }

    // Function to apply a function to all elements of a multi_array
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

    // Function to apply a function to all elements of a multi_array
    // Simple overload
    template <typename Array, typename Functor>
    inline void for_each(Array &A, Functor xform)
    {
        // Dispatch to the proper function
        for_each(boost::type<typename Array::element>(), A.begin(), A.end(), xform);
    }

    // Takes the gradient of a n-dimensional multi_array
    // Todo: Actually implement the gradient calculation
    template <long unsigned int ND, typename... Args>
    inline constexpr boost::multi_array<double, ND> gradient(boost::multi_array<double, ND> inArray, Args... args)
    {

        using arrayIndex = boost::multi_array<double, ND>::index;
        using ndArray = boost::multi_array<ndArrayValue, ND>;

        using ndIndexArray = boost::array<arrayIndex, ND>;

        constexpr std::size_t n = sizeof...(Args);
        std::tuple<Args...> store(args...);
        boost::multi_array<double, ND> my_array;

        ndArrayValue *p = inArray.data();
        ndIndexArray index;
        for (int i = 0; i < inArray.num_elements(); i++)
        {
            index = getIndexArray(inArray, p);
            for (int j = 0; j < n; j++)
            {
                std::cout << index[j] << " ";
            }
            std::cout << " value = " << inArray(index) << "  check = " << *p << std::endl;
            ++p;
        }
        return my_array;
    }

}
#endif