#include "boost/multi_array.hpp"
#include "boost/array.hpp"
#include <cassert>
#include <iostream>

typedef double ndArrayValue;

template <std::size_t ND>
boost::multi_array<ndArrayValue, ND>::index
getIndex(const boost::multi_array<ndArrayValue, ND> &m, const ndArrayValue *requestedElement, const unsigned short int direction)
{
    int offset = requestedElement - m.origin();
    return (offset / m.strides()[direction] % m.shape()[direction] + m.index_bases()[direction]);
}

template <std::size_t ND>
boost::array<typename boost::multi_array<ndArrayValue, ND>::index, ND>
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

template <typename T, typename... Args>
void func(T t, Args... args) // recursive variadic function
{
    std::cout << t << std::endl;

    func(args...);
}

template <typename Array, typename Element, typename Functor>
void for_each(const boost::type<Element> &type_dispatch,
              Array A, Functor &xform)
{
    for_each(type_dispatch, A.begin(), A.end(), xform);
}

template <typename Element, typename Functor>
void for_each(const boost::type<Element> &, Element &Val, Functor &xform)
{
    Val = xform(Val);
}

template <typename Element, typename Iterator, typename Functor>
void for_each(const boost::type<Element> &type_dispatch,
              Iterator begin, Iterator end,
              Functor &xform)
{
    while (begin != end)
    {
        for_each(type_dispatch, *begin, xform);
        ++begin;
    }
}

template <typename Array, typename Functor>
void for_each(Array &A, Functor xform)
{
    // Dispatch to the proper function
    for_each(boost::type<typename Array::element>(), A.begin(), A.end(), xform);
}

template <long unsigned int ND, typename... Args>
constexpr boost::multi_array<double, ND> gradient(boost::multi_array<double, ND> inArray, Args... args)
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
        std::cout << index[0] << " " << index[1] << " " << index[2] << " value = " << inArray(index) << "  check = " << *p << std::endl;
        ++p;
    }
    return my_array;
}

int main()
{
    // Create a 3D array that is 3 x 4 x 2
    typedef boost::multi_array<double, 3>::index index;
    boost::multi_array<double, 3> A(boost::extents[3][4][2]);

    // Assign values to the elements
    int values = 0;
    for (index i = 0; i != 3; ++i)
        for (index j = 0; j != 4; ++j)
            for (index k = 0; k != 2; ++k)
                A[i][j][k] = values++;

    // Verify values
    int verify = 0;
    for (index i = 0; i != 3; ++i)
        for (index j = 0; j != 4; ++j)
            for (index k = 0; k != 2; ++k)
                assert(A[i][j][k] == verify++);

    boost::multi_array<double, 1> dt(boost::extents[3]);
    dt[0] = 1;
    dt[1] = 2;
    dt[2] = 3;
    boost::multi_array<double, 3> my_array = gradient(A, dt, dt, dt);
    return 0;
}