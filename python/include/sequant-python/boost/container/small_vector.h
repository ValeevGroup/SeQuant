#pragma once

#include <pybind11/stl.h>

#include <boost/container/small_vector.hpp>

PYBIND11_NAMESPACE_BEGIN(PYBIND11_NAMESPACE)
PYBIND11_NAMESPACE_BEGIN(detail)

template <typename T, std::size_t N>
struct type_caster<boost::container::small_vector<T, N> >
    : list_caster<boost::container::small_vector<T, N>, T> {};

PYBIND11_NAMESPACE_END(detail)
PYBIND11_NAMESPACE_END(PYBIND11_NAMESPACE)
