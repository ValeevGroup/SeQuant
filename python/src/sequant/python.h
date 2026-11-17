#ifndef SEQUANT2_PYTHON_H
#define SEQUANT2_PYTHON_H

#include <pybind11/complex.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <sequant-python/boost/container/small_vector.h>

namespace py = pybind11;

#include <SeQuant/core/options.hpp>

#include <stdexcept>

namespace sequant::python {

/// @return the named indices that a Python call specifies: deduction, the
/// only form the bindings support (explicit external indices are given as the
/// indices of a ResultExpr, which is not bound yet)
/// @throw std::invalid_argument unless @p deduce is true
inline NamedIndexSpec named_index_spec(bool deduce) {
  if (!deduce)
    throw std::invalid_argument(
        "pass deduce=True to deduce the external indices as those that occur "
        "once in each product; explicit external indices require a ResultExpr");
  return deduce_named_indices;
}

}  // namespace sequant::python

#endif  // SEQUANT2_PYTHON_H
