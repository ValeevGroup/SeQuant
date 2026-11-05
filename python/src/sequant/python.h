#ifndef SEQUANT2_PYTHON_H
#define SEQUANT2_PYTHON_H

#include <pybind11/complex.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <sequant-python/boost/container/small_vector.h>

namespace py = pybind11;

#include <SeQuant/core/options.hpp>

#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace sequant::python {

/// @return the named indices that a Python call specifies: either the
/// explicit external indices, given by their labels, or deduction
/// @throw std::invalid_argument unless exactly one of them is given
inline NamedIndexSpec named_index_spec(
    const std::optional<std::vector<std::wstring>>& external_indices,
    bool deduce) {
  if (deduce == external_indices.has_value())
    throw std::invalid_argument(
        "pass either external_indices (a list of index labels) or "
        "deduce=True");
  if (deduce) return deduce_named_indices;
  container::set<Index> result;
  for (const auto& label : *external_indices) result.emplace(label);
  return result;
}

}  // namespace sequant::python

#endif  // SEQUANT2_PYTHON_H
