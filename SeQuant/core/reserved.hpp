//
// Created by Ajay Melekamburath on 12/24/25.
//

#ifndef SEQUANT_CORE_RESERVED_HPP
#define SEQUANT_CORE_RESERVED_HPP

#include <range/v3/algorithm/contains.hpp>

#include <array>
#include <string>

namespace sequant {

namespace reserved {
/// @brief returns the reserved label for the antisymmetrization operator
inline const std::wstring& antisymm_label() {
  static const std::wstring label = L"A";
  return label;
}

/// @brief returns the reserved label for the symmetrization operator
inline const std::wstring& symm_label() {
  static const std::wstring label = L"S";
  return label;
}

/// @brief overlap/metric tensor label is reserved since it is used by low-level
/// SeQuant machinery. Users can create overlap Tensor using make_overlap()
inline const std::wstring& overlap_label() {
  static const std::wstring label = L"s";
  return label;
}

/// @brief kronecker tensor label is reserved since it is used by low-level
/// SeQuant machinery. Users can create Kronecker Tensor using make_kronecker()
inline const std::wstring& kronecker_label() {
  static const std::wstring label = L"δ";
  return label;
}

/// @brief returns a list of all reserved operator labels
inline const auto& labels() {
  static const std::array reserved{antisymm_label(), symm_label(),
                                   kronecker_label(), overlap_label()};
  return reserved;
}

/// @brief checks if a label is not reserved
inline bool is_nonreserved(const std::wstring& label) {
  return !ranges::contains(labels(), label);
}

}  // namespace reserved

}  // namespace sequant

#endif  // SEQUANT_CORE_RESERVED_HPP
