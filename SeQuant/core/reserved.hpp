//
// Created by Ajay Melekamburath on 12/24/25.
//

#ifndef SEQUANT_CORE_RESERVED_HPP
#define SEQUANT_CORE_RESERVED_HPP

#include <SeQuant/core/expressions/tensor.hpp>

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

/// @brief returns a list of all reserved operator labels
inline const auto& labels() {
  static const std::array reserved{antisymm_label(), symm_label(),
                                   sequant::kronecker_label(),
                                   sequant::overlap_label()};
  return reserved;
}

/// @brief checks if a label is not reserved
inline bool is_nonreserved(const std::wstring& label) {
  return !ranges::contains(labels(), label);
}

}  // namespace reserved

}  // namespace sequant

#endif  // SEQUANT_CORE_RESERVED_HPP
