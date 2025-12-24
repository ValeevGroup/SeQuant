//
// Created by Ajay Melekamburath on 12/24/25.
//

#ifndef SEQUANT_DOMAIN_MBPT_OP_RESERVED_HPP
#define SEQUANT_DOMAIN_MBPT_OP_RESERVED_HPP

#endif  // SEQUANT_DOMAIN_MBPT_OP_RESERVED_HPP

namespace sequant::mbpt {

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
inline auto labels() {
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

}  // namespace sequant::mbpt
