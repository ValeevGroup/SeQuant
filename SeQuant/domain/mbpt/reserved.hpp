//
// Created by Ajay Melekamburath on 12/24/25.
//

#ifndef SEQUANT_DOMAIN_MBPT_OP_RESERVED_HPP
#define SEQUANT_DOMAIN_MBPT_OP_RESERVED_HPP

#endif  // SEQUANT_DOMAIN_MBPT_OP_RESERVED_HPP

namespace sequant::mbpt {

namespace reserved {
/// @brief returns the reserved label for the antisymmetrization operator
inline std::wstring antisymm_label() { return L"A"; }

/// @brief returns the reserved label for the symmetrization operator
inline std::wstring symm_label() { return L"S"; }

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
