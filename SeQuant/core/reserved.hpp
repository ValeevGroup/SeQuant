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
  static const std::wstring label = L"Â";
  return label;
}

/// @brief returns the reserved label for the symmetrization operator
inline const std::wstring& symm_label() {
  static const std::wstring label = L"Ŝ";
  return label;
}

/// @brief returns the reserved label for the transposition operator
inline const std::wstring& transposition_label() {
  static const std::wstring label = L"P̂";
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

/// @brief the (spin-orbital) 1- and multi-body reduced density matrix label
/// is reserved since the extended Wick theorem and mbpt produce it with fixed
/// symmetries, which a Tensor with this label must have. Users can create one
/// using the factories in SeQuant/core/density.hpp
inline const std::wstring& rdm_label() {
  static const std::wstring label = L"γ";
  return label;
}

/// @brief the 1-hole reduced density matrix label; reserved like rdm_label()
inline const std::wstring& hole_rdm_label() {
  static const std::wstring label = L"η";
  return label;
}

/// @brief the density cumulant label; reserved like rdm_label()
inline const std::wstring& cumulant_label() {
  static const std::wstring label = L"κ";
  return label;
}

/// @brief the spin-free reduced density matrix label; reserved like
/// rdm_label()
inline const std::wstring& spinfree_rdm_label() {
  static const std::wstring label = L"Γ";
  return label;
}

/// @brief returns a list of the reserved reference-density tensor labels
inline const auto& density_labels() {
  static const std::array reserved{rdm_label(), hole_rdm_label(),
                                   cumulant_label(), spinfree_rdm_label()};
  return reserved;
}

/// @brief returns a list of all reserved operator labels
inline const auto& labels() {
  static const std::array reserved{
      antisymm_label(),  symm_label(),     transposition_label(),
      kronecker_label(), overlap_label(),  rdm_label(),
      hole_rdm_label(),  cumulant_label(), spinfree_rdm_label()};
  return reserved;
}

/// @brief checks if a label is not reserved
inline bool is_nonreserved(const std::wstring& label) {
  return !ranges::contains(labels(), label);
}

}  // namespace reserved

}  // namespace sequant

#endif  // SEQUANT_CORE_RESERVED_HPP
