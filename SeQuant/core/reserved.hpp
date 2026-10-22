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
inline const std::string& antisymm_label() {
  static const std::string label = "Â";
  return label;
}

/// @brief returns the reserved label for the symmetrization operator
inline const std::string& symm_label() {
  static const std::string label = "Ŝ";
  return label;
}

/// @brief returns the reserved label for the transposition operator
inline const std::string& transposition_label() {
  static const std::string label = "P̂";
  return label;
}

/// @brief overlap/metric tensor label is reserved since it is used by low-level
/// SeQuant machinery. Users can create overlap Tensor using make_overlap()
inline const std::string& overlap_label() {
  static const std::string label = "s";
  return label;
}

/// @brief kronecker tensor label is reserved since it is used by low-level
/// SeQuant machinery. Users can create Kronecker Tensor using make_kronecker()
inline const std::string& kronecker_label() {
  static const std::string label = "δ";
  return label;
}

/// @brief the (spin-orbital) 1- and multi-body reduced density matrix label
/// is reserved since the extended Wick theorem and mbpt produce it with fixed
/// symmetries, which a Tensor with this label must have. Users can create one
/// using the factories in SeQuant/core/density.hpp
inline const std::string& rdm_label() {
  static const std::string label = "γ";
  return label;
}

/// @brief the 1-hole reduced density matrix label; reserved like rdm_label()
inline const std::string& hole_rdm_label() {
  static const std::string label = "η";
  return label;
}

/// @brief the density cumulant label; reserved like rdm_label()
inline const std::string& cumulant_label() {
  static const std::string label = "κ";
  return label;
}

/// @brief the spin-free reduced density matrix label; reserved like
/// rdm_label()
inline const std::string& spinfree_rdm_label() {
  static const std::string label = "Γ";
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
inline bool is_nonreserved(const std::string& label) {
  return !ranges::contains(labels(), label);
}

}  // namespace reserved

}  // namespace sequant

#endif  // SEQUANT_CORE_RESERVED_HPP
