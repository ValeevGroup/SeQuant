//
// Created by Ajay Melekamburath on 12/14/25.
//

#ifndef SEQUANT_DOMAIN_MBPT_OP_REGISTRY_HPP
#define SEQUANT_DOMAIN_MBPT_OP_REGISTRY_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expressions/tensor.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <range/v3/view.hpp>

#include <memory>

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

/// @brief A Registry for MBPT Operators
///
/// A runtime registry that keeps track of MBPT operators by their labels and
/// properties.
class OpRegistry {
 public:
  /// default constructor, creates an empty registry
  OpRegistry()
      : ops_(std::make_shared<container::map<std::wstring, OpClass>>()) {}

  /// constructs an OpRegistry from an existing map of operators and their
  /// classes
  OpRegistry(std::shared_ptr<container::map<std::wstring, OpClass>> ops)
      : ops_(std::move(ops)) {}

  /// copy constructor
  OpRegistry(const OpRegistry& other) : ops_(other.ops_) {}

  /// move constructor
  OpRegistry(OpRegistry&& other) noexcept : ops_(std::move(other.ops_)) {}

  /// @brief Adds a new operator to the registry
  /// @param op the operator label
  /// @param action the class of the operator
  OpRegistry& add(const std::wstring& op, OpClass action);

  OpRegistry& remove(const std::wstring& op);

  /// @brief Checks if the registry contains an operator with the given label
  /// @param op the operator label
  /// @return true if the operator exists, false otherwise
  bool contains(const std::wstring& op) const;

  /// @brief Returns the class of the operator corresponding to the given label
  /// if it exists
  /// @param op the operator label
  /// @return the class of the operator
  OpClass to_class(const std::wstring& op) const;

  /// @brief returns a list of registered operators
  container::svector<std::wstring> ops() const;

  /// @brief clears all registered operators
  void purge() { ops_->clear(); }

 private:
  std::shared_ptr<container::map<std::wstring, OpClass>> ops_;

  /// @brief Validates that the operator label is not reserved and not already
  /// registered
  /// @param op the operator label to validate
  /// @throws sequant::Exception if the label is reserved or already exists
  void validate_op(const std::wstring& op) const;
};
}  // namespace sequant::mbpt

#endif  // SEQUANT_DOMAIN_MBPT_OP_REGISTRY_HPP
