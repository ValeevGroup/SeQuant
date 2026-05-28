//
// Created by Ajay Melekamburath on 12/14/25.
//

#ifndef SEQUANT_DOMAIN_MBPT_OP_REGISTRY_HPP
#define SEQUANT_DOMAIN_MBPT_OP_REGISTRY_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expressions/tensor.hpp>
#include <SeQuant/core/reserved.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <range/v3/view/map.hpp>

#include <memory>

namespace sequant::mbpt {

/// Operator character relative to Fermi vacuum
enum class OpClass { ex, deex, gen };

/// @return the default Hermiticity for an operator of the given OpClass:
/// general operators are matrix elements of (anti-)Hermitian operators and
/// default to Hermitian; (de)excitation operators (cluster amplitudes, etc.)
/// are not Hermitian. This default can be overridden per operator in the
/// OpRegistry (e.g. to keep reference equations that assumed the legacy
/// conjugate-symmetric amplitudes).
inline Hermiticity default_hermiticity(OpClass cls) {
  return cls == OpClass::gen ? Hermiticity::Hermitian
                             : Hermiticity::NonHermitian;
}

/// @brief A Registry for MBPT Operators
///
/// A registry that keeps track of MBPT operators by their labels and
/// properties.
///
/// Copy semantics is shallow (operator map shared via `std::shared_ptr`),
/// allowing multiple mbpt::Context objects to share operator definitions.
/// Use OpRegistry::clone() for deep copies.
class OpRegistry {
 public:
  /// default constructor, creates an empty registry
  OpRegistry()
      : ops_(std::make_shared<container::map<std::wstring, OpClass>>()),
        herm_overrides_(
            std::make_shared<container::map<std::wstring, Hermiticity>>()) {}

  /// constructs an OpRegistry from an existing map of operators and their
  /// classes
  OpRegistry(std::shared_ptr<container::map<std::wstring, OpClass>> ops)
      : ops_(std::move(ops)),
        herm_overrides_(
            std::make_shared<container::map<std::wstring, Hermiticity>>()) {}

  /// copy constructor
  OpRegistry(const OpRegistry& other)
      : ops_(other.ops_), herm_overrides_(other.herm_overrides_) {}

  /// move constructor
  OpRegistry(OpRegistry&& other) noexcept
      : ops_(std::move(other.ops_)),
        herm_overrides_(std::move(other.herm_overrides_)) {}

  /// copy assignment operator
  OpRegistry& operator=(const OpRegistry& other);

  /// move assignment operator
  OpRegistry& operator=(OpRegistry&& other) noexcept;

  /// @brief const iterator to beginning of registry
  [[nodiscard]] decltype(auto) begin() const { return ops_->cbegin(); }

  /// @brief const iterator to end of registry
  [[nodiscard]] decltype(auto) end() const { return ops_->cend(); }

  /// @brief clones this OpRegistry, creates a copy of ops_
  OpRegistry clone() const;

  /// @brief Adds a new operator to the registry
  /// @param op the operator label
  /// @param action the class of the operator
  /// @note the operator's Hermiticity defaults to default_hermiticity(action)
  OpRegistry& add(const std::wstring& op, OpClass action);

  /// @brief Adds a new operator to the registry with an explicit Hermiticity
  /// @param op the operator label
  /// @param action the class of the operator
  /// @param hermiticity the operator's Hermiticity (overrides the default
  ///        implied by @p action)
  OpRegistry& add(const std::wstring& op, OpClass action,
                  Hermiticity hermiticity);

  /// @brief Overrides the Hermiticity of an already-registered operator
  /// @param op the operator label (must already be registered)
  /// @param hermiticity the operator's Hermiticity
  OpRegistry& set_hermiticity(const std::wstring& op, Hermiticity hermiticity);

  /// @brief Removes an operator from the registry
  OpRegistry& remove(const std::wstring& op);

  /// @brief Checks if the registry contains an operator with the given label
  /// @param op the operator label
  /// @return true if the operator exists, false otherwise
  bool contains(const std::wstring& op) const;

  /// @brief Returns the class of the operator corresponding to the given label
  /// if it exists
  /// @param op the operator label
  /// @return the class of the operator
  [[nodiscard]] OpClass to_class(const std::wstring& op) const;

  /// @brief Returns the Hermiticity of the operator with the given label
  /// @param op the operator label
  /// @return the operator's Hermiticity (the per-operator override if set,
  ///         else default_hermiticity(to_class(op)))
  [[nodiscard]] Hermiticity hermiticity(const std::wstring& op) const;

  /// @brief returns a view of registered operator labels
  [[nodiscard]] auto ops() const { return ranges::views::keys(*ops_); }

  /// @brief clears all registered operators (and their Hermiticity overrides)
  void purge() {
    ops_->clear();
    herm_overrides_->clear();
  }

 private:
  std::shared_ptr<container::map<std::wstring, OpClass>> ops_;
  /// sparse per-operator Hermiticity overrides; absence means
  /// default_hermiticity(to_class(op)). Shared (shallow copy) like ops_.
  std::shared_ptr<container::map<std::wstring, Hermiticity>> herm_overrides_;

  /// @brief Validates that the operator label is not reserved and not already
  /// registered
  /// @param op the operator label to validate
  /// @throws std::runtime_error if the label is reserved or already exists
  void validate_op(const std::wstring& op) const;

  /// @brief Equality operator for OpRegistry
  friend bool operator==(const OpRegistry& reg1, const OpRegistry& reg2) {
    return *reg1.ops_ == *reg2.ops_ &&
           *reg1.herm_overrides_ == *reg2.herm_overrides_;
  }
};  // class OpRegistry
}  // namespace sequant::mbpt

#endif  // SEQUANT_DOMAIN_MBPT_OP_REGISTRY_HPP
