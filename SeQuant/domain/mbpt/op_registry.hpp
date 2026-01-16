//
// Created by Ajay Melekamburath on 12/14/25.
//

#ifndef SEQUANT_DOMAIN_MBPT_OP_REGISTRY_HPP
#define SEQUANT_DOMAIN_MBPT_OP_REGISTRY_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expressions/tensor.hpp>
#include <SeQuant/core/reserved.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <memory>

namespace sequant::mbpt {

/// Operator character relative to Fermi vacuum
enum class OpClass { ex, deex, gen };

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
      : ops_(std::make_shared<container::map<std::wstring, OpClass>>()) {}

  /// constructs an OpRegistry from an existing map of operators and their
  /// classes
  OpRegistry(std::shared_ptr<container::map<std::wstring, OpClass>> ops)
      : ops_(std::move(ops)) {}

  /// copy constructor
  OpRegistry(const OpRegistry& other) : ops_(other.ops_) {}

  /// move constructor
  OpRegistry(OpRegistry&& other) noexcept : ops_(std::move(other.ops_)) {}

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
  OpRegistry& add(const std::wstring& op, OpClass action);

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

  /// @brief returns a view of registered operator labels
  [[nodiscard]] auto ops() const { return ranges::views::keys(*ops_); }

  /// @brief clears all registered operators
  void purge() { ops_->clear(); }

 private:
  std::shared_ptr<container::map<std::wstring, OpClass>> ops_;

  /// @brief Validates that the operator label is not reserved and not already
  /// registered
  /// @param op the operator label to validate
  /// @throws std::runtime_error if the label is reserved or already exists
  void validate_op(const std::wstring& op) const;

  /// @brief Equality operator for OpRegistry
  friend bool operator==(const OpRegistry& reg1, const OpRegistry& reg2) {
    return *reg1.ops_ == *reg2.ops_;
  }
};  // class OpRegistry
}  // namespace sequant::mbpt

#endif  // SEQUANT_DOMAIN_MBPT_OP_REGISTRY_HPP
