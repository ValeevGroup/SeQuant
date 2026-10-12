//
// Created by Ajay Melekamburath on 12/14/25.
//

#ifndef SEQUANT_DOMAIN_MBPT_OP_REGISTRY_HPP
#define SEQUANT_DOMAIN_MBPT_OP_REGISTRY_HPP

#include <SeQuant/core/basis.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expressions/abstract_tensor.hpp>
#include <SeQuant/core/expressions/tensor.hpp>
#include <SeQuant/core/reserved.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <range/v3/view/map.hpp>

#include <memory>
#include <optional>
#include <string>
#include <string_view>

namespace sequant::mbpt {

namespace detail {
inline constexpr std::wstring_view pert_superscripts = L"⁰¹²³⁴⁵⁶⁷⁸⁹";

/// @brief decorates a base label with perturbation order as superscript
/// @param base_label the base label to decorate
/// @param pert_order the perturbation order to decorate with
/// @return the decorated label
inline std::wstring decorate_with_pert_order(std::wstring_view base_label,
                                             int pert_order = 0) {
  if (pert_order == 0) return std::wstring(base_label);
  SEQUANT_ASSERT(
      pert_order >= 0 && pert_order <= 9,
      "decorate_with_pert_order: perturbation order out of range [0,9]");

  std::wstring result(base_label);
  result += detail::pert_superscripts[pert_order];
  return result;
}

/// @return @p label without its trailing perturbation-order superscript, if
/// any (the inverse of decorate_with_pert_order)
inline std::wstring_view strip_pert_order(std::wstring_view label) {
  if (!label.empty() &&
      pert_superscripts.find(label.back()) != std::wstring_view::npos)
    label.remove_suffix(1);
  return label;
}
}  // namespace detail

/// Operator character relative to Fermi vacuum
enum class OpClass { Ex, Deex, Gen };

/// @return the default Hermiticity for an operator of the given OpClass:
/// general operators are matrix elements of (anti-)Hermitian operators and
/// default to Hermitian; (de)excitation operators (cluster amplitudes, etc.)
/// are not Hermitian. This default can be overridden per operator in the
/// OpRegistry (e.g. to keep reference equations that assumed the legacy
/// conjugate-symmetric amplitudes).
inline Hermiticity default_hermiticity(OpClass cls) {
  return cls == OpClass::Gen ? Hermiticity::Hermitian
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
            std::make_shared<container::map<std::wstring, Hermiticity>>()),
        basis_grants_(std::make_shared<BasisGrants>()) {}

  /// constructs an OpRegistry from an existing map of operators and their
  /// classes
  OpRegistry(std::shared_ptr<container::map<std::wstring, OpClass>> ops)
      : ops_(std::move(ops)),
        herm_overrides_(
            std::make_shared<container::map<std::wstring, Hermiticity>>()),
        basis_grants_(std::make_shared<BasisGrants>()) {}

  /// copy constructor
  OpRegistry(const OpRegistry& other)
      : ops_(other.ops_),
        herm_overrides_(other.herm_overrides_),
        basis_grants_(other.basis_grants_) {}

  /// move constructor
  OpRegistry(OpRegistry&& other) noexcept
      : ops_(std::move(other.ops_)),
        herm_overrides_(std::move(other.herm_overrides_)),
        basis_grants_(std::move(other.basis_grants_)) {}

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

  /// @brief Grants the basis instance @p instance to the legs of @p op that
  /// OpMaker mints in @p leg_space
  /// @param op the operator label (must be a registered Ex or Deex operator)
  /// @throw Exception if @p op is not registered or is a general operator
  OpRegistry& grant_basis(const std::wstring& op, const IndexSpace& leg_space,
                          IndexBasis::instance_type instance);

  /// @return the registered label that @p op stands for: @p op itself if
  /// registered, else its base label (perturbation-order superscript
  /// stripped) if that is, e.g. `t¹` stands for `t` unless `t¹` is registered;
  /// null if neither is. Never throws
  [[nodiscard]] std::optional<std::wstring> resolve(std::wstring_view op) const;

  /// @return true if the label @p op stands for (see resolve) has a basis
  /// grant on any leg space; never throws
  [[nodiscard]] bool has_basis_grants(const std::wstring& op) const;

  /// @return the basis instance granted to the legs of the label @p op stands
  /// for (see resolve) in exactly @p leg_space, if any; never throws
  [[nodiscard]] IndexBasis::optional_instance basis_grant(
      const std::wstring& op, const IndexSpace& leg_space) const;

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

  /// @brief clears all registered operators (and their Hermiticity overrides
  /// and basis grants)
  void purge() {
    ops_->clear();
    herm_overrides_->clear();
    basis_grants_->clear();
  }

 private:
  std::shared_ptr<container::map<std::wstring, OpClass>> ops_;
  /// sparse per-operator Hermiticity overrides; absence means
  /// default_hermiticity(to_class(op)). Shared (shallow copy) like ops_.
  std::shared_ptr<container::map<std::wstring, Hermiticity>> herm_overrides_;
  /// sparse per-operator, per-leg-space basis grants. Shared like ops_.
  using BasisGrants =
      container::map<std::wstring,
                     container::map<IndexSpace, IndexBasis::instance_type>>;
  std::shared_ptr<BasisGrants> basis_grants_;

  /// @brief Validates that the operator label is not reserved and not already
  /// registered
  /// @param op the operator label to validate
  /// @throws std::runtime_error if the label is reserved or already exists
  void validate_op(const std::wstring& op) const;

  /// @brief Equality operator for OpRegistry
  friend bool operator==(const OpRegistry& reg1, const OpRegistry& reg2) {
    return *reg1.ops_ == *reg2.ops_ &&
           *reg1.herm_overrides_ == *reg2.herm_overrides_ &&
           *reg1.basis_grants_ == *reg2.basis_grants_;
  }
};  // class OpRegistry

/// @return true if the label of @p t, without its adjoint marker, stands for
/// (see OpRegistry::resolve) an Ex or Deex operator of @p reg
bool is_amplitude_tensor(const AbstractTensor& t, const OpRegistry& reg);

}  // namespace sequant::mbpt

#endif  // SEQUANT_DOMAIN_MBPT_OP_REGISTRY_HPP
