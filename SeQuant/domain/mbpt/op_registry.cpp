//
// Created by Ajay Melekamburath on 12/14/25.
//

#include <SeQuant/core/expressions/expr.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/string.hpp>
#include <SeQuant/domain/mbpt/op_registry.hpp>

#include <string>

namespace sequant::mbpt {
void OpRegistry::validate_op(const std::string& op) const {
  if (!reserved::is_nonreserved(op)) {
    throw Exception("mbpt::OpRegistry::add: operator " + op +
                    " uses a reserved label");
  }
  if (this->contains(op)) {
    throw Exception("mbpt::OpRegistry::add: operator " + op +
                    " already exists in registry");
  }
}

OpRegistry& OpRegistry::operator=(const OpRegistry& other) {
  ops_ = other.ops_;
  herm_overrides_ = other.herm_overrides_;
  basis_grants_ = other.basis_grants_;
  return *this;
}

OpRegistry& OpRegistry::operator=(OpRegistry&& other) noexcept {
  ops_ = std::move(other.ops_);
  herm_overrides_ = std::move(other.herm_overrides_);
  basis_grants_ = std::move(other.basis_grants_);
  return *this;
}

OpRegistry& OpRegistry::add(const std::string& op, OpClass action) {
  this->validate_op(op);
  ops_->emplace(op, action);
  return *this;
}

OpRegistry& OpRegistry::add(const std::string& op, OpClass action,
                            Hermiticity hermiticity) {
  this->validate_op(op);
  ops_->emplace(op, action);
  herm_overrides_->insert_or_assign(op, hermiticity);
  return *this;
}

OpRegistry& OpRegistry::set_hermiticity(const std::string& op,
                                        Hermiticity hermiticity) {
  if (!this->contains(op)) {
    throw Exception("mbpt::OpRegistry::set_hermiticity: operator " + op +
                    " does not exist in registry");
  }
  herm_overrides_->insert_or_assign(op, hermiticity);
  return *this;
}

OpRegistry& OpRegistry::grant_basis(const std::string& op,
                                    const IndexSpace& leg_space,
                                    IndexBasis::instance_type instance) {
  if (!this->contains(op)) {
    throw Exception("mbpt::OpRegistry::grant_basis: operator " + op +
                    " does not exist in registry");
  }
  if (this->to_class(op) == OpClass::Gen) {
    throw Exception("mbpt::OpRegistry::grant_basis: operator " + op +
                    " is a general operator; only (de)excitation operators "
                    "carry a basis instance");
  }
  (*basis_grants_)[op].insert_or_assign(leg_space, instance);
  return *this;
}

std::optional<std::string> OpRegistry::resolve(std::string_view op) const {
  if (ops_->contains(std::string(op))) return std::string(op);
  const auto base = detail::strip_pert_order(op);
  if (base.size() != op.size() && ops_->contains(std::string(base)))
    return std::string(base);
  return std::nullopt;
}

bool OpRegistry::has_basis_grants(const std::string& op) const {
  const auto label = resolve(op);
  return label && basis_grants_->contains(*label);
}

IndexBasis::optional_instance OpRegistry::basis_grant(
    const std::string& op, const IndexSpace& leg_space) const {
  const auto label = resolve(op);
  if (!label) return std::nullopt;
  auto it = basis_grants_->find(*label);
  if (it == basis_grants_->end()) return std::nullopt;
  auto leg = it->second.find(leg_space);
  if (leg == it->second.end()) return std::nullopt;
  return leg->second;
}

OpRegistry& OpRegistry::remove(const std::string& op) {
  if (!this->contains(op)) {
    throw Exception("mbpt::OpRegistry::remove: operator " + op +
                    " does not exist in registry");
  }
  ops_->erase(op);
  herm_overrides_->erase(op);
  basis_grants_->erase(op);
  return *this;
}

bool OpRegistry::contains(const std::string& op) const {
  return ops_->contains(op);
}

OpClass OpRegistry::to_class(const std::string& op) const {
  auto it = ops_->find(op);
  if (it == ops_->end()) {
    throw Exception("mbpt::OpRegistry::to_class: operator " + op +
                    " does not exist in registry");
  }
  return it->second;
}

Hermiticity OpRegistry::hermiticity(const std::string& op) const {
  if (auto it = herm_overrides_->find(op); it != herm_overrides_->end())
    return it->second;
  return default_hermiticity(this->to_class(op));
}

OpRegistry OpRegistry::clone() const {
  OpRegistry result(*this);
  result.ops_ = std::make_shared<container::map<std::string, OpClass>>(*ops_);
  result.herm_overrides_ =
      std::make_shared<container::map<std::string, Hermiticity>>(
          *herm_overrides_);
  result.basis_grants_ = std::make_shared<BasisGrants>(*basis_grants_);
  return result;
}

bool is_amplitude_tensor(const AbstractTensor& t, const OpRegistry& reg) {
  const auto label = reg.resolve(strip_adjoint_label(t._label()));
  if (!label) return false;
  const auto cls = reg.to_class(*label);
  return cls == OpClass::Ex || cls == OpClass::Deex;
}

}  // namespace sequant::mbpt
