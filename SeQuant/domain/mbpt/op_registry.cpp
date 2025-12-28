//
// Created by Ajay Melekamburath on 12/14/25.
//

#include <SeQuant/core/expressions/tensor.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/wstring.hpp>
#include <SeQuant/domain/mbpt/op_registry.hpp>

namespace sequant::mbpt {
void OpRegistry::validate_op(const std::wstring& op) const {
  if (!reserved::is_nonreserved(op)) {
    throw std::runtime_error("mbpt::OpRegistry::add: operator " +
                             sequant::to_string(op) + " uses a reserved label");
  }
  if (this->contains(op)) {
    throw std::runtime_error("mbpt::OpRegistry::add: operator " +
                             sequant::to_string(op) +
                             " already exists in registry");
  }
}

OpRegistry& OpRegistry::operator=(const OpRegistry& other) {
  ops_ = other.ops_;
  return *this;
}

OpRegistry& OpRegistry::operator=(OpRegistry&& other) noexcept {
  ops_ = std::move(other.ops_);
  return *this;
}

OpRegistry& OpRegistry::add(const std::wstring& op, OpClass action) {
  this->validate_op(op);
  ops_->emplace(op, action);
  return *this;
}

OpRegistry& OpRegistry::remove(const std::wstring& op) {
  if (!this->contains(op)) {
    throw std::runtime_error("mbpt::OpRegistry::remove: operator " +
                             sequant::to_string(op) +
                             " does not exist in registry");
  }
  ops_->erase(op);
  return *this;
}

bool OpRegistry::contains(const std::wstring& op) const {
  return ops_->contains(op);
}

OpClass OpRegistry::to_class(const std::wstring& op) const {
  auto it = ops_->find(op);
  if (it == ops_->end()) {
    throw std::runtime_error("mbpt::OpRegistry::to_class: operator " +
                             sequant::to_string(op) +
                             " does not exist in registry");
  }
  return it->second;
}

OpRegistry OpRegistry::clone() const {
  OpRegistry result(*this);
  result.ops_ = std::make_shared<container::map<std::wstring, OpClass>>(*ops_);
  return result;
}

}  // namespace sequant::mbpt
