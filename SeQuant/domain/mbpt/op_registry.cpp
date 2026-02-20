//
// Created by Ajay Melekamburath on 12/14/25.
//

#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/string.hpp>
#include <SeQuant/domain/mbpt/op_registry.hpp>

#include <string>

namespace sequant::mbpt {
void OpRegistry::validate_op(const std::wstring& op) const {
  if (!reserved::is_nonreserved(op)) {
    throw Exception("mbpt::OpRegistry::add: operator " + toUtf8(op) +
                    " uses a reserved label");
  }
  if (this->contains(op)) {
    throw Exception("mbpt::OpRegistry::add: operator " + toUtf8(op) +
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
    throw Exception("mbpt::OpRegistry::remove: operator " + toUtf8(op) +
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
    throw Exception("mbpt::OpRegistry::to_class: operator " + toUtf8(op) +
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
