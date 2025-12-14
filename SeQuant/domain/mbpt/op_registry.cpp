//
// Created by Ajay Melekamburath on 12/14/25.
//

#include <SeQuant/core/expressions/tensor.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/wstring.hpp>
#include <SeQuant/domain/mbpt/op_registry.hpp>

namespace sequant::mbpt {

namespace detail {
/// @brief returns a list of reserved operator labels
auto reserved_labels() {
  return std::array<const std::wstring, 4>{antisymm_label(), symm_label(),
                                           overlap_label(), kronecker_label()};
}

/// @brief checks if a label is not reserved
bool is_nonreserved_label(const std::wstring& label) {
  return !ranges::contains(reserved_labels(), label);
}
}  // namespace detail

void OpRegistry::validate_op(const std::wstring& op) const {
  // SEQUANT_ASSERT(detail::is_nonreserved_label(op),
  //                "mbpt::OpRegistry::validate_op: operator " +
  //                sequant::to_string (op) + " uses a reserved label");
  SEQUANT_ASSERT(!this->contains(op), "mbpt::OpRegistry::add: operator " +
                                          sequant::to_string(op) +
                                          " already exists in registry");
}

OpRegistry& OpRegistry::add(const std::wstring& op, OpClass action) {
  this->validate_op(op);
  ops_->emplace(op, action);
  return *this;
}

OpRegistry& OpRegistry::remove(const std::wstring& op) {
  SEQUANT_ASSERT(this->contains(op), "mbpt::OpRegistry::remove: operator " +
                                         sequant::to_string(op) +
                                         " not found in registry");
  ops_->erase(op);
  return *this;
}

bool OpRegistry::contains(const std::wstring& op) const {
  return ops_->contains(op);
}

OpClass OpRegistry::to_class(const std::wstring& op) const {
  SEQUANT_ASSERT(this->contains(op), "mbpt::OpRegistry::action: operator " +
                                         sequant::to_string(op) +
                                         " not found in registry");
  return ops_->at(op);
}

container::svector<std::wstring> OpRegistry::ops() const {
  return ranges::views::keys(*ops_) |
         ranges::to<container::svector<std::wstring>>;
}
}  // namespace sequant::mbpt
