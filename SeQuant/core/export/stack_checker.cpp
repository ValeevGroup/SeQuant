#include <SeQuant/core/export/stack_checker.hpp>

#include <SeQuant/core/export/utils.hpp>
#include <SeQuant/core/utility/exception.hpp>

#include <algorithm>
#include <string>

namespace sequant {

namespace {

std::string describe(const std::variant<Tensor, Variable> &object) {
  return std::visit([](const auto &obj) { return detail::describe(obj); },
                    object);
}

}  // namespace

void MemoryStack::set_batch_indices(std::span<const Index> indices) {
  m_cmp.set_indices(std::vector<Index>(indices.begin(), indices.end()));
}

void MemoryStack::allocate(const Tensor &tensor) {
  m_stack.emplace_back(tensor);
}

void MemoryStack::allocate(const Variable &variable) {
  m_stack.emplace_back(variable);
}

void MemoryStack::release(const Tensor &tensor, std::string_view action) {
  release(Object(tensor), action);
}

void MemoryStack::release(const Variable &variable, std::string_view action) {
  release(Object(variable), action);
}

void MemoryStack::require_empty() const {
  if (!m_stack.empty()) {
    throw Exception("Memory stack is not empty at the end of the expression: " +
                    std::to_string(m_stack.size()) +
                    " object(s) remain allocated, the topmost being " +
                    describe(m_stack.back()));
  }
}

bool MemoryStack::equals(const Object &lhs, const Object &rhs) const {
  if (lhs.index() != rhs.index()) {
    return false;
  }

  if (std::holds_alternative<Tensor>(lhs)) {
    return m_cmp(std::get<Tensor>(lhs), std::get<Tensor>(rhs));
  }

  return std::get<Variable>(lhs) == std::get<Variable>(rhs);
}

void MemoryStack::release(const Object &object, std::string_view action) {
  if (!m_stack.empty() && equals(m_stack.back(), object)) {
    m_stack.pop_back();
    return;
  }

  const std::string what = std::string("Cannot ") + std::string(action) + " " +
                           describe(object) + ": ";

  if (m_stack.empty()) {
    throw Exception(what + "nothing is allocated");
  }

  auto it = std::find_if(m_stack.rbegin(), m_stack.rend(),
                         [&](const auto &o) { return equals(o, object); });
  if (it == m_stack.rend()) {
    throw Exception(what + "it is not allocated");
  }

  throw Exception(what + "violates the stack-like memory order, as " +
                  describe(m_stack.back()) + " has been allocated after it (" +
                  std::to_string(it - m_stack.rbegin()) +
                  " object(s) above it on the stack)");
}

}  // namespace sequant
