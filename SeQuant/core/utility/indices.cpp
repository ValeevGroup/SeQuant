#include <SeQuant/core/utility/indices.hpp>

#include <string>
#include <utility>

namespace sequant {

std::function<bool(const Index&)> unused_label_validator(
    container::set<std::wstring> taken) {
  return [taken = std::move(taken)](const Index& idx) {
    return !taken.contains(std::wstring(idx.label()));
  };
}

std::function<bool(const Index&)> unused_label_validator(const Expr& expr) {
  container::set<std::wstring> labels;
  for (const auto& idx : get_used_indices(expr)) labels.emplace(idx.label());
  return unused_label_validator(std::move(labels));
}

}  // namespace sequant
