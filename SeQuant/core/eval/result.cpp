#include <SeQuant/core/eval/result.hpp>

namespace sequant {

Result::id_t Result::next_id() noexcept {
  static std::atomic<id_t> grand_type_id = 0;
  return ++grand_type_id;
}

bool Result::has_value() const noexcept { return value_.has_value(); }

namespace detail {

std::string annot_label(std::int64_t key) noexcept {
  static container::unordered_map<std::int64_t, std::string> labels;
  static std::size_t next_id = 0;
  if (auto it = labels.find(key); it != labels.end()) return it->second;
  auto name = "i" + std::to_string(next_id++);
  return labels.emplace(key, std::move(name)).first->second;
}

}  // namespace detail

}  // namespace sequant
