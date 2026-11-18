#include <SeQuant/core/eval/result.hpp>

#include <mutex>
#include <string>
#include <string_view>
#include <unordered_map>

namespace sequant {

Result::id_t Result::id_for_name(std::string_view name) {
  static std::mutex mtx;
  static std::unordered_map<std::string, id_t> ids;
  std::scoped_lock lock(mtx);
  return ids.try_emplace(std::string(name), ids.size() + 1).first->second;
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
