#include <SeQuant/core/eval/result.hpp>

namespace sequant {

Result::id_t Result::next_id() noexcept {
  static std::atomic<id_t> grand_type_id = 0;
  return ++grand_type_id;
}

bool Result::has_value() const noexcept { return value_.has_value(); }

}  // namespace sequant
