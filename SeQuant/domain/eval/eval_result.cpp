#include "eval_result.hpp"

namespace sequant {

EvalResult::id_t EvalResult::next_id() noexcept {
  static std::atomic<id_t> grand_type_id = 0;
  return ++grand_type_id;
}

bool EvalResult::has_value() const noexcept { return value_.has_value(); }

}  // namespace sequant
