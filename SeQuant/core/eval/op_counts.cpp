#include <SeQuant/core/eval/op_counts.hpp>

namespace sequant::eval::detail {

// the per-thread diagnostic slots live here, so that the library and its users
// share them (see op_counts.hpp)

double& last_op_flops() noexcept {
  static thread_local double f = 0.0;
  return f;
}

double& last_op_exec() noexcept {
  static thread_local double f = 0.0;
  return f;
}

OrderedOpCounts& ordered_op_counts_slot() noexcept {
  static thread_local OrderedOpCounts c;
  return c;
}

}  // namespace sequant::eval::detail
