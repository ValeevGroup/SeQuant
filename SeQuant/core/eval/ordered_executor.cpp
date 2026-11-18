#include <SeQuant/core/eval/ordered_executor.hpp>

namespace sequant::eval::detail {

// the per-thread diagnostic slots live here, so that the library and its users
// share them (see ordered_executor.hpp)

std::size_t& ordered_last_block_skips_slot() {
  // thread_local, matching ordered_op_counts_slot(): the writer is on the
  // production path (run_ordered_schedule_pre_results, unconditionally), so a
  // process-global would be a data race between two concurrent ordered
  // evaluations, not merely a stale read.
  static thread_local std::size_t n = 0;
  return n;
}

std::size_t& ordered_last_cell_table_size_slot() {
  // thread_local -- see ordered_last_block_skips_slot for why.
  static thread_local std::size_t n = 0;
  return n;
}

OrderedRegistryResidency& ordered_last_registry_residency_slot() {
  // thread_local -- see ordered_last_block_skips_slot for why.
  static thread_local OrderedRegistryResidency r;
  return r;
}

}  // namespace sequant::eval::detail
