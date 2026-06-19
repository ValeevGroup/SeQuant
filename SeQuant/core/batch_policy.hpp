#ifndef SEQUANT_CORE_BATCH_POLICY_HPP
#define SEQUANT_CORE_BATCH_POLICY_HPP

#include <cstddef>
#include <functional>

namespace sequant {

class Index;
class Tensor;

/// One batchability policy shared by the single-term optimizer and the runtime
/// batched evaluator (make_evaluator, Task A3). All predicates default empty.
struct BatchPolicy {
  std::function<bool(Index const&)> is_batchable_index = {};
  std::function<std::size_t(Index const&)> batch_target_size = {};
  std::function<bool(Tensor const&)> is_volatile_leaf = {};
};

}  // namespace sequant

#endif  // SEQUANT_CORE_BATCH_POLICY_HPP
