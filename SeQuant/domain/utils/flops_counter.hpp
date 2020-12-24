#ifndef SEQUANT_UTILS_FLOPS_COUNTER_HPP
#define SEQUANT_UTILS_FLOPS_COUNTER_HPP

#include <SeQuant/core/expr.hpp>
#include <SeQuant/domain/utils/binary_expr.hpp>
#include <SeQuant/domain/utils/eval_expr.hpp>

namespace sequant::utils {

struct flops_counter {
 private:
  using binary_node = binary_expr<eval_expr>::node_ptr;

  const size_t nocc;
  const size_t nvirt;

 public:
  flops_counter(size_t no, size_t nv);

  size_t operator()(const binary_node& expr) const;

  size_t operator()(const binary_node& node, size_t lops, size_t rops) const;

 private:
  size_t flops(size_t oidx_c, size_t vidx_c) const;

  // pair: occupied indices count and virtual indices count
  static inline auto ov_idx_count =
      [](const auto& cont) -> std::pair<size_t, size_t> {
    // occupied indices
    auto oc = ranges::count_if(cont, [](const auto& idx) {
      return idx.space() == IndexSpace::active_occupied;
    });

    return std::pair{oc, cont.size() - oc};
  };

};  // flops counter

}  // namespace sequant::utils

#endif  // SEQUANT_UTILS_FLOPS_COUNTER_HPP
