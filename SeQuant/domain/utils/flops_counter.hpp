#ifndef SEQUANT_UTILS_FLOPS_COUNTER_HPP
#define SEQUANT_UTILS_FLOPS_COUNTER_HPP

#include <SeQuant/core/expr.hpp>
#include <SeQuant/domain/utils/binary_node.hpp>
#include <SeQuant/domain/utils/eval_expr.hpp>

namespace sequant::utils {

struct flops_counter {
 private:
  const size_t nocc;
  const size_t nvirt;

 public:
  flops_counter(size_t no, size_t nv);

  size_t operator()(binary_node<eval_expr> const& expr) const;

  size_t operator()(binary_node<eval_expr> const& expr, size_t lops,
                    size_t rops) const;

 private:
  [[nodiscard]] size_t flops(size_t oidx_c, size_t vidx_c) const;

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
