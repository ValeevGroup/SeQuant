#ifndef SEQUANT_UTILS_FLOPS_COUNTER_HPP
#define SEQUANT_UTILS_FLOPS_COUNTER_HPP

#include <SeQuant/core/expr.hpp>
#include <SeQuant/domain/utils/binary_node.hpp>
#include <SeQuant/domain/utils/eval_expr.hpp>

namespace sequant::utils {

struct FlopsCounter {
 private:
  const size_t nocc;
  const size_t nvirt;

 public:
  FlopsCounter(size_t no, size_t nv);

  size_t operator()(binary_node<eval_expr> const& expr) const;

  size_t operator()(binary_node<eval_expr> const& expr, size_t lops,
                    size_t rops) const;

 private:
  [[nodiscard]] size_t flops(size_t oidx_c, size_t vidx_c) const;

  template <typename Iterable>
  static std::array<size_t, 2> ov_idx_count(Iterable&& container);
};  // flops counter

}  // namespace sequant::utils

#endif  // SEQUANT_UTILS_FLOPS_COUNTER_HPP
