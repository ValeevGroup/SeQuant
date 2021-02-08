#include "flops_counter.hpp"

namespace sequant::utils {

flops_counter::flops_counter(size_t no, size_t nv) : nocc{no}, nvirt{nv} {}

size_t flops_counter::flops(size_t oidx_c, size_t vidx_c) const {
  return (oidx_c == 0 ? 1 : std::pow(nocc, oidx_c)) *
         (vidx_c == 0 ? 1 : std::pow(nvirt, vidx_c));
}

size_t flops_counter::operator()(binary_node<eval_expr> const& expr) const {
  // leaf node
  return 0;
}

size_t flops_counter::operator()(binary_node<eval_expr> const& expr,
                                 size_t lops, size_t rops) const {
  const auto& rtnsr = expr.right()->tensor();

  if (expr->op() == eval_expr::eval_op::Sum) {
    auto [o, v] = flops_counter::ov_idx_count(rtnsr.const_braket());

    return lops + rops + flops_counter::flops(o, v);
  }

  // product of two tensors confirmed
  const auto& ltnsr = expr.left()->tensor();

  auto [o_pre, v_pre] = flops_counter::ov_idx_count(
      ranges::views::concat(ltnsr.const_braket(), rtnsr.const_braket()));

  const auto& tn = expr->tensor();
  auto [o_post, v_post] = flops_counter::ov_idx_count(tn.const_braket());

  auto o_net = o_pre - (o_pre - o_post) / 2;
  auto v_net = v_pre - (v_pre - v_post) / 2;

  return lops + rops + flops_counter::flops(o_net, v_net);
}

}  // namespace sequant::utils
