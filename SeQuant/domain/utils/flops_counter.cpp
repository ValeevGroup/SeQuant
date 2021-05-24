#include "flops_counter.hpp"

namespace sequant::utils {

FlopsCounter::FlopsCounter(size_t no, size_t nv) : nocc{no}, nvirt{nv} {}

size_t FlopsCounter::flops(size_t oidx_c, size_t vidx_c) const {
  return (oidx_c == 0 ? 1 : static_cast<size_t>(std::pow(nocc, oidx_c))) *
         (vidx_c == 0 ? 1 : static_cast<size_t>(std::pow(nvirt, vidx_c)));
}

size_t FlopsCounter::operator()(binary_node<eval_expr> const& expr) const {
  // leaf node
  return 0;
}

size_t FlopsCounter::operator()(binary_node<eval_expr> const& expr,
                                 size_t lops, size_t rops) const {
  const auto& rtnsr = expr.right()->tensor();

  if (expr->op() == eval_expr::eval_op::Sum) {
    auto [o, v] = FlopsCounter::ov_idx_count(rtnsr.const_braket());

    return lops + rops + FlopsCounter::flops(o, v);
  }

  // product of two tensors confirmed
  const auto& ltnsr = expr.left()->tensor();

  auto [o_pre, v_pre] = FlopsCounter::ov_idx_count(
      ranges::views::concat(ltnsr.const_braket(), rtnsr.const_braket()));

  const auto& tn = expr->tensor();
  auto [o_post, v_post] = FlopsCounter::ov_idx_count(tn.const_braket());

  auto o_net = o_pre - (o_pre - o_post) / 2;
  auto v_net = v_pre - (v_pre - v_post) / 2;

  return lops + rops + FlopsCounter::flops(o_net, v_net);
}
template <typename Iterable>
std::array<size_t, 2> FlopsCounter::ov_idx_count(Iterable&& container) {
    size_t oc = ranges::count_if(container, [](const auto& idx) {
      return idx.space() == IndexSpace::active_occupied;
    });
  return {oc, ranges::distance(container) - oc};
}

}  // namespace sequant::utils
