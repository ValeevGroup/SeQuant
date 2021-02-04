#include "flops_counter.hpp"

namespace sequant::utils {

flops_counter::flops_counter(size_t no, size_t nv) : nocc{no}, nvirt{nv} {}

size_t flops_counter::flops(size_t oidx_c, size_t vidx_c) const {
  return (oidx_c == 0 ? 1 : std::pow(nocc, oidx_c)) *
         (vidx_c == 0 ? 1 : std::pow(nvirt, vidx_c));
}

size_t flops_counter::operator()(const binary_node& expr) const {
  // leaf node
  return 0;
}

size_t flops_counter::operator()(const binary_node& node, size_t lops,
                                 size_t rops) const {
  const auto& tr = node->right()->data().tensor();

  if (node->data().op() == eval_expr::eval_op::Sum) {
    auto [o, v] = flops_counter::ov_idx_count(tr.const_braket());
    auto base = flops_counter::flops(o, v);
    auto total = lops + rops + base;
    if (node->left()->data().scalar() != Constant{1}) total += base;
    if (node->right()->data().scalar() != Constant{1}) total += base;

    return total;
  }

  // product of two tensors confirmed
  const auto& tl = node->left()->data().tensor();

  auto [o_pre, v_pre] = flops_counter::ov_idx_count(
      ranges::views::concat(tl.const_braket(), tr.const_braket()));

  const auto& tn = node->data().tensor();
  auto [o_post, v_post] = flops_counter::ov_idx_count(tn.const_braket());

  auto o_net = o_pre - (o_pre - o_post) / 2;
  auto v_net = v_pre - (v_pre - v_post) / 2;

  return lops + rops + flops_counter::flops(o_net, v_net);
}

}  // namespace sequant::utils
