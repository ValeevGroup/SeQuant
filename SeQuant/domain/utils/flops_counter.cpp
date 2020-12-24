#include "flops_counter.hpp"

namespace sequant::utils {

flops_counter::flops_counter(size_t no, size_t nv) : nocc{no}, nvirt{nv} {}

size_t flops_counter::flops(size_t oidx_c, size_t vidx_c) const {
  return (oidx_c == 0 ? 0 : std::pow(nocc, oidx_c)) *
         (vidx_c == 0 ? 0 : std::pow(nvirt, vidx_c));
}

size_t flops_counter::operator()(const binary_node& expr) const { return 0; }

size_t flops_counter::operator()(const binary_node& node, size_t lops,
                                 size_t rops) const {
  // right tensor
  const auto& tr = node->right()->data().seq_expr()->as<Tensor>();
  auto op = node->data().op();

  if (op == eval_expr::eval_op::Scale || op == eval_expr::eval_op::Sum) {
    auto [o, v] = flops_counter::ov_idx_count(tr.const_braket());

    if (op == eval_expr::eval_op::Scale)
      return flops_counter::flops(o, v) + rops;
    else  // sum
      return flops_counter::flops(o, v) + lops + rops;
  }

  // product of two tensors confirmed
  const auto& tl = node->left()->data().seq_expr()->as<Tensor>();

  auto [o_pre, v_pre] = flops_counter::ov_idx_count(
      ranges::views::concat(tl.const_braket(), tr.const_braket()));

  const auto& tn = node->data().seq_expr()->as<Tensor>();
  auto [o_post, v_post] = flops_counter::ov_idx_count(tn.const_braket());

  auto o_net = o_pre - (o_pre - o_post) / 2;
  auto v_net = v_pre - (v_pre - v_post) / 2;

  return lops + rops + flops_counter::flops(o_net, v_net);
}

}  // namespace sequant::utils
