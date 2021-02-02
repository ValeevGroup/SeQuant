#include "optimize.hpp"
#include <iostream>

namespace sequant::factorize {

bin_eval_expr_flops_counter::bin_eval_expr_flops_counter(
    size_t no, size_t nv, const container::set<size_t>& imeds)
    : counter{no, nv}, imed_hashes{imeds} {}

size_t bin_eval_expr_flops_counter::operator()(const binary_node& node) const {
  return 0;
}

size_t bin_eval_expr_flops_counter::operator()(const binary_node& node,
                                               size_t lflops,
                                               size_t rflops) const {
  if (imed_hashes.contains(node->data().hash())) return 0;

  return counter(node, lflops, rflops);
}

sto_result single_term_opt(Product const& flat_prod, size_t nocc, size_t nvirt,
                           container::set<size_t> const& imeds_hash) {
  auto result =
      single_term_opt(flat_prod | ranges::views::transform([](auto const& t) {
                        return utils::eval_expr{t->template as<Tensor>()};
                      }),
                      nocc, nvirt, imeds_hash);
  ranges::for_each(
      result.optimal_seqs, [s = flat_prod.scalar()](auto const& node) {
        node->data().scale(Constant{node->data().scalar().value() * s});
      });
  return result;
}

}  // namespace sequant::factorize
