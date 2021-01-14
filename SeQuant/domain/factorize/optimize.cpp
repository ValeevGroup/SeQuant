#include "optimize.hpp"
#include <iostream>

namespace sequant::factorize {

bin_eval_expr_flops_counter::bin_eval_expr_flops_counter(size_t no, size_t nv,
                                     const container::set<size_t>& imeds)
    : counter{no, nv}, imed_hashes{imeds} {}

size_t bin_eval_expr_flops_counter::operator()(const binary_node& node) const {
  return counter(node);
}

size_t bin_eval_expr_flops_counter::operator()(const binary_node& node, size_t lflops,
                                     size_t rflops) const {
  if (imed_hashes.contains(node->data().hash())) return 0;

  return counter(node, lflops, rflops);
}

}  // namespace sequant::factorize
