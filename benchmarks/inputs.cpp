#include "inputs.hpp"

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/expressions/abstract_tensor.hpp>
#include <SeQuant/core/reserved.hpp>
#include <SeQuant/domain/mbpt/models/cc.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>

#include <map>
#include <utility>
#include <vector>

namespace sequant::bench {

const ExprPtr& cc_residual(std::size_t cc_rank, std::size_t projection) {
  static std::map<std::size_t, std::vector<ExprPtr>> residuals;

  auto it = residuals.find(cc_rank);
  if (it == residuals.end()) {
    it = residuals.emplace(cc_rank, mbpt::CC(cc_rank).t()).first;
  }

  return it->second.at(projection);
}

const ExprPtr& unsimplified_cc_residual(std::size_t cc_rank,
                                        std::size_t projection) {
  static std::map<std::pair<std::size_t, std::size_t>, ExprPtr> cache;

  const auto key = std::make_pair(cc_rank, projection);
  auto it = cache.find(key);
  if (it == cache.end()) {
    // Mirrors the steps closed_shell_spintrace takes before simplifying
    ExprPtr expr = cc_residual(cc_rank, projection)->clone();
    expand(expr);
    if (has_tensor(expr, reserved::antisymm_label())) {
      expr = mbpt::symmetrize_expr(expr);
    }
    expr = mbpt::expand_antisymm(expr);

    it = cache.emplace(key, std::move(expr)).first;
  }

  return it->second;
}

void cc_residual_args(benchmark::Benchmark* benchmark) {
  for (std::size_t n = 1; n <= max_cc_rank; ++n) {
    for (std::size_t p = 1; p <= n; ++p) {
      benchmark->Args({static_cast<int64_t>(n), static_cast<int64_t>(p)});
    }
  }
}

}  // namespace sequant::bench
