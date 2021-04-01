#ifndef SEQUANT_DOMAIN_OPTIMIZE_OPTIMIZE_HPP
#define SEQUANT_DOMAIN_OPTIMIZE_OPTIMIZE_HPP

#include <SeQuant/domain/utils/binarize_expr.hpp>
#include <SeQuant/domain/utils/eval_expr.hpp>
#include <SeQuant/domain/utils/flops_counter.hpp>

namespace sequant::optimize {

utils::binary_node<utils::eval_expr> optimize(ExprPtr const& expr);

/**
 * Omit the first factor from the top level product from given expression.
 * Intended to drop "A" and "S" tensors from CC amplitudes as a preparatory step
 * for evaluation of the amplitudes.
 */
ExprPtr tail_factor(ExprPtr const& expr) noexcept;

/**
 * Function object to perform flops count on binary_expr<eval_expr>.
 */
struct bin_eval_expr_flops_counter {
 private:
  container::set<size_t> const& imed_hashes;

  utils::flops_counter const counter;

 public:
  bin_eval_expr_flops_counter(size_t no, size_t nv,
                              const container::set<size_t>& imeds);

  size_t operator()(utils::binary_node<utils::eval_expr> const& node) const;

  size_t operator()(utils::binary_node<utils::eval_expr> const& node,
                    size_t lflops, size_t rflops) const;
};

/**
 * Result of the single term optimization of a term.
 * Holds operations count.
 *
 * Iterable of one or more binary_expr<eval_expr> root node pointers
 * that lead to the same operations count.
 *
 * ie. degenerate evaluations leading to the minimal operations count are
 * stored as binary tree nodes.
 */
struct sto_result {
  size_t ops;

  container::vector<utils::binary_node<utils::eval_expr>> optimal_seqs;
};

/**
 * Holds the result of the most expensive term scan.
 */
struct met_result {
  // sum of operations count of less expensive terms.
  size_t ops_lets;
  // operations count of most expensive term(s).
  size_t ops_met;

  container::map<ExprPtr, sto_result> mets;
};

sto_result single_term_opt(Product const& prod, size_t nocc, size_t nvirt,
                           container::set<size_t> const& imeds_hash,
                           bool canon = true);

/**
 * @tparam Cont type of @c container.
 *
 * @param container Iterable of ExprPtr to flat Product.
 */
template <typename Cont>
met_result most_expensive(Cont const& iterable, size_t nocc, size_t nvirt,
                          container::set<size_t> const& imed_hashes) {
  using ranges::views::transform;

  auto costs =
      iterable | transform([nocc, nvirt, &imed_hashes](auto const& xpr) {
        return single_term_opt(xpr->template as<Product>(), nocc, nvirt,
                               imed_hashes, true);
      });

  met_result expensive{0, 0, {}};
  auto expr2sto = ranges::views::zip(iterable, costs);
  for (auto&& [xpr, sto] : expr2sto) {
    if (sto.ops == expensive.ops_met) {
      // found another most-expensive term
      expensive.mets.emplace(xpr, std::move(sto));
    } else if (sto.ops > expensive.ops_met) {
      // found a more expensive term than the previous most expensive
      expensive.ops_met = sto.ops;
      expensive.mets.clear();
      expensive.mets.emplace(xpr, std::move(sto));
    } else {
      // found a less expensive
      // do nothing
    }
    expensive.ops_lets += sto.ops;
  }  // for

  expensive.ops_lets -= ranges::distance(expensive.mets) * expensive.ops_met;

  return expensive;
}

///
/// elements in container are ExprPtr to Product
/// container is non-empty
///
template <typename Cont>
container::map<ExprPtr, sto_result> multi_term_opt_hartono(
    Cont const& container, size_t nocc, size_t nvirt) {
  using ranges::none_of;
  using ranges::views::addressof;
  using ranges::views::filter;
  using ranges::views::indirect;
  using ranges::views::keys;

  auto optimized_terms = container::map<ExprPtr, sto_result>{};

  // initial intermediate hash registry is empty
  auto imed_hashes = container::set<size_t>{};

  auto terms_range = container | ranges::to<container::vector<ExprPtr>>;

  while (!terms_range.empty()) {
    // find the most expensive
    met_result expensive = most_expensive(container, nocc, nvirt, imed_hashes);

    // update intermediate hashes
    for (auto&& [xpr, sto] : expensive.mets) {
      // considers only the first of the degenerate evaluation trees
      const auto& opt_tree = *(sto.optimal_seqs.begin());

      // visit tree nodes and copy their hashes to registry

      opt_tree.visit_internal(
          [&imed_hashes](const auto& x) { imed_hashes.emplace(x.hash()); });
    }

    // collect result
    for (auto&& res : expensive.mets) optimized_terms.emplace(std::move(res));

    //
    // update terms_range
    //
    // returns true if xpr is not present in optimized_terms
    auto remaining_filter = [&optimized_terms](auto const& xpr) {
      return none_of(addressof(indirect(keys(optimized_terms))),
                     [&xpr](auto x) { return x == std::addressof(*xpr); });
    };
    //
    terms_range = container | filter(remaining_filter) |
                  ranges::to<container::vector<ExprPtr>>;
  }

  return optimized_terms;
}
}  // namespace sequant::optimize

#endif  // SEQUANT_DOMAIN_OPTIMIZE_OPTIMIZE_HPP
