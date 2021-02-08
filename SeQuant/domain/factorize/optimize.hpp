#ifndef SEQUANT_DOMAIN_FACTORIZE_OPTIMIZE_HPP
#define SEQUANT_DOMAIN_FACTORIZE_OPTIMIZE_HPP

#include <SeQuant/domain/utils/binarize_expr.hpp>
#include <SeQuant/domain/utils/eval_expr.hpp>
#include <SeQuant/domain/utils/flops_counter.hpp>

namespace sequant::factorize {

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
  size_t ops;

  container::map<ExprPtr, sto_result> mets;
};

/**
 * @tparam Cont type of @c container.
 *
 * @param container Iterable of eval_expr objects.
 */
template <typename Cont>
sto_result single_term_opt(Cont const& container, size_t nocc, size_t nvirt,
                           container::set<size_t> const& imeds_hash) {
  using seq_t = utils::eval_seq<size_t>;

  const auto counter = bin_eval_expr_flops_counter{nocc, nvirt, imeds_hash};

  sto_result result{std::numeric_limits<size_t>::max(), {}};

  auto finder = [&result, &container, &counter](const auto& seq) {
    auto tseq = seq.transform([&container](auto x) {
      return utils::eval_expr{*(ranges::begin(container) + x)};
    });

    auto node = tseq.binarize(utils::binarize_eval_expr);

    auto flops = node.evaluate(counter);

    if (flops == result.ops) {
      result.optimal_seqs.emplace_back(std::move(node));
    } else if (flops < result.ops) {
      result.optimal_seqs.clear();
      result.optimal_seqs.emplace_back(std::move(node));
      result.ops = flops;
    } else {
      // flops > optimal flops. do nothing.
    }
  };  // finder

  auto init_seq = ranges::views::iota(size_t{0}) |
                  ranges::views::take(ranges::distance(container)) |
                  ranges::views::transform([](auto x) { return seq_t{x}; }) |
                  ranges::to_vector;

  utils::enumerate_eval_seq<size_t>(init_seq, finder);

  return result;
}

sto_result single_term_opt(Product const& flat_prod, size_t nocc, size_t nvirt,
                           container::set<size_t> const& imeds_hash);

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
                               imed_hashes);
      });

  met_result expensive{0, {}};
  auto zipped = ranges::views::zip(iterable, costs);
  for (auto&& [x, y] : zipped) {
    if (y.ops == expensive.ops) {
      expensive.mets.emplace(x, std::move(y));
    } else if (y.ops > expensive.ops) {
      expensive.ops = y.ops;
      expensive.mets.clear();
      expensive.mets.emplace(x, std::move(y));
    } else {
      // do nothing
    }
  }

  return expensive;
}

///
/// elements in container are ExprPtr to Product
/// container is non-empty
///
template <typename Cont>
auto multi_term_opt_hartono(Cont const& container, size_t nocc, size_t nvirt) {
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
    // std::wcout << "hashes:\n";
    // for (auto x : imed_hashes) std::wcout << "#" << x << "$\n";

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

}  // namespace sequant::factorize

#endif  // SEQUANT_DOMAIN_FACTORIZE_OPTIMIZE_HPP
