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
  using binary_node = utils::binary_expr<utils::eval_expr>::node_ptr;

  container::set<size_t> const& imed_hashes;

  utils::flops_counter const counter;

 public:
  bin_eval_expr_flops_counter(size_t no, size_t nv,
                              const container::set<size_t>& imeds);

  size_t operator()(const binary_node& node) const;

  size_t operator()(const binary_node& node, size_t lflops,
                    size_t rflops) const;
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

  container::vector<utils::binary_expr<utils::eval_expr>::node_ptr>
      optimal_seqs;

  sto_result(sto_result&&) = default;
  sto_result& operator=(sto_result&&) = default;

  sto_result(const sto_result&) = delete;
  sto_result& operator=(const sto_result&) = delete;
};

/**
 * Holds the result of the most expensive term scan.
 */
struct met_result {
  size_t ops;

  container::map<ExprPtr, sto_result> mets;

  met_result(met_result&&) = default;
  met_result& operator=(met_result&&) = default;

  met_result(const met_result&) = delete;
  met_result& operator=(const met_result&) = delete;
};

/**
 * @tparam Cont type of @c container.
 *
 * @param container Iterable of eval_expr objects.
 */
template <typename Cont>
sto_result single_term_opt(Cont const& container, size_t nocc, size_t nvirt,
                           container::set<size_t> const& imeds_hash) {
  using seq_t = utils::eval_sequence<size_t>;

  const auto counter = bin_eval_expr_flops_counter{nocc, nvirt, imeds_hash};

  sto_result result{std::numeric_limits<size_t>::max(), {}};

  auto finder = [&result, &container, &counter](const auto& seq) {
    auto tseq = utils::transform_eval_sequence<size_t, utils::eval_expr>(
        seq, [&container](auto x) {
          return utils::eval_expr{*(ranges::begin(container) + x)};
        });

    auto node =
        utils::binarize_eval_sequence<utils::eval_expr, utils::eval_expr>(
            tseq, utils::binarize_eval_expr);

    auto flops =
        utils::evaluate_binary_expr<utils::eval_expr, size_t>(node, counter);

    if (flops == result.ops) {
      result.optimal_seqs.emplace_back(std::move(node));
    } else if (flops < result.ops) {
      result.optimal_seqs.clear();
      result.optimal_seqs.emplace_back(std::move(node));
      result.ops = flops;
    } else {
      // flops > optimal flops. do nothing.
    }
  };

  auto init_seq = ranges::views::iota(size_t{0}) |
                  ranges::views::take(ranges::distance(container)) |
                  ranges::views::transform([](auto x) { return seq_t{x}; }) |
                  ranges::to<std::vector<seq_t>>;

  utils::enumerate_eval_sequence<size_t>(init_seq, finder);

  return result;
}

/**
 * @tparam Cont type of @c container.
 *
 * @param container Iterable of ExprPtr to flat Product.
 */
template <typename Iter>
met_result most_expensive(Iter const& iterable, size_t nocc, size_t nvirt,
                          container::set<size_t> const& imed_hashes) {
  using ranges::views::transform;

  // prod is a flat Product
  auto evxpr_prod_range = [](Product const& prod) {
    return prod | transform([](ExprPtr const& x) {
             return utils::eval_expr{x->as<Tensor>()};
           });
  };

  auto evxprs = iterable | transform([evxpr_prod_range](const auto& x) {
                  return evxpr_prod_range(x->template as<Product>());
                });

  auto costs = evxprs | transform([nocc, nvirt, &imed_hashes](const auto& x) {
                 return single_term_opt(x, nocc, nvirt, imed_hashes);
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

  return std::move(expensive);
}

// elements in container are ExprPtr to Product
template <typename Cont>
auto multi_term_opt_hartono(Cont const& container, size_t nocc, size_t nvirt) {
  auto imed_hashes = container::set<size_t>{};  // initial intermediate hashes

  // find the most expensive
  auto expensive = most_expensive(container, nocc, nvirt, imed_hashes);
  //
}

}  // namespace sequant::factorize

#endif  // SEQUANT_DOMAIN_FACTORIZE_OPTIMIZE_HPP
