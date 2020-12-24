#ifndef SEQUANT_DOMAIN_FACTORIZE_OPTIMIZE_HPP
#define SEQUANT_DOMAIN_FACTORIZE_OPTIMIZE_HPP

#include <SeQuant/domain/utils/binarize_expr.hpp>
#include <SeQuant/domain/utils/eval_expr.hpp>
#include <SeQuant/domain/utils/flops_counter.hpp>

namespace sequant::factorize {

struct ops_count_result {
  size_t ops;

  container::vector<utils::binary_expr<utils::eval_expr>::node_ptr>
      optimal_seqs;

  ops_count_result(ops_count_result&&) = default;
  ops_count_result& operator=(ops_count_result&&) = default;

  ops_count_result(const ops_count_result&) = delete;
  ops_count_result& operator=(const ops_count_result&) = delete;
};

struct sto_flops_counter {
 private:
  using binary_node = utils::binary_expr<utils::eval_expr>::node_ptr;

  container::set<size_t> const& imed_hashes;

  utils::flops_counter const counter;

 public:
  sto_flops_counter(size_t no, size_t nv, const container::set<size_t>& imeds);

  size_t operator()(const binary_node& node) const;

  size_t operator()(const binary_node& node, size_t lflops,
                    size_t rflops) const;
};

/**
 * @tparam Cont type of @c container.
 *
 * @param container Iterable of eval_expr objects.
 */
template <typename Cont>
ops_count_result single_term_opt(Cont&& container, size_t nocc, size_t nvirt,
                                 container::set<size_t> const& imeds_hash);

// implementations

template <typename Cont>
ops_count_result single_term_opt(Cont&& container, size_t nocc, size_t nvirt,
                                 container::set<size_t> const& imeds_hash) {
  using seq_t = utils::eval_sequence<size_t>;

  const auto counter = sto_flops_counter{nocc, nvirt, imeds_hash};

  ops_count_result result{std::numeric_limits<size_t>::max(), {}};

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

}  // namespace sequant::factorize

#endif  // SEQUANT_DOMAIN_FACTORIZE_OPTIMIZE_HPP
