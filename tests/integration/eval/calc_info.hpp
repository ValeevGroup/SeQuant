//
// Created by Bimal Gaudel on 7/18/21.
//

#ifndef SEQUANT_EVAL_CALC_INFO_HPP
#define SEQUANT_EVAL_CALC_INFO_HPP

#include <data_info.hpp>
#include <options.hpp>

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/cache_manager.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/optimize/optimize.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>

#include <range/v3/view/join.hpp>
#include <range/v3/view/single.hpp>
#include <range/v3/view/tail.hpp>
#include <range/v3/view/transform.hpp>

#include <cstddef>
#include <format>
#include <iostream>

namespace sequant::eval {

struct NoCacheAmplitudeTensor {
  template <typename N>
  [[nodiscard]] bool operator()(N&& node) const noexcept {
    return node->tensor().label() != L"t";
  }
};

struct CalcInfo {
  OptionsEquations const eqn_opts;

  OptionsOptimization const optm_opts;

  OptionsSCF const scf_opts;

  OptionsLog const log_opts;

  DataInfo const fock_eri;

  CalcInfo(OptionsEquations const& equation_options,
           OptionsOptimization const& optmization_options,
           OptionsSCF const& scf_options, OptionsLog const& logging_options,
           DataInfo const& tensor_files);

  [[nodiscard]] container::vector<ExprPtr> exprs() const;

  template <typename ExprT>
  container::vector<container::vector<EvalNode<ExprT>>> nodes(
      const container::vector<sequant::ExprPtr>& exprs) const {
    using namespace ranges::views;
    SEQUANT_ASSERT(exprs.size() == eqn_opts.excit);
    return zip(exprs, iota(size_t{1}, eqn_opts.excit + 1)) |
           transform([this](auto&& pair) {
             return this->template node_<ExprT>(pair.first, pair.second);
           }) |
           ranges::to_vector;
  }

  template <typename ExprT>
  CacheManager<EvalNode<ExprT>> cache_manager_scf(
      container::vector<container::vector<EvalNode<ExprT>>> const& nodes)
      const {
    return optm_opts.reuse_imeds ? cache_manager(ranges::views::join(nodes))
                                 : CacheManager<EvalNode<ExprT>>::empty();
  }

 private:
  /// Omit the first factor from the top level product from given expression.
  /// Intended to drop "A" and "S" tensors from CC amplitudes as a preparatory
  /// step for evaluation of the amplitudes.
  static ExprPtr tail_factor(ExprPtr const& expr) noexcept {
    if (expr->is<Tensor>())
      return expr->clone();
    else if (expr->is<Product>()) {
      auto scalar = expr->as<Product>().scalar();
      if (scalar == 1 && expr->size() == 2) {
        return expr->at(1);
      }
      auto facs = ranges::views::tail(*expr);
      return ex<Product>(
          Product{scalar, ranges::begin(facs), ranges::end(facs)});
    } else {
      // sum
      auto summands = *expr | ranges::views::transform(
                                  [](auto const& x) { return tail_factor(x); });
      return ex<Sum>(Sum{ranges::begin(summands), ranges::end(summands)});
    }
  }

  template <typename ExprT>
  [[nodiscard]] container::vector<EvalNode<ExprT>> node_(ExprPtr const& expr,
                                                         size_t rank) const {
    auto trimmed = tail_factor(expr);

    // Optimize the whole sum at once: reorder() and multi-term factorization
    // both act ACROSS summands, so a per-term optimize() would defeat them.
    // optimize() single-term-optimizes every summand internally and returns the
    // (possibly reordered/factored) sum, whose top-level summands become the
    // eval nodes below.
    //
    // NOTE: multi-term factorization is performed inside optimize(), so it is
    // tied to single_term: it is honored only when single_term is on. When
    // single_term is off, optimize() is not called and multi_term has no effect
    // (the raw expression is binarized as-is). OptimizeOptions has no switch to
    // run multi-term/reorder without single-term optimization, and tying the
    // two keeps this example's option handling simple (see
    // OptionsOptimization).
    OptimizeOptions opts;
    opts.multiterm = optm_opts.multi_term ? MultiTermFactor::Enable
                                          : MultiTermFactor::Disable;
    ExprPtr const optimized =
        optm_opts.single_term ? optimize(trimmed, opts) : trimmed;

    // Each top-level summand of the result is one eval node; a non-Sum result
    // (single-term or atomic equation) is a single summand.
    container::vector<ExprPtr> summands;
    if (optimized->is<Sum>())
      for (auto const& s : *optimized) summands.push_back(s);
    else
      summands.push_back(optimized);

    // SCF reference path: per-summand binarize for energy/residual building;
    // the head's bra/ket layout is consumed by integration helpers that index
    // by slot ordinal and don't depend on conventional layout.
    container::vector<EvalNode<ExprT>> nodes;
    nodes.reserve(summands.size());
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    for (auto const& s : summands) nodes.push_back(binarize<ExprT>(s));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END

    // Print the optimized expressions once (one summand per line) in the
    // deserialize text format, just before returning -- before any evaluation.
    // A header reports the per-equation term counts before/after optimization
    // (multi-term factorization can merge summands), a footer closes the block.
    if (log_opts.print_exprs) {
      const std::size_t nb = expr->is<Sum>() ? expr->size() : 1;
      const std::size_t na = summands.size();
      std::wcout << std::format(
          L"===== R{} terms (# terms before/after opt {}/{}) ======\n", rank,
          nb, na);
      for (auto const& s : summands) std::wcout << serialize(s) << L'\n';
      std::wcout << L"------------\n";
    }

    return nodes;
  }
};

CalcInfo make_calc_info(std::string_view config_file,
                        std::string_view fock_or_eri_file,
                        std::string_view eri_or_fock_file,
                        std::string_view output_file);
}  // namespace sequant::eval

#endif  // SEQUANT_EVAL_CALC_INFO_HPP
