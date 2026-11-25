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
#include <SeQuant/core/optimize/optimize.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>

#include <range/v3/view/join.hpp>
#include <range/v3/view/single.hpp>
#include <range/v3/view/tail.hpp>
#include <range/v3/view/transform.hpp>

#include <cstddef>

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
    using ranges::views::transform;
    auto trimmed = tail_factor(expr);
    container::svector<Index> uocc, occ;
    for (size_t k = 1; k <= rank; ++k) {
      uocc.emplace_back(L"a_" + std::to_wstring(k));
      occ.emplace_back(L"i_" + std::to_wstring(k));
    }
    const Tensor head(L"R", bra(std::move(uocc)), ket(std::move(occ)));
    auto tform_and_save =
        transform([st = optm_opts.single_term, &head](const auto& expr) {
          // every term of a residual has its head, in the layout of the
          // amplitudes: the unoccupied orbitals and then the occupied ones.
          // The external indices of a residual are those of the head, as the
          // equations are derived; binarize throws for a term that lacks one
          return binarize<ExprT>(ResultExpr(head, st ? optimize(expr) : expr));
        }) |
        ranges::to_vector;
    if (trimmed.size() > 0) {
      return *trimmed | tform_and_save;
    } else {  // corner case: trimmed is an atom (i.e. single tensor)
      return ranges::views::single(trimmed) | tform_and_save;
    }
  }
};

CalcInfo make_calc_info(std::string_view config_file,
                        std::string_view fock_or_eri_file,
                        std::string_view eri_or_fock_file,
                        std::string_view output_file);
}  // namespace sequant::eval

#endif  // SEQUANT_EVAL_CALC_INFO_HPP
