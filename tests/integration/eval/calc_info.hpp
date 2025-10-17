//
// Created by Bimal Gaudel on 7/18/21.
//

#ifndef SEQUANT_EVAL_CALC_INFO_HPP
#define SEQUANT_EVAL_CALC_INFO_HPP

#include <data_info.hpp>
#include <options.hpp>

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval_expr.hpp>
#include <SeQuant/core/optimize.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/domain/eval/cache_manager.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>

#include <range/v3/view.hpp>

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
             return node_<ExprT>(pair.first, pair.second);
           }) |
           ranges::to_vector;
  }

  template <typename ExprT>
  CacheManager cache_manager_scf(
      container::vector<container::vector<EvalNode<ExprT>>> const& nodes)
      const {
    return optm_opts.reuse_imeds ? cache_manager(ranges::views::join(nodes))
                                 : CacheManager::empty();
  }

 private:
  template <typename ExprT>
  [[nodiscard]] container::vector<EvalNode<ExprT>> node_(ExprPtr const& expr,
                                                         size_t rank) const {
    using ranges::views::transform;
    auto trimmed = opt::tail_factor(expr);
    auto tform_and_save =
        transform([st = optm_opts.single_term](const auto& expr) {
          return binarize<ExprT>(st ? optimize(expr) : expr);
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
