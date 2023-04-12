//
// Created by Bimal Gaudel on 7/18/21.
//

#ifndef SEQUANT_EVAL_CALC_INFO_HPP
#define SEQUANT_EVAL_CALC_INFO_HPP

#include "examples/eval/data_info.hpp"
#include "examples/eval/options.hpp"

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval_node.hpp>
#include <SeQuant/core/optimize.hpp>
#include <SeQuant/domain/eval/eval.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>

#include <range/v3/view.hpp>

#include <cstddef>
#include <iomanip>

namespace sequant::eval {

struct NoCacheAmplitudeTensor {
  template <typename N>
  [[nodiscard]] bool operator()(N&& node) const noexcept {
    return node->tensor().label() != L"t";
  }
};

struct IndexToSize {
  static const size_t nocc;
  static const size_t nvirt;
  size_t operator()(Index const& idx) const;
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
  container::vector<EvalNode<ExprT>> nodes(
      const container::vector<sequant::ExprPtr>& exprs) const {
    using namespace ranges::views;
    assert(exprs.size() == eqn_opts.excit);
    return zip(exprs, iota(size_t{1}, eqn_opts.excit + 1)) |
           transform([this](auto&& pair) {
             return node_<ExprT>(pair.first, pair.second);
           }) |
           ranges::to<container::vector<EvalNode<ExprT>>>;
  }

  template <typename Data_t, typename ExprT>
  CacheManager<Data_t> cache_manager_scf(
      container::vector<EvalNode<ExprT>> const& nodes) const {
    return optm_opts.reuse_imeds ? cache_manager<Data_t const>(nodes)
                                 : CacheManager<Data_t const>::empty();
  }

 private:
  template <typename ExprT>
  [[nodiscard]] EvalNode<ExprT> node_(ExprPtr const& expr, size_t rank) const {
    auto trimmed = opt::tail_factor(expr);
    return to_eval_node<ExprT>(
        optm_opts.single_term ? optimize(trimmed, IndexToSize{}) : trimmed);
  }
};

CalcInfo make_calc_info(std::string_view config_file,
                        std::string_view fock_or_eri_file,
                        std::string_view eri_or_fock_file,
                        std::string_view output_file);
}  // namespace sequant::eval

#endif  // SEQUANT_EVAL_CALC_INFO_HPP
