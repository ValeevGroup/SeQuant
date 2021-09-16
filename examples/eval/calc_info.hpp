//
// Created by Bimal Gaudel on 7/18/21.
//

#ifndef SEQUANT_EVAL_CALC_INFO_HPP
#define SEQUANT_EVAL_CALC_INFO_HPP

#include "examples/eval/data_info.hpp"
#include "examples/eval/options.hpp"

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval_node.hpp>
#include <SeQuant/core/optimize/optimize.hpp>
#include <SeQuant/domain/eval/cache_manager.hpp>
#include <SeQuant/domain/eval/eval.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>
#include <range/v3/view.hpp>

#include <cstddef>
#include <iomanip>

namespace sequant::eval {

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

  [[nodiscard]] container::vector<EvalNode> nodes(
      container::vector<ExprPtr> const& exprs) const;

  template <typename Data_t>
  CacheManager<Data_t> cache_manager(
      container::vector<EvalNode> const& nodes) const {
    if (optm_opts.reuse_imeds)
      return make_cache_manager<Data_t>(nodes, optm_opts.cache_leaves);
    else if (optm_opts.cache_leaves)
      return make_cache_manager_leaves_only<Data_t>(nodes);
    else
      return CacheManager<Data_t>{{}, {}};
  }

 private:
  [[nodiscard]] EvalNode node_(ExprPtr const& expr, size_t rank) const;
};

CalcInfo make_calc_info(std::string_view config_file,
                        std::string_view fock_or_eri_file,
                        std::string_view eri_or_fock_file,
                        std::string_view output_file);
}  // namespace sequant::eval

#endif  // SEQUANT_EVAL_CALC_INFO_HPP
