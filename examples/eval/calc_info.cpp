//
// Created by Bimal Gaudel on 7/18/21.
//

#include "calc_info.hpp"
#include <SeQuant/domain/eqs/cceqs.hpp>

namespace sequant::eval {

EvalNode CalcInfo::node_(const ExprPtr& expr, size_t rank) const {
  auto trimmed = optimize::tail_factor(expr);
  return optm_opts.single_term ? optimize::optimize(trimmed, false)
                               : to_eval_node(trimmed);
}

container::vector<EvalNode> CalcInfo::nodes(
    const container::vector<sequant::ExprPtr>& exprs) const {
  using namespace ranges::views;
  assert(exprs.size() == eqn_opts.excit);
  return zip(exprs, iota(size_t{1}, eqn_opts.excit + 1)) |
         transform(
             [this](auto&& pair) { return node_(pair.first, pair.second); }) |
         ranges::to<container::vector<EvalNode>>;
}

container::vector<ExprPtr> CalcInfo::exprs() const {
  auto exprs = eqs::cceqvec{eqn_opts.excit, eqn_opts.excit}(true, true, true,
                                                            true, true);
  return exprs | ranges::views::tail |
         ranges::views::transform([this](ExprPtr const& xpr) {
           return eqn_opts.spintrace ? closedshell_cc_spintrace(xpr) : xpr;
         }) |
         ranges::to<container::vector<ExprPtr>>;
}
CalcInfo::CalcInfo(const OptionsEquations& equation_options,
                   const OptionsOptimization& optmization_options,
                   const OptionsSCF& scf_options,
                   const OptionsLog& logging_options,
                   const DataInfo& tensor_files)
    : eqn_opts{equation_options},
      optm_opts{optmization_options},
      scf_opts{scf_options},
      log_opts{logging_options},
      fock_eri{tensor_files} {}

CalcInfo make_calc_info(std::string_view config_file,
                        std::string_view fock_or_eri_file,
                        std::string_view eri_or_fock_file,
                        std::string_view output_file) {
  auto parser = ParseConfigFile{};
  parser.parse(config_file);
  auto const eq_opts = parser.opts_equations();
  auto const optm_opts = parser.opts_optimization();
  auto const scf_opts = parser.opts_scf();
  auto log_opts = parser.opts_log();
  if (!output_file.empty()) log_opts.file = output_file.data();
  auto const data_info = DataInfo{fock_or_eri_file, eri_or_fock_file};
  return CalcInfo{eq_opts, optm_opts, scf_opts, log_opts, data_info};
}

}  // namespace sequant::eval
