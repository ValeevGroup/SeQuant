//
// Created by Bimal Gaudel on 7/18/21.
//

#include "calc_info.hpp"
#include <SeQuant/domain/mbpt/models/cc.hpp>

namespace sequant::eval {

const size_t IndexToSize::nocc = 10;

const size_t IndexToSize::nvirt = 10;

size_t IndexToSize::operator()(Index const& idx) const {
  if (idx.space() == IndexSpace::active_occupied)
    return nocc;
  else if (idx.space() == IndexSpace::active_unoccupied)
    return nvirt;
  else
    throw std::runtime_error("Unsupported IndexSpace type encountered");
}

container::vector<ExprPtr> CalcInfo::exprs() const {
  auto exprs = mbpt::sr::cc{eqn_opts.excit, eqn_opts.excit}.t();
  container::vector<ExprPtr> result{};
  for (auto r = 1; r < exprs.size(); ++r)
    result.emplace_back(
        eqn_opts.spintrace ? closed_shell_CC_spintrace(exprs[r]) : exprs[r]);
  return result;
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
