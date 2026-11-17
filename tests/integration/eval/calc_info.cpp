//
// Created by Bimal Gaudel on 7/18/21.
//

#include "calc_info.hpp"

#include <SeQuant/core/reserved.hpp>
#include <SeQuant/domain/mbpt/models/cc.hpp>

#include <range/v3/algorithm/find_if.hpp>
#include <range/v3/algorithm/sort.hpp>
#include <range/v3/view/concat.hpp>

namespace sequant::eval {

container::vector<ExprPtr> CalcInfo::exprs() const {
  auto exprs = mbpt::CC{eqn_opts.excit}.t();
  container::vector<ExprPtr> result{};
  for (auto r = 1; r < exprs.size(); ++r)
    result.emplace_back(eqn_opts.spintrace
                            ? mbpt::closed_shell_CC_spintrace(exprs[r])
                            : exprs[r]);
  return result;
}

container::vector<container::svector<Index>> CalcInfo::residual_layouts(
    const container::vector<ExprPtr>& exprs) {
  auto is_projector = [](ExprPtr const& x) {
    if (!x->is<Tensor>()) return false;
    auto const& label = x->as<Tensor>().label();
    return label == reserved::antisymm_label() ||
           label == reserved::symm_label();
  };

  container::vector<container::svector<Index>> result;
  for (auto const& expr : exprs) {
    // every term carries the same projector; look at the first one
    auto const& term = expr->is<Sum>() ? expr->at(0) : expr;
    auto it = ranges::find_if(*term, is_projector);
    SEQUANT_ASSERT(it != ranges::end(*term) && "residual without projector");
    auto const& projector = (*it)->as<Tensor>();
    auto kets = projector.ket() | ranges::to<container::svector<Index>>;
    auto bras = projector.bra() | ranges::to<container::svector<Index>>;
    ranges::sort(kets, Index::LabelCompare{});
    ranges::sort(bras, Index::LabelCompare{});
    result.emplace_back(ranges::views::concat(kets, bras) |
                        ranges::to<container::svector<Index>>);
  }
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
