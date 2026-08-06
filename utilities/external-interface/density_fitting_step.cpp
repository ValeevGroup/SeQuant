#include "density_fitting_step.hpp"
#include "processing_data.hpp"
#include "processing_step_factory.hpp"

#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index_space_registry.hpp>
#include <SeQuant/core/utility/string.hpp>
#include <SeQuant/domain/mbpt/rules/df.hpp>

#include <nlohmann/json.hpp>

namespace sequant::util::extint {

SEQUANT_EXTINT_REGISTER_STEP_TYPE(DensityFittingStep, "density_fitting");

std::string DensityFittingStep::kind() const { return "density_fitting"; }

bool DensityFittingStep::accepts_options() const { return true; }

bool DensityFittingStep::requires_options() const { return true; }

void DensityFittingStep::set_options(const nlohmann::json &options) {
  if (!options.is_object()) {
    throw Exception("OutputStep expects a JSON object for its options!");
  }

  for (const auto &[key, value] : options.items()) {
    if (key == "auxiliary_space") {
      if (!value.is_string()) {
        throw Exception("Option '" + key + "' for " + kind() +
                        " requires string argument");
      }

      aux_space_ = get_default_context().index_space_registry()->retrieve(
          value.get<std::string_view>());
    } else if (key == "integral_label") {
      if (!value.is_string()) {
        throw Exception("Option '" + key + "' for " + kind() +
                        " requires string argument");
      }

      two_elec_int_label_ = value.get<std::string>();
    } else if (key == "df_tensor_label") {
      if (!value.is_string()) {
        throw Exception("Option '" + key + "' for " + kind() +
                        " requires string argument");
      }

      df_label_ = value.get<std::string>();
    } else {
      throw Exception("Unknown option key for " + kind() + ": '" + key + "'");
    }
  }

  if (aux_space_ == IndexSpace::null) {
    throw Exception(kind() + " requires the auxiliary_space option to be set");
  }
}

std::size_t DensityFittingStep::process(std::string_view id_prefix,
                                        std::size_t id_start,
                                        ExecutionContext &ctx,
                                        const ExpressionData &data) {
  bool had_effect = false;
  std::vector<ResultExpr> modified;
  for (const ResultExpr &expr : data.expressions) {
    ExprPtr result =
        mbpt::density_fit(expr.expression(), aux_space_,
                          toUtf16(two_elec_int_label_), toUtf16(df_label_));

    had_effect |= result != expr.expression();
    if (expr.produces_tensor()) {
      modified.emplace_back(
          ResultExpr(expr.result_as_tensor(), std::move(result)));
    } else {
      modified.emplace_back(
          ResultExpr(expr.result_as_variable(), std::move(result)));
    }
  }

  if (had_effect) {
    ctx.set_data(id_prefix, id_start,
                 ExpressionData{.expressions = std::move(modified)});
    return 1;
  }

  return 0;
}

bool DensityFittingStep::alias_unchanged_inputs() const { return true; }

}  // namespace sequant::util::extint
