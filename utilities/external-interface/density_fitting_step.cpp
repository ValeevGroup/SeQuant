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
}

std::size_t DensityFittingStep::run(
    std::string_view, ExecutionContext &ctx,
    const std::vector<std::string_view> &inputs) {
  std::size_t new_outputs = 0;

  // TODO: This loop-structure should be abstracted away into a (helper) class
  // so the actual steps only have to implement processing of a single data
  // object (or all at once, if that's what they need)
  for (std::string_view current_input : inputs) {
    for (const auto &current : ctx.get_data(current_input)) {
      const ExpressionData &data =
          convert_data<ExpressionData>(current.data.get());

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

      // TODO: Need proper output group handling. For this, we have to track the
      // relations of every output object to every (input) IDs (also IDs the
      // input objects are associated with but which weren't explicitly provided
      // in the inputs field). Then we can determine which of the input IDs
      // should be created as corresponding output IDs and associated with which
      // objects. This can probably also be implemented in a helper base class
      // so we don't have to do this in every step implementation.
      if (had_effect) {
        ctx.set_data(std::string(current.associated_ids.front()),
                     ExpressionData{.expressions = std::move(modified)});
      } else {
        // TODO: Create alias
      }
    }
  }

  return new_outputs;
}

}  // namespace sequant::util::extint
