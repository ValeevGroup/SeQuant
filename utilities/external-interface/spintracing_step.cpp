#include "spintracing_step.hpp"
#include "processing_data.hpp"
#include "processing_step_factory.hpp"

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>

#include <nlohmann/json.hpp>

#include <string>
#include <string_view>

namespace sequant::util::extint {

SEQUANT_EXTINT_REGISTER_STEP_TYPE(SpintracingStep, "spintracing");

std::string SpintracingStep::kind() const { return "spintracing"; }

bool SpintracingStep::accepts_options() const { return true; }

bool SpintracingStep::requires_options() const { return false; }

void SpintracingStep::set_options(const nlohmann::json &options) {
  if (!options.is_object()) {
    throw Exception(kind() + " expects a JSON object for its options!");
  }

  for (const auto &[key, value] : options.items()) {
    if (key == "algorithm") {
      if (!value.is_string()) {
        throw Exception("Option '" + key + "' for " + kind() +
                        " requires string argument");
      }

      if (value == "rigorous") {
        use_closed_shell_algo_ = false;
      } else if (value == "closed_shell") {
        use_closed_shell_algo_ = true;
      } else {
        throw Exception("Unknown value '" + value.get<std::string>() +
                        "' for option " + key + " of " + kind());
      }
    } else {
      throw Exception("Unknown option key for " + kind() + ": '" + key + "'");
    }
  }
}

std::size_t SpintracingStep::process(std::string_view id_prefix,
                                     std::size_t id_start,
                                     ExecutionContext &ctx,
                                     const ExpressionData &data) {
  std::size_t num_outputs = 0;

  for (const ResultExpr &expr : data.expressions) {
    container::svector<ResultExpr> result;
    if (use_closed_shell_algo_) {
      result = mbpt::closed_shell_spintrace(expr);
    } else {
      result = mbpt::spintrace(expr);
    }

    ExpressionData output;
    output.expressions.insert(output.expressions.end(),
                              std::make_move_iterator(result.begin()),
                              std::make_move_iterator(result.end()));

    ctx.set_data(id_prefix, id_start + num_outputs, std::move(output));

    num_outputs += 1;
  }

  return num_outputs;
}

}  // namespace sequant::util::extint
