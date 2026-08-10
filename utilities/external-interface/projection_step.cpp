#include "projection_step.hpp"
#include "processing_data.hpp"
#include "processing_step_factory.hpp"

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/domain/mbpt/biorthogonalization.hpp>

#include <nlohmann/json.hpp>

#include <string>
#include <string_view>

namespace sequant::util::extint {

SEQUANT_EXTINT_REGISTER_STEP_TYPE(ProjectionStep, "project");

std::string ProjectionStep::kind() const { return "project"; }

bool ProjectionStep::accepts_options() const { return true; }

bool ProjectionStep::requires_options() const { return true; }

void ProjectionStep::set_options(const nlohmann::json &options) {
  if (!options.is_object()) {
    throw Exception(kind() + " expects a JSON object for its options!");
  }

  for (const auto &[key, value] : options.items()) {
    if (key == "method") {
      if (!value.is_string()) {
        throw Exception("Value for " + kind() + " option '" + key +
                        "' must be a string");
      }
      if (value == "biorthogonal") {
        // Since we don't support any other method for now, we don't have to
        // store this option
      } else {
        throw Exception("Invalid value for " + kind() + " option '" + key +
                        "': '" + value.get<std::string>() + "'");
      }
    } else {
      throw Exception("Unknown option key for " + kind() + ": '" + key + "'");
    }
  }
}

std::size_t ProjectionStep::process(std::string_view id_prefix,
                                    std::size_t id_start, ExecutionContext &ctx,
                                    const ExpressionData &data) {
  container::svector<ResultExpr> transformed;

  transformed.insert(transformed.end(), data.expressions.begin(),
                     data.expressions.end());

  mbpt::biorthogonal_transform(transformed);

  ExpressionData data_obj;
  data_obj.expressions.insert(data_obj.expressions.end(),
                              std::make_move_iterator(transformed.begin()),
                              std::make_move_iterator(transformed.end()));
  ctx.set_data(id_prefix, id_start, std::move(data_obj));

  return 1;
}

}  // namespace sequant::util::extint
