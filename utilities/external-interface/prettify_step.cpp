#include "prettify_step.hpp"
#include "processing_data.hpp"
#include "processing_step_factory.hpp"

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/expressions/prettify.hpp>
#include <SeQuant/domain/mbpt/op.hpp>

#include <nlohmann/json.hpp>

#include <iterator>
#include <string>
#include <string_view>
#include <vector>

namespace sequant::util::extint {

SEQUANT_EXTINT_REGISTER_STEP_TYPE(PrettifyStep);

std::string PrettifyStep::kind() const { return "prettify"; }

bool PrettifyStep::accepts_options() const { return true; }

bool PrettifyStep::requires_options() const { return false; }

void PrettifyStep::set_options(const nlohmann::json &options) {
  if (!options.is_object()) {
    throw Exception(kind() + " expects a JSON object for its options!");
  }

  for (const auto &[key, value] : options.items()) {
    if (key == "recanonicalize") {
      if (!value.is_boolean()) {
        throw Exception("Option '" + key + "' for " + kind() +
                        " requires boolean argument");
      }
      recanonicalize_ = value.get<bool>();
    } else {
      throw Exception("Unknown option key for " + kind() + ": '" + key + "'");
    }
  }
}

std::size_t PrettifyStep::process(std::string_view id_prefix,
                                  std::size_t id_start, ExecutionContext &ctx,
                                  const ExpressionData &data) {
  auto opts = mbpt::prettify_options();
  opts.recanonicalize = recanonicalize_;

  std::vector<ResultExpr> outputs;

  for (const ResultExpr &expr : data.expressions) {
    outputs.emplace_back(prettify(expr, opts));
  }

  if (outputs != data.expressions) {
    ExpressionData data_obj;
    data_obj.expressions.insert(data_obj.expressions.end(),
                                std::make_move_iterator(outputs.begin()),
                                std::make_move_iterator(outputs.end()));
    ctx.set_data(id_prefix, id_start, std::move(data_obj));

    return 1;
  }

  return 0;
}

bool PrettifyStep::alias_unchanged_inputs() const { return true; }

}  // namespace sequant::util::extint
