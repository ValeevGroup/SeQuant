#include "canonicalize_step.hpp"
#include "processing_data.hpp"
#include "processing_step_factory.hpp"

#include <SeQuant/core/expr.hpp>

#include <nlohmann/json.hpp>

#include <string>
#include <string_view>

namespace sequant::util::extint {

SEQUANT_EXTINT_REGISTER_STEP_TYPE(CanonicalizeStep, "canonicalize");

std::string CanonicalizeStep::kind() const { return "canonicalize"; }

bool CanonicalizeStep::accepts_options() const { return false; }

bool CanonicalizeStep::requires_options() const { return false; }

void CanonicalizeStep::set_options(const nlohmann::json &) {
  throw Exception(kind() + " doesn't take any options");
}

std::size_t CanonicalizeStep::process(std::string_view id_prefix,
                                      std::size_t id_start,
                                      ExecutionContext &ctx,
                                      const ExpressionData &data) {
  std::vector<ResultExpr> outputs;

  for (const ResultExpr &expr : data.expressions) {
    ResultExpr clone = expr.clone();
    canonicalize(clone);

    outputs.emplace_back(std::move(clone));
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

bool CanonicalizeStep::alias_unchanged_inputs() const { return true; }

}  // namespace sequant::util::extint
