#include "sort_step.hpp"
#include "processing_data.hpp"
#include "processing_step_factory.hpp"

#include <SeQuant/core/expr.hpp>

#include <nlohmann/json.hpp>

#include <iterator>
#include <string>
#include <string_view>
#include <vector>

namespace sequant::util::extint {

SEQUANT_EXTINT_REGISTER_STEP_TYPE(SortStep);

std::string SortStep::kind() const { return "sort"; }

bool SortStep::accepts_options() const { return false; }

bool SortStep::requires_options() const { return false; }

void SortStep::set_options(const nlohmann::json &) {
  throw Exception(kind() + " doesn't take any options");
}

std::size_t SortStep::process(std::string_view id_prefix, std::size_t id_start,
                              ExecutionContext &ctx,
                              const ExpressionData &data) {
  const LexicographicExprCompare cmp;

  std::vector<ResultExpr> outputs;

  for (const ResultExpr &expr : data.expressions) {
    ResultExpr clone = expr.clone();
    sort_expr(clone, cmp);

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

bool SortStep::alias_unchanged_inputs() const { return true; }

}  // namespace sequant::util::extint
