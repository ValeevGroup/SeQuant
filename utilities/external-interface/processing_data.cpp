#include "processing_data.hpp"

#include <SeQuant/core/expr.hpp>

#include <algorithm>
#include <iterator>

namespace sequant::util::extint {

ExpressionData clone(const ExpressionData &data) {
  ExpressionData copy;
  copy.expressions.reserve(data.expressions.size());
  std::ranges::transform(data.expressions, std::back_inserter(copy.expressions),
                         &ResultExpr::clone);
  return copy;
}

ExportTreeData clone(const ExportTreeData &data) { return data; }

}  // namespace sequant::util::extint
