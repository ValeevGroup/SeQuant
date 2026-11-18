#ifndef SEQUANT_EXTERNAL_INTERFACE_FILTERSTEP_HPP
#define SEQUANT_EXTERNAL_INTERFACE_FILTERSTEP_HPP

#include "execution_context.hpp"
#include "processing_step.hpp"

#include <SeQuant/core/expr_fwd.hpp>
#include <SeQuant/core/utility/expr_filter.hpp>
#include <SeQuant/core/utility/expr_matcher.hpp>

#include <nlohmann/json_fwd.hpp>

#include <functional>
#include <map>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

namespace sequant::util::extint {

/// Parses a "contains" rule matching @p spec's "expr" (a serialized
/// expression). @p default_cmp is used if @p spec has no
/// "tensor_equality_mode".
std::unique_ptr<ExpressionFilter::Rule> parse_contains_expr(
    const nlohmann::json &spec,
    TensorComparison default_cmp = TensorComparison::Identity);

/// Parses a "contains" rule matching @p spec's "label" (one or several
/// full-match regular expressions).
std::unique_ptr<ExpressionFilter::Rule> parse_contains_label(
    const nlohmann::json &spec);

/// Parses @p value, a rule or a non-empty array of rules selecting tensors,
/// into a filter that matches if any rule does. A rule is either a label
/// pattern (string) or an object holding "expr" (compared in block mode unless
/// it gives "tensor_equality_mode") or "label". @p name names the option in
/// errors.
std::shared_ptr<ExpressionFilter> parse_tensor_filter(
    const nlohmann::json &value, const std::string &name);

class FilterStep : public ProcessingStep {
 public:
  std::string kind() const override;

  bool accepts_options() const override;
  bool requires_options() const override;
  void set_options(const nlohmann::json &options) override;

  std::size_t run(std::string_view step_id, ExecutionContext &ctx,
                  const std::vector<std::string_view> &inputs = {}) override;

 private:
  std::map<std::string, ExpressionFilter, std::less<>> groups_;
  bool keep_empty_ = true;
};

}  // namespace sequant::util::extint

#endif  // SEQUANT_EXTERNAL_INTERFACE_FILTERSTEP_HPP
