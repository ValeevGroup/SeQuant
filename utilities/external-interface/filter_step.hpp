#ifndef SEQUANT_EXTERNAL_INTERFACE_FILTERSTEP_HPP
#define SEQUANT_EXTERNAL_INTERFACE_FILTERSTEP_HPP

#include "execution_context.hpp"
#include "processing_step.hpp"

#include <SeQuant/core/expr_fwd.hpp>
#include <SeQuant/core/utility/expr_matcher.hpp>

#include <nlohmann/json_fwd.hpp>

#include <functional>
#include <map>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

namespace sequant::util::extint {

class ExpressionFilter {
 public:
  struct Rule {
    bool negate = false;

    virtual ~Rule() = default;

    virtual bool matches(const Expr &expr) const = 0;
  };

  ExpressionFilter() = default;

  bool matches(const Expr &) const;

  void set_require_all(bool require);
  void add_rule(std::unique_ptr<Rule> rule);

 private:
  bool require_all_ = true;
  std::vector<std::unique_ptr<Rule>> rules_;
};

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
