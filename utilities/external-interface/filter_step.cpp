#include "filter_step.hpp"
#include "processing_data.hpp"
#include "processing_step_factory.hpp"

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/io/serialization/serialization.hpp>
#include <SeQuant/core/utility/expr.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <nlohmann/json.hpp>

#include <algorithm>
#include <string>
#include <string_view>
#include <vector>

namespace sequant::util::extint {

bool ExpressionFilter::matches(const Expr &expr) const {
  for (const std::unique_ptr<Rule> &current : rules_) {
    if (current->matches(expr) == !current->negate) {
      if (!require_all_) {
        return true;
      }
    } else if (require_all_) {
      return false;
    }
  }

  return true;
}

void ExpressionFilter::set_require_all(bool require) { require_all_ = require; }

void ExpressionFilter::add_rule(std::unique_ptr<Rule> rule) {
  rules_.emplace_back(std::move(rule));
}

struct ContainsRule : ExpressionFilter::Rule {
  ContainsRule(ExprMatcher matcher) : matcher_(std::move(matcher)) {}

  bool matches(const Expr &expr) const override {
    return matcher_ == expr || std::ranges::find_if(expr, [&](const auto &e) {
                                 return e == matcher_;
                               }) != expr.end();
  }

  ExprMatcher matcher_;
};

SEQUANT_EXTINT_REGISTER_STEP_TYPE(FilterStep);

std::string FilterStep::kind() const { return "filter"; }

bool FilterStep::accepts_options() const { return true; }

bool FilterStep::requires_options() const { return true; }

ExpressionFilter parse_filter(const nlohmann::json &filter) {
  bool require_all = true;
  if (filter.contains("mode")) {
    if (filter.at("mode") == "any") {
      require_all = false;
    } else if (filter.at("mode") != "all") {
      throw Exception("Unknown filter mode '" +
                      filter.at("mode").get<std::string>() + "'");
    }
  }

  ExpressionFilter res;
  res.set_require_all(require_all);

  if (!filter.contains("rules")) {
    throw Exception("Filter is required to have rules");
  }

  if (!filter.at("rules").is_array()) {
    throw Exception("rules is required to be an array");
  }

  for (const nlohmann::json &current : filter.at("rules")) {
    if (!current.contains("type")) {
      throw Exception("Every filter rule must have a type");
    }

    if (current.at("type") == "contains") {
      if (!current.contains("expr")) {
        throw Exception("\"contains\" filter rule requires \"expr\" attribute");
      }

      const io::serialization::DeserializationOptions options{
          .def_perm_symm = Symmetry::Nonsymm,
          .def_braket_symm = BraKetSymmetry::Nonsymm,
          .def_col_symm = ColumnSymmetry::Nonsymm};

      ExprPtr expr = io::serialization::from_string<ExprPtr>(
          current.at("expr").get<std::string>(), options);

      ExprMatcherOptions match_opts;
      if (current.contains("tensor_equality_mode")) {
        const nlohmann::json &mode = current.at("tensor_equality_mode");

        if (!mode.is_string()) {
          throw Exception(
              "\"tensor_equality_mode\" requires a string argument");
        }

        if (mode == "identity") {
          match_opts.tensor_cmp = TensorComparison::Identity;
        } else if (mode == "block") {
          match_opts.tensor_cmp = TensorComparison::Block;
        } else if (mode == "shape") {
          match_opts.tensor_cmp = TensorComparison::Shape;
        } else {
          throw Exception("Unknown tensor equality mode \"" +
                          mode.get<std::string>() + "\"");
        }
      }

      ExprMatcher matcher(std::move(*expr), match_opts);

      auto rule = std::make_unique<ContainsRule>(std::move(matcher));

      if (current.contains("negate")) {
        if (!current.at("negate").is_boolean()) {
          throw Exception("\"negate\" requires a boolean argument");
        }

        rule->negate = current.at("negate").get<bool>();
      }

      res.add_rule(std::move(rule));
    } else {
      throw Exception("Unknown filter rule type \"" +
                      current.at("type").get<std::string>() + "\"");
    }
  }

  return res;
}

void FilterStep::set_options(const nlohmann::json &options) {
  if (!options.is_object()) {
    throw Exception(kind() + " expects a JSON object for its options!");
  }

  for (const auto &[key, value] : options.items()) {
    if (key == "groups") {
      if (!value.is_object()) {
        throw Exception("Option '" + key + "' for " + kind() +
                        " requires object argument");
      }

      for (const auto &[group_name, group_filter] : value.items()) {
        groups_.emplace(group_name, parse_filter(group_filter));
      }
    } else {
      throw Exception("Unknown option key for " + kind() + ": '" + key + "'");
    }
  }

  if (groups_.empty()) {
    throw Exception("Option 'groups' for " + kind() + " is mandatory");
  }
}

/// Filters a single input item's expressions into per-group accumulators,
/// indexed the same way as @p group_names.
std::vector<ExpressionData> filter_into_groups(
    const ExpressionData &data,
    const std::map<std::string, ExpressionFilter, std::less<>> &groups) {
  std::vector<ExpressionData> grouped(groups.size());

  for (const ResultExpr &expr : data.expressions) {
    auto apply_filter = [&](const Expr &term) {
      std::size_t idx = 0;
      for (const auto &[name, filter] : groups) {
        if (filter.matches(term)) {
          ExpressionData &group = grouped.at(idx);

          if (group.expressions.empty()) {
            ResultExpr filtered = expr;
            filtered.expression() = term.clone();
            group.expressions.push_back(std::move(filtered));
          } else {
            group.expressions.back().expression() += term.clone();
          }
        }

        ++idx;
      }
    };

    if (expr.expression().is<Sum>()) {
      for (const ExprPtr &term : expr.expression().as<Sum>().summands()) {
        apply_filter(*term);
      }
    } else {
      apply_filter(*expr.expression());
    }
  }

  return grouped;
}

std::size_t FilterStep::run(std::string_view step_id, ExecutionContext &ctx,
                            const std::vector<std::string_view> &inputs) {
  std::vector<std::string> group_names;
  group_names.reserve(groups_.size());
  for (const auto &[name, filter] : groups_) {
    group_names.push_back(name);
  }

  // Ids of the outputs each group ended up with, across all inputs, for the
  // group's own name alias (e.g. "step_id.res_with_xy") added at the end.
  std::vector<std::vector<std::string>> group_members(groups_.size());
  std::size_t next_slot = 0;

  for (std::string_view current_input : inputs) {
    for (const ExecutionContext::Data<ProcessingData> &current :
         ctx.get_data(current_input)) {
      const ExpressionData &data =
          convert_data<ExpressionData>(current.data.get());

      SEQUANT_ASSERT(!current.associated_ids.empty());

      try {
        std::vector<ExpressionData> grouped = filter_into_groups(data, groups_);

        const std::size_t item_slot_start = next_slot;

        for (std::size_t idx = 0; idx < grouped.size(); ++idx) {
          ExpressionData &group = grouped.at(idx);

          if (!keep_empty_ && group.expressions.empty()) {
            continue;
          }

          ctx.set_data(step_id, next_slot, std::move(group));
          group_members.at(idx).push_back(std::string(step_id) + "." +
                                          std::to_string(next_slot));
          ++next_slot;
        }

        if (next_slot > item_slot_start) {
          // Preserve this input's own alias (e.g. "res2_p2") on whichever
          // group(s) it contributed to, mirroring how
          // OneByOneProcessingStep::run() preserves individual aliases for
          // steps with a single, unambiguous output per input.
          detail::alias_individual_ids(
              ctx, step_id, current.associated_ids,
              detail::make_output_range_id(step_id, item_slot_start,
                                           next_slot - item_slot_start));
        }
      } catch (const std::exception &e) {
        throw Exception("Error in " + kind() + " on input " +
                        std::string(current.associated_ids.front()) + ": " +
                        e.what());
      }
    }
  }

  for (std::size_t idx = 0; idx < group_members.size(); ++idx) {
    if (group_members[idx].empty()) {
      continue;
    }

    ctx.add_data_alias(group_members[idx],
                       std::string(step_id) + "." + group_names[idx]);
  }

  return next_slot;
}

}  // namespace sequant::util::extint
