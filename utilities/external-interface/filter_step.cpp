#include "filter_step.hpp"
#include "processing_data.hpp"
#include "processing_step_factory.hpp"

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/io/serialization/serialization.hpp>
#include <SeQuant/core/utility/expr.hpp>

#include <nlohmann/json.hpp>

#include <algorithm>
#include <ranges>
#include <string>
#include <string_view>

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

      res.add_rule(std::make_unique<ContainsRule>(std::move(matcher)));
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

      groups_.emplace(key, parse_filter(value));
    } else {
      throw Exception("Unknown option key for " + kind() + ": '" + key + "'");
    }
  }
}

std::size_t FilterStep::process(std::string_view id_prefix,
                                std::size_t id_start, ExecutionContext &ctx,
                                const ExpressionData &data) {
  std::vector<ExpressionData> grouped;
  grouped.reserve(groups_.size());

  for (const ResultExpr &expr : data.expressions) {
    auto apply_filter = [&](const Expr &current) {
      std::size_t idx = 0;
      for (const auto &[name, filter] : groups_) {
        if (!filter.matches(current)) {
          ++idx;
          continue;
        }

        if (grouped.size() <= idx) {
          ResultExpr filtered = expr;
          filtered.expression() = current.clone();
          grouped.emplace_back(
              ExpressionData{.expressions = {std::move(filtered)}});
        } else {
          grouped.at(idx).expressions.back().expression() += current.clone();
        }

        ++idx;
      }
    };

    if (expr.expression().is<Sum>()) {
      for (const ExprPtr &current : expr.expression().as<Sum>().summands()) {
        apply_filter(*current);
      }
    } else {
      apply_filter(*expr.expression());
    }
  }

  std::size_t idx = 0;
  std::size_t skipped = 0;
  for (const std::string &name : std::views::keys(groups_)) {
    ExpressionData &data = grouped.at(idx);
    ++idx;

    if (!keep_empty_ && data.expressions.empty()) {
      ++skipped;
      continue;
    }

    ctx.set_data(name, 0, std::move(data));
  }

  return grouped.size() - skipped;
}

}  // namespace sequant::util::extint
