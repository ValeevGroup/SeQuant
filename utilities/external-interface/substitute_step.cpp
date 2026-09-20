#include "substitute_step.hpp"
#include "processing_data.hpp"
#include "processing_step_factory.hpp"

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/io/serialization/serialization.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/expr.hpp>
#include <SeQuant/core/utility/expr_matcher.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <nlohmann/json.hpp>

#include <string>
#include <string_view>

namespace sequant::util::extint {

SEQUANT_EXTINT_REGISTER_STEP_TYPE(SubstituteStep);

std::string SubstituteStep::kind() const { return "substitute"; }

bool SubstituteStep::accepts_options() const { return true; }

bool SubstituteStep::requires_options() const { return true; }

void SubstituteStep::set_options(const nlohmann::json &options) {
  if (!options.is_object()) {
    throw Exception(kind() + " expects a JSON object for its options!");
  }

  const io::serialization::DeserializationOptions parse_options{
      .def_perm_symm = Symmetry::Nonsymm,
      .def_braket_symm = BraKetSymmetry::Nonsymm,
      .def_col_symm = ColumnSymmetry::Nonsymm};

  for (const auto &[key, value] : options.items()) {
    if (key == "substitutions") {
      if (!value.is_object()) {
        throw Exception("Option '" + key + "' for " + kind() +
                        " requires object argument");
      }

      for (const auto &[target, with] : value.items()) {
        if (!with.is_string()) {
          throw Exception("All entries in the '" + key + "' option for " +
                          kind() + " need to be strings");
        }

        substitutions_.emplace_back(
            io::serialization::from_string<ExprPtr>(target, parse_options),
            io::serialization::from_string<ExprPtr>(with.get<std::string>(),
                                                    parse_options));
      }
    } else if (key == "tensor_equality_mode") {
      if (!value.is_string()) {
        throw Exception("Option '" + key + "' for " + kind() +
                        " requires string argument");
      }

      if (value == "identity") {
        tensor_mode_ = TensorComparison::Identity;
      } else if (value == "block") {
        tensor_mode_ = TensorComparison::Block;
      } else if (value == "shape") {
        tensor_mode_ = TensorComparison::Shape;
      } else {
        throw Exception("Unknown value '" + value.get<std::string>() +
                        "' for option '" + key + "' for " + kind());
      }
    } else if (key == "result_relabeling") {
      if (!value.is_object()) {
        throw Exception("Option '" + key + "' for " + kind() +
                        " requires object argument");
      }

      for (const auto &[old_label, new_label] : value.items()) {
        if (!new_label.is_string()) {
          throw Exception("All entries in the '" + key + "' option for " +
                          kind() + " need to be strings");
        }

        result_labels_.emplace(old_label, new_label.get<std::string>());
      }
    } else {
      throw Exception("Unknown option key for " + kind() + ": '" + key + "'");
    }
  }

  if (substitutions_.empty()) {
    throw Exception("Option 'substitutions' for " + kind() + " is mandatory");
  }
}

std::size_t SubstituteStep::process(std::string_view id_prefix,
                                    std::size_t id_start, ExecutionContext &ctx,
                                    const ExpressionData &data) {
  ExpressionData output;
  for (const ResultExpr &expr : data.expressions) {
    output.expressions.emplace_back(expr.clone());

    for (const auto &[target, with] : substitutions_) {
      const ExprMatcher matcher(*target,
                                ExprMatcherOptions{.tensor_cmp = tensor_mode_,
                                                   .cross_comparisons = true});
      replace(output.expressions.back(), matcher, *with);

      auto it = result_labels_.find(toUtf8(output.expressions.back().label()));
      if (it != result_labels_.end()) {
        output.expressions.back().set_label(toUtf16(it->second));
      }
    }
  }

  if (output.expressions == data.expressions) {
    // Nothing changed
    return 0;
  }

  ctx.set_data(id_prefix, id_start, std::move(output));

  return 1;
}

}  // namespace sequant::util::extint
