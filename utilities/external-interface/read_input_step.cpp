#include "read_input_step.hpp"
#include "processing_step_factory.hpp"

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/io/serialization/serialization.hpp>
#include <SeQuant/core/utility/exception.hpp>

#include <nlohmann/json.hpp>

#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

namespace sequant::util::extint {

SEQUANT_EXTINT_REGISTER_STEP_TYPE(ReadInputStep, "read_input");

std::string ReadInputStep::kind() const { return "read_input"; }

bool ReadInputStep::accepts_options() const { return true; }

bool ReadInputStep::requires_options() const { return true; }

void ReadInputStep::set_options(const nlohmann::json &options) {
  if (!options.is_object()) {
    throw Exception(kind() + " expects a JSON object for its options!");
  }

  for (const auto &[key, value] : options.items()) {
    if (key == "file_path") {
      if (value.is_array()) {
        for (const nlohmann::json &current : value) {
          if (!current.is_string()) {
            throw Exception("Array entry for " + kind() + " option '" + key +
                            "' must be strings");
          }

          input_paths_.emplace_back(current.get<std::filesystem::path>());
        }
      } else if (value.is_string()) {
        input_paths_.emplace_back(value.get<std::filesystem::path>());
      } else {
        throw Exception("Invalid data type for " + kind() + " option '" + key +
                        "'");
      }
    } else if (key == "default_symmetry") {
      if (!value.is_string()) {
        throw Exception("Value for " + kind() + " option '" + key +
                        "' must be a string");
      }
      if (value == "none") {
        options_.def_perm_symm = Symmetry::Nonsymm;
        options_.def_col_symm = ColumnSymmetry::Nonsymm;
        options_.def_braket_symm = BraKetSymmetry::Nonsymm;
      } else if (value == "antisymmetric") {
        options_.def_perm_symm = Symmetry::Antisymm;
        options_.def_col_symm = ColumnSymmetry::Symm;
        options_.def_braket_symm = BraKetSymmetry::Nonsymm;
      } else if (value == "symmetric") {
        options_.def_perm_symm = Symmetry::Symm;
        options_.def_col_symm = ColumnSymmetry::Symm;
        options_.def_braket_symm = BraKetSymmetry::Nonsymm;
      } else {
        throw Exception("Invalid value for " + kind() + " option '" + key +
                        "': '" + value.get<std::string>() + "'");
      }
    } else if (key == "simplify") {
      if (!value.is_boolean()) {
        throw Exception("Value for " + kind() + " option '" + key +
                        "' must be a boolean");
      }

      simplify_ = value.get<bool>();
    } else if (key == "canonicalize") {
      if (!value.is_boolean()) {
        throw Exception("Value for " + kind() + " option '" + key +
                        "' must be a boolean");
      }

      canonicalize_ = value.get<bool>();
    } else {
      throw Exception("Unknown option key for " + kind() + ": '" + key + "'");
    }
  }
}

std::size_t ReadInputStep::process(std::string_view id_prefix,
                                   std::size_t id_start,
                                   ExecutionContext &ctx) {
  std::size_t counter = 0;

  for (const std::filesystem::path &current : input_paths_) {
    if (!std::filesystem::exists(current)) {
      throw Exception("Input file '" + current.string() + "' does not exist");
    }

    // Read input file
    std::ifstream in(current);
    const std::string contents(std::istreambuf_iterator<char>(in), {});

    ResultExpr expr =
        io::serialization::from_string<ResultExpr>(contents, options_);

    if (simplify_) {
      simplify(expr);
    } else if (canonicalize_) {
      canonicalize(expr);
    }

    ctx.set_data(id_prefix, id_start + counter++,
                 ExpressionData{.expressions = {std::move(expr)}});
  }

  return counter;
}

}  // namespace sequant::util::extint
