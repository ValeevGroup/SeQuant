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
    } else {
      throw Exception("Unknown option key for " + kind() + ": '" + key + "'");
    }
  }
}

std::size_t ReadInputStep::run(std::string_view step_id, ExecutionContext &ctx,
                               const std::vector<std::string_view> &inputs) {
  if (!inputs.empty()) {
    throw Exception(kind() + " does not take any inputs");
  }

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

    ctx.set_data(std::string(step_id) + "." + std::to_string(counter++),
                 ExpressionData{.expressions = {std::move(expr)}});
  }

  ctx.add_data_alias(
      std::ranges::views::iota(std::size_t(0), counter) |
          std::ranges::views::transform([&step_id](std::size_t num) {
            return std::string(step_id) + "." + std::to_string(num);
          }),
      std::string(step_id));

  return counter;
}

}  // namespace sequant::util::extint
