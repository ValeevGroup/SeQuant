#include "processing_step.hpp"

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/io/serialization/serialization.hpp>

#include <nlohmann/json.hpp>

#include <SeQuant/core/utility/exception.hpp>

#include <filesystem>
#include <string>
#include <fstream>

namespace sequant::util::extint {

ProcessingStep::~ProcessingStep() = default;


std::string ReadInputStep::kind() const { return "ReadInputStep"; }

bool ReadInputStep::accepts_options() const { return true; }

bool ReadInputStep::requires_options() const { return true; }

void ReadInputStep::set_options(const nlohmann::json &options) {
	if (!options.is_object()) {
		throw Exception("ReadInputStep expects a JSON object for its options!");
	}

	for (const auto &[key, value] : options.items()) {
		if (key == "file_path") {
			if (value.is_array()) {
				for (const nlohmann::json &current : value) {
					if (!current.is_string()) {
						throw Exception("Array entry for ReadInputStep option '" + key + "' must be strings");
					}

					input_paths_.emplace_back(current.get<std::filesystem::path>());
				}
			} else if (value.is_string()) {
				input_paths_.emplace_back(value.get<std::filesystem::path>());
			} else {
				throw Exception("Invalid data type for ReadInputStep option '" + key + "'");
			}
		} else {
			throw Exception("Unknown option key for ReadInputStep: '" + key + "'");
		}
	}
}

ProcessingOutput ReadInputStep::run(const ProcessingInput &inp) {
	return run(convert_input<VoidInput>(inp));
}

ExpressionOutput ReadInputStep::run(const VoidInput &) {
	std::vector<ResultExpr> expressions;

	for (const std::filesystem::path &current : input_paths_) {
		if (!std::filesystem::exists(current)) {
			throw Exception("Input file '" + current.string() + "' does not exist");
		}

		// Read input file
		std::ifstream in(current);
		const std::string contents(std::istreambuf_iterator<char>(in), {});

		expressions.emplace_back(io::serialization::from_string<ResultExpr>(contents));
	}

	// TODO: add expression
	return ExpressionOutput{.expressions = std::move(expressions)};
}

}
