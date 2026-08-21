#include "export_step.hpp"
#include "processing_data.hpp"
#include "processing_step_factory.hpp"

#include <SeQuant/core/export/context.hpp>
#include <SeQuant/core/export/export.hpp>
#include <SeQuant/core/export/expression_group.hpp>
#include <SeQuant/core/export/generation_optimizer.hpp>
#include <SeQuant/core/export/generator.hpp>
#include <SeQuant/core/export/itf.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/meta.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <nlohmann/json.hpp>

#include <algorithm>
#include <concepts>
#include <filesystem>
#include <fstream>
#include <functional>
#include <numeric>
#include <string>
#include <string_view>
#include <vector>

namespace sequant::util::extint {

SEQUANT_EXTINT_REGISTER_STEP_TYPE(ExportStep, "export");

std::string ExportStep::kind() const { return "export"; }

bool ExportStep::accepts_options() const { return true; }

bool ExportStep::requires_options() const { return true; }

void ExportStep::set_options(const nlohmann::json &options) {
  if (!options.is_object()) {
    throw Exception(kind() + " expects a JSON object for its options!");
  }

  for (const auto &[key, value] : options.items()) {
    if (key == "language") {
      if (!value.is_string()) {
        throw Exception("Value for " + kind() + " option '" + key +
                        "' must be a string");
      }

      language_ = value.get<std::string>();
    } else if (key == "optimize") {
      if (!value.is_boolean()) {
        throw Exception("Value for " + kind() + " option '" + key +
                        "' must be a boolean");
      }

      optimize_ = value.get<bool>();
    } else if (key == "output") {
      if (!value.is_string()) {
        throw Exception("Value for " + kind() + " option '" + key +
                        "' must be a string");
      }

      filepath_ = value.get<std::string>();
    } else {
      throw Exception("Unknown option key for " + kind() + ": '" + key + "'");
    }

    if (language_.empty()) {
      throw Exception("The 'language' option for " + kind() + " is mandatory");
    }

    if (filepath_.empty()) {
      throw Exception("The 'output' option for " + kind() + " is mandatory");
    }
  }
}

struct ExportTreeDataCompare {
  bool operator()(const ExportTreeData &lhs, const ExportTreeData &rhs) const {
    if (lhs.entries.size() != rhs.entries.size()) {
      return lhs.entries.size() < rhs.entries.size();
    }

    for (std::size_t i = 0; i < lhs.entries.size(); ++i) {
      if (!equal(lhs.entries.at(i), rhs.entries.at(i))) {
        return (*this)(lhs.entries.at(i), rhs.entries.at(i));
      }
    }

    return false;
  }

  bool equal(const ExportTreeData::Entry &lhs,
             const ExportTreeData::Entry &rhs) const {
    if (lhs.symm_contribution_target != rhs.symm_contribution_target) {
      return false;
    }

    if (lhs.tree->is_tensor() != rhs.tree->is_tensor()) {
      return false;
    }

    return lhs.tree->is_tensor()
               ? lhs.tree->as_tensor() == rhs.tree->as_tensor()
               : lhs.tree->as_variable() == rhs.tree->as_variable();
  }

  bool operator()(const ExportTreeData::Entry &lhs,
                  const ExportTreeData::Entry &rhs) const {
    if (lhs.symm_contribution_target != rhs.symm_contribution_target) {
      return lhs.symm_contribution_target < rhs.symm_contribution_target;
    }

    if (lhs.tree->is_tensor() != rhs.tree->is_tensor()) {
      return rhs.tree->is_tensor();
    }

    return lhs.tree->is_tensor()
               ? lhs.tree->as_tensor() < rhs.tree->as_tensor()
               : lhs.tree->as_variable() < rhs.tree->as_variable();
  }
};

std::size_t ExportStep::run(std::string_view, ExecutionContext &exctx,
                            const std::vector<std::string_view> &inputs) {
  std::vector<
      std::reference_wrapper<const ExecutionContext::Data<ProcessingData>>>
      data;
  for (std::string_view current_input : inputs) {
    for (const ExecutionContext::Data<ProcessingData> &current :
         exctx.get_data(current_input)) {
      data.push_back(std::cref(current));
    }
  }

  auto to_export_data = [&](std::size_t idx) {
    return convert_data<ExportTreeData>(data.at(idx).get().data.get());
  };

  std::vector<std::size_t> access_order(data.size());
  std::iota(access_order.begin(), access_order.end(), 0);

  // Sort inputs in partitions of contributions to the same result
  std::ranges::sort(access_order, ExportTreeDataCompare{}, to_export_data);

  std::vector<ExpressionGroup<>> groups;
  for (std::size_t idx : access_order) {
    const ExportTreeData &current_data = to_export_data(idx);

    if (std::ranges::any_of(current_data.entries, [&](const auto &entry) {
          return !ExportTreeDataCompare{}.equal(entry,
                                                current_data.entries.front());
        })) {
      throw Exception(
          "Expected all entries in an ExportTreeData object to contribute to "
          "the same result");
    }

    // TODO: determine group name from metadata
    std::string name = "We'll see";

    auto it =
        std::ranges::find_if(groups, [&name](const ExpressionGroup<> &group) {
          return group.is_named() && group.name() == name;
        });
    if (it == groups.end()) {
      groups.emplace_back(std::move(name));
      it = groups.end() - 1;
    }

    ExpressionGroup<> &group = *it;

    for (const ExportTreeData::Entry &current : current_data.entries) {
      group.add(current.tree);
    }

    // TODO: if current_data was last in its partition (based on above sorting)
    // we need to check if symmetrization is required and if so, generate the
    // code snippet that does it and add it to the group
  }

  auto perform_export = [&](auto &&generator, const auto &genctx)
    requires(
        std::derived_from<std::remove_cvref_t<decltype(genctx)>,
                          ExportContext> &&
        std::derived_from<std::remove_cvref_t<decltype(generator)>,
                          Generator<std::remove_cvref_t<decltype(genctx)>>>)
  {
    // TODO: setup load/store/create/import strategies

	// TODO: handle index batching

    export_groups(groups, generator, genctx);

    return generator.get_generated_code();
  };

  std::string generated;

  if (language_ == "itf") {
    ItfGenerator<ItfContext> generator;
    ItfContext genctx;

    if (optimize_) {
      generated = perform_export(
          GenerationOptimizer<decltype(generator)>(std::move(generator)),
          genctx);
    } else {
      generated = perform_export(generator, genctx);
    }
  } else {
    throw Exception("Unsupported export language '" + language_ + "'");
  }

  SEQUANT_ASSERT(!generated.empty());

  std::ofstream stream(filepath_);
  stream << generated;

  return 0;
}

}  // namespace sequant::util::extint
