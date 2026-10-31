#include "executor.hpp"
#include "processing_step.hpp"
#include "processing_step_factory.hpp"

#include <nlohmann/json.hpp>

#include <spdlog/spdlog.h>

#include <algorithm>
#include <chrono>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <optional>
#include <ranges>
#include <set>
#include <string>
#include <string_view>
#include <vector>

namespace sequant::util::extint {

namespace {

/// @param chain IDs of the steps whose options are currently being resolved
std::optional<nlohmann::json> resolve_options(const nlohmann::json &steps,
                                              const nlohmann::json &step,
                                              std::vector<std::string> &chain) {
  std::optional<nlohmann::json> options;

  if (step.contains("inherit_options_from")) {
    std::string src_id = step.at("inherit_options_from").get<std::string>();

    if (std::ranges::find(chain, src_id) != chain.end()) {
      std::string cycle;
      for (const std::string &id : chain) {
        cycle += "'" + id + "' -> ";
      }
      throw Exception("Cyclic option inheritance: " + cycle + "'" + src_id +
                      "'");
    }

    auto it = std::ranges::find_if(steps, [&](const nlohmann::json &current) {
      return current.contains("id") && current.at("id") == src_id;
    });
    if (it == steps.end()) {
      throw Exception("Can't inherit options from unknown step '" + src_id +
                      "'");
    }

    if (it->at("kind") != step.at("kind")) {
      throw Exception("Step of kind '" + step.at("kind").get<std::string>() +
                      "' can't inherit options from step '" + src_id +
                      "' of kind '" + it->at("kind").get<std::string>() + "'");
    }

    chain.push_back(std::move(src_id));
    options = resolve_options(steps, *it, chain);
    chain.pop_back();
  }

  if (step.contains("options")) {
    if (options) {
      options->merge_patch(step.at("options"));
    } else {
      options = step.at("options");
    }
  }

  return options;
}

/// @returns The index of the step (in @p step_indices) that created the data
/// referenced by @p id (which must not contain selectors), or std::nullopt if
/// there is no such step. IDs always start with the ID of the step that created
/// them, but since step IDs may themselves contain periods, we have to find the
/// longest matching step ID.
std::optional<std::size_t> find_owning_step(
    const std::map<std::string, std::size_t, std::less<>> &step_indices,
    std::string_view id) {
  for (std::size_t len = id.size(); len > 0; --len) {
    if (len < id.size() && id[len] != '.') {
      continue;
    }

    if (auto it = step_indices.find(id.substr(0, len));
        it != step_indices.end()) {
      return it->second;
    }
  }

  return std::nullopt;
}

struct PreparedStep {
  std::string_view kind;
  std::string id;
  std::unique_ptr<ProcessingStep> step;
  std::vector<std::string_view> inputs;
};

}  // namespace

void Executor::execute(const nlohmann::json &steps) {
  if (!steps.is_array()) {
    throw Exception("Steps object must be an array");
  }

  // Data is released as soon as the pipeline no longer needs it, so data from
  // a previous call can't be relied upon anyway
  reset();

  std::vector<PreparedStep> prepared;
  prepared.reserve(steps.size());
  std::map<std::string, std::size_t, std::less<>> step_indices;

  for (const nlohmann::json &step : steps) {
    const std::string_view kind = step.at("kind").get<std::string_view>();

    std::vector<std::string_view> inputs;

    if (step.contains("inputs")) {
      const nlohmann::json &inps = step.at("inputs");
      if (inps.is_string()) {
        inputs.emplace_back(inps.get<std::string_view>());
      } else if (inps.is_array()) {
        for (const auto &current : inps) {
          if (!current.is_string()) {
            throw Exception("Entries in inputs array must be strings");
          }

          inputs.emplace_back(current.get<std::string_view>());
        }
      } else {
        throw Exception("inputs field must be either a string or an array");
      }
    }

    std::unique_ptr<ProcessingStep> proc_step =
        ProcessingStepFactory::instance().instantiate(kind);

    std::vector<std::string> chain;
    if (step.contains("id")) {
      chain.push_back(step.at("id").get<std::string>());
    }

    if (std::optional<nlohmann::json> options =
            resolve_options(steps, step, chain)) {
      if (!proc_step->accepts_options()) {
        throw Exception("Processing step '" + std::string(kind) +
                        "' does not take options but some where given");
      }

      proc_step->set_options(*options);
    } else if (proc_step->requires_options()) {
      throw Exception("Processing step '" + std::string(kind) +
                      "' requires options but none where given");
    }

    if (proc_step->requires_inputs() && inputs.empty()) {
      throw Exception("Processing step '" + std::string(kind) +
                      "' requires at least one input but none were given");
    }

    std::string step_id = step.contains("id")
                              ? step.at("id").get<std::string>()
                              : "step" + std::to_string(prepared.size()) + "." +
                                    proc_step->kind();

    if (!step_indices.emplace(step_id, prepared.size()).second) {
      throw Exception("Duplicate step ID '" + step_id + "'");
    }

    prepared.push_back({.kind = kind,
                        .id = std::move(step_id),
                        .step = std::move(proc_step),
                        .inputs = std::move(inputs)});
  }

  // Determine for every step the index of the last step that accesses any of
  // its data. Afterwards, that data is no longer needed and can be released.
  std::vector<std::size_t> last_use(prepared.size());
  std::iota(last_use.begin(), last_use.end(), std::size_t{0});

  for (std::size_t i = 0; i < prepared.size(); ++i) {
    std::vector<std::string> references = prepared[i].step->referenced_ids();
    references.insert(references.end(), prepared[i].inputs.begin(),
                      prepared[i].inputs.end());

    for (const std::string &ref : references) {
      std::vector<std::string> expanded;
      try {
        expanded = ExecutionContext::expand_id(ref);
      } catch (const std::exception &e) {
        throw Exception("Invalid reference in step '" + prepared[i].id + "' (" +
                        std::string(prepared[i].kind) + "): " + e.what());
      }

      for (const std::string &id : expanded) {
        std::optional<std::size_t> owner = find_owning_step(step_indices, id);

        // References to unknown or future steps will produce an error once
        // the referencing step is executed
        if (owner.has_value() && *owner < i) {
          last_use[*owner] = std::max(last_use[*owner], i);
        }
      }
    }
  }

  // IDs for which no owning step could be determined (reported only once)
  std::set<std::string, std::less<>> unowned_ids;

  for (std::size_t i = 0; i < prepared.size(); ++i) {
    const nlohmann::json &step = steps.at(i);
    const std::string &step_id = prepared[i].id;

    std::cout << "Executing '" << step_id << "' (" << prepared[i].kind
              << ")... ";
    std::cout.flush();

    std::chrono::steady_clock::time_point start =
        std::chrono::steady_clock::now();

    const std::size_t produced_outputs =
        prepared[i].step->run(step_id, context_, prepared[i].inputs);

    std::chrono::steady_clock::duration delta =
        std::chrono::steady_clock::now() - start;
    if (delta > std::chrono::minutes(1)) {
      std::cout << std::chrono::duration_cast<std::chrono::minutes>(delta)
                << std::endl;
    } else if (delta > std::chrono::seconds(1)) {
      std::cout << std::chrono::duration_cast<std::chrono::seconds>(delta)
                << std::endl;
    } else {
      std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(delta)
                << std::endl;
    }

    if (step.contains("outputs")) {
      const nlohmann::json &outputs = step.at("outputs");
      if (!outputs.is_object()) {
        throw Exception("outputs field must be an object");
      }

      for (const auto &[key, val] : outputs.items()) {
        context_.add_data_alias(step_id + "[" + val.get<std::string>() + "]",
                                step_id + "." + key);
      }
    }

    if (produced_outputs > 0) {
      context_.add_data_alias(
          std::ranges::views::iota(std::size_t(0), produced_outputs) |
              std::ranges::views::transform([&step_id](std::size_t num) {
                return std::string(step_id) + "." + std::to_string(num);
              }),
          std::string(step_id));
    }

    // Data is still needed if any of its IDs is referenced by a later step.
    // Note that steps that alias their inputs under their own ID (instead of
    // creating new data) thereby extend the lifetime of these inputs.
    const std::size_t n_released =
        context_.release_data([&](std::string_view id) {
          std::optional<std::size_t> owner = find_owning_step(step_indices, id);
          if (!owner.has_value()) {
            // Be conservative and assume that the data is still needed
            if (unowned_ids.emplace(id).second) {
              spdlog::debug(
                  "Can't determine the step owning ID '{}' - associated data "
                  "will not be released",
                  id);
            }
            return true;
          }
          return last_use[*owner] > i;
        });

    if (n_released > 0) {
      spdlog::debug("Released {} data object(s) after step '{}'", n_released,
                    step_id);
    }
  }
}

void Executor::reset() { context_ = {}; }

}  // namespace sequant::util::extint
