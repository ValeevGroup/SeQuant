#include "executor.hpp"
#include "processing_step.hpp"
#include "processing_step_factory.hpp"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <chrono>
#include <iostream>
#include <optional>
#include <ranges>
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

}  // namespace

void Executor::execute(const nlohmann::json &steps) {
  if (!steps.is_array()) {
    throw Exception("Steps object must be an array");
  }

  std::size_t step_id_counter = 0;

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

    std::string step_id = step.contains("id")
                              ? step.at("id").get<std::string>()
                              : "step" + std::to_string(step_id_counter) + "." +
                                    proc_step->kind();

    if (num_outputs_.find(step_id) != num_outputs_.end()) {
      throw Exception("Duplicate step ID '" + step_id + "'");
    }

    std::cout << "Executing '" << step_id << "' (" << kind << ")... ";
    std::cout.flush();

    std::chrono::steady_clock::time_point start =
        std::chrono::steady_clock::now();

    const std::size_t produced_outputs =
        proc_step->run(step_id, context_, inputs);

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

    num_outputs_.emplace(step_id, produced_outputs);

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

    ++step_id_counter;
  }
}

void Executor::reset() {
  num_outputs_.clear();
  context_ = {};
}

}  // namespace sequant::util::extint
