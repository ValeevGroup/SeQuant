#include "executor.hpp"
#include "processing_step.hpp"
#include "processing_step_factory.hpp"

#include <nlohmann/json.hpp>

namespace sequant::util::extint {

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

    if (step.contains("options")) {
      if (!proc_step->accepts_options()) {
        throw Exception("Processing step '" + std::string(kind) +
                        "' does not take options but some where given");
      }

      proc_step->set_options(step.at("options"));
    } else if (proc_step->requires_options()) {
      throw Exception("Processing step '" + std::string(kind) +
                      "' requires options but none where given");
    }

    std::string step_id = step.contains("id")
                              ? step.at("id").get<std::string>()
                              : "step" + std::to_string(step_id_counter) + "." +
                                    proc_step->kind();

    const std::size_t produced_outputs =
        proc_step->run(step_id, context_, inputs);

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

    ++step_id_counter;
  }
}

void Executor::reset() {
  num_outputs_.clear();
  context_ = {};
}

}  // namespace sequant::util::extint
