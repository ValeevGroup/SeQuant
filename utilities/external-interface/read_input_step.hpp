#ifndef SEQUANT_EXTERNAL_INTERFACE_READINPUTSTEP_HPP
#define SEQUANT_EXTERNAL_INTERFACE_READINPUTSTEP_HPP

#include "processing_step.hpp"

#include <nlohmann/json_fwd.hpp>

#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

namespace sequant::util::extint {

class ReadInputStep : public ProcessingStep {
 public:
  std::string kind() const override;

  bool accepts_options() const override;
  bool requires_options() const override;
  void set_options(const nlohmann::json &options) override;

  std::size_t run(std::string_view step_id, ExecutionContext &ctx,
                  const std::vector<std::string_view> &inputs = {}) override;

 private:
  std::vector<std::filesystem::path> input_paths_;
};

}  // namespace sequant::util::extint

#endif  // SEQUANT_EXTERNAL_INTERFACE_READINPUTSTEP_HPP
