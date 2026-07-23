#ifndef SEQUANT_EXTERNAL_INTERFACE_READINPUTSTEP_HPP
#define SEQUANT_EXTERNAL_INTERFACE_READINPUTSTEP_HPP

#include "processing_step.hpp"

#include <nlohmann/json_fwd.hpp>

#include <filesystem>
#include <string>
#include <vector>

namespace sequant::util::extint {

class ReadInputStep : public ProcessingStep {
 public:
  std::string kind() const override;

  bool accepts_options() const override;
  bool requires_options() const override;
  void set_options(const nlohmann::json &options) override;

  ProcessingOutput run(const ProcessingInput &) override;
  ExpressionOutput run(const VoidInput &);

 private:
  std::vector<std::filesystem::path> input_paths_;
};

}  // namespace sequant::util::extint

#endif  // SEQUANT_EXTERNAL_INTERFACE_READINPUTSTEP_HPP
