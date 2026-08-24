#ifndef SEQUANT_EXTERNAL_INTERFACE_EXPORTSTEP_HPP
#define SEQUANT_EXTERNAL_INTERFACE_EXPORTSTEP_HPP

#include "execution_context.hpp"
#include "processing_step.hpp"

#include <nlohmann/json_fwd.hpp>

#include <filesystem>
#include <functional>
#include <map>
#include <string>
#include <string_view>

namespace sequant::util::extint {

class ExportStep : public ProcessingStep {
 public:
  std::string kind() const override;

  bool accepts_options() const override;
  bool requires_options() const override;
  void set_options(const nlohmann::json &options) override;

  std::size_t run(std::string_view step_id, ExecutionContext &ctx,
                  const std::vector<std::string_view> &inputs = {}) override;

 protected:
  std::string language_;
  bool optimize_ = true;
  std::filesystem::path filepath_;
  std::map<std::string, std::string, std::less<>> group_assoc_;
};

}  // namespace sequant::util::extint

#endif  // SEQUANT_EXTERNAL_INTERFACE_EXPORTSTEP_HPP
