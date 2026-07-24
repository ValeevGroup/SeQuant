#ifndef SEQUANT_EXTERNAL_INTERFACE_OUTPUTSTEP_HPP
#define SEQUANT_EXTERNAL_INTERFACE_OUTPUTSTEP_HPP

#include "processing_step.hpp"

#include <nlohmann/json_fwd.hpp>

namespace sequant::util::extint {

class OutputStep : public ProcessingStep {
 public:
  std::string kind() const override;

  bool accepts_options() const override;
  bool requires_options() const override;
  void set_options(const nlohmann::json &options) override;

  std::size_t run(std::string_view step_id, ExecutionContext &ctx,
                  const std::vector<std::string_view> &inputs = {}) override;

 private:
  bool latex_ = false;
  bool annot_symm_ = true;
};

}  // namespace sequant::util::extint

#endif  // SEQUANT_EXTERNAL_INTERFACE_OUTPUTSTEP_HPP
