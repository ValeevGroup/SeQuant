#ifndef SEQUANT_EXTERNAL_INTERFACE_DENSITYFITTINGSTEP_HPP
#define SEQUANT_EXTERNAL_INTERFACE_DENSITYFITTINGSTEP_HPP

#include "processing_step.hpp"

#include <SeQuant/core/space.hpp>

#include <nlohmann/json_fwd.hpp>

#include <string>
#include <string_view>
#include <vector>

namespace sequant::util::extint {

class DensityFittingStep : public ProcessingStep {
 public:
  std::string kind() const override;

  bool accepts_options() const override;
  bool requires_options() const override;
  void set_options(const nlohmann::json &options) override;

  std::size_t run(std::string_view step_id, ExecutionContext &ctx,
                  const std::vector<std::string_view> &inputs = {}) override;

 private:
  IndexSpace aux_space_ = IndexSpace::null;
  std::string two_elec_int_label_ = "g";
  std::string df_label_ = "DF";
};

}  // namespace sequant::util::extint

#endif  // SEQUANT_EXTERNAL_INTERFACE_DENSITYFITTINGSTEP_HPP
