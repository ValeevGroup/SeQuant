#ifndef SEQUANT_EXTERNAL_INTERFACE_DENSITYFITTINGSTEP_HPP
#define SEQUANT_EXTERNAL_INTERFACE_DENSITYFITTINGSTEP_HPP

#include "execution_context.hpp"
#include "processing_data.hpp"
#include "processing_step.hpp"

#include <SeQuant/core/space.hpp>

#include <nlohmann/json_fwd.hpp>

#include <cstddef>
#include <string>
#include <string_view>

namespace sequant::util::extint {

class DensityFittingStep : public OneByOneProcessingStep<ExpressionData> {
 public:
  std::string kind() const override;

  bool accepts_options() const override;
  bool requires_options() const override;
  void set_options(const nlohmann::json &options) override;

 protected:
  std::size_t process(std::string_view id_prefix, std::size_t id_start,
                      ExecutionContext &ctx,
                      const ExpressionData &data) override;

  bool alias_unchanged_inputs() const override;

 private:
  IndexSpace aux_space_ = IndexSpace::null;
  std::string two_elec_int_label_ = "g";
  std::string df_label_ = "DF";
};

}  // namespace sequant::util::extint

#endif  // SEQUANT_EXTERNAL_INTERFACE_DENSITYFITTINGSTEP_HPP
