#ifndef SEQUANT_EXTERNAL_INTERFACE_VALIDATESTEP_HPP
#define SEQUANT_EXTERNAL_INTERFACE_VALIDATESTEP_HPP

#include "processing_data.hpp"
#include "processing_step.hpp"

#include <nlohmann/json_fwd.hpp>

#include <string>
#include <string_view>
#include <vector>

namespace sequant::util::extint {

class ValidateStep : public OneByOneProcessingStep<ExpressionData> {
 public:
  std::string kind() const override;

  bool accepts_options() const override;
  bool requires_options() const override;
  void set_options(const nlohmann::json &options) override;

  std::size_t process(std::string_view id_prefix, std::size_t id_start,
                      ExecutionContext &ctx,
                      const ExpressionData &data) override;
};

}  // namespace sequant::util::extint

#endif  // SEQUANT_EXTERNAL_INTERFACE_VALIDATESTEP_HPP
