#ifndef SEQUANT_EXTERNAL_INTERFACE_SPINTRACINGSTEP_HPP
#define SEQUANT_EXTERNAL_INTERFACE_SPINTRACINGSTEP_HPP

#include "execution_context.hpp"
#include "processing_data.hpp"
#include "processing_step.hpp"

#include <SeQuant/core/space.hpp>

#include <nlohmann/json_fwd.hpp>

#include <string>
#include <string_view>

namespace sequant::util::extint {

class SpintracingStep : public OneByOneProcessingStep<ExpressionData> {
 public:
  std::string kind() const override;

  bool accepts_options() const override;
  bool requires_options() const override;
  void set_options(const nlohmann::json &options) override;

 protected:
  std::size_t process(std::string_view id_prefix, std::size_t id_start,
                      ExecutionContext &ctx,
                      const ExpressionData &data) override;

 private:
  bool use_closed_shell_algo_ = false;
};

}  // namespace sequant::util::extint

#endif  // SEQUANT_EXTERNAL_INTERFACE_SPINTRACINGSTEP_HPP
