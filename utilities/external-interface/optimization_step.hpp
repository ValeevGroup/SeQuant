#ifndef SEQUANT_EXTERNAL_INTERFACE_OPTIMIZATIONSTEP_HPP
#define SEQUANT_EXTERNAL_INTERFACE_OPTIMIZATIONSTEP_HPP

#include "execution_context.hpp"
#include "processing_data.hpp"
#include "processing_step.hpp"

#include <SeQuant/core/optimize/options.hpp>

#include <nlohmann/json_fwd.hpp>

#include <string>
#include <string_view>

namespace sequant::util::extint {

class OptimizationStep : public OneByOneProcessingStep<ExpressionData> {
 public:
  std::string kind() const override;

  bool accepts_options() const override;
  bool requires_options() const override;
  void set_options(const nlohmann::json &options) override;

 protected:
  OptimizeOptions options_;

  std::size_t process(std::string_view id_prefix, std::size_t id_start,
                      ExecutionContext &ctx,
                      const ExpressionData &data) override;
};

}  // namespace sequant::util::extint

#endif  // SEQUANT_EXTERNAL_INTERFACE_OPTIMIZATIONSTEP_HPP
