#ifndef SEQUANT_EXTERNAL_INTERFACE_OPTIMIZATIONSTEP_HPP
#define SEQUANT_EXTERNAL_INTERFACE_OPTIMIZATIONSTEP_HPP

#include "execution_context.hpp"
#include "filter_step.hpp"
#include "processing_data.hpp"
#include "processing_step.hpp"

#include <SeQuant/core/optimize/options.hpp>

#include <nlohmann/json_fwd.hpp>

#include <deque>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace sequant::util::extint {

class OptimizationStep : public OneByOneProcessingStep<ExpressionData> {
 public:
  std::string kind() const override;

  bool accepts_options() const override;
  bool requires_options() const override;
  void set_options(const nlohmann::json &options) override;

  std::size_t run(std::string_view step_id, ExecutionContext &ctx,
                  const std::vector<std::string_view> &inputs = {}) override;

 protected:
  OptimizeOptions options_;

  std::size_t process(std::string_view id_prefix, std::size_t id_start,
                      ExecutionContext &ctx,
                      const ExpressionData &data) override;

 private:
  struct ReplayOptions {
    std::shared_ptr<const ExpressionFilter> is_volatile;
    std::string intermediate_label = "Persistent";
    std::string persistent_id = "persistent";
    std::string volatile_id = "volatile";
  };

  std::optional<ReplayOptions> replay_;
  /// (persistent, volatile) parts of every input, in processing order
  std::deque<std::pair<ExpressionData, ExpressionData>> pending_;

  void set_replay_options(const nlohmann::json &options);
};

}  // namespace sequant::util::extint

#endif  // SEQUANT_EXTERNAL_INTERFACE_OPTIMIZATIONSTEP_HPP
