#ifndef SEQUANT_EXTERNAL_INTERFACE_INDEX_BATCHING_STEP_HPP
#define SEQUANT_EXTERNAL_INTERFACE_INDEX_BATCHING_STEP_HPP

#include "execution_context.hpp"
#include "processing_data.hpp"
#include "processing_step.hpp"

#include <SeQuant/core/space.hpp>

#include <nlohmann/json_fwd.hpp>

#include <limits>
#include <string>
#include <string_view>

namespace sequant::util::extint {

class IndexBatchingStep : public OneByOneProcessingStep<ExportTreeData> {
 public:
  std::string kind() const override;

  bool accepts_options() const override;
  bool requires_options() const override;
  void set_options(const nlohmann::json &options) override;

 protected:
  std::size_t process(std::string_view id_prefix, std::size_t id_start,
                      ExecutionContext &ctx,
                      const ExportTreeData &data) override;

  bool alias_unchanged_inputs() const override;

 private:
  std::size_t min_unbatched_ = 2;
  std::size_t max_batched_ = std::numeric_limits<std::size_t>::max();
  bool start_with_largest_ = true;
};

}  // namespace sequant::util::extint

#endif  // SEQUANT_EXTERNAL_INTERFACE_INDEX_BATCHING_STEP_HPP
