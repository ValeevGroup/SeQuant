#ifndef SEQUANT_EXTERNAL_INTERFACE_OUTPUTSTEP_HPP
#define SEQUANT_EXTERNAL_INTERFACE_OUTPUTSTEP_HPP

#include "execution_context.hpp"
#include "processing_data.hpp"
#include "processing_step.hpp"

#include <nlohmann/json_fwd.hpp>

#include <cstddef>
#include <string_view>

namespace sequant::util::extint {

class OutputStep : public OneByOneProcessingStep<ExpressionData> {
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
  bool latex_ = false;
  bool annot_symm_ = true;
};

}  // namespace sequant::util::extint

#endif  // SEQUANT_EXTERNAL_INTERFACE_OUTPUTSTEP_HPP
