#ifndef SEQUANT_EXTERNAL_INTERFACE_SUBSTITUTESTEP_HPP
#define SEQUANT_EXTERNAL_INTERFACE_SUBSTITUTESTEP_HPP

#include "execution_context.hpp"
#include "processing_data.hpp"
#include "processing_step.hpp"

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/utility/expr_matcher.hpp>

#include <nlohmann/json_fwd.hpp>

#include <functional>
#include <map>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace sequant::util::extint {

class SubstituteStep : public OneByOneProcessingStep<ExpressionData, false,
                                                     InputAccess::WorkingCopy> {
 public:
  std::string kind() const override;

  bool accepts_options() const override;
  bool requires_options() const override;
  void set_options(const nlohmann::json &options) override;

 protected:
  std::size_t process(std::string_view id_prefix, std::size_t id_start,
                      ExecutionContext &ctx,
                      WorkingCopy<ExpressionData> &input) override;

  // Inputs left unchanged by the substitutions are dropped
  bool needs_original_input() const override;

 private:
  std::vector<std::pair<ExprPtr, ExprPtr>> substitutions_;
  TensorComparison tensor_mode_ = TensorComparison::Block;
  std::map<std::string, std::string, std::less<>> result_labels_;
};

}  // namespace sequant::util::extint

#endif  // SEQUANT_EXTERNAL_INTERFACE_SUBSTITUTESTEP_HPP
