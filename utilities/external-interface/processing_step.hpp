#ifndef SEQUANT_EXTERNAL_INTERFACE_PROCESSINGSTEP_HPP
#define SEQUANT_EXTERNAL_INTERFACE_PROCESSINGSTEP_HPP

#include "execution_context.hpp"

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/meta.hpp>
#include <SeQuant/core/utility/exception.hpp>

#include <nlohmann/json_fwd.hpp>

#include <string>
#include <string_view>
#include <vector>

namespace sequant::util::extint {

class ProcessingStep {
 public:
  ProcessingStep() = default;
  virtual ~ProcessingStep();

  virtual std::string kind() const = 0;

  virtual bool accepts_options() const = 0;
  virtual bool requires_options() const = 0;
  virtual void set_options(const nlohmann::json &options) = 0;

  virtual std::size_t run(std::string_view step_id, ExecutionContext &ctx,
                          const std::vector<std::string_view> &inputs = {}) = 0;
};

}  // namespace sequant::util::extint

#endif  // SEQUANT_EXTERNAL_INTERFACE_PROCESSINGSTEP_HPP
