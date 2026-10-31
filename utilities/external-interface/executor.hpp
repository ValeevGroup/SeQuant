#ifndef SEQUANT_EXTERNAL_INTERFACE_EXECUTOR_HPP
#define SEQUANT_EXTERNAL_INTERFACE_EXECUTOR_HPP

#include "execution_context.hpp"

#include <nlohmann/json_fwd.hpp>

namespace sequant::util::extint {

class Executor {
 public:
  Executor() = default;

  /// Executes the given pipeline of processing steps. Every call starts from
  /// a fresh context, i.e. steps can't reference data produced in previous
  /// calls. Data is released as soon as no later step references it anymore.
  void execute(const nlohmann::json &steps);

  const ExecutionContext &context() const { return context_; }

  void reset();

 private:
  ExecutionContext context_;
};

}  // namespace sequant::util::extint

#endif  // SEQUANT_EXTERNAL_INTERFACE_EXECUTOR_HPP
