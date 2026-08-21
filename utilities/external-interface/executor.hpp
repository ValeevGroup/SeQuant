#ifndef SEQUANT_EXTERNAL_INTERFACE_EXECUTOR_HPP
#define SEQUANT_EXTERNAL_INTERFACE_EXECUTOR_HPP

#include "execution_context.hpp"

#include <nlohmann/json_fwd.hpp>

#include <cstddef>
#include <functional>
#include <map>
#include <string>

namespace sequant::util::extint {

class Executor {
 public:
  Executor() = default;

  void execute(const nlohmann::json &steps);

  void reset();

 private:
  ExecutionContext context_;
  std::map<std::string, std::size_t, std::less<>> num_outputs_;
};

}  // namespace sequant::util::extint

#endif  // SEQUANT_EXTERNAL_INTERFACE_EXECUTOR_HPP
