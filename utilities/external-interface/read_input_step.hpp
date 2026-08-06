#ifndef SEQUANT_EXTERNAL_INTERFACE_READINPUTSTEP_HPP
#define SEQUANT_EXTERNAL_INTERFACE_READINPUTSTEP_HPP

#include "processing_data.hpp"
#include "processing_step.hpp"

#include <SeQuant/core/io/serialization/serialization.hpp>

#include <nlohmann/json_fwd.hpp>

#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

namespace sequant::util::extint {

class ReadInputStep : public OneByOneProcessingStep<void> {
 public:
  std::string kind() const override;

  bool accepts_options() const override;
  bool requires_options() const override;
  void set_options(const nlohmann::json &options) override;

 protected:
  std::size_t process(std::string_view id_prefix, std::size_t id_start,
                      ExecutionContext &ctx) override;

 private:
  std::vector<std::filesystem::path> input_paths_;
  io::serialization::DeserializationOptions options_;
};

}  // namespace sequant::util::extint

#endif  // SEQUANT_EXTERNAL_INTERFACE_READINPUTSTEP_HPP
