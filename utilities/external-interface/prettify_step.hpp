#ifndef SEQUANT_EXTERNAL_INTERFACE_PRETTIFYSTEP_HPP
#define SEQUANT_EXTERNAL_INTERFACE_PRETTIFYSTEP_HPP

#include "execution_context.hpp"
#include "processing_data.hpp"
#include "processing_step.hpp"

#include <nlohmann/json_fwd.hpp>

#include <string>
#include <string_view>

namespace sequant::util::extint {

/// brings the expressions into a readable form whose order does not depend on
/// the hashing algorithm (see sequant::prettify()), using the conventions of
/// mbpt (mbpt::prettify_options()); the parenthesization of the expressions
/// is preserved
class PrettifyStep : public OneByOneProcessingStep<ExpressionData> {
 public:
  std::string kind() const override;

  bool accepts_options() const override;
  bool requires_options() const override;
  void set_options(const nlohmann::json &options) override;

 protected:
  std::size_t process(std::string_view id_prefix, std::size_t id_start,
                      ExecutionContext &ctx,
                      const ExpressionData &data) override;

  bool alias_unchanged_inputs() const override;

 private:
  bool recanonicalize_ = true;
};

}  // namespace sequant::util::extint

#endif  // SEQUANT_EXTERNAL_INTERFACE_PRETTIFYSTEP_HPP
