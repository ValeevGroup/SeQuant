#ifndef SEQUANT_EXTERNAL_INTERFACE_DEPENDENCY_STEP_HPP
#define SEQUANT_EXTERNAL_INTERFACE_DEPENDENCY_STEP_HPP

#include "execution_context.hpp"
#include "processing_step.hpp"

#include <nlohmann/json_fwd.hpp>

#include <cstddef>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace sequant::util::extint {

/// Populates the inter-tree dependencies (ExportTreeData::Entry::dependencies)
/// that the export step uses to order results. Two complementary, additive
/// modes are available and at least one of them must be configured:
/// - autodetect: whenever the result of one tree is used as an input in another
///   tree, the former becomes a dependency of the latter (block-level match);
/// - deps: dependencies given explicitly via external interface IDs.
/// The step annotates its inputs in place and aliases them under its own id so
/// that downstream steps can reference them as usual.
class DependencyStep : public ProcessingStep {
 public:
  std::string kind() const override;

  bool accepts_options() const override;
  bool requires_options() const override;
  void set_options(const nlohmann::json &options) override;

  std::vector<std::string> referenced_ids() const override;

  std::size_t run(std::string_view step_id, ExecutionContext &ctx,
                  const std::vector<std::string_view> &inputs = {}) override;

 private:
  bool autodetect_ = false;
  std::vector<std::pair<std::string, std::vector<std::string>>> deps_;
};

}  // namespace sequant::util::extint

#endif  // SEQUANT_EXTERNAL_INTERFACE_DEPENDENCY_STEP_HPP
