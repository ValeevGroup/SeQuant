#ifndef SEQUANT_EXTERNAL_INTERFACE_COSTREPORTSTEP_HPP
#define SEQUANT_EXTERNAL_INTERFACE_COSTREPORTSTEP_HPP

#include "execution_context.hpp"
#include "processing_step.hpp"

#include <SeQuant/core/utility/expr_filter.hpp>

#include <nlohmann/json_fwd.hpp>

#include <cstddef>
#include <functional>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace sequant::util::extint {

/// Writes a Markdown report on the cost of evaluating its (already optimized)
/// input expressions: largest intermediates, most expensive contractions,
/// peak storage, shape census and, optionally, a simulation of which
/// intermediates a cache would keep. Every summand is binarized as it is, so
/// the factorization reported is the one the inputs carry.
class CostReportStep : public ProcessingStep {
 public:
  std::string kind() const override;

  bool accepts_options() const override;
  bool requires_options() const override;
  void set_options(const nlohmann::json &options) override;

  std::size_t run(std::string_view step_id, ExecutionContext &ctx,
                  const std::vector<std::string_view> &inputs = {}) override;

  struct CacheOptions {
    std::shared_ptr<const ExpressionFilter> is_volatile;
    /// The 'volatile' option as given (serialized JSON), quoted in the report
    std::string volatile_spec;
    std::size_t min_repeats = 2;
    /// Intermediates with more elements than this are not cached; 0 means
    /// unlimited
    double max_footprint = 0.0;
  };

  struct Options {
    /// Path of the report; empty means `<step id>.md`
    std::string output;
    /// Rows in the largest/most-expensive tables
    std::size_t top_n = 20;
    /// Also write each expression's binarized trees to
    /// `<output stem>.<name>.tree.txt` in the directory of `output`
    bool dump_tree = false;
    /// Leave out the SeQuant revision (for reproducible reports)
    bool omit_revision = false;
    /// Index-space label -> character (one code point, UTF-8 encoded) used for
    /// it in shape signatures; spaces not listed are shown as 'X'
    std::map<std::wstring, std::string, std::less<>> space_symbols;
    std::optional<CacheOptions> cache;
  };

 private:
  Options options_;

  /// Throws if a space_symbols key is not a registered, spin-free index-space
  /// label
  void check_space_symbols() const;
};

}  // namespace sequant::util::extint

#endif  // SEQUANT_EXTERNAL_INTERFACE_COSTREPORTSTEP_HPP
