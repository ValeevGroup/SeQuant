
#ifndef SEQUANT_DOMAIN_MBPT_CONTEXT_HPP
#define SEQUANT_DOMAIN_MBPT_CONTEXT_HPP

#include <SeQuant/core/utility/context.hpp>

namespace sequant::mbpt {

/// Whether to use cluster-specific virtuals.
enum class CSV { Yes, No };

/// @brief details of the MBPT formalism
///
/// to be used as a implicit or explicit specification of the formalism
class Context {
 public:
  struct Defaults {
    constexpr static auto csv = CSV::No;
  };

  /// @brief Construct a Formalism object
  explicit Context(CSV csv_formalism = Defaults::csv) noexcept;

  CSV csv() const { return csv_; }

 private:
  CSV csv_ = Defaults::csv;
};

/// old name of Context is a deprecated alias
using Formalism [[deprecated(
    "use sequant::mbpt::Context instead of sequant::mbpt::Formalism")]] =
    Context;

bool operator==(Context const& left, Context const& right);

bool operator!=(Context const& left, Context const& right);

const Context& get_default_formalism();
void set_default_formalism(const Context& ctx);
void reset_default_formalism();
[[nodiscard]] detail::ImplicitContextResetter<Context>
set_scoped_default_formalism(const Context& ctx);

}  // namespace sequant::mbpt

#endif  // SEQUANT_DOMAIN_MBPT_CONTEXT_HPP
