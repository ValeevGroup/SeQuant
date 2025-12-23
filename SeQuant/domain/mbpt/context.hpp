
#ifndef SEQUANT_DOMAIN_MBPT_CONTEXT_HPP
#define SEQUANT_DOMAIN_MBPT_CONTEXT_HPP

#include <SeQuant/domain/mbpt/fwd.hpp>

#include <SeQuant/core/utility/context.hpp>
#include <SeQuant/domain/mbpt/op_registry.hpp>

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
  struct Options {
    /// shared pointer to operator registry
    std::shared_ptr<OpRegistry> op_registry_ptr = nullptr;
    /// optional operator registry, use if no shared pointer is provided
    std::optional<OpRegistry> op_registry = std::nullopt;
    /// whether to use cluster-specific virtuals
    CSV csv = Defaults::csv;
  };

  /// @brief makes default options for mbpt::Context
  static Options make_default_options() {
    return Options{};
  }


  CSV csv() const { return csv_; }

 private:
  CSV csv_ = Defaults::csv;
};

bool operator==(Context const& left, Context const& right);

bool operator!=(Context const& left, Context const& right);

const Context& get_default_mbpt_context();

void set_default_mbpt_context(const Context& ctx);

void reset_default_mbpt_context();
[[nodiscard]] sequant::detail::ImplicitContextResetter<Context>
set_scoped_default_mbpt_context(const Context& ctx);

}  // namespace sequant::mbpt

#endif  // SEQUANT_DOMAIN_MBPT_CONTEXT_HPP
