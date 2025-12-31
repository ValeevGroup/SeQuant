
#ifndef SEQUANT_DOMAIN_MBPT_CONTEXT_HPP
#define SEQUANT_DOMAIN_MBPT_CONTEXT_HPP

#include <SeQuant/domain/mbpt/fwd.hpp>

#include <SeQuant/core/utility/context.hpp>
#include <SeQuant/domain/mbpt/op_registry.hpp>

namespace sequant::mbpt {

/// Whether to use cluster-specific virtuals.
enum class CSV { Yes, No };

// clang-format off
/// @brief Specifies details of the MBPT formalism
///
/// MBPT context contains:
/// - csv: whether to use cluster-specific virtuals
/// - an `OpRegistry`: maps operator labels to their `OpClass`
///
/// @warning Default construction creates a `Context` with null `OpRegistry`.
///          Most MBPT functions require a registry. Can be initialized as:
///          @code
///          set_default_mbpt_context(Context({.op_registry_ptr = make_minimal_registry()}));
///          @endcode
///
/// Can be accessed via `get_default_mbpt_context()`. Functions access this global context to validate operator types and properties.
// clang-format on
class Context {
 public:
  struct Defaults {
    constexpr static auto csv = CSV::No;
  };

  struct Options {
    /// whether to use cluster-specific virtuals
    CSV csv = Defaults::csv;
    /// shared pointer to operator registry
    std::shared_ptr<OpRegistry> op_registry_ptr = nullptr;
    /// optional operator registry, use if no shared pointer is provided
    std::optional<OpRegistry> op_registry = std::nullopt;
  };

  /// @brief makes default options for mbpt::Context
  static Options make_default_options() { return Options{}; }

  /// @brief Construct a Context object, uses default options if none are given
  Context(Options options = make_default_options());

  /// @brief destructor
  ~Context() = default;

  /// @brief move constructor
  Context(Context&&) = default;

  /// @brief copy constructor
  Context(Context const& other) = default;

  /// @brief copy assignment
  Context& operator=(Context const& other) = default;

  /// @brief clones this object and its OpRegistry
  Context clone() const;

  /// @return the value of CSV in this context
  CSV csv() const;

  /// @return a constant pointer to the OpRegistry for this context
  /// @warning can be null when user did not provide one to Context (i.e., it
  /// was default constructed)
  std::shared_ptr<const OpRegistry> op_registry() const;

  /// @return a pointer to the OpRegistry for this context
  /// @warning can be null when user did not provide one to Context (i.e., it
  /// was default constructed)
  std::shared_ptr<OpRegistry> mutable_op_registry() const;

  /// @brief sets the OpRegistry for this context
  Context& set(const OpRegistry& op_registry);

  /// @brief sets the OpRegistry for this context
  Context& set(std::shared_ptr<OpRegistry> op_registry);

  /// @brief sets whether to use cluster-specific virtuals
  Context& set(CSV csv);

 private:
  CSV csv_ = Defaults::csv;
  std::shared_ptr<OpRegistry> op_registry_;
};

bool operator==(Context const& left, Context const& right);

bool operator!=(Context const& left, Context const& right);

const Context& get_default_mbpt_context();

void set_default_mbpt_context(const Context& ctx);

void set_default_mbpt_context(const Context::Options& options);

void reset_default_mbpt_context();

[[nodiscard]] sequant::detail::ImplicitContextResetter<Context>
set_scoped_default_mbpt_context(const Context& ctx);

[[nodiscard]] sequant::detail::ImplicitContextResetter<Context>
set_scoped_default_mbpt_context(const Context::Options& options);

/// predefined operator registries

/// @brief makes a minimal operator registry with only essential operators for
/// MBPT
std::shared_ptr<OpRegistry> make_minimal_registry();

/// @brief make a legacy operator registry with SeQuant's old predefined
/// operators set
std::shared_ptr<OpRegistry> make_legacy_registry();

/// @brief converts an operator label to its OpClass using the default MBPT
/// context
/// @param op the operator label
/// @return the OpClass of the operator
OpClass to_op_class(const std::wstring& op);

}  // namespace sequant::mbpt

#endif  // SEQUANT_DOMAIN_MBPT_CONTEXT_HPP
