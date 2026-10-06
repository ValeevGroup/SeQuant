
#ifndef SEQUANT_DOMAIN_MBPT_CONTEXT_HPP
#define SEQUANT_DOMAIN_MBPT_CONTEXT_HPP

#include <SeQuant/domain/mbpt/fwd.hpp>

#include <SeQuant/core/utility/aggregate.hpp>
#include <SeQuant/core/utility/context.hpp>
#include <SeQuant/domain/mbpt/op_registry.hpp>

namespace sequant::mbpt {

/// Whether to use cluster-specific virtuals.
enum class CSV { Yes, No };

// clang-format off
/// @brief Prefactors of coefficient-summed MBPT operators
///
/// For an operator with `c` creators and `a` annihilators let `m = c! a!` in a
/// spin-orbital basis and `m = c!` in a spin-free basis (where `c = a`). The
/// prefactor depends on the operator's OpClass:
///
/// | operators                                    | Default         | Symmetric |
/// |----------------------------------------------|-----------------|-----------|
/// | OpClass::Ex, OpClass::Deex (e.g. t, λ, R, L) | 1/m             | 1/sqrt(m) |
/// | OpClass::Gen (e.g. h, f, g, θ)               | 1/m             | 1/m       |
/// | projectors Â, Ŝ, P                           | None (Implicit) | 1/sqrt(m) |
///
/// For example, the doubles amplitude operator is
/// `1/4 t^{i1 i2}_{a1 a2} a^{a1 a2}_{i1 i2}` under Default and
/// `1/2 t^{i1 i2}_{a1 a2} a^{a1 a2}_{i1 i2}` under Symmetric.
///
/// The convention is read when an operator is lowered to tensor form.
// clang-format on
enum class NormalizationConvention {
  Default,   ///< factorial prefactors; projectors carry none
  Symmetric  ///< square-root prefactors for amplitudes and projectors
};

// clang-format off
/// @brief Specifies details of the MBPT formalism
///
/// MBPT context contains:
/// - csv: whether to use cluster-specific virtuals
/// - an `OpRegistry`: maps operator labels to their `OpClass`
/// - normalization_convention: prefactors of coefficient-summed operators
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
    constexpr static auto normalization_convention =
        NormalizationConvention::Default;
  };

  struct Options {
    SEQUANT_DESIGNATED_INIT_ONLY;
    /// whether to use cluster-specific virtuals
    CSV csv = Defaults::csv;
    /// shared pointer to operator registry
    std::shared_ptr<OpRegistry> op_registry_ptr = nullptr;
    /// optional operator registry, use if no shared pointer is provided
    std::optional<OpRegistry> op_registry = std::nullopt;
    /// normalization convention for tensor-form generation
    NormalizationConvention normalization_convention =
        Defaults::normalization_convention;
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

  /// @return the normalization convention for tensor-form generation
  NormalizationConvention normalization_convention() const;

  /// @return a constant pointer to the OpRegistry for this context
  /// @note asserts that OpRegistry is not null
  std::shared_ptr<const OpRegistry> op_registry() const;

  /// @return a pointer to the OpRegistry for this context
  /// @note asserts that OpRegistry is not null
  std::shared_ptr<OpRegistry> mutable_op_registry() const;

  /// @brief sets the OpRegistry for this context
  Context& set(const OpRegistry& op_registry);

  /// @brief sets the OpRegistry for this context
  Context& set(std::shared_ptr<OpRegistry> op_registry);

  /// @brief sets whether to use cluster-specific virtuals
  Context& set(CSV csv);

  /// sets the normalization convention for tensor-form generation
  Context& set(NormalizationConvention convention);

 private:
  CSV csv_ = Defaults::csv;
  std::shared_ptr<OpRegistry> op_registry_;
  NormalizationConvention normalization_convention_ =
      Defaults::normalization_convention;

  friend bool operator==(Context const& left, Context const& right);
};

bool operator!=(Context const& left, Context const& right);

const Context& get_default_mbpt_context();

void set_default_mbpt_context(const Context& ctx);

void set_default_mbpt_context(const Context::Options& options);

void reset_default_mbpt_context();

/// @brief changes the default mbpt context until the returned object is
/// destroyed
/// @note the scoped contexts are seen only by the calling thread and by the
/// workers of the parallel primitives (sequant::for_each, etc.) it launches;
/// set_default_mbpt_context() is not seen by this thread until the scope ends
/// @note scopes must end in the reverse order of their creation
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
/// @note returns OpClass::Gen for reserved operator labels
OpClass to_op_class(const std::wstring& op);

/// @brief returns the Hermiticity of an operator label using the default MBPT
/// context
/// @param op the operator label
/// @return the operator's Hermiticity (the registry's per-operator value, or
///         default_hermiticity(to_op_class(op)))
/// @note returns Hermiticity::Hermitian for reserved operator labels (they are
///       OpClass::Gen)
Hermiticity op_hermiticity(const std::wstring& op);

}  // namespace sequant::mbpt

#endif  // SEQUANT_DOMAIN_MBPT_CONTEXT_HPP
