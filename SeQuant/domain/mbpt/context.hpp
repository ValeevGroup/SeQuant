
#ifndef SEQUANT_DOMAIN_MBPT_CONTEXT_HPP
#define SEQUANT_DOMAIN_MBPT_CONTEXT_HPP

#include "SeQuant/core/utility/context.hpp"

namespace sequant::mbpt {

/// @brief details of the MBPT formalism
///
/// to be used as a implicit or explicit specification of the formalism
class Context {
 public:
  /// Whether to employ bra/ket-(anti)symmetrized tensors in definition of
  /// N-body operators
  enum class NBodyInteractionTensorSymm { Yes, No };

  /// Whether to use cluster-specific virtuals.
  enum class CSV { Yes, No };

  struct Defaults {
    constexpr static auto nbody_interaction_tensor_symm =
        NBodyInteractionTensorSymm::Yes;
    constexpr static auto csv = CSV::No;
  };

  /// @brief Construct a Formalism object
  explicit Context(NBodyInteractionTensorSymm nbody_interaction_tensor_symm =
                       Defaults::nbody_interaction_tensor_symm,
                   CSV csv_formalism = Defaults::csv) noexcept;

  /// @brief Construct a Formalism object
  explicit Context(CSV csv_formalism,
                   NBodyInteractionTensorSymm nbody_interaction_tensor_symm =
                       Defaults::nbody_interaction_tensor_symm) noexcept;

  NBodyInteractionTensorSymm nbody_interaction_tensor_symm() const {
    return nbody_interaction_tensor_symm_;
  }
  CSV csv() const { return csv_; }

 private:
  NBodyInteractionTensorSymm nbody_interaction_tensor_symm_ =
      Defaults::nbody_interaction_tensor_symm;
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
detail::ImplicitContextResetter<Context> set_scoped_default_formalism(
    const Context& ctx);

}  // namespace sequant::mbpt

#endif  // SEQUANT_DOMAIN_MBPT_CONTEXT_HPP
