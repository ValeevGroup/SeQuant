
#ifndef SEQUANT_FORMALISM_HPP
#define SEQUANT_FORMALISM_HPP

#include <optional>

namespace sequant::mbpt {

/// @brief the MBPT formalism details
class Formalism {
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
  explicit Formalism(NBodyInteractionTensorSymm nbody_interaction_tensor_symm =
                         Defaults::nbody_interaction_tensor_symm,
                     CSV csv_formalism = Defaults::csv) noexcept;

  /// @brief Construct a Formalism object
  explicit Formalism(CSV csv_formalism,
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

bool operator==(Formalism const& left, Formalism const& right);

bool operator!=(Formalism const& left, Formalism const& right);

const Formalism& get_default_formalism();
void set_default_formalism(const Formalism& ctx);
void reset_default_formalism();

namespace detail {
struct FormalismResetter {
  FormalismResetter() = default;
  FormalismResetter(const Formalism& previous) noexcept;
  ~FormalismResetter() noexcept;

  // FormalismResetter is move-only
  FormalismResetter(const FormalismResetter&) = delete;
  FormalismResetter& operator=(const FormalismResetter&) = delete;

 private:
  std::optional<Formalism> previous_;
};
}  // namespace detail

detail::FormalismResetter set_scoped_default_formalism(const Formalism& ctx);

}  // namespace sequant::mbpt

#endif  // SEQUANT_FORMALISM_HPP
