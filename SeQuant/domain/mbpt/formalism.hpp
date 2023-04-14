
#ifndef SEQUANT_FORMALISM_HPP
#define SEQUANT_FORMALISM_HPP

#include <optional>

namespace sequant::mbpt {

/// Whether to employ antisymmetric operators while considering two-body
/// interactions.
enum class TwoBodyInteraction { Antisymm, Nonsymm };

/// Whether to sum over complete unoccupied orbitals.
enum class SumOverUocc { Complete, Truncated };

/// Whether to use cluster-specific virtual formalism.
enum class CSVFormalism { CSV, NonCSV };

class Formalism {
 public:
  TwoBodyInteraction two_body_interaction() const;

  SumOverUocc sum_over_uocc() const;

  CSVFormalism csv_formalism() const;

  Formalism& set(TwoBodyInteraction);

  Formalism& set(SumOverUocc);

  Formalism& set(CSVFormalism);

  static Formalism make_default();

 private:
  struct Defaults {
    constexpr static auto two_body_interaction_ = TwoBodyInteraction::Antisymm;
    constexpr static auto sum_over_uocc_ = SumOverUocc::Truncated;
    constexpr static auto csv_formalism_ = CSVFormalism::NonCSV;
  };

  TwoBodyInteraction two_body_interaction_;

  SumOverUocc sum_over_uocc_;

  CSVFormalism csv_formalism_;
};

bool operator==(Formalism const& left, Formalism const& right);

bool operator!=(Formalism const& left, Formalism const& right);

namespace detail {

struct FormalismResetter {
  FormalismResetter() = default;

  FormalismResetter(Formalism const& previous) noexcept;

  ~FormalismResetter() noexcept;

  FormalismResetter(FormalismResetter const&) = delete;

  FormalismResetter& operator=(FormalismResetter const&) = delete;

 private:
  std::optional<Formalism> previous_;
};

}  // namespace detail

/// @brief Loads defaults for Formalism @formalism

/// This sets up two-body interaction (antisymmetric or non-symmetric),
/// whether to sum over complete-unoccupieds (complete) or active unoccupieds
/// (truncated), and whether to use cluster-specific-virtuals (eg. PNO)
/// formalism.
void set_default_formalism(Formalism formalism = Formalism::make_default());

/// @brief Get the current Formalism instance in effect.
Formalism const& get_default_formalism();

void reset_default_formalism();

detail::FormalismResetter set_scoped_default_formalism(Formalism const& f);

}  // namespace sequant::mbpt

#endif  // SEQUANT_FORMALISM_HPP
