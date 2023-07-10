#include "formalism.hpp"

namespace sequant::mbpt {

Formalism::Formalism(NBodyInteractionTensorSymm nbody_interaction_tensor_symm,
                     CSV csv) noexcept
    : nbody_interaction_tensor_symm_(nbody_interaction_tensor_symm),
      csv_(csv) {}

Formalism::Formalism(
    CSV csv, NBodyInteractionTensorSymm nbody_interaction_tensor_symm) noexcept
    : Formalism(nbody_interaction_tensor_symm, csv) {}

bool operator==(Formalism const& left, Formalism const& right) {
  return left.nbody_interaction_tensor_symm() ==
             right.nbody_interaction_tensor_symm() &&
         left.csv() == right.csv();
}

bool operator!=(Formalism const& left, Formalism const& right) {
  return !(left == right);
}

namespace detail {

Formalism& default_formalism_instance() {
  static Formalism instance_;
  return instance_;
}

}  // namespace detail

const Formalism& get_default_formalism() {
  return detail::default_formalism_instance();
}

void set_default_formalism(const Formalism& ctx) {
  detail::default_formalism_instance() = ctx;
}

void reset_default_context() {
  detail::default_formalism_instance() = Formalism{};
}

detail::FormalismResetter::FormalismResetter(const Formalism& previous) noexcept
    : previous_{previous} {}

detail::FormalismResetter::~FormalismResetter() noexcept {
  if (previous_) set_default_formalism(*previous_);
}

detail::FormalismResetter set_scoped_default_formalism(const Formalism& f) {
  if (detail::default_formalism_instance() != f) {
    auto prev = detail::default_formalism_instance();
    detail::default_formalism_instance() = f;
    return prev;
  } else {
    return {};
  }
}

}  // namespace sequant::mbpt
