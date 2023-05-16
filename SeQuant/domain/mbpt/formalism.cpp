#include "formalism.hpp"

namespace sequant::mbpt {

TwoBodyInteraction Formalism::two_body_interaction() const {
  return two_body_interaction_;
}

SumOverUocc Formalism::sum_over_uocc() const { return sum_over_uocc_; }

CSVFormalism Formalism::csv_formalism() const { return csv_formalism_; }

Formalism& Formalism::set(TwoBodyInteraction tbi) {
  two_body_interaction_ = tbi;
  return *this;
}

Formalism& Formalism::set(SumOverUocc sou) {
  sum_over_uocc_ = sou;
  return *this;
}

Formalism& Formalism::set(CSVFormalism csvf) {
  csv_formalism_ = csvf;
  return *this;
}

Formalism Formalism::make_default() {
  auto f = Formalism{};

  f.two_body_interaction_ = Defaults::two_body_interaction_;
  f.sum_over_uocc_ = Defaults::sum_over_uocc_;
  f.csv_formalism_ = Defaults::csv_formalism_;

  return f;
}

bool operator==(Formalism const& left, Formalism const& right) {
  return left.two_body_interaction() == right.two_body_interaction() &&
         left.sum_over_uocc() == right.sum_over_uocc() &&
         left.csv_formalism() == right.csv_formalism();
}

bool operator!=(Formalism const& left, Formalism const& right) {
  return !(left == right);
}

Formalism& default_formalism_instance() {
  static Formalism instance_ = Formalism::make_default();
  return instance_;
}

void set_default_formalism(Formalism formalism) {
  default_formalism_instance() = formalism;
}

void reset_default_formalism() {
  default_formalism_instance() = Formalism::make_default();
}

Formalism const& get_default_formalism() {
  return default_formalism_instance();
}
detail::FormalismResetter::FormalismResetter(const Formalism& previous) noexcept
    : previous_{previous} {}

detail::FormalismResetter::~FormalismResetter() noexcept {
  if (previous_) set_default_formalism(*previous_);
}

detail::FormalismResetter set_scoped_default_formalism(const Formalism& f) {
  if (default_formalism_instance() != f) {
    auto prev = default_formalism_instance();
    default_formalism_instance() = f;
    return prev;
  } else {
    return {};
  }
}

}  // namespace sequant::mbpt
