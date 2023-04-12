#include "sequant.hpp"

namespace sequant {

bool operator==(const SeQuant& ctx1, const SeQuant& ctx2) {
  if (&ctx1 == &ctx2)
    return true;
  else
    return ctx1.vacuum() == ctx2.vacuum() && ctx1.metric() == ctx2.metric() &&
           ctx1.braket_symmetry() == ctx2.braket_symmetry() &&
           ctx1.spbasis() == ctx2.spbasis() &&
           ctx1.first_dummy_index_ordinal() == ctx2.first_dummy_index_ordinal();
}

bool operator!=(const SeQuant& ctx1, const SeQuant& ctx2) {
  return !(ctx1 == ctx2);
}

namespace detail {

SeQuant& default_context_instance() {
  static SeQuant instance_;
  return instance_;
}

}  // namespace detail

const SeQuant& get_default_context() {
  return detail::default_context_instance();
}

void set_default_context(const SeQuant& ctx) {
  detail::default_context_instance() = ctx;
}

void reset_default_context() { detail::default_context_instance() = SeQuant{}; }

detail::ContextResetter::ContextResetter(const SeQuant& previous_ctx) noexcept
    : previous_ctx_(previous_ctx) {}

detail::ContextResetter::~ContextResetter() noexcept {
  if (previous_ctx_) set_default_context(*previous_ctx_);
}

detail::ContextResetter set_scoped_default_context(const SeQuant& ctx) {
  if (detail::default_context_instance() != ctx) {
    auto previous_ctx = detail::default_context_instance();
    detail::default_context_instance() = ctx;
    return previous_ctx;
  } else
    return {};
}

SeQuant::SeQuant(Vacuum vac, IndexSpaceMetric m, BraKetSymmetry bks,
                 SPBasis spb, std::size_t fdio)
    : vacuum_(vac),
      metric_(m),
      braket_symmetry_(bks),
      spbasis_(spb),
      first_dummy_index_ordinal_(fdio) {}

Vacuum SeQuant::vacuum() const { return vacuum_; }
IndexSpaceMetric SeQuant::metric() const { return metric_; }
BraKetSymmetry SeQuant::braket_symmetry() const { return braket_symmetry_; }
SPBasis SeQuant::spbasis() const { return spbasis_; }
TwoBodyInteraction SeQuant::two_body_interaction() const {
  return two_body_interaction_;
}

SumOverUocc SeQuant::sum_over_uocc() const { return sum_over_uocc_; }

CSVFormalism SeQuant::csv_formalism() const { return csv_formalism_; }

std::size_t SeQuant::first_dummy_index_ordinal() const {
  return first_dummy_index_ordinal_;
}
SeQuant& SeQuant::set(Vacuum vacuum) {
  vacuum_ = vacuum;
  return *this;
}
SeQuant& SeQuant::set(IndexSpaceMetric metric) {
  metric_ = metric;
  return *this;
}
SeQuant& SeQuant::set(BraKetSymmetry braket_symmetry) {
  braket_symmetry_ = braket_symmetry;
  return *this;
}
SeQuant& SeQuant::set(SPBasis spbasis) {
  spbasis_ = spbasis;
  return *this;
}

SeQuant& SeQuant::set(SumOverUocc sou) {
  sum_over_uocc_ = sou;
  return *this;
}

SeQuant& SeQuant::set(CSVFormalism csvf) {
  csv_formalism_ = csvf;
  return *this;
}

SeQuant& SeQuant::set_first_dummy_index_ordinal(
    std::size_t first_dummy_index_ordinal) {
  first_dummy_index_ordinal_ = first_dummy_index_ordinal;
  return *this;
}

}  // namespace sequant

#ifdef SEQUANT_HAS_MIMALLOC
#include "mimalloc-new-delete.h"
#endif
