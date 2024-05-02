#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/utility/context.hpp>

namespace sequant {

bool operator==(const Context& ctx1, const Context& ctx2) {
  if (&ctx1 == &ctx2)
    return true;
  else
    return ctx1.vacuum() == ctx2.vacuum() && ctx1.metric() == ctx2.metric() &&
           ctx1.braket_symmetry() == ctx2.braket_symmetry() &&
           ctx1.spbasis() == ctx2.spbasis() &&
           ctx1.first_dummy_index_ordinal() ==
               ctx2.first_dummy_index_ordinal() &&
           *ctx1.index_space_registry() == *ctx2.index_space_registry();
}

bool operator!=(const Context& ctx1, const Context& ctx2) {
  return !(ctx1 == ctx2);
}

const Context& get_default_context() {
  return detail::get_implicit_context<Context>();
}

void set_default_context(const Context& ctx) {
  return detail::set_implicit_context(ctx);
}

void reset_default_context() { detail::reset_implicit_context<Context>(); }

[[nodiscard]] detail::ImplicitContextResetter<Context>
set_scoped_default_context(const Context& ctx) {
  return detail::set_scoped_implicit_context(ctx);
}

Context::Context(Vacuum vac, IndexSpaceRegistry isr, IndexSpaceMetric m,
                 BraKetSymmetry bks, SPBasis spb, std::size_t fdio)
    : vacuum_(vac),
      metric_(m),
      braket_symmetry_(bks),
      spbasis_(spb),
      first_dummy_index_ordinal_(fdio),
      idx_space_reg_(std::make_shared<IndexSpaceRegistry>(isr)) {}

Vacuum Context::vacuum() const { return vacuum_; }
std::shared_ptr<IndexSpaceRegistry> Context::index_space_registry() const {
  return idx_space_reg_;
}
IndexSpaceMetric Context::metric() const { return metric_; }
BraKetSymmetry Context::braket_symmetry() const { return braket_symmetry_; }
SPBasis Context::spbasis() const { return spbasis_; }
std::size_t Context::first_dummy_index_ordinal() const {
  return first_dummy_index_ordinal_;
}
Context& Context::set(Vacuum vacuum) {
  vacuum_ = vacuum;
  return *this;
}
Context& Context::set(IndexSpaceRegistry ISR) {
  idx_space_reg_ = std::make_shared<IndexSpaceRegistry>(ISR);
  return *this;
}
Context& Context::set(IndexSpaceMetric metric) {
  metric_ = metric;
  return *this;
}
Context& Context::set(BraKetSymmetry braket_symmetry) {
  braket_symmetry_ = braket_symmetry;
  return *this;
}
Context& Context::set(SPBasis spbasis) {
  spbasis_ = spbasis;
  return *this;
}

Context& Context::set_first_dummy_index_ordinal(
    std::size_t first_dummy_index_ordinal) {
  first_dummy_index_ordinal_ = first_dummy_index_ordinal;
  return *this;
}

}  // namespace sequant

#ifdef SEQUANT_HAS_MIMALLOC
#include <mimalloc-new-delete.h>
#endif
