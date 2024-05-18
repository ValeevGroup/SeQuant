#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/utility/context.hpp>

#include <mutex>

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

static std::recursive_mutex ctx_mtx;  // used to protect the context

const Context& get_default_context(Statistics s) {
  std::scoped_lock lock(ctx_mtx);
  auto& contexts =
      detail::get_implicit_context<std::map<Statistics, Context>>();
  auto it = contexts.find(s);
  /// default for arbitrary statistics is initialized lazily here
  if (it == contexts.end() && s == Statistics::Arbitrary) {
    set_default_context({}, Statistics::Arbitrary);
  }
  it = contexts.find(s);
  // have context for this statistics? else return for arbitrary statistics
  if (it != contexts.end())
    return it->second;
  else
    return get_default_context(Statistics::Arbitrary);
}

void set_default_context(const Context& ctx, Statistics s) {
  std::scoped_lock lock(ctx_mtx);
  auto& contexts =
      detail::implicit_context_instance<std::map<Statistics, Context>>();
  auto it = contexts.find(s);
  if (it != contexts.end()) {
    it->second = ctx;
  } else {
    contexts.emplace(s, ctx);
  }
}

void set_default_context(const std::map<Statistics, Context>& ctxs) {
  for (const auto& [s, ctx] : ctxs) {
    set_default_context(ctx, s);
  }
}

void reset_default_context() {
  std::scoped_lock lock(ctx_mtx);
  detail::reset_implicit_context<std::map<Statistics, Context>>();
}

[[nodiscard]] detail::ImplicitContextResetter<std::map<Statistics, Context>>
set_scoped_default_context(const std::map<Statistics, Context>& ctx) {
  std::scoped_lock lock(ctx_mtx);
  return detail::set_scoped_implicit_context(ctx);
}

[[nodiscard]] detail::ImplicitContextResetter<std::map<Statistics, Context>>
set_scoped_default_context(const Context& ctx) {
  return detail::set_scoped_implicit_context(
      std::map<Statistics, Context>{{Statistics::Arbitrary, ctx}});
}

Context::Context(std::shared_ptr<IndexSpaceRegistry> isr, Vacuum vac,
                 IndexSpaceMetric m, BraKetSymmetry bks, SPBasis spb,
                 std::size_t fdio)
    : idx_space_reg_(std::move(isr)),
      vacuum_(vac),
      metric_(m),
      braket_symmetry_(bks),
      spbasis_(spb),
      first_dummy_index_ordinal_(fdio) {}

Context::Context(IndexSpaceRegistry isr, Vacuum vac, IndexSpaceMetric m,
                 BraKetSymmetry bks, SPBasis spb, std::size_t fdio)
    : idx_space_reg_(std::make_shared<IndexSpaceRegistry>(std::move(isr))),
      vacuum_(vac),
      metric_(m),
      braket_symmetry_(bks),
      spbasis_(spb),
      first_dummy_index_ordinal_(fdio) {}

Context::Context(Vacuum vac, IndexSpaceMetric m, BraKetSymmetry bks,
                 SPBasis spb, std::size_t fdio)
    : idx_space_reg_{},
      vacuum_(vac),
      metric_(m),
      braket_symmetry_(bks),
      spbasis_(spb),
      first_dummy_index_ordinal_(fdio) {}

Vacuum Context::vacuum() const { return vacuum_; }

std::shared_ptr<const IndexSpaceRegistry> Context::index_space_registry()
    const {
  return idx_space_reg_;
}

std::shared_ptr<IndexSpaceRegistry> Context::mutable_index_space_registry()
    const {
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
