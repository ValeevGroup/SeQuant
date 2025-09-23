#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/utility/context.hpp>

#ifdef SEQUANT_CONTEXT_MANIPULATION_THREADSAFE
#include <mutex>
#endif

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
           ctx1.index_space_registry()->spaces() ==
               ctx2.index_space_registry()->spaces() &&
           ctx1.canonicalization_options() == ctx2.canonicalization_options() &&
           ctx1.braket_typesetting() == ctx2.braket_typesetting() &&
           ctx1.braket_slot_typesetting() == ctx2.braket_slot_typesetting() &&
           *ctx1.index_space_registry() == *ctx2.index_space_registry();
}

bool operator!=(const Context& ctx1, const Context& ctx2) {
  return !(ctx1 == ctx2);
}

#ifdef SEQUANT_CONTEXT_MANIPULATION_THREADSAFE
static std::recursive_mutex ctx_mtx;  // used to protect the context
#endif

const Context& get_default_context(Statistics s) {
#ifdef SEQUANT_CONTEXT_MANIPULATION_THREADSAFE
  std::scoped_lock lock(ctx_mtx);
#endif
  auto& contexts =
      detail::get_implicit_context<container::map<Statistics, Context>>();
  auto it = contexts.find(s);
  /// default for arbitrary statistics is initialized lazily here
  if (it == contexts.end() && s == Statistics::Arbitrary) {
    set_default_context(Context{}, Statistics::Arbitrary);
  }
  it = contexts.find(s);
  // have context for this statistics? else return for arbitrary statistics
  if (it != contexts.end())
    return it->second;
  else
    return get_default_context(Statistics::Arbitrary);
}

void set_default_context(Context ctx, Statistics s) {
#ifdef SEQUANT_CONTEXT_MANIPULATION_THREADSAFE
  std::scoped_lock lock(ctx_mtx);
#endif
  auto& contexts =
      detail::implicit_context_instance<container::map<Statistics, Context>>();
  auto it = contexts.find(s);
  if (it != contexts.end()) {
    it->second = std::move(ctx);
  } else {
    contexts.emplace(s, std::move(ctx));
  }
}

void set_default_context(Context::Options ctx_opts, Statistics s) {
  return set_default_context(Context(ctx_opts), s);
}

void set_default_context(const container::map<Statistics, Context>& ctxs) {
  for (const auto& [s, ctx] : ctxs) {
    set_default_context(ctx, s);
  }
}

void reset_default_context() {
#ifdef SEQUANT_CONTEXT_MANIPULATION_THREADSAFE
  std::scoped_lock lock(ctx_mtx);
#endif
  detail::reset_implicit_context<container::map<Statistics, Context>>();
}

[[nodiscard]] detail::ImplicitContextResetter<
    container::map<Statistics, Context>>
set_scoped_default_context(const container::map<Statistics, Context>& ctx) {
#ifdef SEQUANT_CONTEXT_MANIPULATION_THREADSAFE
  std::scoped_lock lock(ctx_mtx);
#endif
  return detail::set_scoped_implicit_context(ctx);
}

[[nodiscard]] detail::ImplicitContextResetter<
    container::map<Statistics, Context>>
set_scoped_default_context(Context ctx) {
  return detail::set_scoped_implicit_context(
      container::map<Statistics, Context>{
          {Statistics::Arbitrary, std::move(ctx)}});
}

[[nodiscard]] detail::ImplicitContextResetter<
    container::map<Statistics, Context>>
set_scoped_default_context(Context::Options ctx_options) {
  return detail::set_scoped_implicit_context(
      container::map<Statistics, Context>{
          {Statistics::Arbitrary, Context(std::move(ctx_options))}});
}

Context::Context(Options options)
    : idx_space_reg_(
          options.index_space_registry_shared_ptr
              ? std::move(options.index_space_registry_shared_ptr)
              : (options.index_space_registry.has_value()
                     ? std::make_shared<IndexSpaceRegistry>(
                           std::move(options.index_space_registry.value()))
                     : nullptr)),
      vacuum_(options.vacuum),
      metric_(options.metric),
      braket_symmetry_(options.braket_symmetry),
      spbasis_(options.spbasis),
      first_dummy_index_ordinal_(options.first_dummy_index_ordinal),
      canonicalization_options_(options.canonicalization_options),
      braket_typesetting_(options.braket_typesetting),
      braket_slot_typesetting_(options.braket_slot_typesetting) {}

Context Context::clone() const {
  Context ctx(*this);
  ctx.idx_space_reg_ =
      std::make_shared<IndexSpaceRegistry>(idx_space_reg_->clone());
  return ctx;
}

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
std::optional<CanonicalizeOptions> Context::canonicalization_options() const {
  return canonicalization_options_;
}

BraKetTypesetting Context::braket_typesetting() const {
  return braket_typesetting_;
}

BraKetSlotTypesetting Context::braket_slot_typesetting() const {
  return braket_slot_typesetting_;
}

Context& Context::set(Vacuum vacuum) {
  vacuum_ = vacuum;
  return *this;
}

Context& Context::set(IndexSpaceRegistry ISR) {
  idx_space_reg_ = std::make_shared<IndexSpaceRegistry>(ISR);
  return *this;
}

Context& Context::set(std::shared_ptr<IndexSpaceRegistry> ISR) {
  idx_space_reg_ = std::move(ISR);
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

Context& Context::set(CanonicalizeOptions copt) {
  canonicalization_options_ = copt;
  return *this;
}

Context& Context::set(BraKetTypesetting bkt) {
  braket_typesetting_ = bkt;
  return *this;
}

Context& Context::set(BraKetSlotTypesetting bkst) {
  braket_slot_typesetting_ = bkst;
  return *this;
}

IndexSpace get_particle_space(const IndexSpace::QuantumNumbers& qn) {
  return get_default_context().index_space_registry()->particle_space(qn);
}

IndexSpace get_hole_space(const IndexSpace::QuantumNumbers& qn) {
  return get_default_context().index_space_registry()->hole_space(qn);
}

IndexSpace get_complete_space(const IndexSpace::QuantumNumbers& qn) {
  return get_default_context().index_space_registry()->complete_space(qn);
}

}  // namespace sequant

#ifdef SEQUANT_HAS_MIMALLOC
#include <mimalloc-new-delete.h>
#endif
