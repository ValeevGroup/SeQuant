#include <SeQuant/core/algorithm.hpp>
#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/reserved.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/core/utility/context.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <atomic>
#include <cstdint>
#include <utility>

#ifdef SEQUANT_CONTEXT_MANIPULATION_THREADSAFE
#include <mutex>
#endif

namespace sequant {

namespace {

// process-wide immutable defaults, shared so that default-constructed
// contexts compare equal
const std::shared_ptr<TensorCanonicalizer>& default_tensor_canonicalizer() {
  static const std::shared_ptr<TensorCanonicalizer> result =
      std::make_shared<DefaultTensorCanonicalizer>();
  return result;
}

const std::shared_ptr<const tensor_index_comparer_t>& default_index_comparer() {
  static const std::shared_ptr<const tensor_index_comparer_t> result =
      std::make_shared<const tensor_index_comparer_t>(
          TensorCanonicalizer::default_index_comparer());
  return result;
}

const std::shared_ptr<const tensor_index_pair_comparer_t>&
default_index_pair_comparer() {
  static const std::shared_ptr<const tensor_index_pair_comparer_t> result =
      std::make_shared<const tensor_index_pair_comparer_t>(
          TensorCanonicalizer::default_index_pair_comparer());
  return result;
}

std::atomic<std::uint64_t> last_context_version{0};

void check_tensor_canonicalizer(
    const std::shared_ptr<TensorCanonicalizer>& canonicalizer) {
  if (!canonicalizer)
    throw Exception("Context: a tensor canonicalizer must not be null");
}

template <typename Comparer>
std::shared_ptr<const Comparer> checked_comparer(
    std::shared_ptr<const Comparer> comparer) {
  SEQUANT_ASSERT(comparer && *comparer);
  return comparer;
}

template <typename Comparer>
std::shared_ptr<const Comparer> make_comparer(Comparer comparer) {
  return checked_comparer(
      std::make_shared<const Comparer>(std::move(comparer)));
}

void check_cardinal_tensor_labels(
    [[maybe_unused]] const container::vector<std::wstring>& labels) {
  SEQUANT_ASSERT(!has_duplicates(labels) &&
                 "cardinal tensor labels must not contain duplicates");
}

}  // namespace

bool default_context_manipulation_threadsafe() {
#ifdef SEQUANT_CONTEXT_MANIPULATION_THREADSAFE
  return true;
#else
  return false;
#endif
}

bool operator==(const Context& ctx1, const Context& ctx2) {
  // a Context need not have a registry
  auto same_registry = [](const auto& r1, const auto& r2) {
    if (!r1 || !r2) return !r1 && !r2;
    return r1->spaces() == r2->spaces() && *r1 == *r2;
  };
  if (&ctx1 == &ctx2)
    return true;
  else
    return ctx1.vacuum() == ctx2.vacuum() && ctx1.metric() == ctx2.metric() &&
           ctx1.assert_strict_braket_symmetry() ==
               ctx2.assert_strict_braket_symmetry() &&
           ctx1.spbasis() == ctx2.spbasis() &&
           ctx1.first_dummy_index_ordinal() ==
               ctx2.first_dummy_index_ordinal() &&
           ctx1.canonicalization_options() == ctx2.canonicalization_options() &&
           ctx1.braket_typesetting() == ctx2.braket_typesetting() &&
           ctx1.braket_slot_typesetting() == ctx2.braket_slot_typesetting() &&
           ctx1.deserialization_symmetry() == ctx2.deserialization_symmetry() &&
           ctx1.deserialization_hermiticity() ==
               ctx2.deserialization_hermiticity() &&
           ctx1.deserialization_column_symmetry() ==
               ctx2.deserialization_column_symmetry() &&
           *ctx1.tensor_canonicalizers_ == *ctx2.tensor_canonicalizers_ &&
           same_registry(ctx1.index_space_registry(),
                         ctx2.index_space_registry());
}

bool operator!=(const Context& ctx1, const Context& ctx2) {
  return !(ctx1 == ctx2);
}

#ifdef SEQUANT_CONTEXT_MANIPULATION_THREADSAFE
static std::recursive_mutex ctx_mtx;  // used to protect the context
// immutable copy of the process-wide contexts for
// get_default_context_snapshot(), and its generation; both replaced, under
// ctx_mtx, by every change of the process-wide contexts
static std::shared_ptr<const container::map<Statistics, Context>>
    published_default_contexts;
static std::atomic<std::uint64_t> default_contexts_generation{1};

static void publish_default_contexts() {
  published_default_contexts = std::make_shared<
      const container::map<Statistics, Context>>(
      detail::implicit_context_instance<container::map<Statistics, Context>>());
  default_contexts_generation.fetch_add(1, std::memory_order_release);
}
#endif

const Context& get_default_context(Statistics s) {
  // a scoped context is thread-local, hence needs no lock
  if (const auto* overlay = detail::implicit_context_overlay<
          container::map<Statistics, Context>>()) {
    auto it = overlay->find(s);
    if (it == overlay->end()) it = overlay->find(Statistics::Arbitrary);
    SEQUANT_ASSERT(it != overlay->end());
    return it->second;
  }
#ifdef SEQUANT_CONTEXT_MANIPULATION_THREADSAFE
  std::scoped_lock lock(ctx_mtx);
#endif
  auto& contexts =
      detail::implicit_context_instance<container::map<Statistics, Context>>();
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

Context get_default_context_snapshot(Statistics s) {
#ifdef SEQUANT_CONTEXT_MANIPULATION_THREADSAFE
  // a scoped context is thread-local, hence needs no lock
  if (detail::implicit_context_overlay<container::map<Statistics, Context>>())
    return get_default_context(s);
  // snapshots are taken concurrently on hot paths, where a lock per read
  // serializes the threads, so this thread holds the published contexts and
  // takes the lock only to fetch them anew after they changed
  struct Cache {
    std::uint64_t generation = 0;
    std::shared_ptr<const container::map<Statistics, Context>> contexts;
  };
  thread_local Cache cache;
  if (cache.generation !=
      default_contexts_generation.load(std::memory_order_acquire)) {
    std::scoped_lock lock(ctx_mtx);
    get_default_context();  // ensures that the arbitrary statistics has one
    cache.contexts = published_default_contexts;
    cache.generation =
        default_contexts_generation.load(std::memory_order_relaxed);
  }
  auto it = cache.contexts->find(s);
  if (it == cache.contexts->end())
    it = cache.contexts->find(Statistics::Arbitrary);
  SEQUANT_ASSERT(it != cache.contexts->end());
  return it->second;
#else
  return get_default_context(s);
#endif
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
#ifdef SEQUANT_CONTEXT_MANIPULATION_THREADSAFE
  publish_default_contexts();
#endif
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
#ifdef SEQUANT_CONTEXT_MANIPULATION_THREADSAFE
  publish_default_contexts();
#endif
}

[[nodiscard]] detail::ImplicitContextResetter<
    container::map<Statistics, Context>>
set_scoped_default_context(container::map<Statistics, Context> ctx) {
  // a scoped context is a thread-local overlay that leaves the process-wide
  // contexts alone, hence needs no lock
  // get_default_context() falls back to the context for arbitrary statistics
  ctx.try_emplace(Statistics::Arbitrary);
  return detail::set_scoped_implicit_context(std::move(ctx));
}

[[nodiscard]] detail::ImplicitContextResetter<
    container::map<Statistics, Context>>
set_scoped_default_context(Context ctx) {
  return set_scoped_default_context(container::map<Statistics, Context>{
      {Statistics::Arbitrary, std::move(ctx)}});
}

[[nodiscard]] detail::ImplicitContextResetter<
    container::map<Statistics, Context>>
set_scoped_default_context(Context::Options ctx_options) {
  return set_scoped_default_context(Context(std::move(ctx_options)));
}

[[nodiscard]] detail::ImplicitContextResetter<
    container::map<Statistics, Context>>
set_scoped_modified_default_context(
    const std::function<void(Context&)>& modify) {
  auto ctxs = [] {
    if (const auto* overlay = detail::implicit_context_overlay<
            container::map<Statistics, Context>>())
      return *overlay;
    get_default_context();  // ensures that the arbitrary statistics has one
#ifdef SEQUANT_CONTEXT_MANIPULATION_THREADSAFE
    std::scoped_lock lock(ctx_mtx);
#endif
    return detail::implicit_context_instance<
        container::map<Statistics, Context>>();
  }();
  for (auto& [s, ctx] : ctxs) modify(ctx);
  return set_scoped_default_context(std::move(ctxs));
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
      assert_strict_braket_symmetry_(options.assert_strict_braket_symmetry),
      spbasis_(options.spbasis),
      first_dummy_index_ordinal_(options.first_dummy_index_ordinal),
      canonicalization_options_(options.canonicalization_options),
      braket_typesetting_(options.braket_typesetting),
      braket_slot_typesetting_(options.braket_slot_typesetting),
      deserialization_symmetry_(options.deserialization_symmetry),
      deserialization_hermiticity_(options.deserialization_hermiticity),
      deserialization_column_symmetry_(
          options.deserialization_column_symmetry) {
  auto tensor_canonicalizers = std::make_shared<TensorCanonicalizers>();
  if (options.tensor_canonicalizers) {
    for (const auto& [label, canonicalizer] : *options.tensor_canonicalizers)
      check_tensor_canonicalizer(canonicalizer);
    tensor_canonicalizers->map = std::move(*options.tensor_canonicalizers);
  }
  tensor_canonicalizers->map.try_emplace(L"", default_tensor_canonicalizer());
  tensor_canonicalizers->index_comparer =
      options.index_comparer ? make_comparer(std::move(*options.index_comparer))
                             : default_index_comparer();
  tensor_canonicalizers->index_pair_comparer =
      options.index_pair_comparer
          ? make_comparer(std::move(*options.index_pair_comparer))
          : default_index_pair_comparer();
  if (options.cardinal_tensor_labels) {
    check_cardinal_tensor_labels(*options.cardinal_tensor_labels);
    tensor_canonicalizers->cardinal_labels =
        std::move(*options.cardinal_tensor_labels);
  } else
    tensor_canonicalizers->cardinal_labels = {reserved::antisymm_label(),
                                              reserved::symm_label(),
                                              reserved::transposition_label()};
  tensor_canonicalizers_ = std::move(tensor_canonicalizers);
  bump_version();
}

Context Context::clone() const {
  Context ctx(*this);
  ctx.idx_space_reg_ =
      std::make_shared<IndexSpaceRegistry>(idx_space_reg_->clone());
  ctx.bump_version();
  return ctx;
}

std::uint64_t Context::version() const { return version_; }

Context::TensorCanonicalizers& Context::mutable_tensor_canonicalizers() {
  auto copy = std::make_shared<TensorCanonicalizers>(*tensor_canonicalizers_);
  auto& result = *copy;
  tensor_canonicalizers_ = std::move(copy);
  return result;
}

void Context::bump_version() {
  version_ = last_context_version.fetch_add(1, std::memory_order_relaxed) + 1;
}

std::uint64_t current_context_version(Statistics s) {
  return get_default_context(s).version();
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

bool Context::assert_strict_braket_symmetry() const {
  return assert_strict_braket_symmetry_;
}

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

Symmetry Context::deserialization_symmetry() const {
  return deserialization_symmetry_;
}

Hermiticity Context::deserialization_hermiticity() const {
  return deserialization_hermiticity_;
}

ColumnSymmetry Context::deserialization_column_symmetry() const {
  return deserialization_column_symmetry_;
}

std::shared_ptr<TensorCanonicalizer> Context::tensor_canonicalizer_ptr(
    std::wstring_view label) const {
  auto result = nondefault_tensor_canonicalizer_ptr(label);
  if (!result) result = nondefault_tensor_canonicalizer_ptr(L"");
  return result;
}

std::shared_ptr<TensorCanonicalizer>
Context::nondefault_tensor_canonicalizer_ptr(std::wstring_view label) const {
  const auto& map = tensor_canonicalizers_->map;
  auto it = map.find(std::wstring{label});
  return it != map.end() ? it->second : nullptr;
}

const TensorCanonicalizer& Context::tensor_canonicalizer(
    std::wstring_view label) const {
  auto ptr = tensor_canonicalizer_ptr(label);
  if (!ptr)
    throw Exception(
        "Context::tensor_canonicalizer: no canonicalizer for this label nor "
        "for the empty label");
  // the map entry keeps *ptr alive
  return *ptr;
}

const tensor_index_comparer_t& Context::index_comparer() const {
  return *tensor_canonicalizers_->index_comparer;
}

const tensor_index_pair_comparer_t& Context::index_pair_comparer() const {
  return *tensor_canonicalizers_->index_pair_comparer;
}

std::shared_ptr<const tensor_index_comparer_t> Context::index_comparer_ptr()
    const {
  return tensor_canonicalizers_->index_comparer;
}

std::shared_ptr<const tensor_index_pair_comparer_t>
Context::index_pair_comparer_ptr() const {
  return tensor_canonicalizers_->index_pair_comparer;
}

const container::vector<std::wstring>& Context::cardinal_tensor_labels() const {
  return tensor_canonicalizers_->cardinal_labels;
}

Context& Context::set(Vacuum vacuum) {
  vacuum_ = vacuum;
  bump_version();
  return *this;
}

Context& Context::set(IndexSpaceRegistry ISR) {
  idx_space_reg_ = std::make_shared<IndexSpaceRegistry>(ISR);
  bump_version();
  return *this;
}

Context& Context::set(std::shared_ptr<IndexSpaceRegistry> ISR) {
  idx_space_reg_ = std::move(ISR);
  bump_version();
  return *this;
}

Context& Context::set(IndexSpaceMetric metric) {
  metric_ = metric;
  bump_version();
  return *this;
}

Context& Context::set(
    AssertStrictBraKetSymmetry assert_strict_braket_symmetry) {
  assert_strict_braket_symmetry_ =
      assert_strict_braket_symmetry == AssertStrictBraKetSymmetry::Yes;
  bump_version();
  return *this;
}

Context& Context::set(SPBasis spbasis) {
  spbasis_ = spbasis;
  bump_version();
  return *this;
}

Context& Context::set_first_dummy_index_ordinal(
    std::size_t first_dummy_index_ordinal) {
  first_dummy_index_ordinal_ = first_dummy_index_ordinal;
  bump_version();
  return *this;
}

Context& Context::set(CanonicalizeOptions copt) {
  canonicalization_options_ = copt;
  bump_version();
  return *this;
}

Context& Context::set(BraKetTypesetting bkt) {
  braket_typesetting_ = bkt;
  bump_version();
  return *this;
}

Context& Context::set(BraKetSlotTypesetting bkst) {
  braket_slot_typesetting_ = bkst;
  bump_version();
  return *this;
}

Context& Context::set(Symmetry symmetry) {
  deserialization_symmetry_ = symmetry;
  bump_version();
  return *this;
}

Context& Context::set(Hermiticity hermiticity) {
  deserialization_hermiticity_ = hermiticity;
  bump_version();
  return *this;
}

Context& Context::set(ColumnSymmetry column_symmetry) {
  deserialization_column_symmetry_ = column_symmetry;
  bump_version();
  return *this;
}

Context& Context::set_tensor_canonicalizer(
    std::wstring_view label,
    std::shared_ptr<TensorCanonicalizer> canonicalizer) {
  check_tensor_canonicalizer(canonicalizer);
  mutable_tensor_canonicalizers().map.insert_or_assign(
      std::wstring{label}, std::move(canonicalizer));
  bump_version();
  return *this;
}

Context& Context::unset_tensor_canonicalizer(std::wstring_view label) {
  mutable_tensor_canonicalizers().map.erase(std::wstring{label});
  bump_version();
  return *this;
}

Context& Context::set_index_comparer(tensor_index_comparer_t comparer) {
  mutable_tensor_canonicalizers().index_comparer =
      make_comparer(std::move(comparer));
  bump_version();
  return *this;
}

Context& Context::set_index_comparer(
    std::shared_ptr<const tensor_index_comparer_t> comparer) {
  mutable_tensor_canonicalizers().index_comparer =
      checked_comparer(std::move(comparer));
  bump_version();
  return *this;
}

Context& Context::set_index_pair_comparer(
    tensor_index_pair_comparer_t comparer) {
  mutable_tensor_canonicalizers().index_pair_comparer =
      make_comparer(std::move(comparer));
  bump_version();
  return *this;
}

Context& Context::set_index_pair_comparer(
    std::shared_ptr<const tensor_index_pair_comparer_t> comparer) {
  mutable_tensor_canonicalizers().index_pair_comparer =
      checked_comparer(std::move(comparer));
  bump_version();
  return *this;
}

Context& Context::set_cardinal_tensor_labels(
    container::vector<std::wstring> labels) {
  check_cardinal_tensor_labels(labels);
  mutable_tensor_canonicalizers().cardinal_labels = std::move(labels);
  bump_version();
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
