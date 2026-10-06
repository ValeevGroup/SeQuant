#include <SeQuant/core/algorithm.hpp>
#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/options.hpp>
#include <SeQuant/core/reserved.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/core/utility/context.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <algorithm>
#include <cstdint>
#include <memory>
#include <mutex>
#include <optional>
#include <string>
#include <utility>
#include <vector>

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

/// identifies a shared object while it lives; the weak reference keeps an
/// object later allocated at the same address from matching
struct ObjectRef {
  const void* ptr = nullptr;
  std::weak_ptr<const void> obj;

  ObjectRef() = default;
  template <typename T>
  explicit ObjectRef(const std::shared_ptr<T>& p) : ptr(p.get()), obj(p) {}

  bool alive() const { return ptr == nullptr || !obj.expired(); }
  bool operator==(const ObjectRef& other) const { return ptr == other.ptr; }
};

/// what canonicalization reads from a Context: its CanonicalizationConfig, its
/// index space registry and its SP basis, with the shared objects referred to
/// weakly and compared by identity, as in
/// operator==(const Context&, const Context&)
struct CanonicalizationKey {
  ObjectRef registry;
  /// determines the symmetry of NormalOperator
  SPBasis spbasis = Context::Defaults::spbasis;
  container::vector<std::pair<std::wstring, ObjectRef>> canonicalizers;
  ObjectRef index_comparer;
  ObjectRef index_pair_comparer;
  container::vector<std::wstring> cardinal_labels;
  std::optional<CanonicalizeOptions> options;

  bool alive() const {
    return registry.alive() && index_comparer.alive() &&
           index_pair_comparer.alive() &&
           std::ranges::all_of(canonicalizers,
                               [](const auto& c) { return c.second.alive(); });
  }

  bool operator==(const CanonicalizationKey& other) const {
    // CanonicalizeOptions::operator== compares only the method
    auto same_options = [](const std::optional<CanonicalizeOptions>& o1,
                           const std::optional<CanonicalizeOptions>& o2) {
      if (!o1 || !o2) return !o1 && !o2;
      return o1->method == o2->method &&
             o1->named_indices == o2->named_indices &&
             o1->ignore_named_index_labels == o2->ignore_named_index_labels;
    };
    return registry == other.registry && spbasis == other.spbasis &&
           canonicalizers == other.canonicalizers &&
           index_comparer == other.index_comparer &&
           index_pair_comparer == other.index_pair_comparer &&
           cardinal_labels == other.cardinal_labels &&
           same_options(options, other.options);
  }
};

/// @return the version of @p key : equal for equal keys of live objects,
/// distinct otherwise, never reused
std::uint64_t canonicalization_version(CanonicalizationKey key) {
  static std::mutex mtx;
  static std::vector<std::pair<CanonicalizationKey, std::uint64_t>> versions;
  static std::uint64_t last_version = 0;
  std::scoped_lock lock(mtx);
  for (const auto& [k, v] : versions)
    if (k.alive() && k == key) return v;
  std::erase_if(versions, [](const auto& e) { return !e.first.alive(); });
  versions.emplace_back(std::move(key), ++last_version);
  return last_version;
}

/// @return @p registry if it is the only owner of its object, else a copy of
/// the object, so that a Context holds the only owners of its registry
std::shared_ptr<const IndexSpaceRegistry> owned(
    std::shared_ptr<const IndexSpaceRegistry> registry) {
  if (!registry || registry.use_count() == 1) return registry;
  return std::make_shared<const IndexSpaceRegistry>(*registry);
}

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
    return *r1 == *r2;
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
           ctx1.braket_typesetting() == ctx2.braket_typesetting() &&
           ctx1.braket_slot_typesetting() == ctx2.braket_slot_typesetting() &&
           ctx1.deserialization_symmetry() == ctx2.deserialization_symmetry() &&
           ctx1.deserialization_hermiticity() ==
               ctx2.deserialization_hermiticity() &&
           ctx1.deserialization_column_symmetry() ==
               ctx2.deserialization_column_symmetry() &&
           *ctx1.canonicalization_config_ == *ctx2.canonicalization_config_ &&
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
  return set_default_context(Context(std::move(ctx_opts)), s);
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
              ? owned(std::move(options.index_space_registry_shared_ptr))
              : (options.index_space_registry.has_value()
                     ? std::make_shared<const IndexSpaceRegistry>(
                           std::move(options.index_space_registry.value()))
                     : nullptr)),
      vacuum_(options.vacuum),
      metric_(options.metric),
      assert_strict_braket_symmetry_(options.assert_strict_braket_symmetry),
      spbasis_(options.spbasis),
      first_dummy_index_ordinal_(options.first_dummy_index_ordinal),
      braket_typesetting_(options.braket_typesetting),
      braket_slot_typesetting_(options.braket_slot_typesetting),
      deserialization_symmetry_(options.deserialization_symmetry),
      deserialization_hermiticity_(options.deserialization_hermiticity),
      deserialization_column_symmetry_(
          options.deserialization_column_symmetry) {
  auto config = std::make_shared<CanonicalizationConfig>();
  if (options.tensor_canonicalizers) {
    for (const auto& [label, canonicalizer] : *options.tensor_canonicalizers)
      check_tensor_canonicalizer(canonicalizer);
    config->tensor_canonicalizers = std::move(*options.tensor_canonicalizers);
  }
  config->tensor_canonicalizers.try_emplace(L"",
                                            default_tensor_canonicalizer());
  config->index_comparer =
      options.index_comparer ? make_comparer(std::move(*options.index_comparer))
                             : default_index_comparer();
  config->index_pair_comparer =
      options.index_pair_comparer
          ? make_comparer(std::move(*options.index_pair_comparer))
          : default_index_pair_comparer();
  if (options.cardinal_tensor_labels) {
    check_cardinal_tensor_labels(*options.cardinal_tensor_labels);
    config->cardinal_labels = std::move(*options.cardinal_tensor_labels);
  } else
    config->cardinal_labels = {reserved::antisymm_label(),
                               reserved::symm_label(),
                               reserved::transposition_label()};
  config->options = std::move(options.canonicalization_options);
  canonicalization_config_ = std::move(config);
  update_version();
}

std::uint64_t Context::version() const { return version_; }

Context::CanonicalizationConfig& Context::mutable_canonicalization_config() {
  auto copy =
      std::make_shared<CanonicalizationConfig>(*canonicalization_config_);
  auto& result = *copy;
  canonicalization_config_ = std::move(copy);
  return result;
}

void Context::update_version() {
  // binds every member, so that a member added to CanonicalizationConfig
  // fails to compile here until the key accounts for it
  const auto& [tensor_canonicalizers, index_comparer, index_pair_comparer,
               cardinal_labels, options] = *canonicalization_config_;
  CanonicalizationKey key{.registry = ObjectRef(idx_space_reg_),
                          .spbasis = spbasis_,
                          .canonicalizers = {},
                          .index_comparer = ObjectRef(index_comparer),
                          .index_pair_comparer = ObjectRef(index_pair_comparer),
                          .cardinal_labels = cardinal_labels,
                          .options = options};
  for (const auto& [label, canonicalizer] : tensor_canonicalizers)
    key.canonicalizers.emplace_back(label, ObjectRef(canonicalizer));
  version_ = canonicalization_version(std::move(key));
}

std::uint64_t current_context_version(Statistics s) {
  return get_default_context(s).version();
}

Vacuum Context::vacuum() const { return vacuum_; }

std::shared_ptr<const IndexSpaceRegistry> Context::index_space_registry()
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
  return canonicalization_config_->options;
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
  const auto& map = canonicalization_config_->tensor_canonicalizers;
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
  return *canonicalization_config_->index_comparer;
}

const tensor_index_pair_comparer_t& Context::index_pair_comparer() const {
  return *canonicalization_config_->index_pair_comparer;
}

std::shared_ptr<const tensor_index_comparer_t> Context::index_comparer_ptr()
    const {
  return canonicalization_config_->index_comparer;
}

std::shared_ptr<const tensor_index_pair_comparer_t>
Context::index_pair_comparer_ptr() const {
  return canonicalization_config_->index_pair_comparer;
}

const container::vector<std::wstring>& Context::cardinal_tensor_labels() const {
  return canonicalization_config_->cardinal_labels;
}

Context& Context::set(Vacuum vacuum) {
  vacuum_ = vacuum;
  return *this;
}

Context& Context::set(IndexSpaceRegistry ISR) {
  idx_space_reg_ = std::make_shared<const IndexSpaceRegistry>(std::move(ISR));
  update_version();
  return *this;
}

Context& Context::set(std::shared_ptr<const IndexSpaceRegistry> ISR) {
  idx_space_reg_ = owned(std::move(ISR));
  update_version();
  return *this;
}

Context& Context::set(IndexSpaceMetric metric) {
  metric_ = metric;
  return *this;
}

Context& Context::set(
    AssertStrictBraKetSymmetry assert_strict_braket_symmetry) {
  assert_strict_braket_symmetry_ =
      assert_strict_braket_symmetry == AssertStrictBraKetSymmetry::Yes;
  return *this;
}

Context& Context::set(SPBasis spbasis) {
  spbasis_ = spbasis;
  update_version();
  return *this;
}

Context& Context::set_first_dummy_index_ordinal(
    std::size_t first_dummy_index_ordinal) {
  first_dummy_index_ordinal_ = first_dummy_index_ordinal;
  return *this;
}

Context& Context::set(CanonicalizeOptions copt) {
  mutable_canonicalization_config().options = copt;
  update_version();
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

Context& Context::set(Symmetry symmetry) {
  deserialization_symmetry_ = symmetry;
  return *this;
}

Context& Context::set(Hermiticity hermiticity) {
  deserialization_hermiticity_ = hermiticity;
  return *this;
}

Context& Context::set(ColumnSymmetry column_symmetry) {
  deserialization_column_symmetry_ = column_symmetry;
  return *this;
}

Context& Context::set_tensor_canonicalizer(
    std::wstring_view label,
    std::shared_ptr<TensorCanonicalizer> canonicalizer) {
  check_tensor_canonicalizer(canonicalizer);
  mutable_canonicalization_config().tensor_canonicalizers.insert_or_assign(
      std::wstring{label}, std::move(canonicalizer));
  update_version();
  return *this;
}

Context& Context::unset_tensor_canonicalizer(std::wstring_view label) {
  mutable_canonicalization_config().tensor_canonicalizers.erase(
      std::wstring{label});
  update_version();
  return *this;
}

Context& Context::set_index_comparer(tensor_index_comparer_t comparer) {
  mutable_canonicalization_config().index_comparer =
      make_comparer(std::move(comparer));
  update_version();
  return *this;
}

Context& Context::set_index_comparer(
    std::shared_ptr<const tensor_index_comparer_t> comparer) {
  mutable_canonicalization_config().index_comparer =
      checked_comparer(std::move(comparer));
  update_version();
  return *this;
}

Context& Context::set_index_pair_comparer(
    tensor_index_pair_comparer_t comparer) {
  mutable_canonicalization_config().index_pair_comparer =
      make_comparer(std::move(comparer));
  update_version();
  return *this;
}

Context& Context::set_index_pair_comparer(
    std::shared_ptr<const tensor_index_pair_comparer_t> comparer) {
  mutable_canonicalization_config().index_pair_comparer =
      checked_comparer(std::move(comparer));
  update_version();
  return *this;
}

Context& Context::set_cardinal_tensor_labels(
    container::vector<std::wstring> labels) {
  check_cardinal_tensor_labels(labels);
  mutable_canonicalization_config().cardinal_labels = std::move(labels);
  update_version();
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
