#ifndef SEQUANT_CORE_CONTEXT_HPP
#define SEQUANT_CORE_CONTEXT_HPP

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/index_basis_registry.hpp>
#include <SeQuant/core/options.hpp>
#include <SeQuant/core/tensor_canonicalizer_fwd.hpp>
#include <SeQuant/core/utility/aggregate.hpp>
#include <SeQuant/core/utility/context.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <atomic>
#include <cstdint>
#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <string_view>

namespace sequant {

// clang-format
/// @brief Specifies SeQuant context, such as vacuum choice, whether index
/// spaces are orthonormal, sizes of index spaces, etc.
///
/// SeQuant context contains the following information:
/// - a IndexBasisRegistry object: contains information about the known
///   IndexSpace objects and their attributes; managed by shared_ptr and
///   shared by the copies of a context. A context owns its registry: it moves
///   or copies in a registry it is given, and adopts a shared_ptr to one only
///   if that is the registry's only owner, so the registry cannot change
///   while any context uses it (short of what
///   Context::set(std::shared_ptr<const IndexBasisRegistry>) warns about).
/// - `vacuum`: the vacuum state used to define normal ordering of
/// `NormalOperator`s
/// - `metric`: whether the plain basis of vector space (ket) modes are
/// orthonormal to their dual (bra) counterparts
///   (`IndexSpaceMetric::Unit`) or not (`IndexSpaceMetric::General`);
///    this affects the value of Wick contractions.
/// - `deserialization_symmetry`, `deserialization_hermiticity`,
///    `deserialization_column_symmetry`: the symmetries given to a
///    *deserialized* tensor that does not specify them; the programmatic
///    Tensor ctors are unaffected (see Tensor::Defaults), hence the
///    `deserialization_` prefix. There is no `braket_symmetry` knob: a
///    tensor's BraKetSymmetry is derived from its Hermiticity and base
///    field.
/// - `spbasis`: whether the bra/ket bases are spinor (`SPBasis::Spinor`) or
/// spin-free (`SPBasis::Spinfree`).
/// - `first_dummy_index_ordinal`: during its operation SeQuant will generate
///    temporary indices with orbitals greater or equal to this; to avoid
///    duplicates user Index objects should have ordinals smaller than this
/// - `canonicalization_options`: if set, this specifies the default options to
///    use for canonicalization of expressions.
/// - `braket_typesetting`: whether `to_latex()` typesets tensor indices of ket
///    (covariant, primal) modes as superscript (`BraKetTypesetting::KetSuper`,
///    default)
//     or as subscript (`BraKetTypesetting::KetSub`); the latter is the
//     traditional tensor convention
/// - `braket_slot_typesetting`: whether `to_latex()` typesets tensor indices
/// using
///   `tensor` LaTeX package (`BraKetSlotTypesetting::TensorPackage`, default)
///   or native typesetting (`BraKetSlotTypesetting::Naive`); the former is
///   preferred for alignment of superscript with subscript slots.
/// - `tensor_canonicalizers`: the TensorCanonicalizer objects applied to
///   tensors, keyed by tensor label; the one keyed by the empty label applies
///   to tensors without a label-specific canonicalizer. The canonicalizer
///   objects are shared by copies of the context and must not be mutated.
///   The empty label maps to a DefaultTensorCanonicalizer unless a map given
///   to the constructor sets it; unset_tensor_canonicalizer() can remove it.
/// - `index_comparer`, `index_pair_comparer`: the objects that order Index
///   objects (and pairs thereof) during tensor canonicalization; default to
///   TensorCanonicalizer::default_index_comparer() and
///   TensorCanonicalizer::default_index_pair_comparer().
/// - `cardinal_tensor_labels`: Tensor labels with lexicographic preference
///   (in order); default to `{reserved::antisymm_label(),
///   reserved::symm_label(), reserved::transposition_label()}`.
///
/// @note canonicalization reads `tensor_canonicalizers`, `index_comparer`,
///   `index_pair_comparer` and `cardinal_tensor_labels` only from the default
///   context for Statistics::Arbitrary; in a context installed for another
///   Statistics they are currently ignored, but keep them identical to those
///   of the Statistics::Arbitrary context, since a future version may consult
///   the statistics-specific context first
// clang-format off
class Context {
 public:
  struct Defaults {
    constexpr static auto vacuum = Vacuum::Physical;
    constexpr static auto metric = IndexSpaceMetric::Unit;
    constexpr static auto assert_strict_braket_symmetry = true;
    constexpr static auto spbasis = SPBasis::Spinor;
    constexpr static auto first_dummy_index_ordinal = 100;
    constexpr static auto braket_typesetting = BraKetTypesetting::ContraSub;
    constexpr static auto braket_slot_typesetting =
        BraKetSlotTypesetting::TensorPackage;
    constexpr static auto deserialization_symmetry = Symmetry::Nonsymm;
    constexpr static auto deserialization_hermiticity =
        Hermiticity::NonHermitian;
    constexpr static auto deserialization_column_symmetry =
        ColumnSymmetry::Nonsymm;
  };

  /// helper for the named-parameter constructor of Context

  // the implicit special members of Options use its deprecated fields
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
  /// see the Context documentation for detailed description
  struct Options {
    SEQUANT_DESIGNATED_INIT_ONLY;
      /// a shared_ptr to an IndexBasisRegistry object; the Context adopts it if it is the only owner of the object (e.g. a temporary or a moved-from shared_ptr), else copies the object; see the warning of Context::set(std::shared_ptr<const IndexBasisRegistry>)
      std::shared_ptr<const IndexBasisRegistry> index_basis_registry_shared_ptr = nullptr;
      /// an IndexBasisRegistry object, moved into the Context; used if index_basis_registry_shared_ptr is null and it is nonnull
      std::optional<IndexBasisRegistry> index_basis_registry = std::nullopt;
      /// @deprecated use index_basis_registry_shared_ptr
      [[deprecated("use index_basis_registry_shared_ptr")]] std::shared_ptr<const IndexBasisRegistry> index_space_registry_shared_ptr = nullptr;
      /// @deprecated use index_basis_registry
      [[deprecated("use index_basis_registry")]] std::optional<IndexBasisRegistry> index_space_registry = std::nullopt;
      /// the Vacuum object
      Vacuum vacuum = Defaults::vacuum;
      /// the IndexSpaceMetric object
      IndexSpaceMetric metric = Defaults::metric;
      /// the flag that controls the strictness of bra-ket checks in
      /// tensor network construction
      bool assert_strict_braket_symmetry = Defaults::assert_strict_braket_symmetry;
      /// the SPBasis object
      SPBasis spbasis = Defaults::spbasis;
      /// the first dummy index ordinal
      std::size_t first_dummy_index_ordinal = Defaults::first_dummy_index_ordinal;
      /// the default canonicalization options
      std::optional<CanonicalizeOptions> canonicalization_options = std::nullopt;
      /// the BraKetTypesetting object
      BraKetTypesetting braket_typesetting = Defaults::braket_typesetting;
      /// the BraKetSlotTypesetting object
      BraKetSlotTypesetting braket_slot_typesetting =
        Defaults::braket_slot_typesetting;
      /// the default bra/ket permutational Symmetry for deserialized tensors
      Symmetry deserialization_symmetry = Defaults::deserialization_symmetry;
      /// the default Hermiticity for deserialized tensors; the braket symmetry
      /// of a deserialized tensor is *derived* from this and its #base_field
      Hermiticity deserialization_hermiticity =
        Defaults::deserialization_hermiticity;
      /// the default ColumnSymmetry (particle-permutation symmetry) for
      /// deserialized tensors
      ColumnSymmetry deserialization_column_symmetry =
        Defaults::deserialization_column_symmetry;
      /// label -> TensorCanonicalizer map, without null entries; the empty label maps to a DefaultTensorCanonicalizer unless given
      std::optional<container::map<std::wstring, std::shared_ptr<TensorCanonicalizer>>> tensor_canonicalizers = std::nullopt;
      /// the Index comparer used by tensor canonicalizers; if not set, TensorCanonicalizer::default_index_comparer()
      std::optional<tensor_index_comparer_t> index_comparer = std::nullopt;
      /// the Index pair comparer used by tensor canonicalizers; if not set, TensorCanonicalizer::default_index_pair_comparer()
      std::optional<tensor_index_pair_comparer_t> index_pair_comparer = std::nullopt;
      /// the cardinal Tensor labels; if not set, the reserved labels
      std::optional<container::vector<std::wstring>> cardinal_tensor_labels = std::nullopt;
  };
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
  static Options make_default_options() { return {}; }

  /// @brief standard named-parameter constructor
  ///
  /// @warning default constructor does not create an IndexBasisRegistry, thus
  /// `this->index_basis_registry()` will return nullptr
  /// Example:
  /// ```cpp
  ///   Context ctx({.vacuum = Vacuum::SingleProduct, .spbasis = SPBasis::Spinfree});
  /// ```
  Context(Options options = make_default_options());

  ~Context() = default;

  /// copy constructor
  /// @param[in] ctx a Context
  /// @note created Context uses the same index basis registry as @p ctx
  Context(const Context& ctx) = default;

  /// copy assignment
  /// @param[in] ctx a Context
  /// @note this object will use the same index basis registry as @p ctx
  /// @return reference to this object
  Context& operator=(const Context& ctx) = default;

  // no move operations, so that an rvalue is copied: a copy is cheap (it
  // shares the registry and the canonicalizer configuration), and a moved-from
  // Context would lack them

  /// @return the version of this context's canonicalization configuration: a
  /// nonzero number that two contexts share if and only if canonicalization
  /// sees the same configuration in both, i.e. they have the same index space
  /// registry, tensor canonicalizers and index comparers (the same objects),
  /// cardinal tensor labels, canonicalization options and SP basis (which
  /// determines the symmetry of NormalOperator); a number is never reused for
  /// another configuration
  /// @note the other settings (vacuum, metric, first dummy index ordinal,
  /// typesetting, deserialization defaults, strict bra-ket checks) do not
  /// affect the version
  /// @note the version does not track in-place mutation of a
  /// TensorCanonicalizer or comparer object that this context refers to
  std::uint64_t version() const;

  /// \return Vacuum of this context
  Vacuum vacuum() const;
  /// @return a constant pointer to the IndexBasisRegistry for this context
  /// @warning can be null when user did not provide one to Context (i.e., it
  /// was default constructed)
  /// @note the registry does not change; to change it, set a modified copy
  std::shared_ptr<const IndexBasisRegistry> index_basis_registry() const;
  /// @deprecated use index_basis_registry()
  [[deprecated("use index_basis_registry()")]] std::shared_ptr<
      const IndexBasisRegistry>
  index_space_registry() const {
    return index_basis_registry();
  }
  /// \return IndexSpaceMetric of this context
  IndexSpaceMetric metric() const;
  /// \return true if strict bra-ket symmetry is asserted;
  /// setting this to false (via `Context::set(AssertStrictBraKetSymmetry::No)`)
  /// allows arbitrary contractions of bra/ket modes as if they were aux indices.
  /// @note disabling strict bra-ket symmetry assertions is a bad idea if
  /// working with complex-valued tensors, or nonunit metric, or in general
  /// @warning if SeQuant was configured with `SEQUANT_ASSERT_BEHAVIOR=IGNORE`
  /// the strict bra-ket contraction checks are turned off even if
  /// this returns true
  bool assert_strict_braket_symmetry() const;
  /// \return SPBasis of this context
  SPBasis spbasis() const;
  /// \return first ordinal of the dummy indices generated by calls to
  /// Index::next_tmp_index when this context is active
  std::size_t first_dummy_index_ordinal() const;
  /// \return canonicalization options to use by default, if nonnull
  std::optional<CanonicalizeOptions> canonicalization_options() const;
  /// \return BraKetTypesetting of this context; if this returns
  /// BraKetTypesetting::ContraSub covariant (ket, creation) and contravariant
  /// (bra, annihilation) indices are typeset in superscript and subscript
  /// LaTeX, respectively.
  BraKetTypesetting braket_typesetting() const;
  /// \return BraKetSlotTypesetting of this context; see BraKetSlotTypesetting
  /// for the meaning of the possible values
  BraKetSlotTypesetting braket_slot_typesetting() const;
  /// \return the default bra/ket permutational Symmetry for deserialized
  /// tensors; programmatic Tensor construction is unaffected by it
  Symmetry deserialization_symmetry() const;
  /// \return the default Hermiticity for deserialized tensors; the braket
  /// symmetry of a deserialized tensor is *derived* from this and its
  /// #base_field. Programmatic Tensor construction is unaffected by it
  Hermiticity deserialization_hermiticity() const;
  /// \return the default ColumnSymmetry (particle-permutation symmetry) for
  /// deserialized tensors; programmatic Tensor construction is unaffected by it
  ColumnSymmetry deserialization_column_symmetry() const;
  /// @param label a Tensor label
  /// @return the TensorCanonicalizer for @p label if any, else the one for
  /// the empty label; null if neither exists
  std::shared_ptr<TensorCanonicalizer> tensor_canonicalizer_ptr(
      std::wstring_view label) const;
  /// @param label a Tensor label
  /// @return the TensorCanonicalizer for exactly @p label, or null
  /// @sa tensor_canonicalizer_ptr
  std::shared_ptr<TensorCanonicalizer> nondefault_tensor_canonicalizer_ptr(
      std::wstring_view label) const;
  /// @param label a Tensor label
  /// @return the TensorCanonicalizer that `tensor_canonicalizer_ptr(label)`
  /// points to
  /// @throw Exception if `tensor_canonicalizer_ptr(label)` is null
  /// @warning the reference is valid only while the map entry that holds it
  /// exists in this context; use tensor_canonicalizer_ptr() if this context
  /// may change
  const TensorCanonicalizer& tensor_canonicalizer(
      std::wstring_view label) const;
  /// \return the object used by tensor canonicalizers to compare Index objects
  const tensor_index_comparer_t& index_comparer() const;
  /// \return the object used by tensor canonicalizers to compare pairs of
  /// Index objects
  const tensor_index_pair_comparer_t& index_pair_comparer() const;
  /// \return the shared object that index_comparer() refers to; passing it to
  /// set_index_comparer() keeps the context equal to its unmodified copies
  std::shared_ptr<const tensor_index_comparer_t> index_comparer_ptr() const;
  /// \return the shared object that index_pair_comparer() refers to; passing
  /// it to set_index_pair_comparer() keeps the context equal to its
  /// unmodified copies
  std::shared_ptr<const tensor_index_pair_comparer_t> index_pair_comparer_ptr()
      const;
  /// \return Tensor labels with lexicographic preference (in order)
  const container::vector<std::wstring>& cardinal_tensor_labels() const;

  /// Sets the Vacuum for this context, convenient for chaining
  /// \param vacuum Vacuum
  /// \return ref to `*this`, for chaining
  Context& set(Vacuum vacuum);
  /// sets the IndexBasisRegistry for this context
  /// \param ISR an IndexBasisRegistry, moved into this context
  /// \return ref to '*this' for chaining
  /// \sa the warning of set(std::shared_ptr<const IndexBasisRegistry>)
  Context& set(IndexBasisRegistry ISR);
  /// sets the IndexBasisRegistry for this context
  /// \param ISR a IndexBasisRegistry shared_ptr; adopted if it is the only
  /// owner of its object (e.g. a temporary or a moved-from shared_ptr), else
  /// the object is copied
  /// \return ref to '*this' for chaining
  /// \warning a registry that is adopted or moved in is the caller's object,
  /// and whatever else the caller kept that reaches it still does: a pointer
  /// or reference to it or into it (e.g. from the non-const
  /// IndexBasisRegistry::retrieve_ptr()), a std::weak_ptr to it, or a
  /// shared_ptr that does not own it (e.g. one with a no-op deleter), which
  /// cannot be told from its only owner. Modifying the registry through any
  /// of these goes unnoticed by the contexts that hold it.
  Context& set(std::shared_ptr<const IndexBasisRegistry> ISR);
  /// Sets the IndexSpaceMetric for this context, convenient for chaining
  /// \param metric IndexSpaceMetric
  /// \return ref to `*this`, for chaining
  Context& set(IndexSpaceMetric metric);
  /// Sets the bra-ket strict assertion flag for this context, convenient for chaining
  /// \param assert_strict_braket_symmetry AssertStrictBraKetSymmetry
  /// \return ref to `*this`, for chaining
  Context& set(AssertStrictBraKetSymmetry assert_strict_braket_symmetry);
  /// Sets the SPBasis for this context, convenient for chaining
  /// \param spbasis SPBasis
  /// \return ref to `*this`, for chaining
  Context& set(SPBasis spbasis);
  /// Sets the first dummy index ordinal for this context, convenient for
  /// chaining \param first_dummy_index_ordinal the first dummy index ordinal
  /// \return ref to `*this`, for chaining
  Context& set_first_dummy_index_ordinal(std::size_t first_dummy_index_ordinal);
  /// Specifies the canonicalization options
  /// \return ref to `*this`, for chaining
  Context& set(CanonicalizeOptions copt);
  /// Sets the BraKetTypesetting for this context, convenient for chaining
  /// \param braket_typeset BraKetTypesetting
  /// \return ref to `*this`, for chaining
  Context& set(BraKetTypesetting braket_typeset);
  /// Sets the BraKetSlotTypesetting for this context, convenient for chaining
  /// \param braket_slot_typeset BraKetSlotTypesetting
  /// \return ref to `*this`, for chaining
  Context& set(BraKetSlotTypesetting braket_slot_typeset);
  /// Sets the default bra/ket permutational Symmetry for deserialized tensors
  /// \return ref to `*this`, for chaining
  Context& set(Symmetry symmetry);
  /// Sets the default Hermiticity for deserialized tensors (the braket symmetry
  /// of a deserialized tensor is derived from this and its base field)
  /// \return ref to `*this`, for chaining
  Context& set(Hermiticity hermiticity);
  /// Sets the default ColumnSymmetry for deserialized tensors
  /// \return ref to `*this`, for chaining
  Context& set(ColumnSymmetry column_symmetry);
  /// Sets the TensorCanonicalizer for Tensor objects labeled @p label ,
  /// replacing the existing one, if any; the empty label applies to Tensor
  /// objects without a label-specific canonicalizer
  /// \param canonicalizer a nonnull TensorCanonicalizer
  /// \throw Exception if @p canonicalizer is null
  /// \return ref to `*this`, for chaining
  /// \warning version() identifies the canonicalizer by the object the
  /// shared_ptr owns, so a shared_ptr that does not own its object (e.g. one
  /// with a no-op deleter) makes the version change when the shared_ptr, not
  /// the object, dies
  Context& set_tensor_canonicalizer(
      std::wstring_view label, std::shared_ptr<TensorCanonicalizer> canonicalizer);
  /// Removes the TensorCanonicalizer for @p label , if any
  /// \return ref to `*this`, for chaining
  Context& unset_tensor_canonicalizer(std::wstring_view label);
  /// Sets the Index comparer used by tensor canonicalizers
  /// \param comparer a nonempty Index comparer
  /// \return ref to `*this`, for chaining
  Context& set_index_comparer(tensor_index_comparer_t comparer);
  /// Sets the Index comparer used by tensor canonicalizers to a shared object
  /// \param comparer a nonnull pointer to a nonempty Index comparer, e.g.
  /// one obtained from index_comparer_ptr()
  /// \return ref to `*this`, for chaining
  /// \warning see the warning of set_tensor_canonicalizer() about a
  /// shared_ptr that does not own its object
  Context& set_index_comparer(
      std::shared_ptr<const tensor_index_comparer_t> comparer);
  /// Sets the Index pair comparer used by tensor canonicalizers
  /// \param comparer a nonempty Index pair comparer
  /// \return ref to `*this`, for chaining
  Context& set_index_pair_comparer(tensor_index_pair_comparer_t comparer);
  /// Sets the Index pair comparer used by tensor canonicalizers to a shared
  /// object
  /// \param comparer a nonnull pointer to a nonempty Index pair comparer,
  /// e.g. one obtained from index_pair_comparer_ptr()
  /// \return ref to `*this`, for chaining
  /// \warning see the warning of set_tensor_canonicalizer() about a
  /// shared_ptr that does not own its object
  Context& set_index_pair_comparer(
      std::shared_ptr<const tensor_index_pair_comparer_t> comparer);
  /// Sets the cardinal Tensor labels
  /// \param labels the complete list of cardinal labels, without duplicates;
  /// the default labels (reserved::antisymm_label(), reserved::symm_label(),
  /// reserved::transposition_label()) are not prepended, so include them where
  /// they should keep their precedence (mbpt::cardinal_tensor_labels() returns
  /// such a complete list)
  /// \return ref to `*this`, for chaining
  Context& set_cardinal_tensor_labels(container::vector<std::wstring> labels);

 private:
  /// the settings that control canonicalization (which also reads the index
  /// space registry and the SP basis, settings not specific to it); these, the
  /// registry and the SP basis define version(), so a setting that
  /// canonicalization comes to read belongs here, and in the key of the
  /// version table (CanonicalizationKey in context.cpp) that mirrors this;
  /// comparers are held by shared_ptr so that equality can be decided by
  /// identity
  struct CanonicalizationConfig {
    container::map<std::wstring, std::shared_ptr<TensorCanonicalizer>>
        tensor_canonicalizers;
    std::shared_ptr<const tensor_index_comparer_t> index_comparer;
    std::shared_ptr<const tensor_index_pair_comparer_t> index_pair_comparer;
    container::vector<std::wstring> cardinal_labels;
    std::optional<CanonicalizeOptions> options;

    /// canonicalizers and comparers compare by identity, labels and options
    /// by value
    bool operator==(const CanonicalizationConfig&) const = default;
  };

  friend bool operator==(const Context& ctx1, const Context& ctx2);

  /// replaces the canonicalization configuration by a copy owned by this
  /// @return the copy, for the caller to modify
  CanonicalizationConfig& mutable_canonicalization_config();

  /// marks the cached version stale, to be recomputed by the next version()
  void invalidate_version();

  /// the cached version(), 0 while stale; atomic because version() may be
  /// called on a context that threads share (e.g. the process-wide default),
  /// copyable so that Context stays copyable
  struct CachedVersion {
    std::atomic<std::uint64_t> value{0};
    CachedVersion() = default;
    CachedVersion(const CachedVersion& other)
        : value(other.value.load(std::memory_order_relaxed)) {}
    CachedVersion& operator=(const CachedVersion& other) {
      value.store(other.value.load(std::memory_order_relaxed),
                  std::memory_order_relaxed);
      return *this;
    }
  };
  mutable CachedVersion version_;

  std::shared_ptr<const IndexBasisRegistry> idx_basis_reg_ = nullptr;
  Vacuum vacuum_ = Defaults::vacuum;
  IndexSpaceMetric metric_ = Defaults::metric;
  bool assert_strict_braket_symmetry_ = Defaults::assert_strict_braket_symmetry;
  SPBasis spbasis_ = Defaults::spbasis;
  std::size_t first_dummy_index_ordinal_ = Defaults::first_dummy_index_ordinal;
  BraKetTypesetting braket_typesetting_ = Defaults::braket_typesetting;
  BraKetSlotTypesetting braket_slot_typesetting_ =
      Defaults::braket_slot_typesetting;
  Symmetry deserialization_symmetry_ = Defaults::deserialization_symmetry;
  Hermiticity deserialization_hermiticity_ =
      Defaults::deserialization_hermiticity;
  ColumnSymmetry deserialization_column_symmetry_ =
      Defaults::deserialization_column_symmetry;
  /// shared by copies of this context, hence never mutated in place
  std::shared_ptr<const CanonicalizationConfig> canonicalization_config_;
};

/// Context object equality comparison
/// \param ctx1
/// \param ctx2
/// \return true if \p ctx1 and \p ctx2 are equal
/// \note index basis registries and cardinal tensor labels are compared by
/// value (contexts without a registry are equal in that respect), the
/// registries including the approximate sizes and fields of their spaces;
/// tensor canonicalizers and index comparers by identity, hence
/// a comparer replaced by a behaviourally identical one compares unequal
/// (re-install a comparer through its shared pointer, e.g.
/// Context::index_comparer_ptr(), to keep contexts equal)
/// \note the versions of the contexts are ignored, and equal contexts may
/// have different ones: Context::version() compares index basis registries as
/// objects, not by value
bool operator==(const Context& ctx1, const Context& ctx2);

/// Context object inequality comparison
/// \param ctx1
/// \param ctx2
/// \return true if \p ctx1 and \p ctx2 are not equal
/// \sa operator==(const Context&, const Context&)
bool operator!=(const Context& ctx1, const Context& ctx2);

/// \name manipulation of implicit context for SeQuant
/// \warning set_default_context(), reset_default_context() and the reads of
/// the process-wide context (by get_default_context(),
/// get_default_context_snapshot() and set_scoped_modified_default_context() on
/// a thread without scoped contexts)
/// are thread-safe only if default_context_manipulation_threadsafe() returns
/// true; set_scoped_default_context(), and the readers on a thread with scoped
/// contexts, touch only thread-local state

/// @{

/// \return whether context manipulation functions are thread-safe
bool default_context_manipulation_threadsafe();

/// @brief the version of the default Context for the given Statistics
/// @param s Statistics
/// @return `get_default_context(s).version()`, i.e. the version of the context
/// in effect on the calling thread (see Context::version())
/// @note canonicalization reads the SP basis, which determines the symmetry
/// of NormalOperator<S>, from the context for the operator's Statistics `S`,
/// and the rest from the context for Statistics::Arbitrary; a cache of
/// canonicalization results is therefore keyed on the versions for all
/// Statistics, not only the one for Statistics::Arbitrary
std::uint64_t current_context_version(Statistics s = Statistics::Arbitrary);

/// @return a value that changes whenever the canonicalization configuration
/// of the effective context of any Statistics changes: the versions of the
/// FermiDirac, BoseEinstein and
/// Arbitrary contexts in effect on the calling thread (see
/// current_context_version()) combined with hash::combine. The canonical mark
/// (see Expr::is_canonical()) is keyed on it.
std::uint64_t current_contexts_version();

/// @brief access default Context for the given Statistics
/// @param s Statistics
/// @return the default context used for Statistics @p s
/// @warning on a thread without scoped contexts the reference is to the
/// process-wide context, which set_default_context() and
/// reset_default_context() replace, on any thread; to hold the context, or
/// anything obtained by reference from it, beyond a brief read use
/// get_default_context_snapshot()
const Context& get_default_context(Statistics s = Statistics::Arbitrary);

/// @brief copy of the default Context for the given Statistics
/// @param s Statistics
/// @return a copy of `get_default_context(s)` that reflects every change of
/// the process-wide context completed before the call; it shares the index
/// space registry and the canonicalizer configuration with its source, hence
/// is cheap, and remains valid when the default context is replaced
/// @note takes the lock that guards the process-wide context only if that
/// changed since the previous snapshot on this thread
Context get_default_context_snapshot(Statistics s = Statistics::Arbitrary);

/// @brief the index basis registry of the default Context for the given
/// Statistics
/// @param s Statistics
/// @return `get_default_context(s).index_basis_registry()`, read like
/// get_default_context_snapshot() reads the context: from the calling thread's
/// innermost scoped context if it has one, else from this thread's copy of the
/// published process-wide contexts, so the lock that guards them is taken only
/// when they changed since the previous read on this thread; costs one
/// shared_ptr copy and no Context copy. Null if that context has no registry
std::shared_ptr<const IndexBasisRegistry> get_default_index_basis_registry(
    Statistics s = Statistics::Arbitrary);

/// @deprecated use get_default_index_basis_registry()
[[deprecated("use get_default_index_basis_registry()")]] inline std::shared_ptr<
    const IndexBasisRegistry>
get_default_index_space_registry(Statistics s = Statistics::Arbitrary) {
  return get_default_index_basis_registry(s);
}

/// @return the entry @p basis is registered as in the registry of
/// get_default_index_basis_registry() (see IndexBasisRegistry::resolve()),
/// which carries a named basis instance's name, approximate size and field;
/// @p basis if there is no such registry
IndexBasis default_registry_resolved(const IndexBasis& basis);

/// @brief sets default Context for the given Statistics
/// @param ctx Context object
/// @param s Statistics
void set_default_context(Context ctx,
                         Statistics s = Statistics::Arbitrary);

/// @brief sets default Context for the given Statistics
/// @param ctx_options Context named-parameter constructor arguments
/// @param s Statistics
void set_default_context(Context::Options ctx_options,
                         Statistics s = Statistics::Arbitrary);

/// @brief sets default Context for several given Statistics
/// @param ctxs a Statistics->Context map
void set_default_context(const container::map<Statistics, Context>& ctxs);

/// @brief resets default Contexts for all statistics to their initial values
void reset_default_context();

/// @brief changes default contexts
/// @param ctx Context objects for one or more statistics
/// @return a move-only ContextResetter object whose destruction will reset the
/// default context to the previous value. Example:
/// ```cpp
/// {
///   auto resetter = set_scoped_default_context({{Statistics::Arbitrary,
///   ctx}});
///   // ctx is now the default context for all statistics
/// } // leaving scope, resetter is destroyed, default context is reset back to
/// the old value
/// ```
/// @note the scoped contexts are seen only by the calling thread and by the
/// workers of the parallel primitives (sequant::for_each, etc.) it launches;
/// set_default_context() is not seen by this thread until the scope ends
/// @note scopes must end in the reverse order of their creation
[[nodiscard]] detail::ImplicitContextResetter<
    container::map<Statistics, Context>>
set_scoped_default_context(container::map<Statistics, Context> ctx);

/// @brief changes default context for arbitrary statistics
/// @note equivalent to `set_scoped_default_context({{Statistics::Arbitrary,
/// ctx}})`
[[nodiscard]] detail::ImplicitContextResetter<
    container::map<Statistics, Context>>
set_scoped_default_context(Context ctx);

/// @brief changes default context for arbitrary statistics
/// @note equivalent to `set_scoped_default_context({{Statistics::Arbitrary,
/// Context{ctx_options}}})`
[[nodiscard]] detail::ImplicitContextResetter<
    container::map<Statistics, Context>>
set_scoped_default_context(Context::Options ctx_options);

/// @brief changes the default contexts for all statistics by a modification
/// @param modify applied to a copy of each current default context
/// @return a move-only ContextResetter object whose destruction will reset the
/// default contexts to the previous values
/// @note see the notes of set_scoped_default_context()
[[nodiscard]] detail::ImplicitContextResetter<
    container::map<Statistics, Context>>
set_scoped_modified_default_context(
    const std::function<void(Context&)>& modify);

/// @brief pins the default contexts in effect for the calling thread
/// @return a move-only ContextResetter object that scopes a copy of the
/// default contexts in effect, if the calling thread has none scoped yet, so
/// that every read of the contexts during its lifetime sees the same ones even
/// if another thread replaces the process-wide contexts meanwhile; a scoped
/// context cannot change under the caller, so then nothing is installed and
/// the object is empty
/// @note the copies have the versions of the contexts in effect, so canonical
/// marks made under one are valid under the other
[[nodiscard]] detail::ImplicitContextResetter<
    container::map<Statistics, Context>>
pin_default_contexts();

///@}

/// \name particle, hole and complete space accessors
/// Syntax sugar for accessing particle and hole spaces in the current context

///@{

/// @brief returns the particle space defined in the current context
/// @param qn QuantumNumbers of the space
/// @return IndexSpace object representing the particle space
[[nodiscard]] IndexSpace get_particle_space(
    const IndexSpace::QuantumNumbers& qn);

/// @brief returns the hole space defined in the current context
/// @param qn QuantumNumbers of the space
/// @return IndexSpace object representing the hole space
[[nodiscard]] IndexSpace get_hole_space(const IndexSpace::QuantumNumbers& qn);

/// @brief returns the complete space defined in the current context
/// @param qn QuantumNumbers of the space
/// @return IndexSpace object representing the complete space
[[nodiscard]] IndexSpace get_complete_space(
    const IndexSpace::QuantumNumbers& qn);

///@}

}  // namespace sequant

#endif
