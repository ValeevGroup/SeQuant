#ifndef SEQUANT_CORE_CONTEXT_HPP
#define SEQUANT_CORE_CONTEXT_HPP

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/index_space_registry.hpp>
#include <SeQuant/core/options.hpp>
#include <SeQuant/core/utility/context.hpp>

namespace sequant {

// clang-format
/// @brief Specifies SeQuant context, such as vacuum choice, whether index
/// spaces are orthonormal, sizes of index spaces, etc.
///
/// SeQuant context contains the following information:
/// - a IndexSpaceRegistry object: contains information about the known
///   IndexSpace objects and their attributes; managed by shared_ptr and
///   can be shared by multiple contexts.
/// - `vacuum`: the vacuum state used to define normal ordering of
/// `NormalOperator`s
/// - `metric`: whether the plain basis of vector space (ket) modes are
/// orthonormal to their dual (bra) counterparts
///   (`IndexSpaceMetric::Unit`) or not (`IndexSpaceMetric::General`);
///    this affects the value of Wick contractions.
/// - `braket_symmetry`: whether the primal (ket) and dual (bra) vector space
/// _bases_
///    are "equivalent" (homogenous to each other; `BraKetSymmetry::symm`),
///    "conjugate to each other" (<a
///    href="https://en.wikipedia.org/wiki/Antilinear_map">conjugate-homogenous</a>
///    to each other; `BraKetSymmetry::conjugate`) or are "nonequivalent" (not
///    homogeneous; `BraKetSymmetry::nonsymm`)
/// - `spbasis`: whether the bra/ket bases are spinor (`SPBasis::spinor`) or
/// spin-free (`SPBasis::spinfree`).
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
// clang-format off
class Context {
 public:
  struct Defaults {
    constexpr static auto vacuum = Vacuum::Physical;
    constexpr static auto metric = IndexSpaceMetric::Unit;
    constexpr static auto braket_symmetry = BraKetSymmetry::conjugate;
    constexpr static auto spbasis = SPBasis::spinor;
    constexpr static auto first_dummy_index_ordinal = 100;
    constexpr static std::optional<CanonicalizeOptions> canonicalization_options = std::nullopt;
    constexpr static auto braket_typesetting = BraKetTypesetting::ContraSub;
    constexpr static auto braket_slot_typesetting =
        BraKetSlotTypesetting::TensorPackage;
  };

  /// helper for the named-parameter constructor of Context

  /// see the Context documentation for detailed description
  struct Options {
      /// a shared_ptr to an IndexSpaceRegistry object
      std::shared_ptr<IndexSpaceRegistry> index_space_registry_shared_ptr = nullptr;
      /// an IndexSpaceRegistry object; used if index_space_registry_shared_ptr is null and it is nonnull
      std::optional<IndexSpaceRegistry> index_space_registry = std::nullopt;
      /// the Vacuum object
      Vacuum vacuum = Defaults::vacuum;
      /// the IndexSpaceMetric object
      IndexSpaceMetric metric = Defaults::metric;
      /// the BraKetSymmetry object
      BraKetSymmetry braket_symmetry = Defaults::braket_symmetry;
      /// the SPBasis object
      SPBasis spbasis = Defaults::spbasis;
      /// the first dummy index ordinal
      std::size_t first_dummy_index_ordinal = Defaults::first_dummy_index_ordinal;
      /// the default canonicalization options
      std::optional<CanonicalizeOptions> canonicalization_options = Defaults::canonicalization_options;
      /// the BraKetTypesetting object
      BraKetTypesetting braket_typesetting = Defaults::braket_typesetting;
      /// the BraKetSlotTypesetting object
      BraKetSlotTypesetting braket_slot_typesetting =
        Defaults::braket_slot_typesetting;
  };

  /// @brief standard named-parameter constructor
  ///
  /// @warning default constructor does not create an IndexSpaceRegistry, thus
  /// `this->index_space_registry()` will return nullptr
  /// Example:
  /// ```cpp
  ///   Context ctx({.vacuum = Vacuum::SingleReference, .spbasis = SPBasis::spinfree});
  /// ```
  Context(Options options = {.index_space_registry_shared_ptr = nullptr, .index_space_registry = std::nullopt, .vacuum = Defaults::vacuum, .metric = Defaults::metric, .braket_symmetry =  Defaults::braket_symmetry, .spbasis = Defaults::spbasis, .first_dummy_index_ordinal = Defaults::first_dummy_index_ordinal, .canonicalization_options = Defaults::canonicalization_options, .braket_typesetting = Defaults::braket_typesetting, .braket_slot_typesetting =
        Defaults::braket_slot_typesetting});

  ~Context() = default;

  /// copy constructor
  /// @param[in] ctx a Context
  /// @warning created Context uses the same index space registry as @p ctx
  /// @sa clone()
  Context(const Context& ctx) = default;

  /// copy assignment
  /// @param[in] ctx a Context
  /// @warning this object will use the same index space registry as @p ctx
  /// @sa clone()
  /// @return reference to this object
  Context& operator=(const Context& ctx) = default;

  /// clones this object AND its index space registry
  /// @note created Context does not use this object's index space registry
  Context clone() const;

  Context(Context&&) = default;

  /// \return Vacuum of this context
  Vacuum vacuum() const;
  /// @return a constant pointer to the IndexSpaceRegistry for this context
  /// @warning can be null when user did not provide one to Context (i.e., it
  /// was default constructed)
  std::shared_ptr<const IndexSpaceRegistry> index_space_registry() const;
  /// @return a pointer to the IndexSpaceRegistry for this context.
  /// @throw std::logic_error if the IndexSpaceRegistry is null
  std::shared_ptr<IndexSpaceRegistry> mutable_index_space_registry() const;
  /// \return IndexSpaceMetric of this context
  IndexSpaceMetric metric() const;
  /// \return BraKetSymmetry of this context
  BraKetSymmetry braket_symmetry() const;
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

  /// Sets the Vacuum for this context, convenient for chaining
  /// \param vacuum Vacuum
  /// \return ref to `*this`, for chaining
  Context& set(Vacuum vacuum);
  /// sets the IndexSpaceRegistry for this context
  /// \param ISR an IndexSpaceRegistry
  /// \return ref to '*this' for chaining
  Context& set(IndexSpaceRegistry ISR);
  /// sets the IndexSpaceRegistry for this context
  /// \param ISR a IndexSpaceRegistry shared_ptr
  /// \return ref to '*this' for chaining
  Context& set(std::shared_ptr<IndexSpaceRegistry> ISR);
  /// Sets the IndexSpaceMetric for this context, convenient for chaining
  /// \param metric IndexSpaceMetric
  /// \return ref to `*this`, for chaining
  Context& set(IndexSpaceMetric metric);
  /// Sets the BraKetSymmetry for this context, convenient for chaining
  /// \param braket_symmetry BraKetSymmetry
  /// \return ref to `*this`, for chaining
  Context& set(BraKetSymmetry braket_symmetry);
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

 private:
  std::shared_ptr<IndexSpaceRegistry> idx_space_reg_ = nullptr;
  Vacuum vacuum_ = Defaults::vacuum;
  IndexSpaceMetric metric_ = Defaults::metric;
  BraKetSymmetry braket_symmetry_ = Defaults::braket_symmetry;
  SPBasis spbasis_ = Defaults::spbasis;
  std::size_t first_dummy_index_ordinal_ = Defaults::first_dummy_index_ordinal;
  std::optional<CanonicalizeOptions> canonicalization_options_ = std::nullopt;
  BraKetTypesetting braket_typesetting_ = Defaults::braket_typesetting;
  BraKetSlotTypesetting braket_slot_typesetting_ =
      Defaults::braket_slot_typesetting;
};

/// Context object equality comparison
/// \param ctx1
/// \param ctx2
/// \return true if \p ctx1 and \p ctx2 are equal
/// \warning does not compare index registries
bool operator==(const Context& ctx1, const Context& ctx2);

/// Context object inequality comparison
/// \param ctx1
/// \param ctx2
/// \return true if \p ctx1 and \p ctx2 are not equal
/// \warning does not compare index registries
bool operator!=(const Context& ctx1, const Context& ctx2);

/// \name manipulation of implicit context for SeQuant
/// \warning all of these are thread-safe only if
/// default_context_manipulation_threadsafe() returns true

/// @{

/// \return whether context manipulation functions are thread-safe
inline constexpr bool default_context_manipulation_threadsafe() {
#ifdef SEQUANT_CONTEXT_MANIPULATION_THREADSAFE
  return true;
#else
  return false;
#endif
}

/// @brief access default Context for the given Statistics
/// @param s Statistics
/// @return the default context used for Statistics @p s
const Context& get_default_context(Statistics s = Statistics::Arbitrary);

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
[[nodiscard]] detail::ImplicitContextResetter<
    container::map<Statistics, Context>>
set_scoped_default_context(const container::map<Statistics, Context>& ctx);

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
