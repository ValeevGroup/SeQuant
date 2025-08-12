#ifndef SEQUANT_CORE_CONTEXT_HPP
#define SEQUANT_CORE_CONTEXT_HPP

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/index_space_registry.hpp>
#include <SeQuant/core/utility/context.hpp>

namespace sequant {

/// Specifies SeQuant context, such as vacuum choice, whether index
/// spaces are orthonormal, sizes of index spaces, etc.
class Context {
 public:
  struct Defaults {
    constexpr static auto vacuum = Vacuum::Physical;
    constexpr static auto metric = IndexSpaceMetric::Unit;
    constexpr static auto braket_symmetry = BraKetSymmetry::conjugate;
    constexpr static auto spbasis = sequant::SPBasis::spinorbital;
    constexpr static auto first_dummy_index_ordinal = 100;
    constexpr static auto braket_typesetting = BraKetTypesetting::ContraSub;
    constexpr static auto braket_slot_typesetting =
        BraKetSlotTypesetting::TensorPackage;
  };

  /// standard full-form constructor

  /// @param isr_sptr a shared_ptr to an IndexSpaceRegistry object
  /// @param vac a Vacuum object
  /// @param m an IndexSpaceMetric object
  /// @param bks a BraKetSymmetry object
  /// @param spb single-particle basis (spin-free or spin-dependent)
  /// @param fdio first dummy index ordinal
  /// @param bkt a BraKetTypesetting object
  /// @param bkst a BraKetSlotTypesetting object
  explicit Context(
      std::shared_ptr<IndexSpaceRegistry> isr_sptr,
      Vacuum vac = Defaults::vacuum, IndexSpaceMetric m = Defaults::metric,
      BraKetSymmetry bks = Defaults::braket_symmetry,
      SPBasis spb = Defaults::spbasis,
      std::size_t fdio = Defaults::first_dummy_index_ordinal,
      BraKetTypesetting bkt = Defaults::braket_typesetting,
      BraKetSlotTypesetting bkst = Defaults::braket_slot_typesetting);

  /// @brief same as the standard ctor, using IndexSpaceRegistry passed by value

  /// @param isr an IndexSpaceRegistry object
  /// @param vac a Vacuum object
  /// @param m an IndexSpaceMetric object
  /// @param bks a BraKetSymmetry object
  /// @param spb single-particle basis (spin-free or spin-dependent)
  /// @param fdio first dummy index ordinal
  /// @param bkt a BraKetTypesetting object
  /// @param bkst a BraKetSlotTypesetting object
  explicit Context(
      IndexSpaceRegistry isr, Vacuum vac = Defaults::vacuum,
      IndexSpaceMetric m = Defaults::metric,
      BraKetSymmetry bks = Defaults::braket_symmetry,
      SPBasis spb = Defaults::spbasis,
      std::size_t fdio = Defaults::first_dummy_index_ordinal,
      BraKetTypesetting bkt = Defaults::braket_typesetting,
      BraKetSlotTypesetting bkst = Defaults::braket_slot_typesetting);

  /// @brief same as the standard ctor, using default-constructed
  /// IndexSpaceRegistry

  /// @param vac a Vacuum object
  /// @param m an IndexSpaceMetric object
  /// @param bks a BraKetSymmetry object
  /// @param spb single-particle basis (spin-free or spin-dependent)
  /// @param fdio first dummy index ordinal
  /// @param bkt a BraKetTypesetting object
  /// @param bkst a BraKetSlotTypesetting object
  explicit Context(
      Vacuum vac, IndexSpaceMetric m = Defaults::metric,
      BraKetSymmetry bks = Defaults::braket_symmetry,
      SPBasis spb = Defaults::spbasis,
      std::size_t fdio = Defaults::first_dummy_index_ordinal,
      BraKetTypesetting bkt = Defaults::braket_typesetting,
      BraKetSlotTypesetting bkst = Defaults::braket_slot_typesetting);

  /// default constructor, equivalent to Context(Vacuum::Physical,
  /// IndexSpaceMetric::Unit, BraKetSymmetry::conjugate,
  /// sequant::SPBasis::spinorbital, 100)
  /// @warning default constructor does not create an IndexSpaceRegistry, thus
  /// `this->index_space_registry()` will return nullptr
  Context() = default;

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
void set_default_context(const Context& ctx,
                         Statistics s = Statistics::Arbitrary);

/// @brief sets default Context for several given Statistics
/// @param ctxs a Statistics->Context map
void set_default_context(const container::map<Statistics, Context>& ctxs);

/// @brief resets default Contexts for all statistics to their initial values
void reset_default_context();

/// @brief changes default contexts
/// @param ctx Context objects for one or more statistics
/// @return a move-only ContextResetter object whose destruction will reset the
/// default context to the previous value
/// @example
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
set_scoped_default_context(const Context& ctx);

/// @brief changes default context for arbitrary statistics
/// @note equivalent to `set_scoped_default_context({{Statistics::Arbitrary,
/// std::move(ctx)}})`
[[nodiscard]] detail::ImplicitContextResetter<
    container::map<Statistics, Context>>
set_scoped_default_context(Context&& ctx);

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
