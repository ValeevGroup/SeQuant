#ifndef SEQUANT_CORE_BASIS_HPP
#define SEQUANT_CORE_BASIS_HPP

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/hash.hpp>
#include <SeQuant/core/space.hpp>

#include <compare>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <utility>

namespace sequant {

/// @brief the basis an Index runs over: an IndexSpace plus an optional basis
/// instance

/// The instance is an opaque integer that distinguishes several bases of the
/// same IndexSpace that meet in one expression, e.g. an eigenbasis and a
/// localized basis of the space, or bases specific to two different tensors.
/// A null instance (the default) is the IndexSpace's own basis; every integer,
/// 0 and negative ones included, is an ordinary instance distinct from null.
///
/// A basis instance may carry the name it is registered under in an
/// IndexBasisRegistry (see IndexBasisRegistry::add(label, space, instance,
/// spec)); a basis
/// obtained from the registry carries it, one built from a space and an
/// instance does not. The name is part of the basis's identity, as a space's
/// label is of the space's: equality, ordering and hashing see it, so a named
/// basis and the unnamed basis of the same space and instance are two bases.
/// The registry keeps an instance number and a name one-to-one within a
/// space (see IndexBasisRegistry::resolve()), so an Index given an instance
/// by number resolves it to the named basis; one given an IndexBasis takes it
/// as given. Such an unnamed basis has no spelling of its own: an index in it
/// prints with the bare instance number, which parses back to the named
/// basis.
///
/// A basis has an extent (the number of functions in it), a metric (whether it
/// is orthonormal) and a scalar field. Unlike the name, these describe the
/// basis but are not part of its identity, so equality, ordering and hashing
/// ignore them. An unnamed basis, whether the space's own basis or an unnamed
/// instance, is the space's own basis up to an orthonormal rotation (e.g. a
/// localized basis): its extent is the dimension of the space, its metric is
/// unit and its field is that of the space's own basis (IndexSpace::field()),
/// and evaluation keys its axes by base_key(), which is the space's key for
/// such an instance, and so sizes, tiles and slices it as the space. A basis
/// that differs in any of these (a truncated or an overcomplete set, a
/// non-orthonormal basis) must be registered under a name, which gives it a
/// key and metadata of its own (see IndexBasisRegistry::add(label, space,
/// instance, spec)).
/// That holds for each subspace of a basis instance: a subspace shares the
/// instance of its space but not its metric (a non-orthonormal basis can have
/// an orthonormal part), so an unregistered subspace instance is orthonormal
/// even if the instance of its space is not.
/// A named basis registered without an extent or a field follows its space's
/// (own_extent() and own_field() are null), and so sees the registry's
/// setters for the space's label; one registered with them, even with the
/// values the space has, does not.
class IndexBasis {
 public:
  using instance_type = std::int32_t;
  using optional_instance = std::optional<instance_type>;

  /// null space, null instance
  IndexBasis() noexcept = default;

  /// the space's own basis (null instance)
  IndexBasis(IndexSpace space) noexcept : space_(std::move(space)) {}

  explicit IndexBasis(IndexSpace space,
                      optional_instance basis_instance) noexcept
      : space_(std::move(space)), basis_instance_(basis_instance) {}

  /// the basis instance @p basis_instance of @p space under the name
  /// @p name, with the metadata of the space's own basis; the entry an
  /// IndexBasisRegistry keeps under that name carries the metadata the
  /// basis was registered with, so compare with it, do not stand in for it
  /// @param name the label @p basis_instance is registered under
  /// @pre @p name is empty or @p basis_instance is non-null
  IndexBasis(IndexSpace space, optional_instance basis_instance,
             std::string name)
      : IndexBasis(std::move(space), basis_instance, std::move(name), {},
                   IndexBasisMetric::Unit, {}) {}

  const IndexSpace& space() const noexcept { return space_; }

  const optional_instance& basis_instance() const noexcept {
    return basis_instance_;
  }

  bool has_basis_instance() const noexcept {
    return basis_instance_.has_value();
  }

  /// @return the name the basis instance is registered under, empty if it
  /// carries none
  const std::string& name() const noexcept { return name_; }

  bool has_name() const noexcept { return !name_.empty(); }

  /// @return name() if the basis carries a name, else space().base_key();
  /// the unnamed instances of a space share its key
  const std::string& base_key() const noexcept {
    return has_name() ? name_ : space_.base_key();
  }

  /// @return the basis instance unless the basis carries a name, which then
  /// stands for it
  optional_instance unnamed_instance() const noexcept {
    return has_name() ? std::nullopt : basis_instance_;
  }

  /// @return the number of functions in the basis: the dimension of its space
  /// (IndexSpace::dimension()) unless registered with an extent of its
  /// own
  std::size_t extent() const noexcept {
    return extent_ ? *extent_ : space_.dimension();
  }

  /// @return whether the basis is orthonormal (IndexBasisMetric::Unit) or not;
  /// unit unless registered otherwise
  IndexBasisMetric metric() const noexcept { return metric_; }

  /// @return the scalar field of the basis: that of its space's own basis
  /// (IndexSpace::field()) unless registered with a field of its own
  Field field() const noexcept { return field_ ? *field_ : space_.field(); }

  /// @return the extent the basis was registered with (or mirrors from the
  /// entry it follows, see IndexBasisRegistry::follow()), null if extent()
  /// follows the dimension of its space (which the registry's
  /// IndexBasisRegistry::extent(label, n) for the space's label then sets)
  const std::optional<std::size_t>& own_extent() const noexcept {
    return extent_;
  }

  /// @return the field the basis was registered with (or mirrors from the
  /// entry it follows, see IndexBasisRegistry::follow()), null if field()
  /// follows that of its space's own basis (which the registry's
  /// IndexBasisRegistry::field(label, f) for the space's label then sets)
  const std::optional<Field>& own_field() const noexcept { return field_; }

  /// @return `";N"` for a non-null instance `N`, else an empty string
  std::string instance_suffix() const;

  /// compares space, instance and name; the extent, metric and field are
  /// ignored
  friend bool operator==(const IndexBasis& b1, const IndexBasis& b2) noexcept {
    return b1.space_ == b2.space_ && b1.basis_instance_ == b2.basis_instance_ &&
           b1.name_ == b2.name_;
  }

  /// orders by space, then by instance (null first), then by name (unnamed
  /// first); the extent, metric and field are ignored
  friend std::strong_ordering operator<=>(const IndexBasis& b1,
                                          const IndexBasis& b2) noexcept {
    if (auto c = b1.space_ <=> b2.space_; c != 0) return c;
    if (auto c = b1.basis_instance_ <=> b2.basis_instance_; c != 0) return c;
    return b1.name_ <=> b2.name_;
  }

  /// @return `hash_value(b.space())` if @p b has no instance; the extent,
  /// metric and field are ignored
  friend std::size_t hash_value(const IndexBasis& b) {
    std::size_t result = hash_value(b.space_);
    if (b.basis_instance_) hash::combine(result, *b.basis_instance_);
    if (b.has_name()) hash::combine(result, b.name_);
    return result;
  }

  /// @name registry keys
  /// An IndexBasisRegistry keys its entries, spaces and named basis instances
  /// alike, by the base key of a label (an Index label with its ordinal
  /// stripped); a name follows the grammar of a space's label, so a basis key
  /// is reduced as a space key is.
  /// @{

  /// @param key a label of a space, of a named basis instance, or of an
  /// Index in either
  /// @return the base key of @p key, as `IndexSpace::reduce_key(key)`
  static std::string_view reduce_key(std::string_view key) {
    return IndexSpace::reduce_key(key);
  }

  /// exception type thrown by a registry lookup when no entry, space or
  /// named basis instance, is registered under a key; an IndexSpace::bad_key,
  /// so that handlers for space lookups catch it
  struct bad_key : IndexSpace::bad_key {
    using IndexSpace::bad_key::bad_key;
  };

  /// @}

 private:
  IndexSpace space_;
  optional_instance basis_instance_;
  std::string name_;
  std::optional<std::size_t> extent_;
  IndexBasisMetric metric_ = IndexBasisMetric::Unit;
  std::optional<Field> field_;

  // the registry builds its named entries with their metadata, and writes
  // the metadata in place
  friend class IndexBasisRegistry;
  /// @param name the label @p basis_instance is registered under
  /// @param extent the number of functions in the basis; null for the
  ///        dimension of @p space
  /// @param metric whether the basis is orthonormal
  /// @param field the scalar field of the basis; null for that of the own
  ///        basis of @p space
  /// @pre @p name is empty or @p basis_instance is non-null; @p name is
  ///      non-empty or every other argument is at its default (an unnamed
  ///      basis is the space's own basis up to an orthonormal rotation)
  IndexBasis(IndexSpace space, optional_instance basis_instance,
             std::string name, std::optional<std::size_t> extent,
             IndexBasisMetric metric, std::optional<Field> field);
  void extent(std::size_t n) noexcept { extent_ = n; }
  void metric(IndexBasisMetric m) noexcept { metric_ = m; }
  void field(Field f) noexcept { field_ = f; }
};

/// @return true if @p b1 and @p b2 are of one space and instance, named or
/// not: the relation the registry keeps one-to-one with a name (see
/// IndexBasisRegistry::resolve())
bool same_instance(const IndexBasis& b1, const IndexBasis& b2) noexcept;

/// @return true if @p basis includes @p subbasis, i.e. its space includes the
/// space of @p subbasis and it is either the space's own basis (which includes
/// every instance of the space) or the same instance, under the same name, as
/// @p subbasis
bool includes(const IndexBasis& basis, const IndexBasis& subbasis);

/// @return true if @p b1 and @p b2 are different basis instances, i.e. both
/// have one and the two differ in instance or name; an identity between
/// functions of such bases is an overlap, not a Kronecker delta
bool different_instances(const IndexBasis& b1, const IndexBasis& b2);

}  // namespace sequant

#endif  // SEQUANT_CORE_BASIS_HPP
