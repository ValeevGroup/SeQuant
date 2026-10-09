#ifndef SEQUANT_CORE_BASIS_HPP
#define SEQUANT_CORE_BASIS_HPP

#include <SeQuant/core/hash.hpp>
#include <SeQuant/core/space.hpp>

#include <compare>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <type_traits>
#include <utility>

namespace sequant {

/// @brief the basis an Index runs over: an IndexSpace plus an optional basis
/// instance

/// The instance is an opaque integer that distinguishes several bases of the
/// same IndexSpace that meet in one expression, e.g. the canonical and the
/// localized orbitals of a perturbation theory, or the pair-specific virtuals
/// of two amplitudes. A null instance (the default) is the IndexSpace's own
/// basis; every integer, 0 and negative ones included, is an ordinary instance
/// distinct from null.
///
/// A basis instance may carry the name it is registered under in an
/// IndexBasisRegistry (see IndexBasisRegistry::add(label, basis)); a basis
/// obtained from the registry carries it, one built from a space and an
/// instance does not. The name is how the basis prints; it is not part of its
/// identity, so equality, ordering and hashing ignore it.
///
/// An unnamed instance spans its space at the space's extent, as a rotation of
/// it does (e.g. the cluster-specific virtuals): evaluation keys its axes by
/// Index::basis_key(), which is the space's key for such an instance, and so
/// sizes, tiles and slices it as the space. A basis of a different extent (a
/// truncated or an overcomplete set, such as the PAOs) must be registered
/// under a name, which gives it a key and an extent of its own.
class IndexBasis {
 public:
  using instance_type = std::int32_t;
  using optional_instance = std::optional<instance_type>;

  /// null space, null instance
  IndexBasis() noexcept = default;

  explicit IndexBasis(IndexSpace space,
                      optional_instance basis_instance = std::nullopt) noexcept
      : space_(std::move(space)), basis_instance_(basis_instance) {}

  /// @param name the label @p basis_instance is registered under
  /// @pre @p name is empty or @p basis_instance is non-null
  IndexBasis(IndexSpace space, optional_instance basis_instance,
             std::wstring name);

  const IndexSpace& space() const noexcept { return space_; }

  const optional_instance& basis_instance() const noexcept {
    return basis_instance_;
  }

  bool has_basis_instance() const noexcept {
    return basis_instance_.has_value();
  }

  /// @return the name the basis instance is registered under, empty if it
  /// carries none
  const std::wstring& name() const noexcept { return name_; }

  bool has_name() const noexcept { return !name_.empty(); }

  /// @return `L";N"` for a non-null instance `N`, else an empty string
  std::wstring instance_suffix() const;

  /// compares space and instance; the name is ignored
  friend bool operator==(const IndexBasis& b1, const IndexBasis& b2) noexcept {
    return b1.space_ == b2.space_ && b1.basis_instance_ == b2.basis_instance_;
  }

  /// orders by space, then by instance (null first); the name is ignored
  friend std::strong_ordering operator<=>(const IndexBasis& b1,
                                          const IndexBasis& b2) noexcept {
    if (auto c = b1.space_ <=> b2.space_; c != 0) return c;
    return b1.basis_instance_ <=> b2.basis_instance_;
  }

  /// @return `hash_value(b.space())` if @p b has no instance; the name is
  /// ignored
  friend std::size_t hash_value(const IndexBasis& b) {
    std::size_t result = hash_value(b.space_);
    if (b.basis_instance_) hash::combine(result, *b.basis_instance_);
    return result;
  }

 private:
  IndexSpace space_;
  optional_instance basis_instance_;
  std::wstring name_;
};

/// @return true if @p basis includes @p subbasis, i.e. its space includes the
/// space of @p subbasis and it is either the space's own basis (which includes
/// every instance of the space) or the same instance as @p subbasis
bool includes(const IndexBasis& basis, const IndexBasis& subbasis);

/// @return true if @p b1 and @p b2 are different basis instances, i.e. both
/// have one and the two differ; an identity between functions of such bases
/// is an overlap, not a Kronecker delta
bool different_instances(const IndexBasis& b1, const IndexBasis& b2);

/// what an Index runs over: an IndexSpace or an IndexBasis. An IndexSpace
/// stands for the space's own basis, so an Index made from one is
/// basis-generic (null basis instance); an IndexBasis names the basis exactly.
/// Functions templated on this default it to IndexSpace, so a braced list in
/// that position still initializes an IndexSpace
template <typename T>
concept space_or_basis = std::is_same_v<std::remove_cvref_t<T>, IndexSpace> ||
                         std::is_same_v<std::remove_cvref_t<T>, IndexBasis>;

}  // namespace sequant

#endif  // SEQUANT_CORE_BASIS_HPP
