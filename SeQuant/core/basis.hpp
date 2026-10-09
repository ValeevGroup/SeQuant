#ifndef SEQUANT_CORE_BASIS_HPP
#define SEQUANT_CORE_BASIS_HPP

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
/// base_key(), which is the space's key for such an instance, and so
/// sizes, tiles and slices it as the space. A basis of a different extent (a
/// truncated or an overcomplete set, such as the PAOs) must be registered
/// under a name, which gives it a key and an extent of its own.
class IndexBasis {
 public:
  using instance_type = std::int32_t;
  using optional_instance = std::optional<instance_type>;

  /// null space, null instance
  IndexBasis() noexcept = default;

  /// the space's own basis (null instance)
  IndexBasis(IndexSpace space) noexcept : space_(std::move(space)) {}

  /// the own basis of the IndexSpace constructed from @p type_label and
  /// @p args, so that a braced IndexSpace initializer also initializes an
  /// IndexBasis
  template <basic_string_convertible S, typename... Args>
    requires(sizeof...(Args) > 0 &&
             std::constructible_from<IndexSpace, S, Args...>)
  IndexBasis(S&& type_label, Args&&... args)
      : space_(std::forward<S>(type_label), std::forward<Args>(args)...) {}

  explicit IndexBasis(IndexSpace space,
                      optional_instance basis_instance) noexcept
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

  /// @return name() if the basis carries a name, else space().base_key();
  /// the unnamed instances of a space share its key
  const std::wstring& base_key() const noexcept {
    return has_name() ? name_ : space_.base_key();
  }

  /// @return the basis instance unless the basis carries a name, which then
  /// stands for it
  optional_instance unnamed_instance() const noexcept {
    return has_name() ? std::nullopt : basis_instance_;
  }

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

}  // namespace sequant

#endif  // SEQUANT_CORE_BASIS_HPP
