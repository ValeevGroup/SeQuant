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
class IndexBasis {
 public:
  using instance_type = std::int32_t;
  using optional_instance = std::optional<instance_type>;

  /// null space, null instance
  IndexBasis() noexcept = default;

  explicit IndexBasis(IndexSpace space,
                      optional_instance basis_instance = std::nullopt) noexcept
      : space_(std::move(space)), basis_instance_(basis_instance) {}

  const IndexSpace& space() const noexcept { return space_; }

  const optional_instance& basis_instance() const noexcept {
    return basis_instance_;
  }

  bool has_basis_instance() const noexcept {
    return basis_instance_.has_value();
  }

  /// @return `L";N"` for a non-null instance `N`, else an empty string
  std::wstring instance_suffix() const;

  friend bool operator==(const IndexBasis&,
                         const IndexBasis&) noexcept = default;

  /// orders by space, then by instance (null first)
  friend std::strong_ordering operator<=>(const IndexBasis&,
                                          const IndexBasis&) noexcept = default;

  /// @return `hash_value(b.space())` if @p b has no instance
  friend std::size_t hash_value(const IndexBasis& b) {
    std::size_t result = hash_value(b.space_);
    if (b.basis_instance_) hash::combine(result, *b.basis_instance_);
    return result;
  }

 private:
  IndexSpace space_;
  optional_instance basis_instance_;
};

/// @return true if @p basis includes @p subbasis, i.e. its space includes the
/// space of @p subbasis and it is either the space's own basis (which includes
/// every instance of the space) or the same instance as @p subbasis
bool includes(const IndexBasis& basis, const IndexBasis& subbasis);

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
