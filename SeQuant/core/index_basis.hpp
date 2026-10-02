#ifndef SEQUANT_CORE_INDEX_BASIS_HPP
#define SEQUANT_CORE_INDEX_BASIS_HPP

#include <SeQuant/core/hash.hpp>
#include <SeQuant/core/space.hpp>

#include <compare>
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
  std::wstring instance_suffix() const {
    return basis_instance_ ? L";" + std::to_wstring(*basis_instance_)
                           : std::wstring{};
  }

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
inline bool includes(const IndexBasis& basis, const IndexBasis& subbasis) {
  return includes(basis.space(), subbasis.space()) &&
         (!basis.has_basis_instance() ||
          basis.basis_instance() == subbasis.basis_instance());
}

}  // namespace sequant

#endif  // SEQUANT_CORE_INDEX_BASIS_HPP
