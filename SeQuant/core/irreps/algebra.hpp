#ifndef SEQUANT_CORE_IRREPS_ALGEBRA_HPP
#define SEQUANT_CORE_IRREPS_ALGEBRA_HPP

#include <SeQuant/core/container.hpp>

#include <algorithm>
#include <compare>
#include <concepts>
#include <format>
#include <initializer_list>
#include <string>
#include <utility>

namespace sequant {

/// Converts the given IRREP into its string representation. This is the
/// customization point for how IRREPs are spelled out (e.g. in generated
/// code). By default, it uses std::format, so the preferred way of
/// customizing it is to provide a std::formatter specialization for the IRREP
/// type. Alternatively, this function template may be fully specialized for
/// a given IRREP type.
template <typename Irrep>
std::string irrep_to_string(const Irrep &irrep) {
  return std::format("{}", irrep);
}

/// A set of IRREPs, i.e. the representation of a quantity spanning several
/// IRREPs (as a direct sum). Besides holding a concrete (possibly empty) list
/// of IRREPs, a set may also be in the special "any" state that represents
/// a quantity of unknown symmetry (spanning all IRREPs of the group). The
/// "any" state is absorbing under unification and under multiplication with
/// a non-empty set. An empty set is the zero of the multiplication.
template <typename Irrep>
class IrrepSet {
 public:
  using value_type = Irrep;
  using container_type = container::svector<Irrep, 2>;
  using const_iterator = typename container_type::const_iterator;

  /// Creates an empty set
  IrrepSet() = default;
  IrrepSet(std::initializer_list<Irrep> irreps) {
    for (const Irrep &current : irreps) {
      insert(current);
    }
  }

  /// @returns A set spanning any (all) IRREPs
  static IrrepSet any() {
    IrrepSet set;
    set.m_any = true;
    return set;
  }

  [[nodiscard]] bool is_any() const { return m_any; }
  /// @returns Whether this set contains no IRREP at all. A set in the "any"
  /// state is not empty.
  [[nodiscard]] bool empty() const { return !m_any && m_irreps.empty(); }
  /// @returns The amount of explicitly stored IRREPs (zero in the "any" state)
  [[nodiscard]] std::size_t size() const { return m_irreps.size(); }

  /// @returns Whether the given IRREP is part of this set. Any IRREP is part
  /// of a set in the "any" state.
  [[nodiscard]] bool contains(const Irrep &irrep) const {
    return m_any || std::ranges::binary_search(m_irreps, irrep);
  }

  /// Adds the given IRREP to this set
  IrrepSet &insert(Irrep irrep) {
    if (m_any) {
      return *this;
    }

    auto it = std::ranges::lower_bound(m_irreps, irrep);
    if (it == m_irreps.end() || *it != irrep) {
      m_irreps.insert(it, std::move(irrep));
    }

    return *this;
  }

  /// Turns this set into the union of itself and the given set
  IrrepSet &unite(const IrrepSet &other) {
    if (other.m_any) {
      *this = any();
    } else {
      for (const Irrep &current : other.m_irreps) {
        insert(current);
      }
    }

    return *this;
  }

  /// Iteration over the explicitly stored IRREPs (in ascending order). A set
  /// in the "any" state doesn't store any IRREPs explicitly.
  const_iterator begin() const { return m_irreps.begin(); }
  const_iterator end() const { return m_irreps.end(); }

  friend bool operator==(const IrrepSet &lhs, const IrrepSet &rhs) {
    return lhs.m_any == rhs.m_any && lhs.m_irreps == rhs.m_irreps;
  }

  friend auto operator<=>(const IrrepSet &lhs, const IrrepSet &rhs) {
    using Ordering = std::compare_three_way_result_t<Irrep>;
    if (lhs.m_any != rhs.m_any) {
      return lhs.m_any ? Ordering(std::strong_ordering::greater)
                       : Ordering(std::strong_ordering::less);
    }

    return std::lexicographical_compare_three_way(
        lhs.m_irreps.begin(), lhs.m_irreps.end(), rhs.m_irreps.begin(),
        rhs.m_irreps.end());
  }

 private:
  container_type m_irreps;
  bool m_any = false;
};

/// @returns The string representation of the given set as "{a,b}" where
/// a and b are stringified via irrep_to_string. A set in the "any" state is
/// represented as "*".
template <typename Irrep>
std::string to_string(const IrrepSet<Irrep> &set) {
  if (set.is_any()) {
    return "*";
  }

  std::string str = "{";
  bool first = true;
  for (const Irrep &current : set) {
    if (!first) {
      str += ",";
    }
    first = false;
    str += irrep_to_string(current);
  }
  str += "}";

  return str;
}

/// Concept for types implementing the multiplication rules of IRREPs.
/// The IRREPs themselves are of the algebra's irrep_type which need not be
/// numeric (they could e.g. be unevaluated symbols). The product of two IRREPs
/// is a set of IRREPs in order to be able to represent the decomposition of
/// direct products for non-abelian groups.
template <typename A>
concept IrrepAlgebra =
    std::copyable<A> && std::totally_ordered<typename A::irrep_type> &&
    std::three_way_comparable<typename A::irrep_type> &&
    std::copyable<typename A::irrep_type> &&
    requires(const A &algebra, const typename A::irrep_type &irrep) {
      {
        algebra.totally_symmetric()
      } -> std::convertible_to<typename A::irrep_type>;
      {
        algebra.product(irrep, irrep)
      } -> std::convertible_to<IrrepSet<typename A::irrep_type>>;
      {
        algebra.conjugate(irrep)
      } -> std::convertible_to<typename A::irrep_type>;
      { irrep_to_string(irrep) } -> std::convertible_to<std::string>;
    };

/// @returns The product of the two given sets, which is the union of the
/// products of all pairs of IRREPs from the two sets
template <IrrepAlgebra A>
IrrepSet<typename A::irrep_type> product(
    const A &algebra, const IrrepSet<typename A::irrep_type> &lhs,
    const IrrepSet<typename A::irrep_type> &rhs) {
  using Set = IrrepSet<typename A::irrep_type>;

  if (lhs.empty() || rhs.empty()) {
    return {};
  }
  if (lhs.is_any() || rhs.is_any()) {
    return Set::any();
  }

  Set result;
  for (const auto &left : lhs) {
    for (const auto &right : rhs) {
      result.unite(algebra.product(left, right));
    }
  }

  return result;
}

/// @returns The set of complex conjugates of the IRREPs in the given set
template <IrrepAlgebra A>
IrrepSet<typename A::irrep_type> conjugate(
    const A &algebra, const IrrepSet<typename A::irrep_type> &set) {
  using Set = IrrepSet<typename A::irrep_type>;

  if (set.is_any()) {
    return Set::any();
  }

  Set result;
  for (const auto &current : set) {
    result.insert(algebra.conjugate(current));
  }

  return result;
}

}  // namespace sequant

#endif  // SEQUANT_CORE_IRREPS_ALGEBRA_HPP
