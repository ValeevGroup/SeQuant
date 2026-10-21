#ifndef SEQUANT_CORE_IRREPS_Z2_ALGEBRA_HPP
#define SEQUANT_CORE_IRREPS_Z2_ALGEBRA_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/irreps/algebra.hpp>

#include <algorithm>
#include <array>
#include <compare>
#include <concepts>
#include <cstdint>
#include <format>
#include <iterator>
#include <string>
#include <string_view>

namespace sequant {

/// Concept for compile-time descriptions of abelian groups in which every
/// element is its own inverse (D2h and its subgroups). The IRREP names are
/// listed in an order such that the index of an IRREP interpreted as a
/// bitmask yields the group multiplication via XOR. The totally symmetric
/// IRREP is the first one.
template <typename G>
concept Z2Group =
    requires {
      { G::names.size() } -> std::convertible_to<std::size_t>;
      { G::names[0] } -> std::convertible_to<std::string_view>;
    } && (std::size(G::names) > 0) && (std::size(G::names) <= 256) &&
    ((std::size(G::names) & (std::size(G::names) - 1)) == 0);

/// Group without any concrete IRREPs (only the totally symmetric one). Using
/// this group, all non-trivial IRREPs have to be symbolic.
struct NoConcreteGroup {
  static constexpr std::array<std::string_view, 1> names{"1"};
};

/// Abelian point groups using Molpro's IRREP ordering
namespace point_groups {

struct C1 {
  static constexpr std::array<std::string_view, 1> names{"A"};
};
struct Ci {
  static constexpr std::array<std::string_view, 2> names{"Ag", "Au"};
};
struct C2 {
  static constexpr std::array<std::string_view, 2> names{"A", "B"};
};
struct Cs {
  static constexpr std::array<std::string_view, 2> names{"A'", "A''"};
};
struct C2v {
  static constexpr std::array<std::string_view, 4> names{"A1", "B1", "B2",
                                                         "A2"};
};
struct C2h {
  static constexpr std::array<std::string_view, 4> names{"Ag", "Au", "Bu",
                                                         "Bg"};
};
struct D2 {
  static constexpr std::array<std::string_view, 4> names{"A", "B3", "B2", "B1"};
};
struct D2h {
  static constexpr std::array<std::string_view, 8> names{
      "Ag", "B3u", "B2u", "B1g", "B1u", "B2g", "B3g", "Au"};
};

}  // namespace point_groups

/// An IRREP of a Z2 group (in which every element is its own inverse) that
/// may contain symbolic (unevaluated) factors. It is represented as the
/// product of a concrete IRREP of Group and a set of symbols, each of which
/// represents an IRREP that is unknown at the time of code generation.
template <Z2Group Group = NoConcreteGroup>
struct Z2Irrep {
  /// Index of the concrete IRREP in Group::names (interpreted as bitmask)
  std::uint8_t concrete = 0;
  /// Sorted list of (unique) symbols
  container::svector<std::string, 2> symbols;

  friend bool operator==(const Z2Irrep &, const Z2Irrep &) = default;

  friend std::strong_ordering operator<=>(const Z2Irrep &lhs,
                                          const Z2Irrep &rhs) {
    if (auto cmp = lhs.concrete <=> rhs.concrete; cmp != 0) {
      return cmp;
    }

    return std::lexicographical_compare_three_way(
        lhs.symbols.begin(), lhs.symbols.end(), rhs.symbols.begin(),
        rhs.symbols.end());
  }
};

/// IRREP algebra for Z2 groups (D2h and its subgroups) that supports symbolic
/// IRREPs. Products are formed by XOR-ing the concrete parts and forming the
/// symmetric difference of the symbols (as every IRREP is its own inverse).
template <Z2Group Group = NoConcreteGroup>
class Z2Algebra {
 public:
  using group_type = Group;
  using irrep_type = Z2Irrep<Group>;

  irrep_type totally_symmetric() const { return {}; }

  IrrepSet<irrep_type> product(const irrep_type &lhs,
                               const irrep_type &rhs) const {
    irrep_type result;
    result.concrete = lhs.concrete ^ rhs.concrete;
    std::ranges::set_symmetric_difference(lhs.symbols, rhs.symbols,
                                          std::back_inserter(result.symbols));

    return {std::move(result)};
  }

  irrep_type conjugate(const irrep_type &irrep) const { return irrep; }
};

/// Creates an IRREP from its name. The name may be a product of several
/// factors separated by '*'. Every factor that is the name of one of Group's
/// IRREPs refers to that concrete IRREP. All other factors are taken to be
/// symbols.
template <Z2Group Group = NoConcreteGroup>
Z2Irrep<Group> make_irrep(std::string_view name) {
  Z2Algebra<Group> algebra;
  Z2Irrep<Group> irrep;

  std::size_t start = 0;
  while (start <= name.size()) {
    std::size_t end = name.find('*', start);
    if (end == std::string_view::npos) {
      end = name.size();
    }
    std::string_view factor_name = name.substr(start, end - start);
    start = end + 1;

    if (factor_name.empty()) {
      continue;
    }

    Z2Irrep<Group> factor;
    auto it = std::ranges::find(Group::names, factor_name);
    if (it != Group::names.end()) {
      factor.concrete =
          static_cast<std::uint8_t>(std::distance(Group::names.begin(), it));
    } else {
      factor.symbols.emplace_back(factor_name);
    }

    irrep = *algebra.product(irrep, factor).begin();
  }

  return irrep;
}

}  // namespace sequant

/// Formats a Z2Irrep as the '*'-separated list of its factors. The concrete
/// factor is omitted if it is the totally symmetric IRREP and there is at
/// least one symbolic factor.
template <sequant::Z2Group Group>
struct std::formatter<sequant::Z2Irrep<Group>> {
  constexpr auto parse(std::format_parse_context &ctx) { return ctx.begin(); }

  auto format(const sequant::Z2Irrep<Group> &irrep,
              std::format_context &ctx) const {
    auto out = ctx.out();
    bool first = true;
    if (irrep.concrete != 0 || irrep.symbols.empty()) {
      out = std::format_to(out, "{}", Group::names[irrep.concrete]);
      first = false;
    }
    for (const std::string &symbol : irrep.symbols) {
      if (!first) {
        *out++ = '*';
      }
      first = false;
      out = std::format_to(out, "{}", symbol);
    }

    return out;
  }
};

static_assert(sequant::IrrepAlgebra<sequant::Z2Algebra<>>);

#endif  // SEQUANT_CORE_IRREPS_Z2_ALGEBRA_HPP
