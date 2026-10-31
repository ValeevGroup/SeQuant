#ifndef SEQUANT_CORE_IRREPS_Z2_ALGEBRA_HPP
#define SEQUANT_CORE_IRREPS_Z2_ALGEBRA_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/irreps/algebra.hpp>
#include <SeQuant/core/utility/exception.hpp>

#include <algorithm>
#include <array>
#include <bit>
#include <compare>
#include <concepts>
#include <cstdint>
#include <format>
#include <iterator>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

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

template <typename Group = NoConcreteGroup>
  requires Z2Group<Group>
class Z2Algebra;

/// An IRREP of a Z2 group (in which every element is its own inverse) that
/// may contain symbolic (unevaluated) factors. It is represented as the
/// product of a concrete IRREP of Group and a set of symbols, each of which
/// represents an IRREP that is unknown at the time of code generation.
template <Z2Group Group = NoConcreteGroup>
class Z2Irrep {
 public:
  using symbol_container = container::svector<std::string, 2>;

  /// Creates the totally symmetric IRREP
  Z2Irrep() = default;

  /// @param concrete Index of the concrete IRREP in Group::names
  /// @param symbols The symbolic factors (in any order). As every IRREP is
  ///        its own inverse, pairs of identical symbols cancel.
  /// @throws Exception if concrete is out of range or if a symbol is empty,
  ///         contains a '*' or is the name of one of Group's IRREPs
  explicit Z2Irrep(std::uint8_t concrete, symbol_container symbols = {})
      : m_concrete(concrete) {
    if (concrete >= Group::names.size()) {
      throw Exception("Concrete IRREP index " + std::to_string(concrete) +
                      " is out of range");
    }

    std::ranges::sort(symbols);
    for (std::string &symbol : symbols) {
      if (symbol.empty() || symbol.find('*') != std::string::npos ||
          std::ranges::find(Group::names, symbol) != Group::names.end()) {
        throw Exception("Invalid IRREP symbol '" + symbol + "'");
      }

      if (!m_symbols.empty() && m_symbols.back() == symbol) {
        m_symbols.pop_back();
      } else {
        m_symbols.push_back(std::move(symbol));
      }
    }
  }

  /// @returns Index of the concrete IRREP in Group::names (interpreted as
  /// bitmask)
  std::uint8_t concrete() const { return m_concrete; }
  /// @returns The sorted list of (unique) symbols
  const symbol_container &symbols() const { return m_symbols; }

  friend bool operator==(const Z2Irrep &, const Z2Irrep &) = default;

  friend std::strong_ordering operator<=>(const Z2Irrep &lhs,
                                          const Z2Irrep &rhs) {
    if (auto cmp = lhs.m_concrete <=> rhs.m_concrete; cmp != 0) {
      return cmp;
    }

    return std::lexicographical_compare_three_way(
        lhs.m_symbols.begin(), lhs.m_symbols.end(), rhs.m_symbols.begin(),
        rhs.m_symbols.end());
  }

 private:
  std::uint8_t m_concrete = 0;
  symbol_container m_symbols;

  friend class Z2Algebra<Group>;
};

/// IRREP algebra for Z2 groups (D2h and its subgroups) that supports symbolic
/// IRREPs. Products are formed by XOR-ing the concrete parts and forming the
/// symmetric difference of the symbols (as every IRREP is its own inverse).
///
/// Optionally, relations between symbolic (and concrete) IRREPs can be
/// specified (e.g. My1*My2 = My3). Without relations, symbols are treated as
/// independent of each other and products of symbols remain unevaluated.
/// With relations, every IRREP is reduced to the representation with the
/// fewest symbols among all representations that are equivalent under the
/// given relations (ties broken by the IRREP ordering). This makes the
/// representation unique and e.g. turns My1*My2 into My3. Note that the cost
/// of this reduction grows exponentially with the number of independent
/// relations.
template <typename Group>
  requires Z2Group<Group>
class Z2Algebra {
 public:
  using group_type = Group;
  using irrep_type = Z2Irrep<Group>;

  /// The maximum number of independent relations
  static constexpr std::size_t max_relations = 16;

  /// Creates an algebra without any relations between symbols
  Z2Algebra() = default;

  /// @param identities Relations of the form lhs = rhs between IRREPs
  /// @throws Exception if the relations imply an identity between different
  ///         concrete IRREPs or if there are more than max_relations
  ///         independent relations
  explicit Z2Algebra(
      const std::vector<std::pair<irrep_type, irrep_type>> &identities) {
    // The span of all relations (in the form x = totally symmetric)
    std::vector<irrep_type> span{totally_symmetric()};

    for (const auto &[lhs, rhs] : identities) {
      irrep_type relation = multiply(lhs, rhs);
      if (std::ranges::find(span, relation) != span.end()) {
        // Trivial or implied by previous relations
        continue;
      }

      if (m_relations.size() >= max_relations) {
        throw Exception("Too many independent IRREP relations (max. " +
                        std::to_string(max_relations) + ")");
      }

      const std::size_t size = span.size();
      for (std::size_t i = 0; i < size; ++i) {
        span.push_back(multiply(span[i], relation));
      }
      m_relations.push_back(std::move(relation));
    }

    for (const irrep_type &current : span) {
      if (current.symbols().empty() && current.concrete() != 0) {
        throw Exception(
            "IRREP relations imply that two different concrete IRREPs are "
            "identical");
      }
    }
  }

  irrep_type totally_symmetric() const { return {}; }

  IrrepSet<irrep_type> product(const irrep_type &lhs,
                               const irrep_type &rhs) const {
    return {reduce(multiply(lhs, rhs))};
  }

  irrep_type conjugate(const irrep_type &irrep) const { return irrep; }

  /// @returns The canonical representation of the given IRREP under this
  /// algebra's relations
  irrep_type reduce(const irrep_type &irrep) const {
    irrep_type best = irrep;
    irrep_type candidate = irrep;

    // Visit all elements of the coset of irrep w.r.t. the span of the
    // relations in Gray code order (only one relation is applied per step)
    const std::size_t n_elements = std::size_t{1} << m_relations.size();
    for (std::size_t i = 1; i < n_elements; ++i) {
      candidate = multiply(candidate, m_relations[std::countr_zero(i)]);

      if (candidate.symbols().size() < best.symbols().size() ||
          (candidate.symbols().size() == best.symbols().size() &&
           candidate < best)) {
        best = candidate;
      }
    }

    return best;
  }

 private:
  /// Independent relations, each of which is equivalent to the totally
  /// symmetric IRREP
  std::vector<irrep_type> m_relations;

  /// @returns The plain product of the given IRREPs (not reduced)
  static irrep_type multiply(const irrep_type &lhs, const irrep_type &rhs) {
    irrep_type result;
    result.m_concrete = lhs.m_concrete ^ rhs.m_concrete;
    std::ranges::set_symmetric_difference(lhs.m_symbols, rhs.m_symbols,
                                          std::back_inserter(result.m_symbols));

    return result;
  }
};

/// Creates an IRREP from its name. The name may be a product of several
/// factors separated by '*'. Every factor that is the name of one of Group's
/// IRREPs refers to that concrete IRREP. All other factors are taken to be
/// symbols.
template <Z2Group Group = NoConcreteGroup>
Z2Irrep<Group> make_irrep(std::string_view name) {
  std::uint8_t concrete = 0;
  typename Z2Irrep<Group>::symbol_container symbols;

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

    auto it = std::ranges::find(Group::names, factor_name);
    if (it != Group::names.end()) {
      concrete ^=
          static_cast<std::uint8_t>(std::distance(Group::names.begin(), it));
    } else {
      symbols.emplace_back(factor_name);
    }
  }

  return Z2Irrep<Group>(concrete, std::move(symbols));
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
    if (irrep.concrete() != 0 || irrep.symbols().empty()) {
      out = std::format_to(out, "{}", Group::names[irrep.concrete()]);
      first = false;
    }
    for (const std::string &symbol : irrep.symbols()) {
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
