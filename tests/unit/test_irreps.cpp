#include <catch2/catch_test_macros.hpp>

#include "catch2_sequant.hpp"

#include <SeQuant/core/irreps/algebra.hpp>
#include <SeQuant/core/utility/exception.hpp>

#include <compare>
#include <format>
#include <string>

namespace {

/// IRREPs of the cyclic group Z3 (which has complex IRREPs)
struct Z3Irrep {
  int value = 0;

  auto operator<=>(const Z3Irrep &) const = default;
};

struct Z3Algebra {
  using irrep_type = Z3Irrep;

  irrep_type totally_symmetric() const { return {0}; }
  sequant::IrrepSet<irrep_type> product(const irrep_type &lhs,
                                        const irrep_type &rhs) const {
    return {Z3Irrep{(lhs.value + rhs.value) % 3}};
  }
  irrep_type conjugate(const irrep_type &irrep) const {
    return {(3 - irrep.value) % 3};
  }
};

}  // namespace

template <>
struct std::formatter<Z3Irrep> {
  constexpr auto parse(std::format_parse_context &ctx) { return ctx.begin(); }
  auto format(const Z3Irrep &irrep, std::format_context &ctx) const {
    return std::format_to(ctx.out(), "G{}", irrep.value);
  }
};

static_assert(sequant::IrrepAlgebra<Z3Algebra>);
TEST_CASE("irreps", "[irreps]") {
  using namespace sequant;

  SECTION("IrrepSet") {
    using Set = IrrepSet<std::string>;

    Set set{"b", "a", "b"};
    REQUIRE(set.size() == 2);
    REQUIRE(to_string(set) == "{a,b}");
    REQUIRE(set.contains("a"));
    REQUIRE(!set.contains("c"));
    REQUIRE(!set.empty());
    REQUIRE(!set.is_any());

    REQUIRE(Set{}.empty());
    REQUIRE(to_string(Set{}) == "{}");

    Set any = Set::any();
    REQUIRE(any.is_any());
    REQUIRE(!any.empty());
    REQUIRE(any.contains("whatever"));
    REQUIRE(to_string(any) == "*");

    Set united = set;
    united.unite(Set{"c"});
    REQUIRE(united == Set{"a", "b", "c"});
    united.unite(any);
    REQUIRE(united == any);
    REQUIRE(Set{"a"}.unite(Set{}) == Set{"a"});

    REQUIRE(Set{"a"} < Set{"b"});
    REQUIRE(Set{"z"} < any);
  }

  SECTION("custom algebra") {
    Z3Algebra alg;
    using Set = IrrepSet<Z3Irrep>;

    REQUIRE(irrep_to_string(Z3Irrep{2}) == "G2");
    REQUIRE(to_string(product(alg, Set{Z3Irrep{1}}, Set{Z3Irrep{1}})) ==
            "{G2}");
    REQUIRE(conjugate(alg, Set{Z3Irrep{1}, Z3Irrep{0}}) ==
            Set{Z3Irrep{0}, Z3Irrep{2}});
  }
}
