#include <catch2/catch_test_macros.hpp>

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/optimize/fusion.hpp>
#include <SeQuant/core/parse_expr.hpp>

#include <array>
#include <memory>
#include <string_view>
#include <vector>

TEST_CASE("TEST_FUSION", "[Fusion]") {
  using sequant::opt::Fusion;
  using namespace sequant;
  std::vector<std::array<std::wstring_view, 3>> fused_terms{
      {
          L"1/2 f{i3;i1}          t{a1,a2;i2,i3}",             // lhs
          L"1/2 f{i3;a3} t{a3;i1} t{a1,a2;i2,i3}",             // rhs
          L"1/2(f{i3;i1} + f{i3;a3} t{a3;i1}) t{a1,a2;i2,i3}"  // fused form
      },

      {L"1/8 g{a1,a2;a3,a4} t{a3,a4;i1,i2}",
       L"1/4 g{a1,a2;a3,a4} t{a3;i1} t{a4;i2}",
       L"g{a1,a2;a3,a4}(1/8 t{a3,a4;i1,i2} + 1/4 t{a3;i1} t{a4;i2})"},

      {L"1/4 g{a1,a2;a3,a4} t{a3;i1} t{a4;i2}",
       L"1/4 g{i3,i4;a3,a4} t{a1;i3} t{a2;i4} t{a3;i1} t{a4;i2}",
       L"1/4(g{a1,a2;a3,a4} + g{i3,i4;a3,a4} t{a1;i3} t{a2;i4}) t{a3;i1} "
       L"t{a4;i2}"},

      {L"1/8 g{i3,i4;a3,a4} t{a1;i3} t{a2;i4} t{a3,a4;i1,i2}",
       L"1/4 g{i3,i4;a3,a4} t{a1;i3} t{a2;i4} t{a3;i1} t{a4;i2}",
       L"g{i3,i4;a3,a4} t{a1;i3} t{a2;i4} "
       L"                      (1/8 t{a3,a4;i1,i2} + 1/4 t{a3;i1} t{a4;i2})"}};

  for (auto&& [l, r, f] : fused_terms) {
    auto const le = parse_expr(l, Symmetry::nonsymm);
    auto const re = parse_expr(r, Symmetry::nonsymm);
    auto const fe = parse_expr(f, Symmetry::nonsymm);
    auto fu = Fusion{le->as<Product>(), re->as<Product>()};
    REQUIRE((fu.left() || fu.right()));

    auto const& fue = fu.left() ? fu.left() : fu.right();

    REQUIRE(fe == fue);
  }
}
