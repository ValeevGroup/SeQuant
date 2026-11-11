#include <catch2/catch_test_macros.hpp>

#include "catch2_sequant.hpp"

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/expressions/expr_algorithms.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/optimize/fusion.hpp>

#include <array>
#include <memory>
#include <string_view>
#include <vector>

TEST_CASE("fusion", "[optimize]") {
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
       L"1/8 g{a1,a2;a3,a4}(t{a3,a4;i1,i2} + 2 t{a3;i1} t{a4;i2})"},

      {L"1/4 g{a1,a2;a3,a4} t{a3;i1} t{a4;i2}",
       L"1/4 g{i3,i4;a3,a4} t{a1;i3} t{a2;i4} t{a3;i1} t{a4;i2}",
       L"1/4(g{a1,a2;a3,a4} + g{i3,i4;a3,a4} t{a1;i3} t{a2;i4}) t{a3;i1} "
       L"t{a4;i2}"},

      {L"1/8 g{i3,i4;a3,a4} t{a1;i3} t{a2;i4} t{a3,a4;i1,i2}",
       L"1/4 g{i3,i4;a3,a4} t{a1;i3} t{a2;i4} t{a3;i1} t{a4;i2}",
       L"1/8 g{i3,i4;a3,a4} t{a1;i3} t{a2;i4} "
       L"                      (t{a3,a4;i1,i2} + 2 t{a3;i1} t{a4;i2})"},

      {L"-1/8 g{i3,i4;a3,a4} t{a1;i3} t{a2;i4} t{a3,a4;i1,i2}",
       L"-1/4 g{i3,i4;a3,a4} t{a1;i3} t{a2;i4} t{a3;i1} t{a4;i2}",
       L"-1/8 g{i3,i4;a3,a4} t{a1;i3} t{a2;i4} "
       L"                      (t{a3,a4;i1,i2} + 2 t{a3;i1} t{a4;i2})"}

  };

  for (auto&& [l, r, f] : fused_terms) {
    auto const le = deserialize(l);
    auto const re = deserialize(r);
    auto const fe = deserialize(f);
    auto fu = Fusion{le->as<Product>(), re->as<Product>()};
    REQUIRE((fu.left() || fu.right()));

    auto const& fue = fu.left() ? fu.left() : fu.right();

    REQUIRE(fe == fue);
  }
}

namespace {

/// \return the binary product \p scalar * \p lhs * \p rhs, kept as written
sequant::ExprPtr binary(sequant::rational scalar, sequant::ExprPtr const& lhs,
                        sequant::ExprPtr const& rhs) {
  using namespace sequant;
  auto result = ex<Product>();
  result->as<Product>().scale(scalar);
  result->as<Product>().append(1, lhs, Product::Flatten::No);
  result->as<Product>().append(1, rhs, Product::Flatten::No);
  return result;
}

/// \return whether \p lhs and \p rhs are equal once expanded and simplified
bool same_value(sequant::ExprPtr const& lhs, sequant::ExprPtr const& rhs) {
  using namespace sequant;
  auto diff =
      ex<Sum>(ExprPtrList{lhs->clone(), ex<Constant>(-1) * rhs->clone()});
  expand(diff);
  simplify(diff);
  return diff->is<Constant>() && diff->as<Constant>().is_zero();
}

}  // namespace

TEST_CASE("fuse_top_factors", "[optimize]") {
  using namespace sequant;
  using opt::fuse_top_factors;
  auto t = [](std::wstring_view s) { return deserialize(s); };
  auto sum_of = [](Sum::summands_type const& s) {
    return ex<Sum>(s.begin(), s.end());
  };

  SECTION("shared factor up to dummy renaming") {
    // t{a1;i3} f{i3;i1} + 2 t{a1;i4} (g{i4,i2;i1,a2} t{a2;i2})
    //   -> t{a1;i3} (f{i3;i1} + 2 g{i3,i2;i1,a2} t{a2;i2})
    Sum::summands_type in{
        binary(1, t(L"t{a1;i3}"), t(L"f{i3;i1}")),
        binary(2, t(L"t{a1;i4}"),
               binary(1, t(L"g{i4,i2;i1,a2}"), t(L"t{a2;i2}")))};
    auto const out = fuse_top_factors(in);
    REQUIRE(out.size() == 1);
    auto const expected = binary(
        1, t(L"t{a1;i3}"),
        ex<Sum>(ExprPtrList{t(L"f{i3;i1}"),
                            binary(2, t(L"g{i3,i2;i1,a2}"), t(L"t{a2;i2}"))}));
    REQUIRE(*out.front() == *expected);
    REQUIRE(same_value(sum_of(out), sum_of(in)));
  }

  SECTION("clashing dummy of the summed factor is renamed") {
    // i3 is a dummy of the second X, and the shared factor's i4 becomes i3
    Sum::summands_type in{
        binary(1, t(L"t{a1;i3}"), t(L"f{i3;i1}")),
        binary(1, t(L"t{a1;i4}"),
               binary(1, t(L"g{i4,i3;i1,a2}"), t(L"t{a2;i3}")))};
    auto const out = fuse_top_factors(in);
    REQUIRE(out.size() == 1);
    REQUIRE(same_value(sum_of(out), sum_of(in)));
  }

  SECTION("shared factor on the right") {
    Sum::summands_type in{
        binary(1, t(L"f{i3;i1}"), t(L"t{a1;i3}")),
        binary(-1, binary(1, t(L"g{i3,i4;i1,a2}"), t(L"t{a2;i4}")),
               t(L"t{a1;i3}"))};
    auto const out = fuse_top_factors(in);
    REQUIRE(out.size() == 1);
    REQUIRE(out.front()->as<Product>().factor(1)->is<Tensor>());
    REQUIRE(same_value(sum_of(out), sum_of(in)));
  }

  SECTION("nested fusion keeps renamed dummies off enclosing indices") {
    // fusing on t{a1;i?} renames the second X's dummy i3; the X's then fuse on
    // t{a2,a?;i2,i?}, which renames the dummy a3; its new name must not be a1,
    // an index of the enclosing product
    Sum::summands_type in{
        binary(1, t(L"t{a1;i3}"),
               binary(1, t(L"t{a2,a3;i2,i4}"), t(L"g{i3,i4;i1,a3}"))),
        binary(-1, t(L"t{a1;i5}"),
               binary(1, t(L"t{a2,a4;i2,i3}"),
                      binary(1, t(L"g{i5,i3;a3,a4}"), t(L"t{a3;i1}"))))};
    auto const out = fuse_top_factors(in);
    REQUIRE(out.size() == 1);
    auto const& inner = out.front()->as<Product>().factor(1);
    REQUIRE(inner->is<Sum>());
    REQUIRE(inner->as<Sum>().size() == 1);  // the X's fused too
    REQUIRE(same_value(sum_of(out), sum_of(in)));
  }

  SECTION("nothing to fuse") {
    Sum::summands_type in{binary(1, t(L"f{a1;a2}"), t(L"t{a2;i1}")),
                          binary(1, t(L"f{i2;i1}"), t(L"t{a1;i2}")),
                          t(L"f{a1;i1}")};
    auto const out = fuse_top_factors(in);
    REQUIRE(out.size() == in.size());
    for (std::size_t i = 0; i < in.size(); ++i) REQUIRE(out[i] == in[i]);
  }

  SECTION("indices with proto indices are not fused") {
    Sum::summands_type in{binary(1, t(L"t{a1<i1>;i3}"), t(L"f{i3;i1}")),
                          binary(1, t(L"t{a1<i1>;i4}"), t(L"f{i4;i1}"))};
    auto const out = fuse_top_factors(in);
    REQUIRE(out.size() == 2);
  }
}
