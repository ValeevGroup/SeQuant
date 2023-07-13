#include "catch.hpp"

#include <SeQuant/core/optimize.hpp>
#include <SeQuant/core/parse_expr.hpp>

sequant::ExprPtr extract(sequant::ExprPtr expr,
                         std::initializer_list<size_t> const& idxs) {
  using namespace sequant;
  ExprPtr result = expr;
  for (auto s : idxs) result = result->at(s);
  return result;
}

TEST_CASE("TEST_OPTIMIZE", "[optimize]") {
  using namespace sequant;

  auto idx2size = [nocc = 4, nvirt = 140](Index const& idx) {
    if (idx.space() == IndexSpace::active_occupied) return nocc;
    if (idx.space() == IndexSpace::active_unoccupied)
      return nvirt;
    else
      throw std::runtime_error("Unsupported IndexSpace type encountered");
  };

  auto single_term_opt = [&idx2size](Product const& prod) {
    return opt::single_term_opt(prod, idx2size);
  };

  auto parse_expr_antisymm = [](auto const& xpr) {
    return parse_expr(xpr, Symmetry::antisymm);
  };

  SECTION("Single term optimization") {
    const auto prod1 = parse_expr_antisymm(
                           L"g_{i3,i4}^{a3,a4}"     // T1
                           " * t_{a1,a2}^{i3,i4}"   // T2
                           " * t_{a3,a4}^{i1,i2}")  // T3
                           ->as<Product>();
    //
    // Cost of evaluation prod1:
    //
    // ((T1 * T2) * T3)  : 2 * O^2 * V^4  best if nocc > nvirt
    //
    // this is the one we want to find
    // ((T1 * T3) * T2)  : 2 * O^4 * V^2  best if nvirt > nocc
    //
    // (T1 * (T2 * T3))  : 2 * O^4 * V^4  worst sequence of evaluation
    //

    const auto res1 = single_term_opt(prod1);

    REQUIRE(extract(res1, {0, 0}) == prod1.at(0));
    REQUIRE(extract(res1, {0, 1}) == prod1.at(2));
    REQUIRE(extract(res1, {1}) == prod1.at(1));

    const auto prod2 = parse_expr_antisymm(
                           L"   g_{i3,i4}^{a3,a4}"
                           L" * t_{a3,a4}^{i1,i2}"
                           L" * t_{a1}^{i3}"
                           L" * t_{a2}^{i4}")
                           ->as<Product>();

    const auto res2 = single_term_opt(prod2);

    REQUIRE(extract(res2, {0, 0, 0}) == prod2.at(0));
    REQUIRE(extract(res2, {0, 0, 1}) == prod2.at(1));
    REQUIRE(extract(res2, {0, 1}) == prod2.at(2));
    REQUIRE(extract(res2, {1}) == prod2.at(3));

    const auto prod3 = parse_expr_antisymm(
                           L""                   //
                           " g_{i3,i4}^{a3,a4}"  //
                           " t_{a1}^{i3}"        //
                           " t_{a2}^{i4}"        //
                           " t_{a3,a4}^{i1,i2}"  //
                           )
                           ->as<Product>();
    auto res3 = single_term_opt(prod3);

    REQUIRE(extract(res3, {0, 0, 0}) == prod3.at(0));
    REQUIRE(extract(res3, {0, 0, 1}) == prod3.at(3));
    REQUIRE(extract(res3, {0, 1}) == prod3.at(1));
    REQUIRE(extract(res3, {1}) == prod3.at(2));

    //
    // single-term optimization when a dot product occurs in the tensor network
    // ========================

    auto prod4 =
        parse_expr_antisymm(L"1/4 λ{i1;a1} g{i2,i3;a2,a3} t{a2,a3;i2,i3}")
            ->as<Product>();
    auto res4 = single_term_opt(prod4);

    REQUIRE(extract(res4, {0}) == prod4.at(0));
    REQUIRE(extract(res4, {1, 0}) == prod4.at(1));
    REQUIRE(extract(res4, {1, 1}) == prod4.at(2));

    auto prod5 =
        parse_expr_antisymm(L"x{i1,i2;a3,a4} y{a1,a2;i1,i2} z{a3,a4;a1,a2}")
            ->as<Product>();
    auto res5 = single_term_opt(prod5);
    REQUIRE(extract(res5, {0, 0}) == prod5.at(0));
    REQUIRE(extract(res5, {0, 1}) == prod5.at(2));
    REQUIRE(extract(res5, {1}) == prod5.at(1));
  }
}