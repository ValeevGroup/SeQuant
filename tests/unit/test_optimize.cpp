#include <catch2/catch_test_macros.hpp>

#include "catch2_sequant.hpp"

#include <SeQuant/core/algorithm.hpp>
#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/optimize.hpp>
#include <SeQuant/core/parse.hpp>
#include <SeQuant/core/space.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>

#include <cstddef>
#include <initializer_list>
#include <memory>
#include <stdexcept>

#include <range/v3/all.hpp>

sequant::ExprPtr extract(sequant::ExprPtr expr,
                         std::initializer_list<size_t> const& idxs) {
  using namespace sequant;
  ExprPtr result = expr;
  for (auto s : idxs) result = result->at(s);
  return result;
}

TEST_CASE("optimize", "[optimize]") {
  using namespace sequant;

  // for optimization tests, set index space sizes
  {
    auto reg = get_default_context().mutable_index_space_registry();
    auto occ = reg->retrieve_ptr(L"i");
    auto uocc = reg->retrieve_ptr(L"a");
    auto aux = reg->retrieve_ptr(L"x");
    assert(occ);
    assert(uocc);
    assert(aux);
    occ->approximate_size(10);
    uocc->approximate_size(100);
    aux->approximate_size(4);
    assert(uocc->approximate_size() == 100);
  }

  auto single_term_opt = [](Product const& prod) {
    return opt::single_term_opt(prod, [](Index const& ix) {
      auto lbl = to_string(ix.label());
      auto sz = ix.space().approximate_size();
      return ix.space().approximate_size();
    });
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

    //
    // single-term optimization when sequant::Variables appear in a product
    //
    auto prod6 = parse_expr(
                     L"α * β * γ * "
                     "g_{i3,i4}^{a3,a4}"      // T1
                     " * t_{a1,a2}^{i3,i4}"   // T2
                     " * t_{a3,a4}^{i1,i2}")  // T3
                     ->as<Product>();
    auto res6 = single_term_opt(prod6);

    // this is the one we want to find
    // α * β * γ * ((T1 * T3) * T2)  : 2 * O^4 * V^2  best if nvirt > nocc
    REQUIRE(extract(res6, {0}) == prod6.at(0));
    REQUIRE(extract(res6, {1}) == prod6.at(1));
    REQUIRE(extract(res6, {2}) == prod6.at(2));
    REQUIRE(extract(res6, {3, 0}) == prod6.at(3));
    REQUIRE(extract(res6, {3, 1}) == prod6.at(5));
    REQUIRE(extract(res6, {4}) == prod6.at(4));

    //
    // single-term optimization including tensors with auxiliary indices
    //
    auto prod7 = parse_expr(
                     L"DF{a_1;a_3;x_1} "  // T1
                     "DF{a_2;i_1;x_1} "   // T2
                     "t{a_3;i_2}"         // T3
                     )
                     ->as<Product>();
    auto res7 = single_term_opt(prod7);

    // this is the one we want to find
    // (T1 T3) T2: V^2 O^1 A^1 + V^2 O^2 A^1 best if nvirt > nocc and nvirt >
    // nact
    REQUIRE(extract(res7, {0, 0}) == prod7.at(0));
    REQUIRE(extract(res7, {0, 1}) == prod7.at(2));
    REQUIRE(extract(res7, {1}) == prod7.at(1));

    auto prod8 = parse_expr(
                     L"T1{i_1;i_2;x_1,x_2,x_3,x_4} T2{i_2;i_1;x_5,x_6,x_7,x_8} "
                     L"T3{i_3;;x_1,x_2,x_3,x_4} T4{i_4;;x_5,x_6,x_7,x_8}")
                     ->as<Product>();
    auto res8 = single_term_opt(prod8);

    // this is the one we want to find
    // (T1 T3)(T2 T4)
    REQUIRE(extract(res8, {0, 0}) == prod8.at(0));
    REQUIRE(extract(res8, {0, 1}) == prod8.at(2));
    REQUIRE(extract(res8, {1, 0}) == prod8.at(1));
    REQUIRE(extract(res8, {1, 1}) == prod8.at(3));
  }

  SECTION("Ensure single-value sums/products are not discarded") {
    auto sum = ex<Sum>();
    sum->as<Sum>().append(ex<Product>(ExprPtrList{parse_expr(L"f{a_1;i_1}")}));
    REQUIRE(sum->as<Sum>().summand(0).as<Product>().factors().size() == 1);
    auto optimized = optimize(sum);
    REQUIRE(optimized->is<Sum>());
    REQUIRE(optimized->as<Sum>().summands().size() == 1);
    REQUIRE(sum->as<Sum>().summand(0).as<Product>().factors().size() == 1);
  }
}
