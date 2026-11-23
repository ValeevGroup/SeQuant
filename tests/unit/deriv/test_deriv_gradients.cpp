#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/domain/mbpt/op.hpp>
#include <SeQuant/domain/mbpt/rotation_gradient.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "catch2_sequant.hpp"
#include "contexts.hpp"
#include "judge.hpp"

#include <string>
#include <string_view>

using namespace sequant;

TEST_CASE("deriv/gradients/rotation", "[deriv][gradients]") {
  using deriv::test::require_same;
  using mbpt::rotation_gradient;
  auto T = [](std::wstring_view s) { return deserialize(s); };

  SECTION("the Hartree-Fock energy: the Fock operator") {
    const auto E = T(L"h{i_1;i_1}:N-N-N + 1/2 g{i_1,i_2;i_1,i_2}:A-N-S");
    const Index i(L"i_5"), a(L"a_5");
    require_same(rotation_gradient(E, {L"h", L"g"}, i, a),
                 T(L"h{a_5;i_5}:N-N-N + g{a_5,i_1;i_5,i_1}:A-N-S + "
                   L"h{i_5;a_5}:N-N-N + g{i_5,i_1;a_5,i_1}:A-N-S"),
                 {i, a}, {});
    // a virtual p: no slot of the energy runs over it
    REQUIRE(rotation_gradient(E, {L"h", L"g"}, a, i)->is_zero());
    // only the leaves named rotate
    require_same(rotation_gradient(E, {L"h"}, i, a),
                 T(L"h{a_5;i_5}:N-N-N + h{i_5;a_5}:N-N-N"), {i, a}, {});
  }

  SECTION("a Hadamard index and an auxiliary slot") {
    const Index i(L"i_5"), a(L"a_5");
    require_same(
        rotation_gradient(
            T(L"B{i_1;a_1;x_1}:N-N-N C{a_1;i_1;x_1}:N-N-N D{;;x_1}:N-N-N"),
            {L"B"}, i, a),
        T(L"B{a_5;a_1;x_1}:N-N-N C{a_1;i_5;x_1}:N-N-N D{;;x_1}:N-N-N"), {i, a},
        {});
  }

  SECTION("rejections") {
    const auto E = T(L"f{i_1;a_1}:N-N-N t{a_1;i_1}:N-N-N");
    const Index i(L"i_5"), a(L"a_5");
    const auto X = [&](const ExprPtr& e, const Index& p, const Index& q,
                       std::string text) {
      REQUIRE_THROWS_WITH(rotation_gradient(e, {L"f"}, p, q),
                          "mbpt::rotation_gradient: " + text);
    };
    X(nullptr, i, a, "null argument");
    X(E, i, i, "p and q are one index");
    X(E, Index(L"i_1"), a, "the function uses the label of p or q");
    X(E, Index(L"p_5"), a,
      "a slot of a leaf to rotate overlaps the space of p without including "
      "it");
    X(T(L"f{i_1;a_1}:N-N-N"), i, a, "a slot of a leaf to rotate is free");
    X(pow(T(L"f{i_1;a_1}:N-N-N"), 2) * T(L"t{a_1;i_1}:N-N-N"), i, a,
      "a leaf to rotate stands under a TensorPower");
    X(mbpt::op::H(1), i, a, "the function holds operators");
  }
}
