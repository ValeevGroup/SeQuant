#include <SeQuant/core/context.hpp>
#include <SeQuant/core/deriv/deriv.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/expressions/result_expr.hpp>
#include <SeQuant/core/expressions/tensor_power.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/rational.hpp>
#include <SeQuant/domain/mbpt/op.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "catch2_sequant.hpp"
#include "contexts.hpp"
#include "judge.hpp"

#include <string>
#include <string_view>
#include <utility>

using namespace sequant;

using sequant::deriv::test::real_spaces;

TEST_CASE("deriv/rules/scalars", "[deriv][rules]") {
  const auto diff = [](const ExprPtr& E, const ExprPtr& x) {
    return sequant::deriv::diff(E, {x});
  };
  auto T = [](std::wstring_view s) { return deserialize(s); };
  const auto x = ex<Variable>(L"x");
  const auto x_star = adjoint(x);

  REQUIRE(diff(ex<Constant>(3), x)->is_zero());
  REQUIRE(diff(T(L"y z"), x)->is_zero());
  REQUIRE_THAT(diff(x, x), EquivalentTo(L"1"));
  REQUIRE_THAT(diff(T(L"3 x^(2) y + x"), x), EquivalentTo(L"6 x y + 1"));
  REQUIRE_THAT(diff(T(L"x y x"), x), EquivalentTo(L"2 x y"));
  REQUIRE_THAT(diff(T(L"x^(1/2)"), x), EquivalentTo(L"1/2 x^(-1/2)"));
  REQUIRE_THAT(diff(T(L"x^(-1)"), x), EquivalentTo(L"-1 x^(-2)"));
  // x and x^* are independent
  REQUIRE_THAT(diff(x * x_star, x), EquivalentTo(x_star));
  REQUIRE_THAT(diff(x * x_star, x_star), EquivalentTo(x));
  // a conjugated power is a power of x^*
  auto p = ex<Power>(x, rational(3));
  p->adjoint();
  auto q = ex<Power>(x, rational(2));
  q->adjoint();
  REQUIRE(diff(p, x)->is_zero());
  REQUIRE_THAT(diff(p, x_star), EquivalentTo(ex<Constant>(3) * q));
  REQUIRE_THAT(diff(q, x_star), EquivalentTo(ex<Constant>(2) * x_star));
  // a power of x^* is one too
  const auto r = ex<Power>(x_star, rational(3));
  REQUIRE(diff(r, x)->is_zero());
  REQUIRE_THAT(diff(r, x_star),
               EquivalentTo(ex<Constant>(3) * ex<Power>(x_star, rational(2))));
  // a nested product with its own scalar
  REQUIRE_THAT(
      diff(ex<Product>(2,
                       ExprPtrList{x, ex<Product>(3, ExprPtrList{x, T(L"y")})}),
           x),
      EquivalentTo(L"12 x y"));
  // a deeper nested sum
  REQUIRE_THAT(diff(x * (T(L"y") + x * (x + T(L"z"))), x),
               EquivalentTo(L"y + 3 x^(2) + 2 x z"));
  // tensors are constants
  REQUIRE_THAT(diff(T(L"x^(3) f{i_1;a_1}:N-N-N t{a_1;i_1}:N-N-N"), x),
               EquivalentTo(L"3 x^(2) f{i_1;a_1}:N-N-N t{a_1;i_1}:N-N-N"));
  // an input that is not expanded
  REQUIRE_THAT(diff((x + T(L"y")) * (x + T(L"z")), x),
               EquivalentTo(L"2 x + y + z"));
  // any order, left to right; no variable is the function itself
  const auto y = ex<Variable>(L"y");
  REQUIRE_THAT(sequant::deriv::diff(T(L"x^(3) y^(2)"), {x, x, y}),
               EquivalentTo(L"12 x y"));
  REQUIRE_THAT(sequant::deriv::diff(T(L"x^(3) y^(2)"), {}),
               EquivalentTo(L"x^(3) y^(2)"));
}

TEST_CASE("deriv/rules/tensors", "[deriv][rules]") {
  using sequant::deriv::test::check_rule;
  const auto diff = [](const ExprPtr& E, const ExprPtr& x) {
    return sequant::deriv::diff(E, {x});
  };
  auto T = [](std::wstring_view s) { return deserialize(s); };
  const auto R = [](std::wstring_view f, std::wstring_view x,
                    std::wstring_view expected) { check_rule(f, x, expected); };

  // a contracted leaf: the rest of the term under x's names
  R(L"f{i_1;a_1}:N-N-N t{a_1;i_1}:N-N-N", L"t{a_2;i_2}:N-N-N",
    L"f{i_2;a_2}:N-N-N");
  // a leaf by itself: both slots are external, the Kronecker deltas stay
  R(L"t{a_1;i_1}:N-N-N", L"t{a_2;i_2}:N-N-N",
    L"δ{a_1;a_2}:N-S-S δ{i_2;i_1}:N-S-S");
  // the product rule
  R(L"1/2 g{i_1,i_2;a_1,a_2}:N-N-N t{a_1;i_1}:N-N-N t{a_2;i_2}:N-N-N",
    L"t{a_3;i_3}:N-N-N",
    L"1/2 g{i_3,i_2;a_3,a_2}:N-N-N t{a_2;i_2}:N-N-N + "
    L"1/2 g{i_1,i_3;a_1,a_3}:N-N-N t{a_1;i_1}:N-N-N");
  // index names: x spelled with names E sums over
  R(L"g{i_1,i_2;a_1,a_2}:N-N-N t{a_1;i_1}:N-N-N t{a_2;i_2}:N-N-N",
    L"t{a_1;i_1}:N-N-N",
    L"g{i_1,i_3;a_1,a_3}:N-N-N t{a_3;i_3}:N-N-N + "
    L"g{i_3,i_1;a_3,a_1}:N-N-N t{a_3;i_3}:N-N-N");
  // ... and E holds the names make_tmp_index mints next
  {
    const auto occupied = Index(L"i_1").basis();
    Index::reset_tmp_index();
    const auto i_a = Index::make_tmp_index(occupied);
    const auto i_b = Index::make_tmp_index(occupied);
    const auto loop = [](const Index& i) {
      return ex<Tensor>(L"C", bra{i}, ket{}) * ex<Tensor>(L"Q", bra{}, ket{i});
    };
    Index::reset_tmp_index();
    check_rule(T(L"w{i_1;a_1}:N-N-N t{a_1;i_1}:N-N-N") * loop(i_a) * loop(i_b),
               T(L"t{a_1;i_1}:N-N-N"),
               T(L"w{i_1;a_1}:N-N-N") * loop(i_a) * loop(i_b));
  }
  // another block of the same symbol is a constant
  R(L"f{i_1;a_1}:N-N-N t{a_1;i_1}:N-N-N + w{a_1;i_1}:N-N-N t{i_1;a_1}:N-N-N",
    L"t{a_2;i_2}:N-N-N", L"f{i_2;a_2}:N-N-N");
  // a leaf of the same symbol with other ranks is another object, whatever
  // symmetries it declares
  R(L"f{i_1;a_1}:N-N-N t{a_1;i_1}:N-N-N + "
    L"1/4 g{i_1,i_2;a_1,a_2}:A-N-S t{a_1,a_2;i_1,i_2}:A-N-S",
    L"t{a_3;i_3}:N-N-N", L"f{i_3;a_3}:N-N-N");
  // E does not hold x
  REQUIRE(diff(T(L"f{i_1;a_1}:N-N-N u{a_1;i_1}:N-N-N"), T(L"t{a_2;i_2}:N-N-N"))
              ->is_zero());
  // E with its own Â: its external indices are not renamed
  R(L"Â{i_1;a_1}:A-N-S v{a_1,i_2;i_1,a_2}:N-N-N t{a_2;i_2}:N-N-N",
    L"t{a_3;i_3}:N-N-N", L"Â{i_1;a_1}:A-N-S v{a_1,i_3;i_1,a_3}:N-N-N");
  // trace: a leaf that holds one index in two slots
  R(L"h{p_1;p_1}:N-N-N", L"h{p_2;p_3}:N-N-N", L"δ{p_3;p_2}:N-S-S");
  // ... by a variable whose slots run over disjoint spaces: no element of it
  // is on the diagonal
  R(L"h{p_1;p_1}:N-N-N", L"h{i_3;a_4}:N-N-N", L"0");
  // a free index that appears only as a protoindex is external: it keeps its
  // name
  R(L"w{;a_5<i_2>}:N-N-N v{a_5<i_2>;}:N-N-N f{i_1;a_1}:N-N-N t{a_1;i_1}:N-N-N",
    L"t{a_3;i_3}:N-N-N",
    L"w{;a_5<i_2>}:N-N-N v{a_5<i_2>;}:N-N-N f{i_3;a_3}:N-N-N");
  // unexpanded input
  check_rule((T(L"f{i_1;a_1}:N-N-N") + T(L"w{i_1;a_1}:N-N-N")) *
                 T(L"t{a_1;i_1}:N-N-N"),
             T(L"t{a_2;i_2}:N-N-N"), T(L"f{i_2;a_2}:N-N-N + w{i_2;a_2}:N-N-N"));
  // a scalar factor passes through
  R(L"x^(2) f{i_1;a_1}:N-N-N t{a_1;i_1}:N-N-N", L"t{a_2;i_2}:N-N-N",
    L"x^(2) f{i_2;a_2}:N-N-N");
  // an aux slot
  check_rule(ex<Tensor>(L"B", bra{Index(L"a_1")}, ket{Index(L"i_1")},
                        aux{Index(L"p_1")}) *
                 ex<Tensor>(L"C", bra{Index(L"i_1")}, ket{Index(L"a_1")},
                            aux{Index(L"p_1")}),
             ex<Tensor>(L"B", bra{Index(L"a_2")}, ket{Index(L"i_2")},
                        aux{Index(L"p_2")}),
             ex<Tensor>(L"C", bra{Index(L"i_2")}, ket{Index(L"a_2")},
                        aux{Index(L"p_2")}));
  // protoindices that are the variable's own slots
  R(L"f{i_1;a_1<i_1>}:N-N-N t{a_1<i_1>;i_1}:N-N-N", L"t{a_2<i_2>;i_2}:N-N-N",
    L"f{i_2;a_2<i_2>}:N-N-N");
  // a leaf and x on different pairs: Kronecker deltas on every label, never an
  // overlap
  check_rule(T(L"t{a_1<i_1>;i_1}:N-N-N"), T(L"t{a_2<i_2>;i_2}:N-N-N"),
             make_kronecker(Index(L"a_1", {L"i_1"}), Index(L"a_2", {L"i_2"})) *
                 make_kronecker(Index(L"i_2"), Index(L"i_1")));
  // labels that are not slots: a per-pair rotation; its pair has no order, so
  // both orders appear
  R(L"u{;i_1,i_2}:N-N-N M{a_3<i_1,i_2>;a_1<i_1,i_2>}:N-N-N "
    L"k{a_1<i_1,i_2>;a_3<i_1,i_2>}:N-N-N v{i_1,i_2;}:N-N-N",
    L"k{a_2<i_5,i_6>;a_4<i_5,i_6>}:N-N-N",
    L"u{;i_5,i_6}:N-N-N M{a_4<i_5,i_6>;a_2<i_5,i_6>}:N-N-N v{i_5,i_6;}:N-N-N + "
    L"u{;i_6,i_5}:N-N-N M{a_4<i_5,i_6>;a_2<i_5,i_6>}:N-N-N v{i_6,i_5;}:N-N-N");
  // a single label: one term
  R(L"u{;i_1}:N-N-N M{a_3<i_1>;a_1<i_1>}:N-N-N k{a_1<i_1>;a_3<i_1>}:N-N-N "
    L"v{i_1;}:N-N-N",
    L"k{a_2<i_5>;a_4<i_5>}:N-N-N",
    L"u{;i_5}:N-N-N M{a_4<i_5>;a_2<i_5>}:N-N-N v{i_5;}:N-N-N");
  // overlapping spaces: the leaf runs over the union space, the variable over
  // a block: renamed
  R(L"W{i_1;p_1}:N-N-N k{p_1;i_1}:N-N-N", L"k{a_2;i_2}:N-N-N",
    L"W{i_2;a_2}:N-N-N");
  // the reverse: the leaf's summed indices run over less than the variable's,
  // so the deltas stay
  R(L"W{i_1;a_1}:N-N-N k{a_1;i_1}:N-N-N", L"k{p_2;p_3}:N-N-N",
    L"δ{a_1;p_2}:N-S-S δ{p_3;i_1}:N-S-S W{i_1;a_1}:N-N-N");
  // a block of a pair density from amplitudes over the full space
  R(L"v{p_3;a_1}:N-N-N u{a_1,p_1;i_1,i_2}:N-N-N w{i_1,i_2;p_3,p_1}:N-N-N",
    L"u{p_5,p_6;i_5,i_6}:N-N-N",
    L"δ{a_1;p_5}:N-S-S v{p_3;a_1}:N-N-N w{i_5,i_6;p_3,p_6}:N-N-N");
  // the derivative by the adjoint spelling (complex field: t and t⁺ are
  // independent)
  R(L"f{i_1;a_1}:N-N-N t{a_1;i_1}:N-N-N + w{a_1;i_1}:N-N-N t⁺{i_1;a_1}:N-N-N",
    L"t{a_2;i_2}:N-N-N", L"f{i_2;a_2}:N-N-N");
  R(L"f{i_1;a_1}:N-N-N t{a_1;i_1}:N-N-N + w{a_1;i_1}:N-N-N t⁺{i_1;a_1}:N-N-N",
    L"t⁺{i_2;a_2}:N-N-N", L"w{a_2;i_2}:N-N-N");
  // a Kronecker delta follows a rename that a later label of the leaf decides
  R(L"h{a_1;a_1}:N-N-N", L"h{p_2;a_3}:N-N-N", L"δ{a_3;p_2}:N-S-S");
  R(L"W{i_1;a_1<i_1>}:N-N-N k{a_1<i_1>;i_1}:N-N-N", L"k{e_2<i_2>;i_2}:N-N-N",
    L"W{i_2;a_1<i_2>}:N-N-N δ{a_1<i_2>;e_2<i_2>}:N-S-S");
  const auto R_nonzero = [&](std::wstring_view f, std::wstring_view x) {
    check_rule(f, x, L"");
    REQUIRE_FALSE(diff(T(f), T(x))->is_zero());
  };
  R_nonzero(
      L"u{;i_1}:N-N-N W{;a_1<i_1>}:N-N-N k{a_1<i_1>;}:N-N-N "
      L"v{i_1;}:N-N-N",
      L"k{e_2<i_2>;}:N-N-N");
  R_nonzero(
      L"u{;i_1,i_2}:N-N-N M{a_3<i_1,i_2>;a_1<i_1,i_2>}:N-N-N "
      L"k{a_1<i_1,i_2>;a_3<i_1,i_2>}:N-N-N v{i_1,i_2;}:N-N-N",
      L"k{e_2<i_5,i_6>;a_4<i_5,i_6>}:N-N-N");
  R_nonzero(L"W{;i_1}:N-N-N h{a_1<i_1>;a_1<i_1>}:N-N-N V{i_1;}:N-N-N",
            L"h{e_2<i_2>;a_3<i_2>}:N-N-N");
  // two slots with different pair labels: the orders are taken per list
  R_nonzero(
      L"Y{;i_3,i_4}:N-N-N W{a_3<i_3>;a_4<i_4>}:N-N-N "
      L"s{a_4<i_4>;a_3<i_3>}:N-N-N Z{i_3,i_4;}:N-N-N",
      L"s{a_5<i_5>;a_6<i_6>}:N-N-N");
  R_nonzero(
      L"u{;i_1,i_2,i_3,i_4}:N-N-N M{a_3<i_3,i_4>;a_1<i_1,i_2>}:N-N-N "
      L"k{a_1<i_1,i_2>;a_3<i_3,i_4>}:N-N-N v{i_1,i_2,i_3,i_4;}:N-N-N",
      L"k{a_2<i_5,i_6>;a_4<i_7,i_8>}:N-N-N");
}

TEST_CASE("deriv/rules/symmetry", "[deriv][rules]") {
  using sequant::deriv::test::check_rule;
  const auto diff = [](const ExprPtr& E, const ExprPtr& x) {
    return sequant::deriv::diff(E, {x});
  };
  auto T = [](std::wstring_view s) { return deserialize(s); };
  const auto R = [](std::wstring_view f, std::wstring_view x,
                    std::wstring_view expected) { check_rule(f, x, expected); };

  // antisymmetric: 4 Â on the gradient; Â g = g here
  R(L"1/4 g{i_1,i_2;a_1,a_2}:A-C-S t{a_1,a_2;i_1,i_2}:A-N-S",
    L"t{a_3,a_4;i_3,i_4}:A-N-S",
    L"Â{a_3,a_4;i_3,i_4}:A-N-S g{i_3,i_4;a_3,a_4}:A-C-S");
  // a leaf by itself
  R(L"t{a_1,a_2;i_1,i_2}:A-N-S", L"t{a_3,a_4;i_3,i_4}:A-N-S",
    L"4 Â{a_3,a_4;i_3,i_4}:A-N-S δ{a_1;a_3}:N-S-S δ{a_2;a_4}:N-S-S "
    L"δ{i_3;i_1}:N-S-S δ{i_4;i_2}:N-S-S");
  // column symmetry: 2 Ŝ
  R(L"g{i_1,i_2;a_1,a_2}:N-C-S t{a_1,a_2;i_1,i_2}:N-N-S",
    L"t{a_3,a_4;i_3,i_4}:N-N-S",
    L"2 Ŝ{a_3,a_4;i_3,i_4}:N-N-S g{i_3,i_4;a_3,a_4}:N-C-S");
  // mixed spaces: no operation keeps the block; the leaf is spelled the other
  // way in E
  R(L"W{i_1,a_1;a_2,i_2}:N-N-N g{a_2,i_2;i_1,a_1}:A-C-S",
    L"g{i_3,a_3;i_4,a_4}:A-C-S", L"-1 W{i_4,a_4;a_3,i_3}:N-N-N");
  // the ket slots alone share a basis: a one-sided Â
  R(L"W{i_1,i_2;i_3,a_1}:N-N-N g{i_3,a_1;i_1,i_2}:A-C-S",
    L"g{i_4,a_4;i_5,i_6}:A-C-S",
    L"2 Â{;i_5,i_6}:A-N-S W{i_5,i_6;i_4,a_4}:N-N-N");
  // complex Hermitian over two spaces: the other block is the conjugate, a
  // constant
  R(L"1/4 W{i_1,i_2;a_1,a_2}:A-N-S g{a_1,a_2;i_1,i_2}:A-C-S + 1/4 "
    L"V{a_1,a_2;i_1,i_2}:A-N-S g{i_1,i_2;a_1,a_2}:A-C-S",
    L"g{a_3,a_4;i_3,i_4}:A-C-S",
    L"Â{a_3,a_4;i_3,i_4}:A-N-S W{i_3,i_4;a_3,a_4}:A-N-S");
  // complex Hermitian over one space: every element is its own variable
  R(L"G{p_2;p_1}:N-N-N h{p_1;p_2}:N-C-N", L"h{p_3;p_4}:N-C-N",
    L"G{p_4;p_3}:N-N-N");
  // a CSV amplitude
  R(L"1/4 g{i_1,i_2;a_1<i_1,i_2>,a_2<i_1,i_2>}:A-C-S "
    L"t{a_1<i_1,i_2>,a_2<i_1,i_2>;i_1,i_2}:A-N-S",
    L"t{a_3<i_3,i_4>,a_4<i_3,i_4>;i_3,i_4}:A-N-S",
    L"Â{a_3<i_3,i_4>,a_4<i_3,i_4>;i_3,i_4}:A-N-S "
    L"g{i_3,i_4;a_3<i_3,i_4>,a_4<i_3,i_4>}:A-C-S");
  // pair labels that are not slots under the variable's projector: labels of
  // its free indices, so free (the judge alone decides, and the result is not
  // zero)
  const std::pair<std::wstring_view, std::wstring_view> pair_labels[]{
      {L"u{;i_1,i_2}:N-N-N "
       L"M{a_3<i_1,i_2>,a_4<i_1,i_2>;a_1<i_1,i_2>,a_2<i_1,i_2>}:N-N-N "
       L"U{a_1<i_1,i_2>,a_2<i_1,i_2>;a_3<i_1,i_2>,a_4<i_1,i_2>}:A-N-S "
       L"v{i_1,i_2;}:N-N-N",
       L"U{a_5<i_5,i_6>,a_6<i_5,i_6>;a_7<i_5,i_6>,a_8<i_5,i_6>}:A-N-S"},
      {L"u{;i_1,i_2}:N-N-N "
       L"M{a_3<i_1,i_2>,a_4<i_1,i_2>;a_1<i_1,i_2>,a_2<i_1,i_2>}:N-N-N "
       L"U{a_1<i_1,i_2>,a_2<i_1,i_2>;a_3<i_1,i_2>,a_4<i_1,i_2>}:A-S-S "
       L"v{i_1,i_2;}:N-N-N",
       L"U{a_5<i_5,i_6>,a_6<i_5,i_6>;a_7<i_5,i_6>,a_8<i_5,i_6>}:A-S-S"},
      {L"u{;i_1}:N-N-N M{a_3;a_1<i_1>,a_2<i_1>}:N-N-N "
       L"U{a_1<i_1>,a_2<i_1>;a_3}:A-N-S v{i_1;}:N-N-N",
       L"U{a_5<i_5>,a_6<i_5>;a_7}:A-N-S"},
      {L"u{;i_1,i_2}:N-N-N "
       L"M{a_3<i_2>,a_4<i_2>;a_1<i_1>,a_2<i_1>}:N-N-N "
       L"U{a_1<i_1>,a_2<i_1>;a_3<i_2>,a_4<i_2>}:N-N-S v{i_1,i_2;}:N-N-N",
       L"U{a_5<i_5>,a_6<i_5>;a_7<i_6>,a_8<i_6>}:N-N-S"},
      {L"u{;i_1,i_2}:N-N-N "
       L"M{a_3<i_2>,a_4<i_2>;a_1<i_1>,a_2<i_1>}:N-N-N "
       L"U{a_1<i_1>,a_2<i_1>;a_3<i_2>,a_4<i_2>}:A-S-S v{i_1,i_2;}:N-N-N",
       L"U{a_5<i_5>,a_6<i_5>;a_7<i_6>,a_8<i_6>}:A-S-S"}};
  for (const auto& [f, x] : pair_labels) {
    check_rule(f, x, L"");
    REQUIRE_FALSE(diff(T(f), T(x))->is_zero());
  }
  // ... differentiated again by a leaf that holds the pair labels
  for (const auto& [k, x2] : {std::pair{0, L"u{;i_9,i_10}:N-N-N"},
                              std::pair{3, L"v{i_9,i_10;}:N-N-N"}}) {
    const auto& [f, x] = pair_labels[k];
    const auto first = diff(T(f), T(x));
    check_rule(first, T(x2), nullptr);
    REQUIRE_FALSE(diff(first, T(x2))->is_zero());
  }
  // a group that is part of a side: the two occupied bra slots and the two
  // unoccupied ket slots of a rank-3 tensor (the judge alone decides, and the
  // result is not zero)
  check_rule(
      T(L"W{a_2,a_3,i_3;i_1,i_2,a_1}:N-N-N w{i_1,i_2,a_1;a_2,a_3,i_3}:A-N-S"),
      T(L"w{i_4,i_5,a_4;a_5,a_6,i_6}:A-N-S"), nullptr);
  REQUIRE_FALSE(diff(T(L"W{a_2,a_3,i_3;i_1,i_2,a_1}:N-N-N "
                       L"w{i_1,i_2,a_1;a_2,a_3,i_3}:A-N-S"),
                     T(L"w{i_4,i_5,a_4;a_5,a_6,i_6}:A-N-S"))
                    ->is_zero());
  // a leaf over the union space by a variable whose related slots are alike:
  // the operations that keep the block are all its arrangements
  R(L"1/4 W{p_3,p_4;p_1,p_2}:N-N-N g{p_1,p_2;p_3,p_4}:A-N-S",
    L"g{i_5,i_6;a_5,a_6}:A-N-S",
    L"Â{i_5,i_6;a_5,a_6}:A-N-S W{a_5,a_6;i_5,i_6}:N-N-N");
  R(L"1/4 W{i_1,i_2;p_1,p_2}:N-N-N t{p_1,p_2;i_1,i_2}:A-N-S",
    L"t{a_3,a_4;i_3,i_4}:A-N-S",
    L"Â{a_3,a_4;i_3,i_4}:A-N-S W{i_3,i_4;a_3,a_4}:N-N-N");
  // E with its own Â, x with a projector: both stay as factors
  check_rule(T(L"Â{i_1,i_2;a_1,a_2}:A-N-S g{a_1,a_2;i_1,i_2}:A-C-S"),
             T(L"g{a_3,a_4;i_3,i_4}:A-C-S"), nullptr);
  // second derivatives commute
  {
    const auto E =
        T(L"1/2 w{i_1,i_2,i_3;a_1,a_2,a_3}:N-N-N t{a_1,a_2;i_1,i_2}:A-N-S "
          L"t{a_3;i_3}:N-N-N");
    const auto x2 = T(L"t{a_51,a_52;i_51,i_52}:A-N-S"),
               x1 = T(L"t{a_71;i_71}:N-N-N");
    REQUIRE_THAT(diff(diff(E, x2), x1), EquivalentTo(diff(diff(E, x1), x2)));
    REQUIRE_THAT(sequant::deriv::diff(E, {x2, x1}),
                 EquivalentTo(diff(diff(E, x2), x1)));
    REQUIRE_THAT(diff(diff(E, x2), x1),
                 EquivalentTo(L"2 Â{a_51,a_52;i_51,i_52}:A-N-S "
                              L"w{i_51,i_52,i_71;a_51,a_52,a_71}:N-N-N"));
  }
}

TEST_CASE("deriv/rules/bra-ket-symmetric", "[deriv][rules]") {
  namespace t = sequant::deriv::test;
  const auto R = [](std::wstring_view f, std::wstring_view x,
                    std::wstring_view expected, const t::Data& data) {
    t::check_rule(f, x, expected, data);
  };
  // a Symm tensor's data is exchange-symmetric
  const t::Data sym{.real = {L"h", L"g"}};

  // bra-ket symmetric over one space: G + Gᵀ, whatever the field and whatever
  // G is
  R(L"G{p_2;p_1}:N-N-N h{p_1;p_2}:N-S-N", L"h{p_3;p_4}:N-S-N",
    L"G{p_4;p_3}:N-N-N + G{p_3;p_4}:N-N-N", sym);
  // the same spaces in another order: core keeps both spellings, E holds the
  // other one
  R(L"W{i_1,a_1;a_2,i_2}:N-N-N g{a_2,i_2;i_1,a_1}:N-S-N",
    L"g{i_3,a_3;a_4,i_4}:N-S-N", L"W{i_3,a_3;a_4,i_4}:N-N-N", sym);
  // a leaf over the union space holds both spellings
  R(L"W{p_3,p_4;p_1,p_2}:N-N-N g{p_1,p_2;p_3,p_4}:N-S-N",
    L"g{i_3,a_3;a_4,i_4}:N-S-N",
    L"W{a_4,i_4;i_3,a_3}:N-N-N + W{i_3,a_3;a_4,i_4}:N-N-N", sym);
  // unequal bra and ket ranks: the exchange is no operation, and a leaf
  // written either way holds the variable's elements, whichever way the
  // variable is written
  for (const auto x : {L"g{;a_5}:N-S-N", L"g{a_5;}:N-S-N"}) {
    R(L"W{a_1;}:N-N-N g{;a_1}:N-S-N", x, L"W{a_5;}:N-N-N", sym);
    R(L"W{a_1;}:N-N-N g{;a_1}:N-S-N + V{;a_1}:N-N-N g{a_1;}:N-S-N", x,
      L"W{a_5;}:N-N-N + V{;a_5}:N-N-N", sym);
  }
  for (const auto x : {L"g{i_5;a_5,a_6}:N-S-N", L"g{a_5,a_6;i_5}:N-S-N"}) {
    R(L"W{a_1,a_2;i_1}:N-N-N g{i_1;a_1,a_2}:N-S-N", x, L"W{a_5,a_6;i_5}:N-N-N",
      sym);
    R(L"W{a_1,a_2;i_1}:N-N-N g{i_1;a_1,a_2}:N-S-N + "
      L"V{i_1;a_1,a_2}:N-N-N g{a_1,a_2;i_1}:N-S-N",
      x, L"W{a_5,a_6;i_5}:N-N-N + V{i_5;a_5,a_6}:N-N-N", sym);
  }
  for (const auto x : {L"g{i_5;a_5,a_6}:A-S-S", L"g{a_5,a_6;i_5}:A-S-S"})
    R(L"W{a_1,a_2;i_1}:N-N-N g{i_1;a_1,a_2}:A-S-S + "
      L"V{i_1;a_1,a_2}:N-N-N g{a_1,a_2;i_1}:A-S-S",
      x,
      L"2 Â{;a_5,a_6}:A-N-S W{a_5,a_6;i_5}:N-N-N + "
      L"2 Â{a_5,a_6;}:A-N-S V{i_5;a_5,a_6}:N-N-N",
      sym);
  // ... and with an aux slot
  R(L"W{a_1;;p_3}:N-N-N B{;a_1;p_3}:N-S-N", L"B{;a_5;p_6}:N-S-N",
    L"W{a_5;;p_6}:N-N-N", {});
}

TEST_CASE("deriv/rules/real", "[deriv][rules]") {
  namespace t = sequant::deriv::test;
  const auto R = [](std::wstring_view f, std::wstring_view x,
                    std::wstring_view expected, const t::Data& data) {
    t::check_rule(f, x, expected, data);
  };
  // every space real
  const auto resetter = real_spaces();
  const t::Data real{.real = {L"t", L"f", L"w", L"g", L"G", L"K", L"h", L"W",
                              L"V", L"v", L"s", L"u", L"Y", L"Z"}};

  // a real t: t⁺ holds the same elements
  R(L"f{i_1;a_1}:N-N-N t{a_1;i_1}:N-N-N + w{a_1;i_1}:N-N-N t⁺{i_1;a_1}:N-N-N",
    L"t{a_2;i_2}:N-N-N", L"f{i_2;a_2}:N-N-N + w{a_2;i_2}:N-N-N", real);
  R(L"t⁺{i_1;a_1}:N-N-N g{a_1,i_2;i_1,a_2}:N-N-N t{a_2;i_2}:N-N-N",
    L"t{a_3;i_3}:N-N-N",
    L"t⁺{i_1;a_1}:N-N-N g{a_1,i_3;i_1,a_3}:N-N-N + "
    L"g{a_3,i_2;i_3,a_2}:N-N-N t{a_2;i_2}:N-N-N",
    real);
  // each term carries the projector in the layout of its own spelling
  R(L"1/4 g{a_1,a_2;i_1,i_2}:A-N-S t⁺{i_1,i_2;a_1,a_2}:A-N-S + 1/4 "
    L"w{i_1,i_2;a_1,a_2}:A-N-S t{a_1,a_2;i_1,i_2}:A-N-S",
    L"t{a_3,a_4;i_3,i_4}:A-N-S",
    L"Â{a_3,a_4;i_3,i_4}:A-N-S w{i_3,i_4;a_3,a_4}:A-N-S + "
    L"Â{i_3,i_4;a_3,a_4}:A-N-S g{a_3,a_4;i_3,i_4}:A-N-S",
    real);
  // bra-ket symmetric over one space: G + Gᵀ
  R(L"G{p_2;p_1}:N-N-N h{p_1;p_2}:N-S-N", L"h{p_3;p_4}:N-S-N",
    L"G{p_4;p_3}:N-N-N + G{p_3;p_4}:N-N-N", real);
  // bra-ket symmetric over two spaces: one spelling, no second term
  R(L"G{i_1;a_1}:N-N-N h{a_1;i_1}:N-S-N + K{a_1;i_1}:N-N-N h{i_1;a_1}:N-S-N",
    L"h{a_2;i_2}:N-S-N", L"G{i_2;a_2}:N-N-N + K{a_2;i_2}:N-N-N", real);
  // bra-ket symmetric over one space, with a projector
  R(L"1/4 W{p_3,p_4;p_1,p_2}:N-N-N g{p_1,p_2;p_3,p_4}:A-S-S",
    L"g{p_5,p_6;p_7,p_8}:A-S-S",
    L"Â{p_5,p_6;p_7,p_8}:A-N-S W{p_7,p_8;p_5,p_6}:N-N-N + "
    L"Â{p_7,p_8;p_5,p_6}:A-N-S W{p_5,p_6;p_7,p_8}:N-N-N",
    real);
  // a leaf over the union space: it holds both spellings of a Hermitian
  // variable over two spaces, and one of a variable without bra-ket symmetry
  R(L"G{p_2;p_1}:N-N-N h{p_1;p_2}:N-C-N", L"h{i_3;a_4}:N-C-N",
    L"G{a_4;i_3}:N-N-N + G{i_3;a_4}:N-N-N", real);
  R(L"W{p_2;p_1}:N-N-N t{p_1;p_2}:N-N-N", L"t{a_3;i_3}:N-N-N",
    L"W{i_3;a_3}:N-N-N", real);
  // a bra-ket symmetric variable with a union-space slot, the leaf spelled
  // otherwise
  R(L"W{a_3,i_4;p_1,a_2}:N-N-N g{p_1,a_2;a_3,i_4}:N-S-S",
    L"g{p_5,a_6;a_7,i_8}:N-S-S",
    L"W{a_6,i_4;i_8,a_7}:N-N-N δ{p_5;i_4}:N-S-S + W{a_7,i_8;p_5,a_6}:N-N-N",
    real);
  // Hermitian over real spaces: the adjoint spelling holds the same elements
  R(L"G{p_2;p_1}:N-N-N h{p_1;p_2}:N-C-N", L"h{p_3;p_4}:N-C-N",
    L"G{p_4;p_3}:N-N-N + G{p_3;p_4}:N-N-N", real);
  R(L"W{i_1,a_1;a_2,i_2}:N-N-N g{a_2,i_2;i_1,a_1}:N-C-N + "
    L"V{a_1,i_1;i_2,a_2}:N-N-N g{i_2,a_2;a_1,i_1}:N-C-N",
    L"g{i_3,a_3;a_4,i_4}:N-C-N",
    L"V{a_4,i_4;i_3,a_3}:N-N-N + W{i_3,a_3;a_4,i_4}:N-N-N", real);
  // bra-ket symmetric with unequal bra and ket ranks: a leaf written either
  // way holds the variable's elements
  R(L"W{a_1,a_2;i_1}:N-N-N g{i_1;a_1,a_2}:N-S-N", L"g{i_5;a_5,a_6}:N-S-N",
    L"W{a_5,a_6;i_5}:N-N-N", real);
  R(L"W{a_1,a_2;i_1}:N-N-N g{i_1;a_1,a_2}:N-S-N + "
    L"V{i_1;a_1,a_2}:N-N-N g{a_1,a_2;i_1}:N-S-N",
    L"g{i_5;a_5,a_6}:N-S-N", L"W{a_5,a_6;i_5}:N-N-N + V{i_5;a_5,a_6}:N-N-N",
    real);
  // a CSV amplitude with t⁺ leaves
  R(L"1/4 g{a_1<i_1,i_2>,a_2<i_1,i_2>;i_1,i_2}:A-N-S "
    L"t⁺{i_1,i_2;a_1<i_1,i_2>,a_2<i_1,i_2>}:A-N-S + "
    L"1/4 w{i_1,i_2;a_1<i_1,i_2>,a_2<i_1,i_2>}:A-N-S "
    L"t{a_1<i_1,i_2>,a_2<i_1,i_2>;i_1,i_2}:A-N-S",
    L"t{a_3<i_3,i_4>,a_4<i_3,i_4>;i_3,i_4}:A-N-S",
    L"Â{a_3<i_3,i_4>,a_4<i_3,i_4>;i_3,i_4}:A-N-S "
    L"w{i_3,i_4;a_3<i_3,i_4>,a_4<i_3,i_4>}:A-N-S + "
    L"Â{i_3,i_4;a_3<i_3,i_4>,a_4<i_3,i_4>}:A-N-S "
    L"g{a_3<i_3,i_4>,a_4<i_3,i_4>;i_3,i_4}:A-N-S",
    real);
  // a second derivative by a real variable
  {
    const auto E = deserialize(
        L"t⁺{i_1;a_1}:N-N-N g{a_1,i_2;i_1,a_2}:N-N-N t{a_2;i_2}:N-N-N");
    const auto x1 = deserialize(L"t{a_3;i_3}:N-N-N"),
               x2 = deserialize(L"t{a_5;i_5}:N-N-N");
    const auto first = sequant::deriv::diff(E, {x1});
    t::check_rule(first, x2,
                  deserialize(L"g{a_3,i_5;i_3,a_5}:N-N-N + "
                              L"g{a_5,i_3;i_5,a_3}:N-N-N"),
                  real);
    REQUIRE_THAT(sequant::deriv::diff(E, {x1, x2}),
                 EquivalentTo(sequant::deriv::diff(first, {x2})));
  }
  // a scalar in the term is left as it is
  R(L"ε w{a_1;i_1}:N-N-N t⁺{i_1;a_1}:N-N-N + ε^(2) f{i_1;a_1}:N-N-N "
    L"t{a_1;i_1}:N-N-N",
    L"t{a_2;i_2}:N-N-N", L"ε w{a_2;i_2}:N-N-N + ε^(2) f{i_2;a_2}:N-N-N", real);
  // a residual with its own Â
  R(L"Â{i_1;a_1}:A-N-S w{a_1,a_2;i_1,i_2}:N-N-N t⁺{i_2;a_2}:N-N-N + "
    L"Â{i_1;a_1}:A-N-S v{a_1,i_2;i_1,a_2}:N-N-N t{a_2;i_2}:N-N-N",
    L"t{a_3;i_3}:N-N-N",
    L"Â{i_1;a_1}:A-N-S v{a_1,i_3;i_1,a_3}:N-N-N + "
    L"Â{i_1;a_1}:A-N-S w{a_1,a_3;i_1,i_3}:N-N-N",
    real);
  // no bra-ket symmetry: the adjoint is another label, so no slots of the
  // variable are related (the judge alone decides, and the result is not zero)
  for (const auto& [f, x] :
       {std::pair{L"Y{;i_3,i_4}:N-N-N W{a_3<i_3>;a_4<i_4>}:N-N-N "
                  L"s{a_4<i_4>;a_3<i_3>}:N-N-N Z{i_3,i_4;}:N-N-N",
                  L"s{a_5<i_5>;a_6<i_6>}:N-N-N"},
        std::pair{L"W{i_1,i_2;p_1,p_2}:N-N-N u{p_1,p_2;i_1,i_2}:N-N-N",
                  L"u{p_5,p_6;i_5,i_6}:N-N-N"},
        std::pair{
            L"W{i_2,a_2;i_1,a_1<i_1>}:N-N-N w{i_1,a_1<i_1>;i_2,a_2}:N-N-N",
            L"w{i_3,a_3<i_3>;i_4,a_4}:N-N-N"}}) {
    R(f, x, L"", real);
    REQUIRE_FALSE(
        sequant::deriv::diff(deserialize(f), {deserialize(x)})->is_zero());
  }
}

TEST_CASE("deriv/rules/spellings", "[deriv][rules]") {
  namespace t = sequant::deriv::test;
  const auto diff = [](const ExprPtr& E, const ExprPtr& x) {
    return sequant::deriv::diff(E, {x});
  };
  auto T = [](std::wstring_view s) { return deserialize(s); };
  const auto R = [](std::wstring_view f, std::wstring_view x,
                    std::wstring_view expected, const t::Data& data = {}) {
    t::check_rule(f, x, expected, data);
  };
  // the judge alone decides; the result is a sum of n terms
  const auto R_terms = [&](std::wstring_view f, std::wstring_view x,
                           std::size_t n, const t::Data& data = {}) {
    t::check_rule(f, x, L"", data);
    const auto r = diff(T(f), T(x));
    REQUIRE(r->is<Sum>());
    REQUIRE(r->as<Sum>().size() == n);
  };

  // a leaf over a union space meets several spellings of the variable: a
  // term per spelling
  R(L"G{p_2;p_1}:N-N-N h{p_1;p_2}:N-S-N", L"h{i_3;a_4}:N-S-N",
    L"G{i_3;a_4}:N-N-N + G{a_4;i_3}:N-N-N");
  R(L"1/4 W{p_3,p_4;p_1,p_2}:N-N-N g{p_1,p_2;p_3,p_4}:A-N-S",
    L"g{i_5,a_5;i_6,a_6}:A-N-S",
    L"1/4 W{i_6,a_6;i_5,a_5}:N-N-N - 1/4 W{a_6,i_6;i_5,a_5}:N-N-N - "
    L"1/4 W{i_6,a_6;a_5,i_5}:N-N-N + 1/4 W{a_6,i_6;a_5,i_5}:N-N-N");
  R(L"W{i_1,i_2;p_1,p_2}:N-N-N g{p_1,p_2;i_1,i_2}:N-N-S",
    L"g{a_3,i_5;i_3,i_4}:N-N-S",
    L"W{i_4,i_3;i_5,a_3}:N-N-N + W{i_3,i_4;a_3,i_5}:N-N-N");
  R(L"W{p_3,p_4;p_1,p_2}:N-N-N g{p_1,p_2;p_3,p_4}:N-C-S",
    L"g{i_3,a_3;a_4,i_4}:N-C-S",
    L"W{a_4,i_4;i_3,a_3}:N-N-N + W{i_4,a_4;a_3,i_3}:N-N-N");
  R(L"W{p_3,p_4;p_1,p_2}:N-N-N g{p_1,p_2;p_3,p_4}:N-S-S",
    L"g{i_3,a_3;a_4,i_4}:N-S-S",
    L"W{a_4,i_4;i_3,a_3}:N-N-N + W{a_3,i_3;i_4,a_4}:N-N-N + "
    L"W{i_3,a_3;a_4,i_4}:N-N-N + W{i_4,a_4;a_3,i_3}:N-N-N");
  // a leaf that meets one spelling, not the one the variable is written in
  R(L"W{p_2,a_2;p_1,a_1}:N-N-N g{p_1,a_1;p_2,a_2}:A-N-S",
    L"g{i_5,a_5;i_6,a_6}:A-N-S", L"W{i_6,a_6;i_5,a_5}:N-N-N");
  R(L"G{a_1;p_1}:N-N-N h{p_1;a_1}:N-S-N", L"h{i_3;a_3}:N-S-N",
    L"G{a_3;i_3}:N-N-N");
  // variables whose operations relate slots that are not alike
  for (const auto x :
       {L"t{a_1<i_1>,a_2<i_2>;i_1,i_2}:A-N-S", L"h{p_1;i_1}:N-S-N",
        L"t{p_1,i_1;a_1,a_2}:A-N-S", L"t{a_1<i_1>,a_2<i_1>;i_1,i_2}:A-N-S"})
    R(L"f{i_1;a_1}:N-N-N t{a_1;i_1}:N-N-N", x, L"0");
  R(L"1/4 W{a_1,a_2;p_1,p_2}:N-N-N t{p_1,p_2;a_1,a_2}:A-N-S",
    L"t{p_3,i_3;a_3,a_4}:A-N-S",
    L"1/2 Â{;a_3,a_4}:A-N-S W{a_3,a_4;p_3,i_3}:N-N-N - "
    L"1/2 Â{;a_3,a_4}:A-N-S W{a_3,a_4;i_3,p_3}:N-N-N");
  // a bra-ket symmetric variable with a union-space slot, the leaf spelled
  // otherwise (over real spaces in deriv/rules/real)
  R(L"W{a_3,i_4;p_1,a_2}:N-N-N g{p_1,a_2;a_3,i_4}:N-S-S",
    L"g{p_5,a_6;a_7,i_8}:N-S-S",
    L"W{a_6,i_4;i_8,a_7}:N-N-N δ{p_5;i_4}:N-S-S + W{a_7,i_8;p_5,a_6}:N-N-N",
    {.real = {L"g", L"W"}});
  // a bra-ket exchange with like slots on both sides
  R(L"W{a_4,a_5;i_1,i_2}:N-N-N g{i_1,i_2;a_4,a_5}:A-S-S",
    L"g{a_12,a_13;i_14,i_15}:A-S-S",
    L"4 Â{i_14,i_15;a_12,a_13}:A-N-S W{a_12,a_13;i_14,i_15}:N-N-N");
  // a three-slot cycle
  R_terms(L"W{p_4;p_1,p_2,p_3}:N-N-N g{p_1,p_2,p_3;p_4}:A-N-S",
          L"g{i_5,a_6,p_7;a_8}:A-N-S", 6);
  // a union-space leaf with free indices
  R(L"h{p_1;p_2}:N-S-N", L"h{i_3;a_4}:N-S-N",
    L"δ{p_1;i_3}:N-S-S δ{a_4;p_2}:N-S-S + δ{p_1;a_4}:N-S-S δ{i_3;p_2}:N-S-S");
  // a variable with a union-space slot: first and second derivative
  const t::Data h_data{.real = {L"h"}};
  R(L"G{p_2;p_1}:N-N-N h{p_1;p_2}:N-S-N", L"h{p_3;i_4}:N-S-N",
    L"G{i_4;p_3}:N-N-N + G{p_3;i_4}:N-N-N", h_data);
  {
    const auto E = T(L"1/2 h{p_1;p_2}:N-S-N h{p_2;p_1}:N-S-N");
    const auto x1 = T(L"h{p_3;i_3}:N-S-N"), x2 = T(L"h{p_5;i_5}:N-S-N");
    t::check_rule(E, x1, T(L"h{i_3;p_3}:N-S-N + h{p_3;i_3}:N-S-N"), h_data);
    const auto first = diff(E, x1);
    t::check_rule(first, x2, nullptr, h_data);
    REQUIRE_FALSE(sequant::deriv::diff(E, {x1, x2})->is_zero());
  }
}

TEST_CASE("deriv/rules/mixed-field", "[deriv][rules]") {
  namespace t = sequant::deriv::test;
  // the occupied space real, the others complex
  const auto resetter = real_spaces({L"i"});
  const t::Data data{.real = {L"h"}};
  const auto E =
      L"f{i_2;i_1}:N-N-N h{i_1;i_2}:N-N-N + w{i_1;i_2}:N-N-N h⁺{i_1;i_2}:N-N-N "
      L"+ g{i_1;a_1}:N-N-N t{a_1;i_1}:N-N-N + v{a_1;i_1}:N-N-N "
      L"t⁺{i_1;a_1}:N-N-N";
  // over real spaces only: the adjoint spelling holds the same elements
  t::check_rule(E, L"h{i_3;i_4}:N-N-N", L"f{i_4;i_3}:N-N-N + w{i_4;i_3}:N-N-N",
                data);
  // over a complex space: t and t⁺ are independent
  t::check_rule(E, L"t{a_2;i_2}:N-N-N", L"g{i_2;a_2}:N-N-N", data);
}

TEST_CASE("deriv/rules/rejections", "[deriv][rules]") {
  using sequant::deriv::test::check_rule;
  const auto diff = [](const ExprPtr& E, const ExprPtr& x) {
    return sequant::deriv::diff(E, {x});
  };
  auto T = [](std::wstring_view s) { return deserialize(s); };
  const auto X = [&](const ExprPtr& E, const ExprPtr& x, std::string text) {
    REQUIRE_THROWS_WITH(diff(E, x),
                        Catch::Matchers::ContainsSubstring(std::move(text)));
  };

  const auto E = T(L"f{i_1;a_1}:N-N-N t{a_1;i_1}:N-N-N");
  X(E, T(L"x + y"),
    "deriv::diff: the variable is not a Variable or a tensor leaf");
  X(E, T(L"Â{i_1;a_1}:A-N-S"),
    "the variable is not a Variable or a tensor leaf");
  X(E, make_kronecker(Index(L"a_2"), Index(L"a_3")),
    "the variable is not a Variable or a tensor leaf");
  // (core itself refuses an index repeated within the bra or within the ket)
  X(E, T(L"h{p_1;p_1}:N-N-N"), "the variable repeats an index");
  X(E, T(L"x{a_1<i_1,i_2>;i_1}:N-N-N"),
    "an index of the variable mixes protoindices that are its slots with "
    "others");
  X(E, T(L"x{a_1<i_1,i_2>;a_2<i_2,i_3>}:N-N-N"),
    "the protoindex lists of the variable's slots share some but not all "
    "labels");
  X(E, T(L"t{a_1,a_2;i_1,i_2}:S-N-S"),
    "a variable with Symmetry::Symm is not supported");
  X(E, T(L"w{i_1,i_2,a_1,a_2;i_3,i_4,a_3,a_4}:A-N-S"),
    "the variable has two groups of like slots");
  X(mbpt::op::H(1), T(L"t{a_2;i_2}:N-N-N"), "the function holds operators");
  // a hyperindex, by the scalar rule (by a tensor: deriv/rules/hadamard)
  REQUIRE_THAT(
      diff(T(L"x^(2) B{a_1;i_1;i_3}:N-N-N C{i_1;a_1;i_3}:N-N-N D{;;i_3}:N-N-N"),
           ex<Variable>(L"x")),
      EquivalentTo(
          L"2 x B{a_1;i_1;i_3}:N-N-N C{i_1;a_1;i_3}:N-N-N D{;;i_3}:N-N-N"));
  X(nullptr, ex<Variable>(L"x"), "deriv::diff: null argument");
  X(E, nullptr, "deriv::diff: null argument");
  REQUIRE_THROWS_WITH(
      sequant::deriv::at_zero(nullptr, {ex<Variable>(L"x")}),
      Catch::Matchers::ContainsSubstring("deriv::at_zero: null argument"));
  REQUIRE_THROWS_WITH(
      sequant::deriv::at_zero(E, {nullptr}),
      Catch::Matchers::ContainsSubstring("deriv::at_zero: null argument"));
  X(T(L"t{a_1,a_2;i_1,i_2}:N-N-S"), T(L"t{a_3,a_4;i_3,i_4}:A-N-S"),
    "a leaf declares other symmetries than the variable");
  // the leaf's own symmetry spells it in another slot order than the variable's
  for (const auto& [f, v] :
       {std::pair{L"W{a_1;i_1}:N-N-N f{i_1;a_1}:N-N-N + "
                  L"V{a_1;i_1}:N-N-N f{i_1;a_1}:N-S-N",
                  L"f{i_3;a_3}:N-N-N"},
        std::pair{L"W{a_1;i_1}:N-N-N h{i_1;a_1}:N-S-N", L"h{i_3;a_3}:N-N-N"},
        std::pair{L"W{a_1;i_1}:N-N-N h{i_1;a_1}:N-N-N", L"h{a_3;i_3}:N-S-N"},
        std::pair{L"W{i_1;a_1}:N-N-N h{a_1;i_1}:N-S-N", L"h{i_3;a_3}:N-N-N"},
        std::pair{L"g{a_1,i_1;a_2,i_2}:N-N-N", L"g{a_3,i_3;a_4,i_4}:A-N-S"},
        std::pair{L"g{a_1,i_1;a_2,i_2}:A-N-S", L"g{a_3,i_3;a_4,i_4}:N-N-N"}})
    X(T(f), T(v), "a leaf declares other symmetries than the variable");
  // the leaf's pair labels are not its slots, the variable's are
  X(T(L"g{i_1,i_2;a_1<i_3,i_4>,a_2<i_3,i_4>}:N-N-N "
      L"t{a_1<i_3,i_4>,a_2<i_3,i_4>;i_1,i_2}:N-N-N"),
    T(L"t{a_5<i_5,i_6>,a_6<i_5,i_6>;i_5,i_6}:N-N-N"),
    "a leaf's protoindices do not follow its slots as the variable's do");
  // the variable names one pair where the leaf has two
  X(T(L"u{;i_1,i_2,i_3,i_4}:N-N-N M{a_3<i_3,i_4>;a_1<i_1,i_2>}:N-N-N "
      L"k{a_1<i_1,i_2>;a_3<i_3,i_4>}:N-N-N v{i_1,i_2,i_3,i_4;}:N-N-N"),
    T(L"k{a_2<i_5,i_6>;a_4<i_5,i_6>}:N-N-N"),
    "a leaf's protoindices do not follow its slots as the variable's do");
  // ket slots in different protoindex lists are not alike: the exchange of
  // the variable's ket slots is a spelling the leaf meets
  X(T(L"W{i_1,i_2;a_1<i_1>,a_2<i_1,i_2>}:N-N-N "
      L"t{a_1<i_1>,a_2<i_1,i_2>;i_1,i_2}:A-N-S"),
    T(L"t{a_5<i_5>,a_6<i_5,i_6>;i_5,i_6}:A-N-S"),
    "a leaf's protoindices do not follow its slots as the variable's do");
  X(T(L"Â{i_1;a_1}:A-N-S v{a_1,i_2;i_1,a_2}:N-N-N t{a_2;i_2}:N-N-N"),
    T(L"t{a_1;i_1}:N-N-N"),
    "the variable names an external index of the function");
  // the same spelling twice: the first result has x's indices as external
  const auto x = T(L"t{a_2;i_2}:N-N-N");
  X(diff(T(L"g{i_1,i_2;a_1,a_2}:N-N-N t{a_1;i_1}:N-N-N t{a_2;i_2}:N-N-N"), x),
    x, "the variable names an external index of the function");
  // a free index that appears only as a protoindex is external, too
  X(T(L"w{;a_5<i_2>}:N-N-N v{a_5<i_2>;}:N-N-N f{i_1;a_1}:N-N-N "
      L"t{a_1;i_1}:N-N-N"),
    x, "the variable names an external index of the function");
  // a pair label of a leaf of the variable that is not summed, or summed over
  // less than the variable's label runs over
  for (const auto& [f, v] :
       {std::pair{L"Y{;i_1}:N-N-N W{;a_1<i_1>}:N-N-N k{a_1<i_1>;}:N-N-N",
                  L"k{e_2<i_2>;}:N-N-N"},
        std::pair{L"h{a_1<i_1>;a_1<i_1>}:N-N-N", L"h{e_2<i_2>;a_3<i_4>}:N-N-N"},
        std::pair{L"k{a_1<i_1,i_2>;a_3<i_1,i_2>}:N-N-N "
                  L"Z{a_3<i_1,i_2>;a_1<i_1,i_2>}:N-N-N",
                  L"k{a_2<i_5,i_6>;a_4<i_7,i_8>}:N-N-N"},
        std::pair{L"Y{;i_1}:N-N-N W{;a_1<i_1>}:N-N-N k{a_1<i_1>;}:N-N-N "
                  L"Z{i_1;}:N-N-N",
                  L"k{a_2<m_2>;}:N-N-N"}})
    X(T(f), T(v),
      "a protoindex label of a leaf of the variable must be summed over the "
      "space of the variable's label");
  // protoindex patterns the symmetry of the variable would not keep
  const Index o1(L"i_1"), o2(L"i_2");
  const Index a1(Index(L"a_1"), {o1, o2}, false),
      a2(Index(L"a_2"), {o1, o2}, false);
  X(E, ex<Tensor>(L"t", bra{a1, a2}, ket{o1, o2}, Symmetry::Antisymm),
    "ordered protoindices are not supported");
  // ordered lists that share some labels
  X(E,
    ex<Tensor>(L"x", bra{a1},
               ket{Index(Index(L"a_2"), {o2, Index(L"i_3")}, false)}),
    "ordered protoindices are not supported");
  X(E, T(L"r{a_1;i_1,i_2,i_3}:N-N-S"),
    "a column-symmetric variable with unequal ranks is not supported");
  // with at most one slot per side there is nothing to permute
  REQUIRE_THAT(sequant::deriv::diff(T(L"e{;a_1}:N-N-S + e{i_1;}:N-N-S"),
                                    {T(L"e{;a_2}:N-N-S")}),
               EquivalentTo(L"δ{a_2;a_1}"));
  // the layout is checked against every spelling of the variable the leaf
  // meets: here the one with the bra slots exchanged
  X(T(L"W{i_1,i_2;a_1<i_1>,a_2<i_2>}:N-N-N "
      L"t{a_1<i_1>,a_2<i_2>;i_1,i_2}:A-N-S"),
    T(L"t{a_5<i_5>,a_6<i_6>;i_5,i_6}:A-N-S"),
    "a leaf's protoindices do not follow its slots as the variable's do");

  // what is accepted: a TensorPower of another tensor is a constant
  check_rule(ex<TensorPower>(T(L"D{i_1;a_1}:N-N-N"), rational(-1)) *
                 T(L"r{a_1;i_1}:N-N-N"),
             T(L"r{a_2;i_2}:N-N-N"),
             ex<TensorPower>(T(L"D{i_2;a_2}:N-N-N"), rational(-1)));
  // columns related only where both their bra and their ket spaces overlap
  check_rule(T(L"W{a_1,i_2;p_1,i_1}:N-N-N w{p_1,i_1;a_1,i_2}:N-N-S"),
             T(L"w{p_2,i_3;a_2,i_4}:N-N-S"), nullptr);
  REQUIRE_FALSE(diff(T(L"W{a_1,i_2;p_1,i_1}:N-N-N w{p_1,i_1;a_1,i_2}:N-N-S"),
                     T(L"w{p_2,i_3;a_2,i_4}:N-N-S"))
                    ->is_zero());
  // the slots of the variable carry different pair labels, the leaf's one
  check_rule(T(L"W{;i_1}:N-N-N h{a_1<i_1>;a_1<i_1>}:N-N-N V{i_1;}:N-N-N"),
             T(L"h{e_2<i_2>;a_3<i_4>}:N-N-N"), nullptr);
  // neither the order of the variable's pair labels nor the name of a summed
  // index decides whether a pair label is renamed
  const auto R_nonzero = [&](const ExprPtr& f, const ExprPtr& x,
                             const ExprPtr& expected) {
    check_rule(f, x, expected);
    REQUIRE_FALSE(diff(f, x)->is_zero());
  };
  R_nonzero(T(L"W{;i_1}:N-N-N h{a_1<i_1>;a_1<i_1>}:N-N-N V{i_1;}:N-N-N"),
            T(L"h{e_2<i_2>;a_3<m_4>}:N-N-N"), nullptr);
  R_nonzero(T(L"W{;i_1}:N-N-N h{a_1<i_1>;a_1<i_1>}:N-N-N V{i_1;}:N-N-N"),
            T(L"h{e_2<m_2>;a_3<i_4>}:N-N-N"),
            T(L"W{;i_4}:N-N-N V{i_4;}:N-N-N δ{a_3<i_4>;e_2<m_2>} δ{i_4;m_2}"));
  R_nonzero(T(L"W{;i_2}:N-N-N h{a_1<i_2>;a_1<i_2>}:N-N-N V{i_2;}:N-N-N"),
            T(L"h{e_2<i_2>;a_3<m_4>}:N-N-N"), nullptr);
  // a protoindex label of an external index may stand in more than two slots,
  // as in this first derivative, which is differentiated again
  {
    const auto F =
        T(L"λ{i_1;a_5<i_1>}:N-N-N s{a_5<i_1>;a_6<i_1,i_2>}:N-N-N "
          L"t{a_6<i_1,i_2>;i_1,i_2}:N-N-N u{i_2;}:N-N-N");
    const auto x1 = T(L"s{a_51<i_51>;a_52<i_52,i_53>}:N-N-N");
    const auto first = diff(F, x1);
    for (const auto& x2 :
         {T(L"u{i_71;}:N-N-N"), T(L"λ{i_71;a_71<i_71>}:N-N-N")}) {
      R_nonzero(first, x2, nullptr);
      REQUIRE_THAT(sequant::deriv::diff(F, {x1, x2}),
                   EquivalentTo(diff(first, x2)));
    }
  }
  // bra and ket with different pair labels, bra-ket symmetric: the exchanged
  // element is spelled the other way, the adjoint spelling's term
  check_rule(T(L"Y{;i_3,i_4}:N-N-N W{a_3<i_3>;a_4<i_4>}:N-N-N "
               L"h{a_4<i_4>;a_3<i_3>}:N-S-N Z{i_3,i_4;}:N-N-N"),
             T(L"h{a_5<i_5>;a_6<i_6>}:N-S-N"), nullptr,
             sequant::deriv::test::Data{.real = {L"h"}});
  REQUIRE_FALSE(diff(T(L"Y{;i_3,i_4}:N-N-N W{a_3<i_3>;a_4<i_4>}:N-N-N "
                       L"h{a_4<i_4>;a_3<i_3>}:N-S-N Z{i_3,i_4;}:N-N-N"),
                     T(L"h{a_5<i_5>;a_6<i_6>}:N-S-N"))
                    ->is_zero());
}

TEST_CASE("deriv/rules/powers-and-definitions", "[deriv][rules]") {
  using sequant::deriv::diff;
  using sequant::deriv::test::compare_tangents;
  using sequant::deriv::test::require_same;
  auto T = [](std::wstring_view s) { return deserialize(s); };
  const auto P = [&](std::wstring_view s, int e) {
    return pow(T(s), rational(e));
  };
  // an elementwise derivative is a Hadamard product: its free indices stand
  // in two slots, which core's strict bra-ket check rejects
  auto ctx = get_default_context();
  ctx.set(AssertStrictBraKetSymmetry::No);
  const auto scope = set_scoped_default_context(ctx);
  // the head d{i;a} of a definition, nonsymmetric
  const auto d = [](std::wstring_view i, std::wstring_view a) {
    return Tensor(L"d", bra{i}, ket{a}, Symmetry::Nonsymm);
  };
  // d{i;a} = e{;a} - e{i;}, a broadcast sum
  const container::svector<ResultExpr> defs{
      ResultExpr(d(L"i_9", L"a_9"), T(L"e{;a_9}:N-N-S - e{i_9;}:N-N-S"))};

  SECTION("power rule") {
    // elementwise: d(X^e)/dX = e X^(e-1), at the variable's element; checked
    // against the evaluator, as check_rule() does, but for the slot counts of
    // a Hadamard product
    const auto R = [](const ExprPtr& f, const ExprPtr& x,
                      const ExprPtr& expected) {
      const auto r = diff(f, {x});
      if (expected) REQUIRE_THAT(r, EquivalentTo(expected));
      const auto c = compare_tangents(f, r, x, {}, {.count = 256, .seed = 1});
      INFO(c.difference);
      REQUIRE(c.difference.empty());
      REQUIRE(c.nonzero != 0);
    };
    R(P(L"X{i_1;a_1}:N-N-N", -1) * T(L"r{a_1;i_1}:N-N-N"),
      T(L"X{i_2;a_2}:N-N-N"),
      ex<Constant>(-1) * P(L"X{i_2;a_2}:N-N-N", -2) * T(L"r{a_2;i_2}:N-N-N"));
    R(P(L"X{i_1;a_1}:N-N-N", 3) * T(L"r{a_1;i_1}:N-N-N"),
      T(L"X{i_2;a_2}:N-N-N"),
      ex<Constant>(3) * P(L"X{i_2;a_2}:N-N-N", 2) * T(L"r{a_2;i_2}:N-N-N"));
    // an antisymmetric base: the spellings of the variable are summed by Â
    R(P(L"t{a_1,a_2;i_1,i_2}:A-N-S", 2) * T(L"r{i_1,i_2;a_1,a_2}:N-N-N"),
      T(L"t{a_3,a_4;i_3,i_4}:A-N-S"), nullptr);
  }

  SECTION("chain rule through a definition") {
    // a leaf defined alone: as the expression substituted for it
    const auto x = T(L"X{a_1;i_1}:N-N-N");
    require_same(
        diff(x * T(L"d{i_1;a_1}:N-N-N"), {T(L"e{;a_2}:N-N-S")}, defs),
        diff(x * T(L"e{;a_1}:N-N-S - e{i_1;}:N-N-S"), {T(L"e{;a_2}:N-N-S")}),
        {Index(L"a_2")}, {});
    // under a TensorPower: e d^(e-1) times the derivative of the definition
    REQUIRE_THAT(
        diff(x * P(L"d{i_1;a_1}:N-N-N", -1), {T(L"e{;a_2}:N-N-S")}, defs),
        EquivalentTo(ex<Constant>(-1) * T(L"X{a_2;i_1}:N-N-N") *
                     P(L"d{i_1;a_2}:N-N-N", -2)));
    REQUIRE_THAT(
        diff(x * P(L"d{i_1;a_1}:N-N-N", -1), {T(L"e{i_2;}:N-N-S")}, defs),
        EquivalentTo(T(L"X{a_1;i_2}:N-N-N") * P(L"d{i_2;a_1}:N-N-N", -2)));
    // the leaf's own variable is not reached through the definition
    REQUIRE_THAT(
        diff(x * P(L"d{i_1;a_1}:N-N-N", -1), {T(L"X{a_2;i_2}:N-N-N")}, defs),
        EquivalentTo(P(L"d{i_2;a_2}:N-N-N", -1)));
    // without the definition, d is a constant
    REQUIRE(
        diff(x * P(L"d{i_1;a_1}:N-N-N", -1), {T(L"e{;a_2}:N-N-S")})->is_zero());
    // definitions nest: c{i;a} = s d{i;a}, by a scalar variable through both
    const container::svector<ResultExpr> nested{
        defs.front(),
        ResultExpr(Tensor(L"c", bra{L"i_8"}, ket{L"a_8"}, Symmetry::Nonsymm),
                   ex<Variable>(L"s") * T(L"d{i_8;a_8}:N-N-N"))};
    REQUIRE_THAT(
        diff(x * P(L"c{i_1;a_1}:N-N-N", 2), {T(L"e{;a_2}:N-N-S")}, nested),
        EquivalentTo(ex<Constant>(2) * T(L"X{a_2;i_1}:N-N-N") *
                     T(L"c{i_1;a_2}:N-N-N") * ex<Variable>(L"s")));
    REQUIRE_THAT(diff(x * T(L"c{i_1;a_1}:N-N-N"), {ex<Variable>(L"s")}, nested),
                 EquivalentTo(x * T(L"d{i_1;a_1}:N-N-N")));
  }

  SECTION("definitions written like the function") {
    const auto E = T(L"X{a_1;i_1}:N-N-N") * P(L"d{i_1;a_1}:N-N-N", -1);
    const auto expected =
        ex<Constant>(-1) * T(L"X{a_2;i_1}:N-N-N") * P(L"d{i_1;a_2}:N-N-N", -2);
    // the head named like the variable: its names are dummies
    REQUIRE_THAT(diff(E, {T(L"e{;a_2}:N-N-S")},
                      {ResultExpr(d(L"i_2", L"a_2"),
                                  T(L"e{;a_2}:N-N-S - e{i_2;}:N-N-S"))}),
                 EquivalentTo(expected));
    // an expression that is not expanded
    auto twice = ex<Product>(ExprPtrList{});
    twice->as<Product>().append(2, T(L"e{;a_9}:N-N-S - e{i_9;}:N-N-S"),
                                Product::Flatten::No);
    REQUIRE_THAT(
        diff(E, {T(L"e{;a_2}:N-N-S")}, {ResultExpr(d(L"i_9", L"a_9"), twice)}),
        EquivalentTo(ex<Constant>(2) * expected));
    // a summed index of the expression whose protoindex is a slot of the head
    REQUIRE_THAT(
        diff(T(L"X{a_1;i_1}:N-N-N d{i_1;a_1}:N-N-N"), {ex<Variable>(L"s")},
             {ResultExpr(d(L"i_9", L"a_9"),
                         T(L"s h{i_9;a_8<i_9>}:N-N-N k{a_8<i_9>;a_9}:N-N-N"))}),
        EquivalentTo(T(L"X{a_1;i_1}:N-N-N h{i_1;a_8<i_1>}:N-N-N "
                       L"k{a_8<i_1>;a_1}:N-N-N")));
  }

  SECTION("rejected definitions") {
    REQUIRE_THROWS_WITH(
        diff(T(L"X{a_1<i_1>;i_1}:N-N-N Z{i_1;}:N-N-N") *
                 P(L"d{i_1;a_1<i_1>}:N-N-N", -1),
             {T(L"e{;a_2<i_2>}:N-N-S")}, defs),
        "deriv::diff: a leaf with a definition has protoindices");
    const auto E = T(L"X{a_1;i_1}:N-N-N") * P(L"d{i_1;a_1}:N-N-N", -1);
    const auto e = T(L"e{;a_2}:N-N-S");
    const auto rejects = [&](const container::svector<ResultExpr>& defs,
                             std::string text) {
      REQUIRE_THROWS_WITH(diff(E, {e}, defs), "deriv::diff: " + text);
    };
    rejects({ResultExpr(d(L"i_9", L"a_9"), T(L"e{;a_9}:N-N-S")),
             ResultExpr(d(L"i_8", L"a_8"), T(L"e{i_8;}:N-N-S"))},
            "two definitions have one head");
    rejects(
        {ResultExpr(d(L"i_9", L"a_9"), T(L"c{i_9;a_9}:N-N-N")),
         ResultExpr(Tensor(L"c", bra{L"i_8"}, ket{L"a_8"}, Symmetry::Nonsymm),
                    T(L"d{i_8;a_8}:N-N-N"))},
        "a definition depends on itself");
    rejects({ResultExpr(Tensor(L"e", bra{}, ket{L"a_9"}, Symmetry::Nonsymm),
                        T(L"d{i_9;a_9}:N-N-N"))},
            "a definition's head has the label of the variable");
    rejects({ResultExpr(Variable(L"z"), T(L"e{;a_9}:N-N-S"))},
            "a definition's head is not a labeled tensor");
    rejects(
        {ResultExpr(Tensor(L"d", bra{L"i_9"}, ket{L"i_9"}, Symmetry::Nonsymm),
                    T(L"e{i_9;}:N-N-S"))},
        "a definition's head repeats an index");
    rejects({ResultExpr(d(L"i_9", L"a_9"), T(L"e{;a_8}:N-N-S"))},
            "a definition's expression has a free index its head lacks");
  }
}

TEST_CASE("deriv/rules/hadamard", "[deriv][rules]") {
  using sequant::deriv::diff;
  using sequant::deriv::test::group_order;
  using sequant::deriv::test::require_same;
  auto T = [](std::wstring_view s) { return deserialize(s); };
  auto ctx = get_default_context();
  ctx.set(AssertStrictBraKetSymmetry::No);
  const auto scope = set_scoped_default_context(ctx);
  // f homogeneous of @p degree in the variable x, every index summed:
  // Euler's identity, r = diff(f, x) contracted with x is degree times f
  // times the number of operations of x; the evaluator reads an index in
  // three slots as free, so the arrays are compared with their free labels
  // given (none)
  const auto R = [&](const ExprPtr& f, const ExprPtr& x,
                     const ExprPtr& expected, int degree = 1) {
    const auto r = diff(f, {x});
    if (expected) REQUIRE_THAT(r, EquivalentTo(expected));
    require_same(r * x,
                 ex<Constant>(degree * static_cast<int>(group_order(x))) * f,
                 {}, {});
  };
  // an index in three slots is summed: a Hadamard product under a contraction
  R(T(L"B{a_1;i_1;i_3}:N-N-N C{i_1;a_1;i_3}:N-N-N D{;;i_3}:N-N-N"),
    T(L"B{a_2;i_2;i_4}:N-N-N"), T(L"C{i_2;a_2;i_4}:N-N-N D{;;i_4}:N-N-N"));
  // the shape of a perturbative correction, L ⊙ D⁻² · R, every index summed
  // over three slots (nonsymmetric: the evaluator applies no symmetrizer to
  // a Hadamard product)
  R(T(L"L{i_1,i_2;a_1,a_2}:N-N-N R{a_1,a_2;i_1,i_2}:N-N-N") *
        pow(T(L"D{i_1,i_2;a_1,a_2}:N-N-N"), -2),
    T(L"R{a_3,a_4;i_3,i_4}:N-N-N"),
    T(L"L{i_3,i_4;a_3,a_4}:N-N-N") * pow(T(L"D{i_3,i_4;a_3,a_4}:N-N-N"), -2));
  R(T(L"L{i_1,i_2;a_1,a_2}:N-N-N R{a_1,a_2;i_1,i_2}:N-N-N") *
        pow(T(L"D{i_1,i_2;a_1,a_2}:N-N-N"), -2),
    T(L"D{i_3,i_4;a_3,a_4}:N-N-N"),
    ex<Constant>(-2) * T(L"L{i_3,i_4;a_3,a_4}:N-N-N R{a_3,a_4;i_3,i_4}:N-N-N") *
        pow(T(L"D{i_3,i_4;a_3,a_4}:N-N-N"), -3),
    -2);
  // antisymmetric: the spellings of the variable are summed by Â
  REQUIRE_THAT(diff(T(L"L{i_1,i_2;a_1,a_2}:A-S-S R{a_1,a_2;i_1,i_2}:A-S-S") *
                        pow(T(L"D{i_1,i_2;a_1,a_2}:N-N-N"), -2),
                    {T(L"R{a_3,a_4;i_3,i_4}:A-S-S")}),
               EquivalentTo(ex<Constant>(4) * T(L"Â{a_3,a_4;i_3,i_4}:A-N-S") *
                            T(L"L{i_3,i_4;a_3,a_4}:A-S-S") *
                            pow(T(L"D{i_3,i_4;a_3,a_4}:N-N-N"), -2)));
  // a residual under Â whose indices stand in two other slots each: they are
  // its external indices, kept with Kronecker deltas
  REQUIRE_THAT(diff(T(L"Â{i_1,i_2;a_1,a_2}:A-N-S t{a_1,a_2;i_1,i_2}:A-S-S") *
                        pow(T(L"D{i_1,i_2;a_1,a_2}:N-N-N"), -1),
                    {T(L"t{a_3,a_4;i_3,i_4}:A-S-S")}),
               EquivalentTo(ex<Constant>(4) * T(L"Â{a_3,a_4;i_3,i_4}:A-N-S") *
                            T(L"Â{i_1,i_2;a_1,a_2}:A-N-S") *
                            pow(T(L"D{i_1,i_2;a_1,a_2}:N-N-N"), -1) *
                            T(L"δ{i_1;i_3}:N-S-S δ{i_2;i_4}:N-S-S "
                              L"δ{a_3;a_1}:N-S-S δ{a_4;a_2}:N-S-S")));
  // the indices of the variable stay free in a second derivative, though the
  // first holds them in two slots
  const auto f = pow(T(L"X{i_1;a_1}:N-N-N"), -1) * T(L"r{a_1;i_1}:N-N-N");
  REQUIRE_THAT(diff(f, {T(L"X{i_2;a_2}:N-N-N"), T(L"X{i_3;a_3}:N-N-N")}),
               EquivalentTo(ex<Constant>(2) * pow(T(L"X{i_2;a_2}:N-N-N"), -3) *
                            T(L"r{a_2;i_2}:N-N-N") * T(L"δ{i_2;i_3}:N-S-S") *
                            T(L"δ{a_3;a_2}:N-S-S")));
  REQUIRE_THAT(diff(f, {T(L"r{a_2;i_2}:N-N-N"), T(L"X{i_3;a_3}:N-N-N")}),
               EquivalentTo(ex<Constant>(-1) * pow(T(L"X{i_2;a_2}:N-N-N"), -2) *
                            T(L"δ{i_2;i_3}:N-S-S") * T(L"δ{a_3;a_2}:N-S-S")));
}

TEST_CASE("deriv/rules/open-shell", "[deriv][rules]") {
  using sequant::deriv::test::check_rule;
  const auto diff = [](const ExprPtr& E, const ExprPtr& x) {
    return sequant::deriv::diff(E, {x});
  };
  auto T = [](std::wstring_view s) { return deserialize(s); };
  // spin blocks of one label declare different symmetries; a block whose
  // spins are not the variable's is a constant (the judge alone decides, and
  // the result is not zero)
  const auto same = T(L"t{a↑3,a↑4;i↑3,i↑4}:A-N-S"),
             opposite = T(L"t{a↑3,a↓4;i↑3,i↓4}:N-N-S");
  const auto R_nonzero = [&](const ExprPtr& f, const ExprPtr& x) {
    check_rule(f, x, nullptr);
    REQUIRE_FALSE(diff(f, x)->is_zero());
  };
  for (const auto& x : {same, opposite})
    R_nonzero(T(L"1/4 g{i↑1,i↑2;a↑1,a↑2}:A-N-S t{a↑1,a↑2;i↑1,i↑2}:A-N-S + "
                L"g{i↑1,i↓1;a↑1,a↓1}:N-N-S t{a↑1,a↓1;i↑1,i↓1}:N-N-S"),
              x);
  // the spin-orbital CCSD energy, open-shell spin-traced
  auto context = get_default_context();
  context.set(CanonicalizeOptions{.method = CanonicalizationMethod::Complete});
  const auto scope = set_scoped_default_context(context);
  const auto energy = mbpt::open_shell_CC_spintrace(deserialize(
      L"f{i_1;a_1} t{a_1;i_1} + 1/4 g{i_1,i_2;a_1,a_2} t{a_1,a_2;i_1,i_2} + "
      L"1/2 g{i_1,i_2;a_1,a_2} t{a_1;i_1} t{a_2;i_2}",
      {.def_perm_symm = Symmetry::Antisymm}));
  REQUIRE(energy.size() == 1);
  for (const auto& x : {same, opposite}) R_nonzero(energy[0], x);
}

TEST_CASE("deriv/rules/at-zero", "[deriv][rules]") {
  using sequant::deriv::at_zero;
  const auto diff = [](const ExprPtr& E, const ExprPtr& x) {
    return sequant::deriv::diff(E, {x});
  };
  auto T = [](std::wstring_view s) { return deserialize(s); };
  const auto x = ex<Variable>(L"x"), eps = ex<Variable>(L"ε");
  const auto t1 = T(L"t{a_2;i_2}:N-N-N");
  REQUIRE_THAT(at_zero(T(L"2 x + x y + y"), {x}), EquivalentTo(L"y"));
  REQUIRE_THAT(at_zero(T(L"x^(2) y + y"), {x}), EquivalentTo(L"y"));
  REQUIRE_THAT(at_zero(adjoint(x) * T(L"y") + T(L"y"), {x}),
               EquivalentTo(L"y"));
  REQUIRE(at_zero(T(L"x y"), {x})->is_zero());
  REQUIRE_THAT(at_zero(T(L"x y + z"), {x, ex<Variable>(L"z")}),
               EquivalentTo(ex<Constant>(0)));
  REQUIRE(at_zero(T(L"x y + z"), {x, ex<Variable>(L"z")})->is_zero());
  REQUIRE(at_zero(ex<Constant>(0), {x})->is_zero());
  // E without the variable: E itself, not a Sum of one term
  {
    const auto E = T(L"f{i_1;a_1}:N-N-N u{a_1;i_1}:N-N-N");
    const auto r = at_zero(E, {x, t1});
    REQUIRE_THAT(r, EquivalentTo(E));
    REQUIRE(r->is<Product>());
  }
  // a tensor: its block and its adjoint spelling vanish, another block stays
  REQUIRE_THAT(
      at_zero(T(L"f{i_1;a_1}:N-N-N t{a_1;i_1}:N-N-N + f{i_1;a_1}:N-N-N "
                L"u{a_1;i_1}:N-N-N"),
              {t1}),
      EquivalentTo(L"f{i_1;a_1}:N-N-N u{a_1;i_1}:N-N-N"));
  REQUIRE_THAT(at_zero(T(L"w{a_1;i_1}:N-N-N t⁺{i_1;a_1}:N-N-N + "
                         L"f{i_1;a_1}:N-N-N u{a_1;i_1}:N-N-N"),
                       {t1}),
               EquivalentTo(L"f{i_1;a_1}:N-N-N u{a_1;i_1}:N-N-N"));
  REQUIRE_THAT(at_zero(T(L"w{a_1;i_1}:N-N-N t{i_1;a_1}:N-N-N"), {t1}),
               EquivalentTo(L"w{a_1;i_1}:N-N-N t{i_1;a_1}:N-N-N"));
  // a leaf that lies inside the variable vanishes; one that also holds
  // elements outside it does not
  REQUIRE(
      at_zero(T(L"W{i_1;a_1}:N-N-N k{a_1;i_1}:N-N-N"), {T(L"k{p_2;i_2}:N-N-N")})
          ->is_zero());
  REQUIRE_THROWS_WITH(
      at_zero(T(L"W{i_1;p_1}:N-N-N k{p_1;i_1}:N-N-N"),
              {T(L"k{a_2;i_2}:N-N-N")}),
      Catch::Matchers::ContainsSubstring(
          "deriv::at_zero: a leaf holds elements of a variable and lies "
          "inside none of the variables"));
  // the same for each spelling of the variable: a leaf inside one vanishes;
  // one that meets some and lies inside none throws
  REQUIRE(
      at_zero(T(L"G{i_2;p_1}:N-N-N h{p_1;i_2}:N-S-N"), {T(L"h{p_3;i_3}:N-S-N")})
          ->is_zero());
  REQUIRE_THAT(
      at_zero(T(L"W{i_3,a_3;a_1,i_1}:N-N-N g{a_1,i_1;i_3,a_3}:A-N-S + y"),
              {T(L"g{i_5,a_5;i_6,a_6}:A-N-S")}),
      EquivalentTo(L"y"));
  // several variables: a leaf inside one of them vanishes; one that lies
  // inside their union only throws
  {
    const container::svector<ExprPtr> xs{T(L"k{m_2;i_2}:N-N-N"),
                                         T(L"k{e_3;i_3}:N-N-N")};
    REQUIRE_THAT(at_zero(T(L"W{i_1;a_1}:N-N-N k{a_1;i_1}:N-N-N + y"), xs),
                 EquivalentTo(L"y"));
    REQUIRE_THROWS_WITH(
        at_zero(T(L"W{i_1;p_1}:N-N-N k{p_1;i_1}:N-N-N"), xs),
        Catch::Matchers::ContainsSubstring(
            "deriv::at_zero: a leaf holds elements of a variable and lies "
            "inside none of the variables"));
  }
  // bra-ket symmetric with unequal bra and ket ranks: a leaf written either
  // way vanishes
  for (const auto f : {L"W{a_1,a_2;i_1}:N-N-N g{i_1;a_1,a_2}:N-S-N + y",
                       L"V{i_1;a_1,a_2}:N-N-N g{a_1,a_2;i_1}:N-S-N + y"})
    REQUIRE_THAT(at_zero(T(f), {T(L"g{i_5;a_5,a_6}:N-S-N")}),
                 EquivalentTo(L"y"));
  for (const auto& [f, v] :
       {std::pair{L"G{p_2;p_1}:N-N-N h{p_1;p_2}:N-S-N", L"h{p_3;i_3}:N-S-N"},
        std::pair{L"W{p_2,a_2;p_1,a_1}:N-N-N g{p_1,a_1;p_2,a_2}:A-N-S",
                  L"g{i_5,a_5;i_6,a_6}:A-N-S"},
        std::pair{L"G{a_1;p_1}:N-N-N h{p_1;a_1}:N-S-N", L"h{i_3;a_3}:N-S-N"}})
    REQUIRE_THROWS_WITH(
        at_zero(T(f), {T(v)}),
        Catch::Matchers::ContainsSubstring(
            "deriv::at_zero: a leaf holds elements of a variable and lies "
            "inside none of the variables"));
  // a leaf whose protoindices do not follow its slots as the variable's do,
  // as in diff: only its elements with i_2 = i_3 are the variable's
  {
    const auto E =
        T(L"W{i_3;a_2<i_2>}:N-N-N V{i_2;}:N-N-N t{a_2<i_2>;i_3}:N-N-N");
    const auto x = T(L"t{a_1<i_1>;i_1}:N-N-N");
    const std::string text =
        "a leaf's protoindices do not follow its slots as the variable's do";
    REQUIRE_THROWS_WITH(
        diff(E, x), Catch::Matchers::ContainsSubstring("deriv::diff: " + text));
    REQUIRE_THROWS_WITH(at_zero(E, {x}), Catch::Matchers::ContainsSubstring(
                                             "deriv::at_zero: " + text));
  }
  // a leaf of the variable's block that declares other symmetries, as in
  // diff
  REQUIRE_THROWS_WITH(
      at_zero(T(L"W{i_1,i_2;a_1,a_2}:N-N-N t{a_1,a_2;i_1,i_2}:N-N-S"),
              {T(L"t{a_3,a_4;i_3,i_4}:A-N-S")}),
      Catch::Matchers::ContainsSubstring(
          "deriv::at_zero: a leaf declares other symmetries than the "
          "variable"));
  // the response pattern: differentiate in the strength, then set it to zero
  REQUIRE_THAT(at_zero(diff(T(L"ε^(2) y + 3 ε z + w"), eps), {eps}),
               EquivalentTo(L"3 z"));
  // the checks on the variable alone, as in diff
  for (const auto& [v, text] :
       {std::pair{T(L"Â{i_1;a_1}:A-N-S"),
                  "the variable is not a Variable or a tensor leaf"},
        std::pair{make_kronecker(Index(L"a_2"), Index(L"a_3")),
                  "the variable is not a Variable or a tensor leaf"},
        std::pair{T(L"h{p_1;p_1}:N-N-N"), "the variable repeats an index"},
        std::pair{T(L"x{a_1<i_1,i_2>;i_1}:N-N-N"),
                  "an index of the variable mixes protoindices that are its "
                  "slots with others"},
        std::pair{T(L"t{a_1,a_2;i_1,i_2}:S-N-S"),
                  "a variable with Symmetry::Symm is not supported"}})
    REQUIRE_THROWS_WITH(at_zero(T(L"f{i_1;a_1}:N-N-N t{a_1;i_1}:N-N-N"), {v}),
                        Catch::Matchers::ContainsSubstring(
                            std::string("deriv::at_zero: ") + text));
  // a variable under a non-positive power
  REQUIRE_THROWS_WITH(
      at_zero(T(L"x^(-1) y"), {x}),
      Catch::Matchers::ContainsSubstring(
          "deriv::at_zero: a variable stands under a non-positive power"));
  REQUIRE_THROWS_WITH(
      at_zero(ex<TensorPower>(T(L"t{a_1;i_1}:N-N-N"), rational(-1)) *
                  T(L"r{i_1;a_1}:N-N-N"),
              {t1}),
      Catch::Matchers::ContainsSubstring(
          "a variable stands under a non-positive power"));
}
