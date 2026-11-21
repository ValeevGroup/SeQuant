#include <SeQuant/core/deriv/deriv.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/rational.hpp>

#include <catch2/catch_test_macros.hpp>
#include "catch2_sequant.hpp"
#include "judge.hpp"

#include <string_view>

using namespace sequant;

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
