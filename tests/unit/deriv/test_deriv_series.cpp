#include <SeQuant/core/context.hpp>
#include <SeQuant/core/deriv/deriv.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/domain/mbpt/op.hpp>
#include <SeQuant/domain/mbpt/series.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "contexts.hpp"
#include "evaluator.hpp"
#include "judge.hpp"

#include <cstddef>
#include <string>
#include <string_view>

using namespace sequant;

namespace {

using deriv::test::require_same;
using deriv::test::with_batching_spaces;

/// the tensor leaves of @p e labelled @p label
container::svector<ExprPtr> leaves_of(const ExprPtr& e,
                                      std::wstring_view label) {
  container::svector<ExprPtr> result;
  e->visit(
      [&](const ExprPtr& leaf) {
        if (leaf->is<Tensor>() && leaf->as<Tensor>().label() == label)
          result.push_back(leaf);
      },
      /*atoms_only=*/true);
  return result;
}

}  // namespace

TEST_CASE("deriv/series", "[deriv][series]") {
  using Catch::Matchers::StartsWith;
  using deriv::at_zero;
  using deriv::diff;
  const auto T = [](std::wstring_view s) { return deserialize(s); };
  const auto eps = ex<Variable>(L"ε");

  SECTION("scalar strength") {
    const auto E = T(L"f{i_1;a_1}:N-C-S κₒ{a_1;i_1}:N-N-N");
    const auto S = mbpt::series(E, L"f", eps, 2);
    require_same(S, T(L"(f{i_1;a_1}:N-C-S + ε f¹{i_1;a_1}:N-C-S + "
                      L"1/2 ε^(2) f²{i_1;a_1}:N-C-S) κₒ{a_1;i_1}:N-N-N"));
    const auto& f = E->at(0)->as<Tensor>();
    for (const auto* label : {L"f¹", L"f²"}) {
      const auto leaves = leaves_of(S, label);
      REQUIRE(leaves.size() == 1);
      const auto& c = leaves.front()->as<Tensor>();
      REQUIRE(c.symmetry() == f.symmetry());
      REQUIRE(c.braket_symmetry() == f.braket_symmetry());
      REQUIRE(c.column_symmetry() == f.column_symmetry());
      REQUIRE(c.bra() == f.bra());
      REQUIRE(c.ket() == f.ket());
      REQUIRE(c.aux().empty());
    }
  }

  SECTION("Taylor property") {
    const auto E =
        T(L"1/2 g{i_1,i_2;a_1,a_2}:A-C-S t{a_1;i_1}:N-N-N t{a_2;i_2}:N-N-N");
    const auto S = mbpt::series(E, L"t", eps, 2);
    require_same(at_zero(S, {eps}), E);
    require_same(at_zero(diff(S, {eps}), {eps}),
                 T(L"g{i_1,i_2;a_1,a_2}:A-C-S t¹{a_1;i_1}:N-N-N "
                   L"t{a_2;i_2}:N-N-N"));
    require_same(at_zero(diff(S, {eps, eps}), {eps}),
                 T(L"g{i_1,i_2;a_1,a_2}:A-C-S t¹{a_1;i_1}:N-N-N "
                   L"t¹{a_2;i_2}:N-N-N + g{i_1,i_2;a_1,a_2}:A-C-S "
                   L"t²{a_1;i_1}:N-N-N t{a_2;i_2}:N-N-N"));
  }

  SECTION("third order") {
    const auto E = T(L"f{i_1;a_1}:N-C-S κₒ{a_1;i_1}:N-N-N");
    require_same(
        at_zero(diff(mbpt::series(E, L"f", eps, 3), {eps, eps, eps}), {eps}),
        T(L"f³{i_1;a_1}:N-C-S κₒ{a_1;i_1}:N-N-N"));
    const auto scope = with_batching_spaces();
    const auto strength = T(L"ε{;;z_1}:N-N-N");
    const Index z3(L"z_3"), z4(L"z_4"), z5(L"z_5");
    const auto S = mbpt::series(T(L"h{p_1;p_2}:N-C-S D{p_2;p_1}:N-N-N"), L"h",
                                strength, 3);
    std::wstring expected;
    for (const auto* z : {L"z_3,z_4,z_5", L"z_3,z_5,z_4", L"z_4,z_3,z_5",
                          L"z_4,z_5,z_3", L"z_5,z_3,z_4", L"z_5,z_4,z_3"})
      expected += (expected.empty() ? L"" : L" + ") +
                  std::wstring(L"1/6 h³{p_1;p_2;") + z +
                  L"}:N-C-S D{p_2;p_1}:N-N-N";
    require_same(at_zero(diff(S, {T(L"ε{;;z_3}:N-N-N"), T(L"ε{;;z_4}:N-N-N"),
                                  T(L"ε{;;z_5}:N-N-N")}),
                         {strength}),
                 T(expected), {z3, z4, z5});
  }

  SECTION("label match") {
    const auto E =
        T(L"t{a_1;i_1}:N-N-N f{i_1;a_1}:N-C-S + "
          L"1/4 t{a_1,a_2;i_1,i_2}:A-N-S g{i_1,i_2;a_1,a_2}:A-C-S + "
          L"t²{a_1;i_1}:N-N-N w{i_1;a_1}:N-N-N");
    const auto S = mbpt::series(E, L"t", eps, 1);
    REQUIRE(leaves_of(S, L"t¹").size() == 2);
    REQUIRE(leaves_of(S, L"t²¹").empty());
    REQUIRE(leaves_of(S, L"f¹").empty());
    require_same(at_zero(S, {eps}), E);
    require_same(at_zero(diff(S, {eps}), {eps}),
                 T(L"t¹{a_1;i_1}:N-N-N f{i_1;a_1}:N-C-S + "
                   L"1/4 t¹{a_1,a_2;i_1,i_2}:A-N-S g{i_1,i_2;a_1,a_2}:A-C-S"));
  }

  SECTION("every block") {
    const auto E =
        T(L"f{i_1;i_1}:N-C-S + f{i_1;a_1}:N-C-S κₒ{a_1;i_1}:N-N-N + "
          L"f{a_1;a_2}:N-C-S X{a_2;a_1}:N-N-N");
    require_same(at_zero(diff(mbpt::series(E, L"f", eps, 1), {eps}), {eps}),
                 T(L"f¹{i_1;i_1}:N-C-S + f¹{i_1;a_1}:N-C-S κₒ{a_1;i_1}:N-N-N "
                   L"+ f¹{a_1;a_2}:N-C-S X{a_2;a_1}:N-N-N"));
  }

  SECTION("adjoint") {
    const auto kappa_adjoint = adjoint(T(L"κₒ{a_1;i_1}:N-N-N"));
    const auto E = kappa_adjoint * T(L"f{a_1;i_1}:N-C-S");
    const auto S = mbpt::series(E, L"κₒ", eps, 1);
    const auto zeroth = leaves_of(S, kappa_adjoint->as<Tensor>().label());
    REQUIRE(zeroth.size() == 1);
    REQUIRE(*zeroth.front() == *kappa_adjoint);
    const auto coefficient = adjoint(T(L"κₒ¹{a_1;i_1}:N-N-N"));
    const auto leaves =
        leaves_of(S, coefficient->as<Tensor>().label());  // κₒ¹⁺
    REQUIRE(leaves.size() == 1);
    REQUIRE(*leaves.front() == *coefficient);
    require_same(at_zero(diff(S, {eps}), {eps}),
                 coefficient * T(L"f{a_1;i_1}:N-C-S"));
  }

  SECTION("strength tensor") {
    const auto scope = with_batching_spaces();
    const auto strength = T(L"ε{;;z_1}:N-N-N");
    const auto e3 = T(L"ε{;;z_3}:N-N-N"), e4 = T(L"ε{;;z_4}:N-N-N");
    const Index z3(L"z_3"), z4(L"z_4");
    const auto S = mbpt::series(T(L"h{p_1;p_2}:N-C-S D{p_2;p_1}:N-N-N"), L"h",
                                strength, 2);
    require_same(at_zero(diff(S, {e3}), {strength}),
                 T(L"h¹{p_1;p_2;z_3}:N-C-S D{p_2;p_1}:N-N-N"), {z3});
    require_same(at_zero(diff(S, {e3, e4}), {strength}),
                 T(L"1/2 h²{p_1;p_2;z_3,z_4}:N-C-S D{p_2;p_1}:N-N-N + "
                   L"1/2 h²{p_1;p_2;z_4,z_3}:N-C-S D{p_2;p_1}:N-N-N"),
                 {z3, z4});
  }

  SECTION("two leaves, strength tensor") {
    const auto scope = with_batching_spaces();
    const auto strength = T(L"ε{;;z_1}:N-N-N");
    const auto e3 = T(L"ε{;;z_3}:N-N-N"), e4 = T(L"ε{;;z_4}:N-N-N");
    const Index z3(L"z_3"), z4(L"z_4");
    const auto S = mbpt::series(T(L"h{p_1;p_2}:N-C-S h{p_2;p_1}:N-C-S"), L"h",
                                strength, 1);
    require_same(at_zero(diff(S, {e3, e4}), {strength}),
                 T(L"h¹{p_1;p_2;z_3}:N-C-S h¹{p_2;p_1;z_4}:N-C-S + "
                   L"h¹{p_1;p_2;z_4}:N-C-S h¹{p_2;p_1;z_3}:N-C-S"),
                 {z3, z4});
  }

  SECTION("two series in one term") {
    const auto scope = with_batching_spaces();
    const auto strength = T(L"ε{;;z_1}:N-N-N");
    const auto e3 = T(L"ε{;;z_3}:N-N-N"), e4 = T(L"ε{;;z_4}:N-N-N");
    const Index z3(L"z_3"), z4(L"z_4");
    const auto h12 = T(L"h{p_1;p_2}:N-C-S"), h21 = T(L"h{p_2;p_1}:N-C-S");
    const auto S = mbpt::series(h12, L"h", strength, 1) *
                   mbpt::series(h21, L"h", strength, 1);
    require_same(S, mbpt::series(h12 * h21, L"h", strength, 1));
    require_same(at_zero(diff(S, {e3, e4}), {strength}),
                 T(L"h¹{p_1;p_2;z_3}:N-C-S h¹{p_2;p_1;z_4}:N-C-S + "
                   L"h¹{p_1;p_2;z_4}:N-C-S h¹{p_2;p_1;z_3}:N-C-S"),
                 {z3, z4});
  }

  SECTION("a summed index named like the variable's") {
    const auto scope = with_batching_spaces();
    const auto strength = T(L"ε{;;z_1}:N-N-N");
    const Index z1(L"z_1");
    Index::reset_tmp_index();
    const auto S = mbpt::series(
        T(L"h{p_1;p_2}:N-C-S h{p_3;p_3}:N-C-S B{p_2;p_1;z_1}:N-N-N "
          L"C{;;z_1}:N-N-N"),
        L"h", strength, 1);
    require_same(at_zero(diff(S, {strength}), {strength}),
                 T(L"h¹{p_1;p_2;z_1}:N-C-S h{p_3;p_3}:N-C-S "
                   L"B{p_2;p_1;z_2}:N-N-N C{;;z_2}:N-N-N + "
                   L"h{p_1;p_2}:N-C-S h¹{p_3;p_3;z_1}:N-C-S "
                   L"B{p_2;p_1;z_2}:N-N-N C{;;z_2}:N-N-N"),
                 {z1});
  }

  SECTION("a name in use") {
    const auto scope = with_batching_spaces();
    const auto strength = T(L"ε{;;z_1}:N-N-N");
    const auto e3 = T(L"ε{;;z_3}:N-N-N");
    const Index z3(L"z_3");
    // the first temporary index of z, then the counter back at it
    Index::reset_tmp_index();
    const auto z_tmp = Index::make_tmp_index(z3.basis());
    Index::reset_tmp_index();
    const auto B = ex<Tensor>(L"B", bra{Index(L"p_2")}, ket{Index(L"p_1")},
                              aux{z_tmp}, Symmetry::Nonsymm,
                              BraKetSymmetry::Nonsymm, ColumnSymmetry::Symm);
    const auto S = mbpt::series(T(L"h{p_1;p_2}:N-C-S") * B, L"h", strength, 1);
    require_same(at_zero(diff(S, {e3}), {strength}),
                 T(L"h¹{p_1;p_2;z_3}:N-C-S") * B, {z3, z_tmp});
  }

  SECTION("a leaf with an auxiliary index") {
    const auto scope = with_batching_spaces();
    const auto strength = T(L"ε{;;z_1}:N-N-N");
    const auto e3 = T(L"ε{;;z_3}:N-N-N");
    const Index z3(L"z_3"), z9(L"z_9");
    const auto S = mbpt::series(T(L"B{p_1;p_2;z_9}:N-N-N D{p_2;p_1}:N-N-N"),
                                L"B", strength, 1);
    const auto leaves = leaves_of(S, L"B¹");
    REQUIRE(leaves.size() == 1);
    const auto& aux = leaves.front()->as<Tensor>().aux();
    REQUIRE(aux.size() == 2);
    REQUIRE(aux[0] == z9);
    REQUIRE(aux[1].space() == z9.space());
    require_same(at_zero(diff(S, {e3}), {strength}),
                 T(L"B¹{p_1;p_2;z_9,z_3}:N-N-N D{p_2;p_1}:N-N-N"), {z3, z9});
  }

  SECTION("rejections") {
    const auto scope = with_batching_spaces();
    const auto E = T(L"h{p_1;p_2}:N-C-S D{p_2;p_1}:N-N-N");
    for (const auto& strength :
         {T(L"ε{p_1;p_2}:N-N-N"), T(L"ε{;;z_1,z_2}:N-N-N"), ex<Constant>(2)}) {
      INFO(deriv::test::text_of(strength));
      REQUIRE_THROWS_WITH(mbpt::series(E, L"h", strength, 1),
                          StartsWith("mbpt::series: "));
    }
    REQUIRE_THROWS_WITH(mbpt::series(E, L"h", eps, 10),
                        StartsWith("mbpt::series: "));
    REQUIRE_NOTHROW(mbpt::series(E, L"h", eps, 9));
    const auto rejects = [&](const ExprPtr& e, std::wstring_view label,
                             const ExprPtr& strength, std::size_t order = 1) {
      REQUIRE_THROWS_WITH(mbpt::series(e, label, strength, order),
                          StartsWith("mbpt::series: "));
    };
    rejects(nullptr, L"h", eps);
    rejects(E, L"h", nullptr);
    rejects(E, L"", eps);
    rejects(E, L"h⁺", eps);
    rejects(E, L"h", adjoint(eps));
    rejects(E, L"h", T(L"ε{;;z_1<p_1>}:N-N-N"));
    // a strength labelled like the leaves or a coefficient
    rejects(E, L"h", T(L"h{;;z_1}:N-N-N"));
    rejects(E * T(L"ε{;;z_2}:N-N-N"), L"ε", T(L"ε{;;z_1}:N-N-N"));
    rejects(E, L"h", T(L"h¹{;;z_1}:N-N-N"));
    rejects(E, L"h", ex<Variable>(L"h"));
    rejects(E, L"h", ex<Variable>(L"h²"), 2);
    // input series() would return unchanged: operators, a Variable of the label
    REQUIRE_THROWS_WITH(
        mbpt::series(mbpt::op::t(1), L"t", eps, 1),
        "mbpt::series: the expression holds operators; expand its tensor form");
    rejects(ex<Variable>(L"h") * E, L"h", eps);
    rejects(ex<Power>(ex<Variable>(L"h"), 2) * E, L"h", eps);
    REQUIRE_THROWS_WITH(
        mbpt::series(
            ex<TensorPower>(T(L"t{a_1;i_1}:N-N-N"), 2) * T(L"f{i_1;a_1}:N-C-S"),
            L"t", eps, 1),
        "mbpt::series: a tensor to expand stands under a power");
    // a coefficient label already in use: two meanings of h¹ (h²⁺)
    rejects(E + T(L"h¹{p_1;p_2}:N-C-S W{p_2;p_1}:N-N-N"), L"h", eps);
    const auto h2_adjoint = adjoint(T(L"h²{p_1;p_2}:N-N-N"));
    rejects(E + h2_adjoint * T(L"W{p_1;p_2}:N-N-N"), L"h", eps, 2);
    REQUIRE_NOTHROW(
        mbpt::series(E + h2_adjoint * T(L"W{p_1;p_2}:N-N-N"), L"h", eps, 1));
  }
}
