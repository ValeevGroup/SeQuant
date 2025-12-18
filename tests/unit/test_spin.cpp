//
// Created by Nakul Teke on 12/20/19.
//

#include <catch2/catch_test_macros.hpp>
#include "catch2_sequant.hpp"

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/hash.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/latex.hpp>
#include <SeQuant/core/parse.hpp>
#include <SeQuant/core/rational.hpp>
#include <SeQuant/core/space.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>

#include <cassert>
#include <cstddef>
#include <initializer_list>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

#include <range/v3/all.hpp>

TEST_CASE("spin", "[spin]") {
  using namespace sequant;
  using namespace sequant::mbpt;

  // P.S. ref outputs produced with complete canonicalization
  auto ctx = get_default_context();
  ctx.set(CanonicalizeOptions{.method = CanonicalizationMethod::Complete});
  auto _ = set_scoped_default_context(ctx);

  TensorCanonicalizer::register_instance(
      std::make_shared<DefaultTensorCanonicalizer>());

  auto reset_idx_tags = [](ExprPtr& expr) {
    if (expr->is<AbstractTensor>())
      ranges::for_each(expr->as<AbstractTensor>()._slots(),
                       [](const Index& idx) { idx.reset_tag(); });
  };

  SECTION("protoindices supported") {
    auto isr = get_default_context().index_space_registry();
    Index i1(L"i_1");
    Index a1(L"a_1", {i1});

    const auto expr =
        ex<Tensor>(L"t", bra{i1}, ket{a1}) * ex<Tensor>(L"F", bra{a1}, ket{i1});
    REQUIRE_NOTHROW(spintrace(expr));
    {  // assume spin-free spaces
      auto expr_st = spintrace(expr);
      simplify(expr_st);
      REQUIRE_THAT(expr_st, EquivalentTo("2 t{i1;a1<i1>} F{a1<i1>;i1}"));
    }
    {  // assume spin-dependent spaces
      auto expr_st = spintrace(expr, {}, /* assume_spin_free_spaces */ false);
      simplify(expr_st);
      REQUIRE_THAT(expr_st, EquivalentTo("t{i↓1;a↓1<i↓1>} F{a↓1<i↓1>;i↓1} + "
                                         "t{i↑1;a↑1<i↑1>} F{a↑1<i↑1>;i↑1}"));
    }
  }

  SECTION("ASCII label") {
    IndexSpace pup(L"p↑", 0b011, mbpt::Spin::alpha);
    IndexSpace pdown(L"p↓", 0b011, mbpt::Spin::beta);
    IndexSpace alphaup(L"α↑", 0b110, mbpt::Spin::alpha);

    auto p1 = Index(L"p↑_1", pup);
    auto p2 = Index(L"p↓_2", pdown);
    auto p3 = Index(L"p↑_3", pup);
    auto p4 = Index(L"p↓_4", pdown);
    auto alpha1 = Index(L"α↑_1", alphaup);

    SEQUANT_PRAGMA_CLANG(diagnostic push)
    SEQUANT_PRAGMA_CLANG(diagnostic ignored "-Wdeprecated-declarations")
    SEQUANT_PRAGMA_GCC(diagnostic push)
    SEQUANT_PRAGMA_GCC(diagnostic ignored "-Wdeprecated-declarations")
    REQUIRE(p1.ascii_label() == "pa_1");
    REQUIRE(p2.ascii_label() == "pb_2");
    REQUIRE(p3.ascii_label() == "pa_3");
    REQUIRE(p4.ascii_label() == "pb_4");
    REQUIRE(alpha1.ascii_label() == "alphaa_1");
    SEQUANT_PRAGMA_GCC(diagnostic pop)
    SEQUANT_PRAGMA_CLANG(diagnostic pop)
  }

  SECTION("Index: add/remove spin") {
    auto i = Index(L"i", {L"i", 0b01, mbpt::Spin::any});
    auto i1 = Index(L"i_1", {L"i", 0b01, mbpt::Spin::any});
    auto p = Index(L"p", {L"p", 0b11, mbpt::Spin::any});
    auto p1 = Index(L"p_1", {L"p", 0b11, mbpt::Spin::any});
    auto p1_a = Index(L"p↑_1", {L"p↑", 0b11, mbpt::Spin::alpha});
    auto p2 = Index(L"p_2", {L"p", 0b11, mbpt::Spin::any});
    auto p2_b = Index(L"p↓_2", {L"p↓", 0b11, mbpt::Spin::beta});

    auto p_i = Index({L"p", 0b11, mbpt::Spin::any}, {i});
    auto p1_i = Index({L"p", 0b11, mbpt::Spin::any}, 1, {i});
    auto p_i1 = Index({L"p", 0b11, mbpt::Spin::any}, {i1});
    auto p1_i1 = Index({L"p", 0b11, mbpt::Spin::any}, 1, {i1});

    // make_spinalpha
    {
      // plain
      REQUIRE_NOTHROW(make_spinalpha(p));
      REQUIRE(make_spinalpha(p).label() == L"p↑");
      IndexSpace p_a(L"p↑", 0b11, mbpt::Spin::alpha);
      REQUIRE(make_spinalpha(p).space() == p_a);
      REQUIRE_NOTHROW(make_spinalpha(p1));
      REQUIRE(make_spinalpha(p1) == p1_a);
      // idempotent
      REQUIRE_NOTHROW(make_spinalpha(p1_a));
      REQUIRE(make_spinalpha(p1_a) == p1_a);
      // can flip spin
      REQUIRE_NOTHROW(make_spinalpha(p2_b));
      REQUIRE(make_spinalpha(p2_b) == make_spinalpha(p2));

      // proto
      REQUIRE_NOTHROW(make_spinalpha(p_i));
      REQUIRE(make_spinalpha(p_i).label() == L"p↑");
      REQUIRE(make_spinalpha(p_i).full_label() == L"p↑<i↑>");
      REQUIRE(make_spinalpha(p_i).to_latex() == L"{p↑^{{i↑}}}");
      REQUIRE_NOTHROW(make_spinalpha(p1_i));
      REQUIRE(make_spinalpha(p1_i).label() == L"p↑_1");
      REQUIRE(make_spinalpha(p1_i).full_label() == L"p↑_1<i↑>");
      REQUIRE(make_spinalpha(p1_i).to_latex() == L"{p↑_1^{{i↑}}}");
      REQUIRE_NOTHROW(make_spinalpha(p_i1));
      REQUIRE(make_spinalpha(p_i1).label() == L"p↑");
      REQUIRE(make_spinalpha(p_i1).full_label() == L"p↑<i↑_1>");
      REQUIRE(make_spinalpha(p_i1).to_latex() == L"{p↑^{{i↑_1}}}");
      REQUIRE_NOTHROW(make_spinalpha(p1_i1));
      REQUIRE(make_spinalpha(p1_i1).label() == L"p↑_1");
      REQUIRE(make_spinalpha(p1_i1).full_label() == L"p↑_1<i↑_1>");
      REQUIRE(make_spinalpha(p1_i1).to_latex() == L"{p↑_1^{{i↑_1}}}");
    }

    // make_spinbeta
    {
      REQUIRE_NOTHROW(make_spinbeta(p1));
      REQUIRE(make_spinbeta(p2) == p2_b);
      // idempotent
      REQUIRE_NOTHROW(make_spinbeta(p2_b));
      REQUIRE(make_spinbeta(p2_b) == p2_b);
      // can flip spin
      REQUIRE_NOTHROW(make_spinbeta(p1_a));
      REQUIRE(make_spinbeta(p1_a) == make_spinbeta(p1));

      // proto
      // N.B. only test spin flip
      REQUIRE_NOTHROW(make_spinbeta(make_spinalpha(p1_i1)));
      REQUIRE(make_spinbeta(make_spinalpha(p1_i1)) == make_spinbeta((p1_i1)));
      REQUIRE(make_spinbeta(make_spinalpha(p1_i1)) ==
              make_spinbeta(make_spinbeta((p1_i1))));
    }

    // make spinnull
    {
      // plain
      REQUIRE_NOTHROW(make_spinfree(p1_a));
      REQUIRE(make_spinfree(p1_a) == p1);
      REQUIRE_NOTHROW(make_spinfree(p2_b));
      REQUIRE(make_spinfree(p2_b) == p2);
      REQUIRE_NOTHROW(make_spinfree(p1));
      REQUIRE(make_spinfree(p1) == p1);
      // idempotent
      REQUIRE_NOTHROW(make_spinfree(p2));
      REQUIRE(make_spinfree(p2) == p2);

      // proto
      REQUIRE_NOTHROW(make_spinfree(make_spinalpha(p1_i1)));
      REQUIRE(make_spinfree(make_spinalpha(p1_i1)) == p1_i1);
      REQUIRE(make_spinfree(make_spinalpha(p1_i1)) == make_spinfree(p1_i1));
    }
  }

  SECTION("Tensor: can_expand, ms_conserving_columns, remove_spin") {
    auto p1 = Index(L"p↑_1");
    auto p2 = Index(L"p↓_2");
    auto p3 = Index(L"p↑_3");
    auto p4 = Index(L"p↓_4");

    auto input = ex<Tensor>(L"t", bra{p1, p2}, ket{p3, p4});
    REQUIRE(can_expand(input->as<Tensor>()) == true);
    REQUIRE(ms_conserving_columns(input->as<Tensor>()) == true);

    auto spin_swap_tensor = swap_spin(input->as<Tensor>());
    REQUIRE_THAT(spin_swap_tensor, EquivalentTo("t{p↓1,p↑2;p↓3,p↑4}"));

    auto result = remove_spin(input);
    for (auto& i : result->as<Tensor>().const_braket_indices())
      REQUIRE(i.space().base_key() == L"p");

    input = ex<Tensor>(L"t", bra{p1, p3}, ket{p2, p4});
    REQUIRE_THAT(swap_spin(input), EquivalentTo("t{p↓1,p↓3;p↑2,p↑4}"));
    REQUIRE(can_expand(input->as<Tensor>()) == false);
    REQUIRE(ms_conserving_columns(input->as<Tensor>()) == false);
  }

  SECTION("Tensor: expand_antisymm") {
    // 1-body
    auto input = ex<Tensor>(L"t", bra{L"a_1"}, ket{L"i_1"});
    auto result = expand_antisymm(input->as<Tensor>());
    REQUIRE(input->as<Tensor>() == result->as<Tensor>());
    REQUIRE(!result->is<Sum>());
    REQUIRE_THAT(result, EquivalentTo("t{a1;i1}"));

    // 1-body
    input = ex<Tensor>(L"t", bra{L"a_1"}, ket{L"i_1"}, Symmetry::Antisymm);
    result = expand_antisymm(input->as<Tensor>());
    REQUIRE(input->as<Tensor>().symmetry() == Symmetry::Antisymm);
    REQUIRE(result->as<Tensor>().symmetry() == Symmetry::Nonsymm);
    REQUIRE_THAT(result, EquivalentTo("t{a1;i1}"));

    // 2-body
    input = ex<Tensor>(L"g", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"},
                       Symmetry::Antisymm);
    result = expand_antisymm(input->as<Tensor>());
    REQUIRE_THAT(result, EquivalentTo("g{i1,i2;a1,a2} - g{i2,i1;a1,a2}"));

    // 3-body
    input = ex<Tensor>(L"t", bra{L"a_1", L"a_2", L"a_3"},
                       ket{L"i_1", L"i_2", L"i_3"}, Symmetry::Antisymm);
    result = expand_antisymm(input->as<Tensor>());
    REQUIRE_THAT(result,
                 EquivalentTo("t{a1,a2,a3;i1,i2,i3} - t{a1,a3,a2;i1,i2,i3} - "
                              "t{a2,a1,a3;i1,i2,i3} + t{a2,a3,a1;i1,i2,i3} + "
                              "t{a3,a1,a2;i1,i2,i3} - t{a3,a2,a1;i1,i2,i3}"));
  }

  SECTION("Constant") {
    auto exprPtr = ex<Constant>(rational{1, 4});
    auto result = spintrace(exprPtr);
    REQUIRE(result->is<Constant>());
    REQUIRE(result->is_atom());
    REQUIRE_THAT(result, EquivalentTo("1/4"));
    REQUIRE_THAT(swap_spin(exprPtr), EquivalentTo("1/4"));
  }
  SECTION("Variable") {
    auto exprPtr = ex<Variable>(L"Var");
    auto result = spintrace(exprPtr);
    REQUIRE(result->is<Variable>());
    REQUIRE(result->is_atom());
    REQUIRE(result->as<Variable>().label() == L"Var");
    REQUIRE(to_latex(swap_spin(exprPtr->clone())) == to_latex(exprPtr));
  }

  SECTION("Tensor") {
    {
      const auto expr = ex<Constant>(rational{1, 4}) *
                        ex<Tensor>(L"g", bra{L"p_1", L"p_2"},
                                   ket{L"p_3", L"p_4"}, Symmetry::Antisymm);
      auto result = spintrace(expr, {{L"p_1", L"p_3"}, {L"p_2", L"p_4"}});
      REQUIRE_THAT(result,
                   EquivalentTo("-1/2 g{p1,p2;p4,p3} + g{p1,p2;p3,p4}"));
    }
    {
      // Note the provided external index pairings which is different from the
      // way this tensor is written down Also, the prefactor of -1 is important
      // as this forces the code to take a path where it traces a product which
      // ends up as a single tensor (due to an additional factor of -1 induces
      // by the chosen pairing of external indices)
      ExprPtr expr = parse_expr(L"- g{p1,p2;p3,p4}:A");
      ExprPtr result = spintrace(expr, {{L"p_1", L"p_4"}, {L"p_2", L"p_3"}});

      REQUIRE_THAT(result,
                   EquivalentTo(L"4 g{p1,p2;p4,p3} - 2 g{p1,p2;p3,p4}"));
    }
    {  // non-closed-shell -> keep spin labels
      const auto expr = ex<Constant>(rational{1, 4}) *
                        ex<Tensor>(L"g", bra{L"p_1", L"p_2"},
                                   ket{L"p_3", L"p_4"}, Symmetry::Antisymm);
      auto result = spintrace(expr, {{L"p_1", L"p_3"}, {L"p_2", L"p_4"}},
                              /* spinfree = */ false);
      REQUIRE_THAT(
          result,
          EquivalentTo(
              "1/4 g{p↑_1,p↑_2;p↑_3,p↑_4}:N-C-S - 1/4 "
              "g{p↑_2,p↑_1;p↑_3,p↑_4}:N-C-S + 1/4 g{p↑_1,p↓_2;p↑_3,p↓_4}:N-C-S "
              "+ 1/4 g{p↑_2,p↓_1;p↑_4,p↓_3}:N-C-S + 1/4 "
              "g{p↓_1,p↓_2;p↓_3,p↓_4}:N-C-S - 1/4 "
              "g{p↓_2,p↓_1;p↓_3,p↓_4}:N-C-S"));
    }
  }

  SECTION("Product") {
    const auto expr = ex<Tensor>(L"f", bra{L"i_1"}, ket{L"a_1"}) *
                      ex<Tensor>(L"t", bra{L"a_1"}, ket{L"i_1"});
    auto result = spintrace(expr);
    canonicalize(result);
    REQUIRE_THAT(result, EquivalentTo("2 f{i1;a1} t{a1;i1}"));
  }

  SECTION("Scaled Product") {
    {
      // 1/2 * g * t1 * t1
      const auto expr = ex<Constant>(rational{1, 2}) *
                        ex<Tensor>(L"g", bra{L"i_1", L"i_2"},
                                   ket{L"a_1", L"a_2"}, Symmetry::Antisymm) *
                        ex<Tensor>(L"t", bra{L"a_1"}, ket{L"i_1"}) *
                        ex<Tensor>(L"t", bra{L"a_2"}, ket{L"i_2"});
      auto result = spintrace(expr);
      canonicalize(result);
      REQUIRE_THAT(result,
                   EquivalentTo("- g{i1,i2;a1,a2} t{a1;i2} t{a2;i1} "
                                "+ 2 g{i1,i2;a1,a2} t{a1;i1} t{a2;i2}"));
    }
  }

  SECTION("Scaled Product with variable") {
    ExprPtr expr = parse_expr(L"1/2 Var g{i1,i2;a1,a2}:A t{a1;i1} t{a2;i2}");
    auto result = spintrace(expr);
    canonicalize(result);
    REQUIRE_THAT(
        result,
        EquivalentTo(
            L"2 Var * g{i_1,i_2;a_1,a_2}:N * t{a_1;i_1}:N * t{a_2;i_2}:N"
            " - 1 Var * g{i_1,i_2;a_1,a_2}:N * t{a_1;i_2}:N * t{a_2;i_1}:N"));
  }

  SECTION("Tensor times variable") {
    ResultExpr expr = parse_result_expr(
        L"R2{a1,a2;i1,i2}:A = 1/4 A{i1,i2;a1,a2}:A INTkx{a1,a2;i1,i2}:A H");
    auto results = closed_shell_spintrace(expr);
    REQUIRE_THAT(
        results.at(0),
        EquivalentTo(L"R2{a_1,a_2;i_1,i_2}:N = -1 H * S{i_1,i_2;a_1,a_2}:N "
                     L"* INTkx{a_1,a_2;i_2,i_1}:N + 2 H * "
                     L"S{i_1,i_2;a_1,a_2}:N * INTkx{a_1,a_2;i_1,i_2}:N"));
  }

  SECTION("Sum") {
    // f * t1 + 1/2 * g * t1 * t1 + 1/4 * g * t2
    const auto ex1 = ex<Tensor>(L"f", bra{L"i_1"}, ket{L"a_1"}) *
                     ex<Tensor>(L"t", bra{L"a_1"}, ket{L"i_1"});
    const auto ex2 = ex<Constant>(rational{1, 2}) *
                     ex<Tensor>(L"g", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"},
                                Symmetry::Antisymm) *
                     ex<Tensor>(L"t", bra{L"a_1"}, ket{L"i_1"}) *
                     ex<Tensor>(L"t", bra{L"a_2"}, ket{L"i_2"});
    const auto ex3 = ex<Constant>(rational{1, 4}) *
                     ex<Tensor>(L"g", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"},
                                Symmetry::Antisymm) *
                     ex<Tensor>(L"t", bra{L"a_1", L"a_2"}, ket{L"i_1", L"i_2"},
                                Symmetry::Antisymm);

    auto expr = ex1 + ex2 + ex3;
    auto result = ex<Constant>(rational{1, 2}) * spintrace(expr);
    expand(result);
    rapid_simplify(result);
    canonicalize(result);
    REQUIRE(result->is<Sum>());
    REQUIRE(result->size() == 5);
    REQUIRE_THAT(result, EquivalentTo("-1/2 g{i1,i2;a1,a2} t{a1,a2;i2,i1} "
                                      "+ g{i1,i2;a1,a2} t{a1,a2;i1,i2} "
                                      "+ f{i1;a1} t{a1;i1} "
                                      "- 1/2 g{i1,i2;a1,a2} t{a1;i2} t{a2;i1} "
                                      "+ g{i1,i2;a1,a2} t{a1;i1} t{a2;i2}"));
  }  // Sum

  SECTION("Expand Antisymmetrizer"){// 0-body
                                    {auto input = ex<Constant>(1);
  auto result = expand_A_op(input);
  REQUIRE(result->size() == 0);
  REQUIRE(result->is_atom());

  input = ex<Constant>(1) *
          ex<Tensor>(L"A", bra{L"i_1"}, ket{L"a_1"}, Symmetry::Antisymm);
  result = expand_A_op(input);
  REQUIRE(result->size() == 0);
  REQUIRE(result->is_atom());
}

// 1-body
{
  auto input = ex<Tensor>(L"A", bra{L"i_1"}, ket{L"a_1"}, Symmetry::Antisymm) *
               ex<Tensor>(L"t", bra{L"a_1"}, ket{L"i_1"}, Symmetry::Antisymm);
  auto result = expand_A_op(input);
  REQUIRE(result->size() == 1);
  REQUIRE(!result->is<Sum>());
}

// 2-body
{
  auto input = ex<Constant>(rational{1, 4}) *
               ex<Tensor>(L"g", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"},
                          Symmetry::Antisymm);
  auto result = expand_A_op(input);
  REQUIRE_THAT(result, EquivalentTo("1/4 g{i1,i2;a1,a2}:A"));

  input = ex<Constant>(rational{1, 4}) *
          ex<Tensor>(L"A", bra{L"a_1", L"a_2"}, ket{L"i_1", L"i_2"},
                     Symmetry::Antisymm) *
          ex<Tensor>(L"g", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"},
                     Symmetry::Antisymm);
  result = expand_A_op(input);
  REQUIRE_THAT(result, SimplifiesTo("1/4 g{i1,i2;a1,a2}:A "
                                    "- 1/4 g{i1,i2;a2,a1}:A "
                                    "- 1/4 g{i2,i1;a1,a2}:A "
                                    "+ 1/4 g{i2,i1;a2,a1}:A"));

  // 1/4 * A * g * t1 * t1
  input = ex<Constant>(rational{1, 4}) *
          ex<Tensor>(L"A", bra{L"a_1", L"a_2"}, ket{L"i_1", L"i_2"},
                     Symmetry::Antisymm) *
          ex<Tensor>(L"g", bra{L"a_1", L"a_2"}, ket{L"a_3", L"a_4"},
                     Symmetry::Antisymm) *
          ex<Tensor>(L"t", bra{L"a_3"}, ket{L"i_1"}) *
          ex<Tensor>(L"t", bra{L"a_4"}, ket{L"i_2"});
  result = expand_A_op(input);
  REQUIRE(result->is<Sum>());
  REQUIRE(result->size() == 4);
  REQUIRE_THAT(result,
               SimplifiesTo("1/4 g{a1,a2;a3,a4}:A t{a3;i1} t{a4;i2} "
                            "- 1/4 g{a2,a1;a3,a4}:A t{a3;i1} t{a4;i2} "
                            "- 1/4 g{a1,a2;a3,a4}:A t{a3;i2} t{a4;i1} "
                            "+ 1/4 g{a2,a1;a3,a4}:A t{a3;i2} t{a4;i1}"));

  // 1/4 * A * g * t1 * t1 * t1 * t1
  input = ex<Constant>(rational{1, 4}) *
          ex<Tensor>(L"A", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"},
                     Symmetry::Antisymm) *
          ex<Tensor>(L"g", bra{L"i_3", L"i_4"}, ket{L"a_3", L"a_4"},
                     Symmetry::Antisymm) *
          ex<Tensor>(L"t", bra{L"a_3"}, ket{L"i_1"}) *
          ex<Tensor>(L"t", bra{L"a_4"}, ket{L"i_2"}) *
          ex<Tensor>(L"t", bra{L"a_1"}, ket{L"i_3"}) *
          ex<Tensor>(L"t", bra{L"a_2"}, ket{L"i_4"});
  result = expand_A_op(input);
  REQUIRE_THAT(
      result,
      SimplifiesTo(
          "1/4 g{i3,i4;a3,a4}:A t{a3;i1} t{a4;i2} t{a1;i3} t{a2;i4} "
          "- 1/4 g{i3,i4;a3,a4}:A t{a3;i2} t{a4;i1} t{a1;i3} t{a2;i4} "
          "- 1/4 g{i3,i4;a3,a4}:A t{a3;i1} t{a4;i2} t{a2;i3} t{a1;i4} "
          "+ 1/4 g{i3,i4;a3,a4}:A t{a3;i2} t{a4;i1} t{a2;i3} t{a1;i4}"));
}

// 3-body
{
  auto input = ex<Tensor>(L"t", bra{L"a_1", L"a_2", L"a_3"},
                          ket{L"i_1", L"i_2", L"i_3"}, Symmetry::Antisymm);
  auto result = expand_A_op(input);
  REQUIRE_THAT(result, SimplifiesTo("t{a1,a2,a3;i1,i2,i3}:A"));

  input = ex<Tensor>(L"A", bra{L"i_1", L"i_2", L"i_3"},
                     ket{L"a_1", L"a_2", L"a_3"}, Symmetry::Antisymm) *
          ex<Tensor>(L"t", bra{L"a_1", L"a_2", L"a_3"},
                     ket{L"i_1", L"i_2", L"i_3"}, Symmetry::Antisymm);
  result = expand_A_op(input);
  REQUIRE(result->is<Sum>());
  REQUIRE(result->size() == 36);
}

{  // 4-body
  const auto input =
      ex<Tensor>(L"A", bra{L"i_1", L"i_2", L"i_3", L"i_4"},
                 ket{L"a_1", L"a_2", L"a_3", L"a_4"}, Symmetry::Antisymm) *
      ex<Tensor>(L"t", bra{L"a_1", L"a_2", L"a_3", L"a_4"},
                 ket{L"i_1", L"i_2", L"i_3", L"i_4"}, Symmetry::Antisymm);
  auto asm_input = expand_A_op(input);
  REQUIRE(asm_input->size() == 576);
  REQUIRE(asm_input->is<Sum>());
}

#ifndef SEQUANT_SKIP_LONG_TESTS
{  // 5-body
  const auto input =
      ex<Tensor>(L"A", bra{L"i_1", L"i_2", L"i_3", L"i_4", L"i_5"},
                 ket{L"a_1", L"a_2", L"a_3", L"a_4", L"a_5"},
                 Symmetry::Antisymm) *
      ex<Tensor>(L"t", bra{L"a_1", L"a_2", L"a_3", L"a_4", L"a_5"},
                 ket{L"i_1", L"i_2", L"i_3", L"i_4", L"i_5"},
                 Symmetry::Antisymm);
  auto asm_input = expand_A_op(input);
  REQUIRE(asm_input->size() == 14400);
  REQUIRE(asm_input->is<Sum>());
}
#endif
}

SECTION("Expand Symmetrizer") {
  {  // 2-body
    const auto input = ex<Tensor>(L"S", bra{L"i_1", L"i_2"},
                                  ket{L"a_1", L"a_2"}, Symmetry::Nonsymm) *
                       ex<Tensor>(L"t", bra{L"a_1", L"a_2"},
                                  ket{L"i_1", L"i_2"}, Symmetry::Antisymm);
    auto result = S_maps(input);
    REQUIRE(result->size() == 2);
    REQUIRE(result->is<Sum>());
    REQUIRE_THAT(result, SimplifiesTo("t{a1,a2;i1,i2}:A + t{a2,a1;i2,i1}:A"));
  }

  {  // 3-body
    const auto input =
        ex<Tensor>(L"S", bra{L"i_1", L"i_2", L"i_3"},
                   ket{L"a_1", L"a_2", L"a_3"}, Symmetry::Nonsymm) *
        ex<Tensor>(L"t", bra{L"a_1", L"a_2", L"a_3"},
                   ket{L"i_1", L"i_2", L"i_3"}, Symmetry::Antisymm);
    auto result = S_maps(input);
    REQUIRE_THAT(result, SimplifiesTo("t{a1,a2,a3;i1,i2,i3}:A "
                                      "+ t{a1,a3,a2;i1,i3,i2}:A "
                                      "+ t{a2,a1,a3;i2,i1,i3}:A "
                                      "+ t{a2,a3,a1;i2,i3,i1}:A "
                                      "+ t{a3,a1,a2;i3,i1,i2}:A "
                                      "+ t{a3,a2,a1;i3,i2,i1}:A"));
  }

  {  // 4-body
    const auto input =
        ex<Tensor>(L"S", bra{L"i_1", L"i_2", L"i_3", L"i_4"},
                   ket{L"a_1", L"a_2", L"a_3", L"a_4"}, Symmetry::Nonsymm) *
        ex<Tensor>(L"t", bra{L"a_1", L"a_2", L"a_3", L"a_4"},
                   ket{L"i_1", L"i_2", L"i_3", L"i_4"}, Symmetry::Antisymm);
    auto result = S_maps(input);
    REQUIRE_THAT(result, SimplifiesTo("t{a1,a2,a3,a4;i1,i2,i3,i4}:A "
                                      "+ t{a1,a2,a4,a3;i1,i2,i4,i3}:A "
                                      "+ t{a1,a3,a2,a4;i1,i3,i2,i4}:A "
                                      "+ t{a1,a3,a4,a2;i1,i3,i4,i2}:A "
                                      "+ t{a1,a4,a2,a3;i1,i4,i2,i3}:A "
                                      "+ t{a1,a4,a3,a2;i1,i4,i3,i2}:A "
                                      "+ t{a2,a1,a3,a4;i2,i1,i3,i4}:A "
                                      "+ t{a2,a1,a4,a3;i2,i1,i4,i3}:A "
                                      "+ t{a2,a3,a1,a4;i2,i3,i1,i4}:A "
                                      "+ t{a2,a3,a4,a1;i2,i3,i4,i1}:A "
                                      "+ t{a2,a4,a1,a3;i2,i4,i1,i3}:A "
                                      "+ t{a2,a4,a3,a1;i2,i4,i3,i1}:A "
                                      "+ t{a3,a1,a2,a4;i3,i1,i2,i4}:A "
                                      "+ t{a3,a1,a4,a2;i3,i1,i4,i2}:A "
                                      "+ t{a3,a2,a1,a4;i3,i2,i1,i4}:A "
                                      "+ t{a3,a2,a4,a1;i3,i2,i4,i1}:A "
                                      "+ t{a3,a4,a1,a2;i3,i4,i1,i2}:A "
                                      "+ t{a3,a4,a2,a1;i3,i4,i2,i1}:A "
                                      "+ t{a4,a1,a2,a3;i4,i1,i2,i3}:A "
                                      "+ t{a4,a1,a3,a2;i4,i1,i3,i2}:A "
                                      "+ t{a4,a2,a1,a3;i4,i2,i1,i3}:A "
                                      "+ t{a4,a2,a3,a1;i4,i2,i3,i1}:A "
                                      "+ t{a4,a3,a1,a2;i4,i3,i1,i2}:A "
                                      "+ t{a4,a3,a2,a1;i4,i3,i2,i1}:A"));
  }

  {
    const auto input =
        ex<Constant>(4) *
        ex<Tensor>(L"S", bra{L"i_1", L"i_2", L"i_3"},
                   ket{L"a_1", L"a_2", L"a_3"}, Symmetry::Nonsymm) *
        ex<Tensor>(L"g", bra{L"i_4", L"i_5"}, ket{L"a_4", L"a_5"},
                   Symmetry::Nonsymm) *
        ex<Tensor>(L"t", bra{L"a_3"}, ket{L"i_4"}) *
        ex<Tensor>(L"t", bra{L"a_5"}, ket{L"i_1"}) *
        ex<Tensor>(L"t", bra{L"a_4"}, ket{L"i_2"}) *
        ex<Tensor>(L"t", bra{L"a_1", L"a_2"}, ket{L"i_5", L"i_3"});
    auto result = S_maps(input);
    REQUIRE(result->is<Sum>());
    REQUIRE(result->size() == 6);
    result->canonicalize();
    rapid_simplify(result);
    REQUIRE_THAT(
        result,
        EquivalentTo(
            "4 g{i4,i5;a4,a5} t{a2;i4} t{a4;i3} t{a5;i1} t{a1,a3;i5,i2} + "
            "4 g{i4,i5;a4,a5} t{a1;i5} t{a4;i2} t{a5;i3} t{a2,a3;i4,i1} + "
            "4 g{i4,i5;a4,a5} t{a3;i4} t{a4;i1} t{a5;i2} t{a1,a2;i3,i5} + "
            "4 g{i4,i5;a4,a5} t{a1;i5} t{a4;i3} t{a5;i2} t{a2,a3;i1,i4} + "
            "4 g{i4,i5;a4,a5} t{a2;i4} t{a4;i1} t{a5;i3} t{a1,a3;i2,i5} + "
            "4 g{i4,i5;a4,a5} t{a3;i4} t{a4;i2} t{a5;i1} t{a1,a2;i5,i3}"));
  }
}

SECTION("partial expansion + S_maps = full expansion") {
  auto input = ex<Constant>(rational{1, 4}) *
               ex<Tensor>(L"A", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"},
                          Symmetry::Antisymm) *
               ex<Tensor>(L"t", bra{L"a_1", L"a_2"}, ket{L"i_1", L"i_2"},
                          Symmetry::Antisymm);
  auto result = symmetrize_expr(input);
  REQUIRE_THAT(
      result, SimplifiesTo(
                  "1/4 S{i_1,i_2;a_1,a_2}:N-C-S * t{a_1,a_2;i_1,i_2}:A-C-S "
                  "- 1/4 S{i_1,i_2;a_1,a_2}:N-C-S * t{a_1,a_2;i_2,i_1}:A-C-S"));
  //(canonicalized: 1/2 S{i_1,i_2;a_1,a_2}:N-C-S * t{a_1,a_2;i_1,i_2}:A-C-S)

  result = S_maps(result);
  REQUIRE_THAT(
      result,
      SimplifiesTo(
          "1/4 t{a_1,a_2;i_1,i_2}:A-C-S + 1/4 t{a_2,a_1;i_2,i_1}:A-C-S "
          "- 1/4 t{a_1,a_2;i_2,i_1}:A-C-S - 1/4 t{a_2,a_1;i_1,i_2}:A-C-S"));
  // (canonicalized: t{a_1,a_2;i_1,i_2}:A-C-S)

  result = expand_A_op(input);
  REQUIRE_THAT(
      result,
      SimplifiesTo(
          "1/4 t{a_1,a_2;i_1,i_2}:A-C-S - 1/4 t{a_1,a_2;i_2,i_1}:A-C-S "
          "- 1/4 t{a_2,a_1;i_1,i_2}:A-C-S + 1/4 t{a_2,a_1;i_2,i_1}:A-C-S"));
  // (canonicalized: t{a_1,a_2;i_1,i_2}:A-C-S)
}

SECTION("partial spintracing + S_maps = full spintracing") {
  auto input = ex<Constant>(rational{1, 4}) *
               ex<Tensor>(L"A", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"},
                          Symmetry::Antisymm) *
               ex<Tensor>(L"t", bra{L"a_1", L"a_2"}, ket{L"i_1", L"i_2"},
                          Symmetry::Antisymm);
  auto result =
      closed_shell_spintrace(input, {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}});
  REQUIRE_THAT(
      result,
      EquivalentTo("-1 S{i_1,i_2;a_1,a_2}:N-C-S * t{a_1,a_2;i_2,i_1}:N-C-S "
                   "+ 2 S{i_1,i_2;a_1,a_2}:N-C-S * t{a_1,a_2;i_1,i_2}:N-C-S"));
  result = S_maps(result);
  REQUIRE_THAT(
      result, EquivalentTo(
                  "-1 t{a_1,a_2;i_2,i_1}:N-C-S - 1 t{a_2,a_1;i_1,i_2}:N-C-S "
                  "+ 2 t{a_1,a_2;i_1,i_2}:N-C-S + 2 t{a_2,a_1;i_2,i_1}:N-C-S"));
  // (canonicalized: -2 t{a_1,a_2;i_2,i_1}:N-C-S + 4 t{a_1,a_2;i_1,i_2}:N-C-S)

  result =
      closed_shell_spintrace(input, {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}}, true);
  REQUIRE_THAT(
      result, EquivalentTo(
                  "-1 t{a_1,a_2;i_2,i_1}:N-C-S - 1 t{a_2,a_1;i_1,i_2}:N-C-S "
                  "+ 2 t{a_1,a_2;i_1,i_2}:N-C-S + 2 t{a_2,a_1;i_2,i_1}:N-C-S"));
  //(canonicalized: -2 t{a_1,a_2;i_2,i_1}:N-C-S + 4 t{a_1,a_2;i_1,i_2}:N-C-S)
}

SECTION("Symmetrize expression") {
  {
    // g * t1 + g * t1
    auto input = ex<Tensor>(L"g", bra{L"a_1", L"a_2"}, ket{L"i_1", L"a_3"},
                            Symmetry::Symm) *
                     ex<Tensor>(L"t", bra{L"a_3"}, ket{L"i_2"}) +
                 ex<Tensor>(L"g", bra{L"a_2", L"a_1"}, ket{L"i_2", L"a_3"},
                            Symmetry::Symm) *
                     ex<Tensor>(L"t", bra{L"a_3"}, ket{L"i_1"});
    auto result =
        factorize_S(input, {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}}, true);
    REQUIRE_THAT(result,
                 EquivalentTo("S{i1,i2;a1,a2} g{a1,a2;i2,a3}:S t{a3;i1}"));
  }

  {
    // g * t1 * t1 * t1 + g * t1 * t1 * t1
    auto input = ex<Tensor>(L"g", bra{L"i_3", L"i_4"}, ket{L"i_1", L"a_3"},
                            Symmetry::Symm) *
                     ex<Tensor>(L"t", bra{L"a_1"}, ket{L"i_3"}) *
                     ex<Tensor>(L"t", bra{L"a_2"}, ket{L"i_4"}) *
                     ex<Tensor>(L"t", bra{L"a_3"}, ket{L"i_2"}) +
                 ex<Tensor>(L"g", bra{L"i_3", L"i_4"}, ket{L"i_2", L"a_3"},
                            Symmetry::Symm) *
                     ex<Tensor>(L"t", bra{L"a_2"}, ket{L"i_3"}) *
                     ex<Tensor>(L"t", bra{L"a_1"}, ket{L"i_4"}) *
                     ex<Tensor>(L"t", bra{L"a_3"}, ket{L"i_1"});
    auto result = factorize_S(input, {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}});
    REQUIRE_THAT(
        result,
        EquivalentTo(
            "S{i3,i4;a2,a1} g{i1,i2;i4,a3}:S t{a1;i1} t{a2;i2} t{a3;i3}"));
  }

  {
    // g * t1 * t1 * t2 + g * t1 * t1 * t2
    auto input =
        ex<Constant>(2) *
            ex<Tensor>(L"g", bra{L"i_3", L"i_4"}, ket{L"a_3", L"a_4"},
                       Symmetry::Symm) *
            ex<Tensor>(L"t", bra{L"a_3"}, ket{L"i_3"}) *
            ex<Tensor>(L"t", bra{L"a_2"}, ket{L"i_4"}) *
            ex<Tensor>(L"t", bra{L"a_1", L"a_4"}, ket{L"i_1", L"i_2"}) +
        ex<Constant>(2) *
            ex<Tensor>(L"g", bra{L"i_3", L"i_4"}, ket{L"a_3", L"a_4"},
                       Symmetry::Symm) *
            ex<Tensor>(L"t", bra{L"a_3"}, ket{L"i_3"}) *
            ex<Tensor>(L"t", bra{L"a_1"}, ket{L"i_4"}) *
            ex<Tensor>(L"t", bra{L"a_2", L"a_4"}, ket{L"i_2", L"i_1"});
    auto result =
        factorize_S(input, {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}}, true);
    REQUIRE(result->is<Sum>() == false);
    REQUIRE_THAT(result, EquivalentTo("2 S{i1,i2;a1,a2} g{i3,i4;a3,a4}:S "
                                      "t{a4;i3} t{a2;i4} t{a1,a3;i1,i2}"));
  }
}

SECTION("Swap bra kets") {
  // Constant
  {
    auto input = ex<Constant>(rational{1, 2});
    auto result = swap_bra_ket(input);
    REQUIRE_THAT(result, EquivalentTo("1/2"));
  }

  // Tensor
  {
    auto input = ex<Tensor>(L"g", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"},
                            Symmetry::Nonsymm);
    auto result = swap_bra_ket(input);
    REQUIRE_THAT(result, EquivalentTo("g{a1,a2;i1,i2}"));
  }

  // Product
  {
    auto input = ex<Tensor>(L"g", bra{L"a_5", L"a_6"}, ket{L"i_5", L"i_6"},
                            Symmetry::Nonsymm) *
                 ex<Tensor>(L"t", bra{L"i_2"}, ket{L"a_6"});
    auto result = swap_bra_ket(input);
    REQUIRE_THAT(result, EquivalentTo("g{i5,i6;a5,a6} t{a6;i2}"));
  }

  // Sum
  {
    auto input = ex<Tensor>(L"f", bra{L"i_1"}, ket{L"i_5"}) +
                 ex<Tensor>(L"g", bra{L"a_5", L"a_6"}, ket{L"i_5", L"i_6"},
                            Symmetry::Nonsymm) *
                     ex<Tensor>(L"t", bra{L"i_2"}, ket{L"a_6"});
    auto result = swap_bra_ket(input);
    // TODO: This should be EquivalentTo but canonicalization currently doesn't
    // permit expressions that break co/contra variance of index contractions.
    REQUIRE_THAT(result, SimplifiesTo("f{i5;i1} + g{i5,i6;a5,a6} t{a6;i2}"));
  }
}

SECTION("Closed-shell spintrace CCD") {
  // Energy expression
  {
    {  // standard = v1
      const auto input = ex<Sum>(ExprPtrList{parse_expr(
          L"1/4 g{i_1,i_2;a_1,a_2} t{a_1,a_2;i_1,i_2}", Symmetry::Antisymm)});
      auto result = closed_shell_CC_spintrace_v1(input);
      REQUIRE_THAT(result,
                   EquivalentTo(L"- g{i_1,i_2;a_1,a_2} t{a_1,a_2;i_2,i_1} + "
                                L"2 g{i_1,i_2;a_1,a_2} t{a_1,a_2;i_1,i_2}"));
    }
    {  // compact = v2
      const auto input = ex<Sum>(ExprPtrList{parse_expr(
          L"1/4 g{i_1,i_2;a_1,a_2} t{a_1,a_2;i_1,i_2}", Symmetry::Antisymm)});

      auto result = closed_shell_CC_spintrace_v2(input);
      REQUIRE_THAT(result,
                   EquivalentTo(L"- g{i_1,i_2;a_1,a_2} t{a_1,a_2;i_2,i_1} + "
                                L"2 g{i_1,i_2;a_1,a_2} t{a_1,a_2;i_1,i_2}"));
    }
    {  // CSV (aka PNO) for regular cs
      const auto pno_ccd_energy_so = parse_expr(
          L"1/4 g{a1<i1,i2>, a2<i1,i2>;i1, i2}:A-C t{i1,i2;a1<i1,i2>, "
          L"a2<i1,i2>}:A");

      // why???
      const auto pno_ccd_energy_so_as_sum =
          ex<Sum>(ExprPtrList{pno_ccd_energy_so});
      auto pno_ccd_energy_sf =
          closed_shell_CC_spintrace_v1(pno_ccd_energy_so_as_sum);
      REQUIRE_THAT(pno_ccd_energy_sf,
                   EquivalentTo("2 g{a1<i1,i2>,a2<i1,i2>;i1,i2}:N-C "
                                "t{i1,i2;a1<i1,i2>,a2<i1,i2>}:N-C - "
                                "g{a1<i1,i2>,a2<i1,i2>;i1,i2}:N-C "
                                "t{i1,i2;a2<i1,i2>,a1<i1,i2>}:N-C"));
    }
    {  // CSV (aka PNO) for more compact equations
      const auto pno_ccd_energy_so = parse_expr(
          L"1/4 g{a1<i1,i2>, a2<i1,i2>;i1, i2}:A-C t{i1,i2;a1<i1,i2>, "
          L"a2<i1,i2>}:A");

      // why???
      const auto pno_ccd_energy_so_as_sum =
          ex<Sum>(ExprPtrList{pno_ccd_energy_so});
      auto pno_ccd_energy_sf =
          closed_shell_CC_spintrace_v2(pno_ccd_energy_so_as_sum);
      REQUIRE_THAT(pno_ccd_energy_sf,
                   EquivalentTo("2 g{a1<i1,i2>,a2<i1,i2>;i1,i2}:N-C "
                                "t{i1,i2;a1<i1,i2>,a2<i1,i2>}:N-C - "
                                "g{a1<i1,i2>,a2<i1,i2>;i1,i2}:N-C "
                                "t{i1,i2;a2<i1,i2>,a1<i1,i2>}:N-C"));
    }
  }
}

SECTION("Closed-shell spintrace CCSD") {
  // These terms from CCSD R1 equations
  {
    // A * f
    const auto input = ex<Tensor>(L"A", bra{L"i_1"}, ket{L"a_1"}) *
                       ex<Tensor>(L"f", bra{L"a_1"}, ket{L"i_1"});
    auto result =
        ex<Constant>(rational{1, 2}) * spintrace(input, {{L"i_1", L"a_1"}});
    expand(result);
    rapid_simplify(result);
    canonicalize(result);
    REQUIRE_THAT(result, EquivalentTo("f{a1;i1}"));
  }

  {
    // - A * f * t1
    const auto input = ex<Constant>(-1) *
                       ex<Tensor>(L"A", bra{L"i_1"}, ket{L"a_1"}) *
                       ex<Tensor>(L"f", bra{L"i_2"}, ket{L"i_1"}) *
                       ex<Tensor>(L"t", bra{L"a_1"}, ket{L"i_2"});
    auto result =
        ex<Constant>(rational{1, 2}) * spintrace(input, {{L"i_1", L"a_1"}});
    expand(result);
    REQUIRE_THAT(result, EquivalentTo("- f{i2;i1} t{a1;i2}"));
  }

  {
    // A * f * t1
    const auto input = ex<Tensor>(L"A", bra{L"i_1"}, ket{L"a_1"}) *
                       ex<Tensor>(L"f", bra{L"a_1"}, ket{L"a_2"}) *
                       ex<Tensor>(L"t", bra{L"a_2"}, ket{L"i_1"});
    auto result =
        ex<Constant>(rational{1, 2}) * spintrace(input, {{L"i_1", L"a_1"}});
    expand(result);
    REQUIRE_THAT(result, EquivalentTo("f{a1;a2} t{a2;i1}"));
  }

  {
    // -1/2 * A * g * t2
    const auto input = ex<Constant>(rational{-1, 2}) *
                       ex<Tensor>(L"A", bra{L"i_1"}, ket{L"a_1"}) *
                       ex<Tensor>(L"g", bra{L"i_2", L"i_3"},
                                  ket{L"i_1", L"a_2"}, Symmetry::Antisymm) *
                       ex<Tensor>(L"t", bra{L"a_1", L"a_2"},
                                  ket{L"i_2", L"i_3"}, Symmetry::Antisymm);
    auto result =
        ex<Constant>(rational{1, 2}) * spintrace(input, {{L"i_1", L"a_1"}});
    expand(result);
    REQUIRE_THAT(
        result,
        EquivalentTo(
            "g{i2,i3;i1,a2} t{a1,a2;i3,i2} - 2 g{i2,i3;a2,i1} t{a1,a2;i3,i2}"));
  }

  {
    // -1/2 * A * g * t2
    const auto input = ex<Constant>(rational{-1, 2}) *
                       ex<Tensor>(L"A", bra{L"i_1"}, ket{L"a_1"}) *
                       ex<Tensor>(L"g", bra{L"i_2", L"a_1"},
                                  ket{L"a_2", L"a_3"}, Symmetry::Antisymm) *
                       ex<Tensor>(L"t", bra{L"a_2", L"a_3"},
                                  ket{L"i_1", L"i_2"}, Symmetry::Antisymm);
    auto result =
        ex<Constant>(rational{1, 2}) * spintrace(input, {{L"i_1", L"a_1"}});
    expand(result);

    REQUIRE_THAT(result, EquivalentTo("- g{a1,i2;a3,a2} t{a2,a3;i1,i2} + 2 "
                                      "g{a1,i2;a3,a2} t{a2,a3;i2,i1}"));
  }

  {
    // A * f * t2
    const auto input =
        ex<Tensor>(L"A", bra{L"i_1"}, ket{L"a_1"}) *
        ex<Tensor>(L"f", bra{L"i_2"}, ket{L"a_2"}, Symmetry::Antisymm) *
        ex<Tensor>(L"t", bra{L"a_1", L"a_2"}, ket{L"i_1", L"i_2"},
                   Symmetry::Antisymm);
    auto result =
        ex<Constant>(rational{1, 2}) * spintrace(input, {{L"i_1", L"a_1"}});
    expand(result);
    REQUIRE_THAT(
        result,
        EquivalentTo("-f{i2;a2} t{a1,a2;i2,i1} + 2 f{i2;a2} t{a1,a2;i1,i2}"));
  }

  {
    // A * g * t1 * t1
    const auto input = ex<Tensor>(L"A", bra{L"i_1"}, ket{L"a_1"}) *
                       ex<Tensor>(L"g", bra{L"i_2", L"a_1"},
                                  ket{L"a_2", L"a_3"}, Symmetry::Antisymm) *
                       ex<Tensor>(L"t", bra{L"a_2"}, ket{L"i_2"}) *
                       ex<Tensor>(L"t", bra{L"a_3"}, ket{L"i_1"});
    auto result =
        ex<Constant>(rational{1, 2}) * spintrace(input, {{L"i_1", L"a_1"}});
    expand(result);
    REQUIRE_THAT(result, EquivalentTo("2 g{a1,i2;a3,a2} t{a2;i2} t{a3;i1} - "
                                      "g{a1,i2;a3,a2} t{a2;i1} t{a3;i2}"));
  }

  {
    // A * g * t2 * t2
    const auto input = ex<Tensor>(L"A", bra{L"i_1"}, ket{L"a_1"}) *
                       ex<Tensor>(L"g", bra{L"i_2", L"i_3"},
                                  ket{L"i_1", L"a_2"}, Symmetry::Antisymm) *
                       ex<Tensor>(L"t", bra{L"a_2"}, ket{L"i_2"}) *
                       ex<Tensor>(L"t", bra{L"a_1"}, ket{L"i_3"});
    auto result =
        ex<Constant>(rational{1, 2}) * spintrace(input, {{L"i_1", L"a_1"}});
    expand(result);
    REQUIRE_THAT(result, EquivalentTo("g{i2,i3;i1,a2} t{a1;i3} t{a2;i2} - 2 "
                                      "g{i2,i3;a2,i1} t{a1;i3} t{a2;i2}"));
  }

  {
    // A * f * t1 * t1
    const auto input = ex<Constant>(-1) *
                       ex<Tensor>(L"A", bra{L"i_1"}, ket{L"a_1"}) *
                       ex<Tensor>(L"f", bra{L"i_2"}, ket{L"a_2"}) *
                       ex<Tensor>(L"t", bra{L"a_2"}, ket{L"i_1"}) *
                       ex<Tensor>(L"t", bra{L"a_1"}, ket{L"i_2"});
    auto result =
        ex<Constant>(rational{1, 2}) * spintrace(input, {{L"i_1", L"a_1"}});
    expand(result);
    REQUIRE_THAT(result, EquivalentTo("-f{i2;a2} t{a1;i2} t{a2;i1}"));
  }

  {
    // -1/2 * A * g * t1 * t2
    const auto input = ex<Constant>(rational{-1, 2}) *
                       ex<Tensor>(L"A", bra{L"i_1"}, ket{L"a_1"}) *
                       ex<Tensor>(L"g", bra{L"i_2", L"i_3"},
                                  ket{L"a_2", L"a_3"}, Symmetry::Antisymm) *
                       ex<Tensor>(L"t", bra{L"a_1"}, ket{L"i_2"}) *
                       ex<Tensor>(L"t", bra{L"a_2", L"a_3"},
                                  ket{L"i_1", L"i_3"}, Symmetry::Antisymm);
    auto result =
        ex<Constant>(rational{1, 2}) * spintrace(input, {{L"i_1", L"a_1"}});
    expand(result);
    REQUIRE_THAT(result,
                 EquivalentTo("g{i2,i3;a2,a3} t{a1;i3} t{a2,a3;i1,i2} - 2 "
                              "g{i2,i3;a2,a3} t{a1;i3} t{a2,a3;i2,i1}"));
  }

  {
    // -1/2 * A * g * t1 * t2
    const auto input = ex<Constant>(rational{-1, 2}) *
                       ex<Tensor>(L"A", bra{L"i_1"}, ket{L"a_1"}) *
                       ex<Tensor>(L"g", bra{L"i_2", L"i_3"},
                                  ket{L"a_2", L"a_3"}, Symmetry::Antisymm) *
                       ex<Tensor>(L"t", bra{L"a_2"}, ket{L"i_1"}) *
                       ex<Tensor>(L"t", bra{L"a_1", L"a_3"},
                                  ket{L"i_2", L"i_3"}, Symmetry::Antisymm);
    auto result =
        ex<Constant>(rational{1, 2}) * spintrace(input, {{L"i_1", L"a_1"}});
    expand(result);
    REQUIRE_THAT(result,
                 EquivalentTo("g{i2,i3;a2,a3} t{a2;i1} t{a1,a3;i3,i2} - 2 "
                              "g{i2,i3;a2,a3} t{a3;i1} t{a1,a2;i3,i2}"));
  }

  {
    // A * g * t1 * t2
    const auto input = ex<Constant>(1) *
                       ex<Tensor>(L"A", bra{L"i_1"}, ket{L"a_1"}) *
                       ex<Tensor>(L"g", bra{L"i_2", L"i_3"},
                                  ket{L"a_2", L"a_3"}, Symmetry::Antisymm) *
                       ex<Tensor>(L"t", bra{L"a_2"}, ket{L"i_2"}) *
                       ex<Tensor>(L"t", bra{L"a_1", L"a_3"},
                                  ket{L"i_1", L"i_3"}, Symmetry::Antisymm);
    auto result =
        ex<Constant>(rational{1, 2}) * spintrace(input, {{L"i_1", L"a_1"}});
    expand(result);
    REQUIRE_THAT(result,
                 EquivalentTo("-2 g{i2,i3;a2,a3} t{a3;i2} t{a1,a2;i1,i3} + 4 "
                              "g{i2,i3;a2,a3} t{a2;i2} t{a1,a3;i1,i3} + "
                              "g{i2,i3;a2,a3} t{a3;i2} t{a1,a2;i3,i1} - 2 "
                              "g{i2,i3;a2,a3} t{a2;i2} t{a1,a3;i3,i1}"));
  }

  {
    // - A * g * t1 * t1 * t1
    auto input = ex<Constant>(-1) *
                 ex<Tensor>(L"g", bra{L"i_2", L"i_3"}, ket{L"a_2", L"a_3"},
                            Symmetry::Antisymm) *
                 ex<Tensor>(L"t", bra{L"a_2"}, ket{L"i_2"}) *
                 ex<Tensor>(L"t", bra{L"a_3"}, ket{L"i_1"}) *
                 ex<Tensor>(L"t", bra{L"a_1"}, ket{L"i_3"});
    auto result =
        ex<Constant>(rational{1, 2}) * spintrace(input, {{L"i_1", L"a_1"}});
    expand(result);
    REQUIRE_THAT(result,
                 EquivalentTo("-2 g{i2,i3;a2,a3} t{a1;i3} t{a2;i2} t{a3;i1} + "
                              "g{i2,i3;a2,a3} t{a1;i3} t{a2;i1} t{a3;i2}"));
  }
}  // CCSD R1

SECTION("Closed-shell spintrace CCSDT terms") {
  SECTION("A3 * f * t3, , spintracing with partial-expansion") {
    auto input = ex<Constant>(rational{1, 12}) *
                 ex<Tensor>(L"A", bra{L"i_1", L"i_2", L"i_3"},
                            ket{L"a_1", L"a_2", L"a_3"}, Symmetry::Antisymm) *
                 ex<Tensor>(L"f", bra{L"i_4"}, ket{L"i_1"}) *
                 ex<Tensor>(L"t", bra{L"a_1", L"a_2", L"a_3"},
                            ket{L"i_2", L"i_3", L"i_4"}, Symmetry::Antisymm);

    auto result = expand_A_op(input);
    REQUIRE(result->size() == 36);
    result = expand_antisymm(result);
    result = closed_shell_spintrace(
        input, {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}, {L"i_3", L"a_3"}});
    simplify(result);
    REQUIRE(result->size() == 4);
    REQUIRE_THAT(
        result,
        EquivalentTo("2 S{i1,i2,i3;a1,a2,a3} f{i4;i3} t{a1,a2,a3;i4,i1,i2} - 4 "
                     "S{i1,i2,i3;a1,a2,a3} f{i4;i3} t{a1,a2,a3;i1,i4,i2} + 4 "
                     "S{i1,i2,i3;a1,a2,a3} f{i4;i3} t{a1,a2,a3;i1,i2,i4} - 2 "
                     "S{i1,i2,i3;a1,a2,a3} f{i4;i3} t{a1,a2,a3;i2,i1,i4}"));
  }

  SECTION("ppl term: A3 * g * t3, spintracing with direct full-expansion") {
    auto input = ex<Constant>(rational{1, 24}) *
                 ex<Tensor>(L"A", bra{L"i_1", L"i_2", L"i_3"},
                            ket{L"a_1", L"a_2", L"a_3"}, Symmetry::Antisymm) *
                 ex<Tensor>(L"g", bra{L"a_1", L"a_2"}, ket{L"a_4", L"a_5"},
                            Symmetry::Antisymm) *
                 ex<Tensor>(L"t", bra{L"a_3", L"a_4", L"a_5"},
                            ket{L"i_1", L"i_2", L"i_3"}, Symmetry::Antisymm);

    auto result_1 = expand_A_op(input);
    REQUIRE(result_1->size() == 36);
    result_1 = expand_antisymm(result_1);
    result_1 = closed_shell_spintrace(
        input, {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}, {L"i_3", L"a_3"}}, true);
    simplify(result_1);
    REQUIRE(result_1->size() == 18);  // 18 raw terms
    REQUIRE_THAT(
        result_1,
        EquivalentTo(
            "  8 g{a_1,a_2;a_4,a_5}:N-C-S * t{a_3,a_4,a_5;i_3,i_1,i_2}:N-C-S + "
            "2 "
            "g{a_1,a_2;a_4,a_5}:N-C-S * t{a_3,a_4,a_5;i_2,i_3,i_1}:N-C-S - 4 "
            "g{a_1,a_3;a_4,a_5}:N-C-S * t{a_2,a_4,a_5;i_3,i_1,i_2}:N-C-S - 4 "
            "g{a_2,a_3;a_4,a_5}:N-C-S * t{a_1,a_4,a_5;i_1,i_3,i_2}:N-C-S - 4 "
            "g{a_1,a_2;a_4,a_5}:N-C-S * t{a_3,a_4,a_5;i_2,i_1,i_3}:N-C-S - 4 "
            "g{a_2,a_3;a_4,a_5}:N-C-S * t{a_1,a_4,a_5;i_2,i_1,i_3}:N-C-S + 2 "
            "g{a_2,a_3;a_4,a_5}:N-C-S * t{a_1,a_4,a_5;i_3,i_1,i_2}:N-C-S - 4 "
            "g{a_1,a_2;a_4,a_5}:N-C-S * t{a_3,a_4,a_5;i_3,i_2,i_1}:N-C-S + 2 "
            "g{a_2,a_3;a_4,a_5}:N-C-S * t{a_1,a_4,a_5;i_2,i_3,i_1}:N-C-S - 4 "
            "g{a_1,a_2;a_4,a_5}:N-C-S * t{a_3,a_4,a_5;i_1,i_3,i_2}:N-C-S - 4 "
            "g{a_1,a_3;a_4,a_5}:N-C-S * t{a_2,a_4,a_5;i_1,i_2,i_3}:N-C-S + 8 "
            "g{a_2,a_3;a_4,a_5}:N-C-S * t{a_1,a_4,a_5;i_1,i_2,i_3}:N-C-S + 8 "
            "g{a_1,a_3;a_4,a_5}:N-C-S * t{a_2,a_4,a_5;i_2,i_1,i_3}:N-C-S + 2 "
            "g{a_1,a_3;a_4,a_5}:N-C-S * t{a_2,a_4,a_5;i_3,i_2,i_1}:N-C-S - 4 "
            "g{a_1,a_3;a_4,a_5}:N-C-S * t{a_2,a_4,a_5;i_2,i_3,i_1}:N-C-S - 4 "
            "g{a_2,a_3;a_4,a_5}:N-C-S * t{a_1,a_4,a_5;i_3,i_2,i_1}:N-C-S + 2 "
            "g{a_1,a_2;a_4,a_5}:N-C-S * t{a_3,a_4,a_5;i_1,i_2,i_3}:N-C-S + 2 "
            "g{a_1,a_3;a_4,a_5}:N-C-S * t{a_2,a_4,a_5;i_1,i_3,i_2}:N-C-S"));

    // the new efficient method, spintracing with partial expansion, then
    // expanding by S_map ( this method is used in
    // closed_shell_CC_spintrace_v2)
    auto result_2 = closed_shell_spintrace(
        input, {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}, {L"i_3", L"a_3"}});
    simplify(result_2);
    result_2 = S_maps(result_2);
    simplify(result_2);

    REQUIRE(result_2->size() == 18);  // 18 raw terms
    REQUIRE_THAT(
        result_2,
        EquivalentTo(
            "8 g{a_1,a_2;a_4,a_5}:N-C-S * t{a_3,a_4,a_5;i_3,i_1,i_2}:N-C-S + 2 "
            "g{a_1,a_2;a_4,a_5}:N-C-S * t{a_3,a_4,a_5;i_2,i_3,i_1}:N-C-S - 4 "
            "g{a_1,a_3;a_4,a_5}:N-C-S * t{a_2,a_4,a_5;i_3,i_1,i_2}:N-C-S - 4 "
            "g{a_2,a_3;a_4,a_5}:N-C-S * t{a_1,a_4,a_5;i_1,i_3,i_2}:N-C-S - 4 "
            "g{a_1,a_2;a_4,a_5}:N-C-S * t{a_3,a_4,a_5;i_2,i_1,i_3}:N-C-S - 4 "
            "g{a_2,a_3;a_4,a_5}:N-C-S * t{a_1,a_4,a_5;i_2,i_1,i_3}:N-C-S + 2 "
            "g{a_2,a_3;a_4,a_5}:N-C-S * t{a_1,a_4,a_5;i_3,i_1,i_2}:N-C-S - 4 "
            "g{a_1,a_2;a_4,a_5}:N-C-S * t{a_3,a_4,a_5;i_3,i_2,i_1}:N-C-S + 2 "
            "g{a_2,a_3;a_4,a_5}:N-C-S * t{a_1,a_4,a_5;i_2,i_3,i_1}:N-C-S - 4 "
            "g{a_1,a_2;a_4,a_5}:N-C-S * t{a_3,a_4,a_5;i_1,i_3,i_2}:N-C-S - 4 "
            "g{a_1,a_3;a_4,a_5}:N-C-S * t{a_2,a_4,a_5;i_1,i_2,i_3}:N-C-S + 8 "
            "g{a_2,a_3;a_4,a_5}:N-C-S * t{a_1,a_4,a_5;i_1,i_2,i_3}:N-C-S + 8 "
            "g{a_1,a_3;a_4,a_5}:N-C-S * t{a_2,a_4,a_5;i_2,i_1,i_3}:N-C-S + 2 "
            "g{a_1,a_3;a_4,a_5}:N-C-S * t{a_2,a_4,a_5;i_3,i_2,i_1}:N-C-S - 4 "
            "g{a_1,a_3;a_4,a_5}:N-C-S * t{a_2,a_4,a_5;i_2,i_3,i_1}:N-C-S - 4 "
            "g{a_2,a_3;a_4,a_5}:N-C-S * t{a_1,a_4,a_5;i_3,i_2,i_1}:N-C-S + 2 "
            "g{a_1,a_2;a_4,a_5}:N-C-S * t{a_3,a_4,a_5;i_1,i_2,i_3}:N-C-S + 2 "
            "g{a_1,a_3;a_4,a_5}:N-C-S * t{a_2,a_4,a_5;i_1,i_3,i_2}:N-C-S"));
  }

  SECTION("ppl term in optimal") {  // results in 1 term
    const auto input = ex<Sum>(ExprPtrList{
        parse_expr(L"1/24 A{i_1,i_2,i_3;a_1,a_2,a_3} * "
                   L"g{a_1,a_2;a_4,a_5} * t{a_3,a_4,a_5;i_1,i_2,i_3}",
                   Symmetry::Antisymm)});

    auto result = closed_shell_CC_spintrace_v2(input);

    REQUIRE_THAT(
        result,
        EquivalentTo(
            L"1/2 S{i_1,i_2,i_3;a_1,a_2,a_3}:N-C-S * "
            "g{a_1,a_2;a_4,a_5}:N-C-S * t{a_3,a_4,a_5;i_3,i_1,i_2}:N-C-S"));
  }

  SECTION("ppl term in regular_cs") {  // results in 4 terms
    const auto input = ex<Sum>(ExprPtrList{
        parse_expr(L"1/24 A{i_1,i_2,i_3;a_1,a_2,a_3} * g{a_1,a_2;a_4,a_5} * "
                   "t{a_3,a_4,a_5;i_1,i_2,i_3}",
                   Symmetry::Antisymm)});

    auto result = closed_shell_CC_spintrace_v1(input);
    REQUIRE(result->size() == 4);
    REQUIRE_THAT(
        result,
        EquivalentTo(
            L"-1/5 S{i_1,i_2,i_3;a_1,a_2,a_3}:N-C-S * g{a_1,a_2;a_4,a_5}:N-C-S "
            L"* "
            "t{a_3,a_4,a_5;i_1,i_2,i_3}:N-C-S + 1/2 "
            "S{i_1,i_2,i_3;a_1,a_2,a_3}:N-C-S * "
            "g{a_1,a_2;a_4,a_5}:N-C-S * t{a_3,a_4,a_5;i_3,i_1,i_2}:N-C-S - "
            "1/10 "
            "S{i_1,i_2,i_3;a_1,a_2,a_3}:N-C-S * g{a_1,a_2;a_4,a_5}:N-C-S * "
            "t{a_3,a_4,a_5;i_3,i_2,i_1}:N-C-S - 1/5 "
            "S{i_1,i_2,i_3;a_1,a_2,a_3}:N-C-S * "
            "g{a_1,a_2;a_4,a_5}:N-C-S * t{a_3,a_4,a_5;i_2,i_1,i_3}:N-C-S"));
  }

  SECTION("f * t3") {
    auto input = ex<Tensor>(L"f", bra{L"i_4"}, ket{L"i_1"}) *
                 ex<Tensor>(L"t", bra{L"a_1", L"a_2", L"a_3"},
                            ket{L"i_2", L"i_3", L"i_4"}, Symmetry::Antisymm);

    auto result = expand_A_op(input);
    REQUIRE(result->size() == 2);
    result = expand_antisymm(result);
    result = closed_shell_spintrace(
        input, {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}, {L"i_3", L"a_3"}});
    REQUIRE(result->size() == 6);
    simplify(result);
    REQUIRE(result->size() == 6);
  }

  SECTION("A * g * t3") {
    auto input = ex<Constant>(rational{-1, 4}) *
                 ex<Tensor>(L"A", bra{L"i_1", L"i_2", L"i_3"},
                            ket{L"a_1", L"a_2", L"a_3"}, Symmetry::Antisymm) *
                 ex<Tensor>(L"g", bra{L"i_4", L"a_1"}, ket{L"i_1", L"a_4"},
                            Symmetry::Antisymm) *
                 ex<Tensor>(L"t", bra{L"a_2", L"a_3", L"a_4"},
                            ket{L"i_2", L"i_3", L"i_4"}, Symmetry::Antisymm);
    auto result = expand_A_op(input);
    REQUIRE(result->size() == 36);
    result = expand_antisymm(result);
    result = closed_shell_spintrace(
        input, {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}, {L"i_3", L"a_3"}});
    REQUIRE(result->size() == 20);
  }

  SECTION("g * t3") {
    auto input = ex<Tensor>(L"g", bra{L"i_4", L"a_1"}, ket{L"i_1", L"a_4"},
                            Symmetry::Antisymm) *
                 ex<Tensor>(L"t", bra{L"a_2", L"a_3", L"a_4"},
                            ket{L"i_2", L"i_3", L"i_4"}, Symmetry::Antisymm);
    auto result = expand_A_op(input);
    REQUIRE(result->size() == 2);
    result = expand_antisymm(result);
    result = closed_shell_spintrace(
        input, {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}, {L"i_3", L"a_3"}});
    REQUIRE(result->size() == 12);
    simplify(result);
    REQUIRE(result->size() == 12);
  }
}

SECTION("Merge P operators") {
  auto P1 = Tensor(L"P", bra{L"i_1", L"i_2"}, ket{}, Symmetry::Symm);
  auto P2 = Tensor(L"P", bra{}, ket{L"a_1", L"a_2"}, Symmetry::Symm);
  auto P3 =
      Tensor(L"P", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"}, Symmetry::Symm);
  auto P4 = Tensor(L"P", bra{}, ket{}, Symmetry::Symm);
  auto P12 = merge_tensors(P1, P2);
  auto P34 = merge_tensors(P3, P4);
  REQUIRE_THAT(P12, EquivalentTo("P{i1,i2;a1,a2}:S"));
  REQUIRE_THAT(P34, EquivalentTo("P{i1,i2;a1,a2}:S"));
}

SECTION("Permutation operators") {
  auto A_12 = ex<Tensor>(L"A", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"},
                         Symmetry::Antisymm);
  auto A_23 = ex<Tensor>(L"A", bra{L"i_2", L"i_3"}, ket{L"a_2", L"a_3"},
                         Symmetry::Antisymm);
  auto A2 = ex<Tensor>(L"A", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"},
                       Symmetry::Antisymm);
  auto A3 = ex<Tensor>(L"A", bra{L"i_1", L"i_2", L"i_3"},
                       ket{L"a_1", L"a_2", L"a_3"}, Symmetry::Antisymm);
  auto A4 = ex<Tensor>(L"A", bra{L"i_1", L"i_2", L"i_3", L"i_4"},
                       ket{L"a_1", L"a_2", L"a_3", L"a_4"}, Symmetry::Antisymm);
  auto A5 = ex<Tensor>(L"A", bra{L"i_1", L"i_2", L"i_3", L"i_4", L"i_5"},
                       ket{L"a_1", L"a_2", L"a_3", L"a_4", L"a_5"},
                       Symmetry::Antisymm);

  auto Avec2 = open_shell_A_op(A2->as<Tensor>());
  auto P3 = open_shell_P_op_vector(A3->as<Tensor>());
  auto Avec3 = open_shell_A_op(A3->as<Tensor>());
  REQUIRE(P3[0]->size() == 0);
  REQUIRE(P3[1]->size() == 9);
  REQUIRE(P3[2]->size() == 9);
  REQUIRE(P3[3]->size() == 0);

  auto P4 = open_shell_P_op_vector(A4->as<Tensor>());
  auto Avec4 = open_shell_A_op(A4->as<Tensor>());
  REQUIRE(P4[0]->size() == 0);
  REQUIRE(P4[1]->size() == 16);
  REQUIRE(P4[2]->size() == 36);
  REQUIRE(P4[3]->size() == 16);
  REQUIRE(P4[4]->size() == 0);

  auto P5 = open_shell_P_op_vector(A5->as<Tensor>());
  auto Avec5 = open_shell_A_op(A5->as<Tensor>());
  REQUIRE(P5[0]->size() == 0);
  REQUIRE(P5[1]->size() == 25);
  REQUIRE(P5[2]->size() == 100);
  REQUIRE(P5[3]->size() == 100);
  REQUIRE(P5[4]->size() == 25);
  REQUIRE(P5[5]->size() == 0);
}

SECTION("Relation in spin P operators") {
  auto input = ex<Tensor>(L"g", bra{L"i_4", L"a_1"}, ket{L"i_1", L"i_2"},
                          Symmetry::Antisymm) *
               ex<Tensor>(L"t", bra{L"a_2", L"a_3"}, ket{L"i_3", L"i_4"},
                          Symmetry::Antisymm);

  auto P13_b = ex<Tensor>(L"P", bra{}, ket{L"a_1", L"a_3"}, Symmetry::Nonsymm);
  auto P13_k = ex<Tensor>(L"P", bra{L"i_1", L"i_3"}, ket{}, Symmetry::Nonsymm);
  auto P12_b = ex<Tensor>(L"P", bra{}, ket{L"a_1", L"a_2"}, Symmetry::Nonsymm);
  auto P12_k = ex<Tensor>(L"P", bra{L"i_1", L"i_2"}, ket{}, Symmetry::Nonsymm);

  auto P23_b = ex<Tensor>(L"P", bra{}, ket{L"a_2", L"a_3"}, Symmetry::Nonsymm);
  auto P23_k = ex<Tensor>(L"P", bra{L"i_2", L"i_3"}, ket{}, Symmetry::Nonsymm);

  auto P4_1313 = ex<Tensor>(L"P", bra{L"i_1", L"i_3"}, ket{L"a_1", L"a_3"},
                            Symmetry::Nonsymm);
  auto P4_1323 = ex<Tensor>(L"P", bra{L"i_1", L"i_3"}, ket{L"a_2", L"a_3"},
                            Symmetry::Nonsymm);
  auto P4_2313 = ex<Tensor>(L"P", bra{L"i_2", L"i_3"}, ket{L"a_1", L"a_3"},
                            Symmetry::Nonsymm);
  auto P4_2323 = ex<Tensor>(L"P", bra{L"i_2", L"i_3"}, ket{L"a_2", L"a_3"},
                            Symmetry::Nonsymm);

  auto P4_1212 = ex<Tensor>(L"P", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"},
                            Symmetry::Nonsymm);
  auto P4_1213 = ex<Tensor>(L"P", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_3"},
                            Symmetry::Nonsymm);
  auto P4_1312 = ex<Tensor>(L"P", bra{L"i_1", L"i_3"}, ket{L"a_1", L"a_2"},
                            Symmetry::Nonsymm);

  auto p_aab = ex<Constant>(1) - P13_b - P23_b - P13_k - P23_k + P4_1313 +
               P4_1323 + P4_2313 + P4_2323;

  auto p_abb = ex<Constant>(1) - P13_b - P12_b - P13_k - P12_k + P4_1212 +
               P4_1213 + P4_1312 + P4_1313;

  auto p6_input = p_aab * input;
  expand(p6_input);
  auto p6_result = expand_P_op(p6_input);
  p6_result->visit(reset_idx_tags);
  simplify(p6_result);

  auto A_12 = ex<Tensor>(L"A", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"},
                         Symmetry::Antisymm);
  auto A_23 = ex<Tensor>(L"A", bra{L"i_2", L"i_3"}, ket{L"a_2", L"a_3"},
                         Symmetry::Antisymm);
  auto A3 = ex<Tensor>(L"A", bra{L"i_1", L"i_2", L"i_3"},
                       ket{L"a_1", L"a_2", L"a_3"}, Symmetry::Antisymm);

  p6_result = A_12 * p6_result;
  expand(p6_result);
  canonicalize(p6_result);
  p6_result = expand_A_op(p6_result);
  p6_result->visit(reset_idx_tags);
  simplify(p6_result);

  auto p7_input = p_abb * input;
  expand(p7_input);
  auto p7_result = expand_P_op(p7_input);
  p7_result->visit(reset_idx_tags);
  simplify(p7_result);

  p7_result = A_23 * p7_result;
  expand(p7_result);
  p7_result = expand_A_op(p7_result);
  p7_result->visit(reset_idx_tags);
  simplify(p7_result);

  auto expanded_A = A3 * input;
  expanded_A = expand_A_op(expanded_A);
  expanded_A->visit(reset_idx_tags);
  simplify(expanded_A);
  REQUIRE(p6_result == p7_result);
  REQUIRE(p6_result == expanded_A);
}

SECTION("Expand P operator pair-wise") {
  auto P1 = Tensor(L"P", bra{L"i_1", L"i_2"}, ket{});
  auto P2 = Tensor(L"P", bra{L"i_1", L"i_2", L"i_3", L"i_4"}, ket{});
  auto P3 = Tensor(L"P", bra{}, ket{L"a_1", L"a_2"});
  auto P4 = Tensor(L"P", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"});
  auto P5 = Tensor(L"P", bra{L"i_1", L"i_2", L"i_3", L"i_4"},
                   ket{L"a_1", L"a_2", L"a_3", L"a_4"});

  // g* t3
  auto input = ex<Tensor>(L"g", bra{L"i_4", L"a_1"}, ket{L"i_1", L"a_4"},
                          Symmetry::Antisymm) *
               ex<Tensor>(L"t", bra{L"a_2", L"a_3", L"a_4"},
                          ket{L"i_2", L"i_3", L"i_4"}, Symmetry::Antisymm);

  size_t n_p = 0;
  for (auto& P : {P1, P2, P3, P4, P5}) {
    auto term = ex<Tensor>(P) * input;
    expand(term);
    auto result = expand_P_op(term);
    switch (n_p) {
      case 0:
        REQUIRE_THAT(result,
                     EquivalentTo("g{i4,a1;i2,a4}:A t{a2,a3,a4;i1,i3,i4}:A"));
        break;
      case 1:
        REQUIRE_THAT(result,
                     EquivalentTo("g{i3,a1;i2,a4}:A t{a2,a3,a4;i1,i4,i3}:A"));
        break;
      case 2:
        REQUIRE_THAT(result,
                     EquivalentTo("g{i4,a2;i1,a4}:A t{a1,a3,a4;i2,i3,i4}:A"));
        break;
      case 3:
        REQUIRE_THAT(result,
                     EquivalentTo("g{i4,a2;i2,a4}:A t{a1,a3,a4;i1,i3,i4}:A"));
        break;
      case 4:
        REQUIRE_THAT(result,
                     EquivalentTo("g{i3,a2;i2,a3}:A t{a1,a4,a3;i1,i4,i3}:A"));
        break;
      default:
        break;
    }
    ++n_p;
  }

  auto input2 = ex<Tensor>(L"g", bra{L"a_1", L"a_2"}, ket{L"i_1", L"i_2"},
                           Symmetry::Antisymm);
  auto term = ex<Tensor>(P2) * input2;
  expand(term);
  auto result = expand_P_op(term);
  REQUIRE_THAT(result, EquivalentTo("-g{a1,a2;i1,i2}:A"));
}

SECTION("Open-shell spin-tracing") {
  const auto i1A = Index(L"i↑_1");
  const auto i2A = Index(L"i↑_2");
  const auto i3A = Index(L"i↑_3");
  const auto i4A = Index(L"i↑_4");
  const auto i5A = Index(L"i↑_5");
  const auto i1B = Index(L"i↓_1");
  const auto i2B = Index(L"i↓_2");
  const auto i3B = Index(L"i↓_3");

  const auto a1A = Index(L"a↑_1");
  const auto a2A = Index(L"a↑_2");
  const auto a3A = Index(L"a↑_3");
  const auto a1B = Index(L"a↓_1");
  const auto a2B = Index(L"a↓_2");
  const auto a3B = Index(L"a↓_3");

  // Tensor canonicalize
  {
    auto t3 = ex<Tensor>(Tensor(L"t", bra{a3A, a2B, a2A}, ket{i1A, i2B, i3A}));
    auto f = ex<Tensor>(Tensor(L"f", bra{a1A}, ket{a2A}));
    auto ft3 = f * t3;
    ft3->canonicalize();
    REQUIRE_THAT(ft3, EquivalentTo("f{a↑1;a↑2} t{a↑3,a↑2,a↓2;i↑1,i↑3,i↓2}"));
  }

  //  g
  {
    auto input = ex<Constant>(rational{1, 4}) *
                 ex<Tensor>(L"g", bra{L"a_1", L"a_2"}, ket{L"i_1", L"i_2"},
                            Symmetry::Antisymm);
    auto result =
        open_shell_spintrace(input, {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}});
    REQUIRE(result.size() == 3);
    REQUIRE_THAT(result[0], EquivalentTo("1/4 g{a↑1,a↑2;i↑1,i↑2}:A"));
    REQUIRE_THAT(result[1], EquivalentTo("1/4 g{a↑1,a↓2;i↑1,i↓2}"));
    REQUIRE_THAT(result[2], EquivalentTo("1/4 g{a↓1,a↓2;i↓1,i↓2}:A"));
  }

  // f_oo * t2
  {
    auto input = ex<Constant>(rational{1, 2}) *
                 ex<Tensor>(L"f", bra{L"i_3"}, ket{L"i_1"}) *
                 ex<Tensor>(L"t", bra{L"a_1", L"a_2"}, ket{L"i_2", L"i_3"},
                            Symmetry::Antisymm);

    auto result =
        open_shell_spintrace(input, {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}});
    REQUIRE(result.size() == 3);
    REQUIRE_THAT(result[0],
                 EquivalentTo("1/2 f{i↑3;i↑1} t{a↑1,a↑2;i↑2,i↑3}:A"));
    REQUIRE_THAT(result[1], EquivalentTo("-1/2 f{i↑2;i↑1} t{a↑1,a↓2;i↑2,i↓2}"));
    REQUIRE_THAT(result[2],
                 EquivalentTo("1/2 f{i↓3;i↓1} t{a↓1,a↓2;i↓2,i↓3}:A"));
  }

  // g * t1
  {
    auto input = ex<Constant>(rational{1, 2}) *
                 ex<Tensor>(L"g", bra{L"i_3", L"a_1"}, ket{L"i_1", L"i_2"},
                            Symmetry::Antisymm) *
                 ex<Tensor>(L"t", bra{L"a_2"}, ket{L"i_3"}, Symmetry::Nonsymm);
    auto result =
        open_shell_spintrace(input, {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}});
    REQUIRE(result.size() == 3);
    REQUIRE(
        toUtf8(to_latex(result[0])) ==
        toUtf8(L"{{{\\frac{1}{2}}}{\\bar{g}^{{i↑_1}{i↑_2}}_{{i↑_3}{a↑_1}}}{t^{{"
               L"i↑_3}}_{{a↑_2}}}}"));
    REQUIRE(to_latex(result[1]) ==
            L"{{{-\\frac{1}{2}}}{g^{{i↑_1}{i↓_2}}_{{a↑_1}{i↓_1}}}{t^{{i↓_1}}_"
            L"{{a↓_2}}}}");
    REQUIRE(to_latex(result[2]) ==
            L"{{{\\frac{1}{2}}}{\\bar{g}^{{i↓_1}{i↓_2}}_{{i↓_3}{a↓_1}}}{t^{{"
            L"i↓_3}}_{{a↓_2}}}}");
  }

  // f * t3
  {
    auto input = ex<Constant>(rational{1, 12}) *
                 ex<Tensor>(L"f", bra{L"a_1"}, ket{L"a_4"}) *
                 ex<Tensor>(L"t", bra{L"a_2", L"a_3", L"a_4"},
                            ket{L"i_1", L"i_2", L"i_3"}, Symmetry::Antisymm);
    auto result = open_shell_spintrace(
        input, {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}, {L"i_3", L"a_3"}});
    REQUIRE(result.size() == 4);
    REQUIRE_THAT(result[0],
                 EquivalentTo("1/12 f{a↑1;a↑4} t{a↑2,a↑3,a↑4;i↑1,i↑2,i↑3}:A"));
    REQUIRE_THAT(result[1],
                 EquivalentTo("-1/12 f{a↑1;a↑3} t{a↑2,a↑3,a↓3;i↑_1,i↑2,i↓3} + "
                              "1/12 f{a↑1;a↑3} t{a↑2,a↑3,a↓3;i↑2,i↑1,i↓3}"));
    REQUIRE_THAT(result[2],
                 EquivalentTo("-1/12 f{a↑1;a↑2} t{a↑2,a↓2,a↓3;i↑1,i↓3,i↓2} + "
                              "1/12 f{a↑1;a↑2} t{a↑2,a↓2,a↓3;i↑1,i↓2,i↓3}"));
    REQUIRE_THAT(result[3],
                 EquivalentTo("1/12 f{a↓1;a↓4} t{a↓2,a↓3,a↓4;i↓1,i↓2,i↓3}:A"));
  }

  // aab: g*t3 (CCSDT R3 4)
  {
    auto A2_aab =
        Tensor(L"A", bra{i1A, i2A}, ket{a1A, a2A}, Symmetry::Antisymm);
    auto A2_abb =
        Tensor(L"A", bra{i2B, i3B}, ket{a2B, a3B}, Symmetry::Antisymm);

    auto g = Tensor(L"g", bra{i3A, i4A}, ket{i1A, i2A}, Symmetry::Antisymm);
    auto t3 =
        Tensor(L"t", bra{a1A, a2A, a3B}, ket{i3A, i4A, i3B}, Symmetry::Nonsymm);

    auto input = ex<Constant>(rational{1, 12}) * ex<Tensor>(A2_aab) *
                 ex<Tensor>(g) * ex<Tensor>(t3);
    auto result = expand_A_op(input);
    result->visit(reset_idx_tags);
    REQUIRE_THAT(
        result,
        EquivalentTo("1/3 g{i↑3,i↑4;i↑1,i↑2}:A t{a↑1,a↑2,a↓3;i↑3,i↑4,i↓3}:N"));

    g = Tensor(L"g", bra{i4A, i5A}, ket{i1A, i2A}, Symmetry::Antisymm);
    t3 =
        Tensor(L"t", bra{a1A, a2A, a3B}, ket{i4A, i5A, i3B}, Symmetry::Nonsymm);

    input = ex<Constant>(rational{1, 12}) * ex<Tensor>(A2_aab) * ex<Tensor>(g) *
            ex<Tensor>(t3);
    result = expand_A_op(input);
    result->visit(reset_idx_tags);
    REQUIRE_THAT(
        result,
        EquivalentTo("1/3 g{i↑3,i↑4;i↑1,i↑2}:A t{a↑_1,a↑2,a↓3;i↑3,i↑4,i↓3}:N"));
  }

  // CCSDT R3 10 aaa, bbb
  {
    auto input = ex<Constant>(rational{1, 8}) *
                 ex<Tensor>(L"g", bra{L"i_4", L"i_5"}, ket{L"a_4", L"a_5"},
                            Symmetry::Antisymm) *
                 ex<Tensor>(L"t", bra{L"a_1", L"a_4"}, ket{L"i_1", L"i_2"},
                            Symmetry::Antisymm) *
                 ex<Tensor>(L"t", bra{L"a_2", L"a_3", L"a_5"},
                            ket{L"i_3", L"i_4", L"i_5"}, Symmetry::Antisymm);

    auto result = open_shell_spintrace(
        input, {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}, {L"i_3", L"a_3"}});
    REQUIRE(result[0]->size() == 3);
    auto A3_aaa = Tensor(L"A", bra{i1A, i2A, i3A}, ket{a1A, a2A, a3A},
                         Symmetry::Antisymm);
    auto result2 = ex<Tensor>(A3_aaa) * result[0];
    expand(result2);
    result2 = expand_A_op(result2);
    result2->visit(reset_idx_tags);
    canonicalize(result2);
    rapid_simplify(result2);
    REQUIRE(result2->size() == 27);

    auto A3_bbb = Tensor(L"A", bra{i1B, i2B, i3B}, ket{a1B, a2B, a3B},
                         Symmetry::Antisymm);
    auto result3 = ex<Tensor>(A3_bbb) * result[3];
    expand(result3);
    result3 = expand_A_op(result3);
    result3->visit(reset_idx_tags);
    canonicalize(result3);
    rapid_simplify(result3);
    REQUIRE(result3->size() == 27);
  }

  // CCSDT R3 10 aab
  {
    auto input =
        ex<Constant>(rational{1, 8}) *
            ex<Tensor>(L"P", bra{L"i_1", L"i_3"}, ket{L"a_1", L"a_3"},
                       Symmetry::Nonsymm) *
            ex<Tensor>(L"g", bra{L"i_4", L"i_5"}, ket{L"a_4", L"a_5"},
                       Symmetry::Antisymm) *
            ex<Tensor>(L"t", bra{L"a_1", L"a_4"}, ket{L"i_1", L"i_2"},
                       Symmetry::Antisymm) *
            ex<Tensor>(L"t", bra{L"a_2", L"a_3", L"a_5"},
                       ket{L"i_3", L"i_4", L"i_5"}, Symmetry::Antisymm) +
        ex<Constant>(rational{1, 8}) *
            ex<Tensor>(L"P", bra{L"i_2", L"i_3"}, ket{L"a_2", L"a_3"},
                       Symmetry::Nonsymm) *
            ex<Tensor>(L"g", bra{L"i_4", L"i_5"}, ket{L"a_4", L"a_5"},
                       Symmetry::Antisymm) *
            ex<Tensor>(L"t", bra{L"a_1", L"a_4"}, ket{L"i_1", L"i_2"},
                       Symmetry::Antisymm) *
            ex<Tensor>(L"t", bra{L"a_2", L"a_3", L"a_5"},
                       ket{L"i_3", L"i_4", L"i_5"}, Symmetry::Antisymm);

    input = expand_P_op(input);
    input->visit(reset_idx_tags);
    auto result = open_shell_spintrace(
        input, {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}, {L"i_3", L"a_3"}});

    auto result_aab = ex<Tensor>(Tensor(L"A", bra{i1A, i2A}, ket{a1A, a2A},
                                        Symmetry::Antisymm)) *
                      result[1];
    expand(result_aab);
    result_aab = expand_A_op(result_aab);
    result_aab->visit(reset_idx_tags);
    canonicalize(result_aab);
    rapid_simplify(result_aab);
    REQUIRE(result_aab->size() == 18);

    auto input2 = ex<Constant>(rational{1, 8}) *
                  ex<Tensor>(L"A", bra{L"i_1", L"i_2", L"i_3"},
                             ket{L"a_1", L"a_2", L"a_3"}, Symmetry::Antisymm) *
                  ex<Tensor>(L"g", bra{L"i_4", L"i_5"}, ket{L"a_4", L"a_5"},
                             Symmetry::Antisymm) *
                  ex<Tensor>(L"t", bra{L"a_1", L"a_4"}, ket{L"i_1", L"i_2"},
                             Symmetry::Antisymm) *
                  ex<Tensor>(L"t", bra{L"a_2", L"a_3", L"a_5"},
                             ket{L"i_3", L"i_4", L"i_5"}, Symmetry::Antisymm);

    auto result2 = open_shell_spintrace(
        input2, {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}, {L"i_3", L"a_3"}});
    REQUIRE(result2[1]->size() == 24);
  }
}

SECTION("ResultExpr") {
  auto ctx = get_default_context();
  ctx.set(mbpt::make_mr_spaces());
  auto resetter = set_scoped_default_context(ctx);

  const std::vector<std::wstring> inputs = {
      L"R = 1/4",
      L"R = Var",
      L"R{p1,p2;p3,p4}:A = g{p1,p2;p3,p4}:A",
      L"R = f{i1;a1} t{a1;i1}",
      L"R = 1/2 Var g{i1,i2;a1,a2}:A t{a1;i1} t{a2;i2}",
      L"R{a1,a2;i1,u1}:A = g{a1,a2;i1,u1}:A",
      L"R{a1,u1;i1,i2}:A = g{a1,u1;i1,i2}:A",
      L"R{a1,u1;i1,u2}:A = g{a1,u1;i1,u2}:A",
      L"R{a1,u1;i1,u2}:A = f{a1;i1}:A γ{u1;u2}:A + g{a1,u1;i1,u3}:A γ{u3;u2}:A",
  };
  const std::vector<std::vector<std::wstring>> expected_outputs = {
      {L"R = 1/4"},
      {L"R = Var"},
      {L"R{p1,p2;p3,p4} = 4 g{p1,p2;p3,p4} - 2 g{p2,p1;p3,p4}"},
      {L"R = 2 f{i1;a1} t{a1;i1}"},
      {L"R = 2 Var * g{i_1,i_2;a_1,a_2} * t{a_1;i_1} * t{a_2;i_2}"
       " - 1 Var * g{i_1,i_2;a_1,a_2} * t{a_1;i_2} * t{a_2;i_1}"},
      {L"R{a1,a2;i1,u1} = 4 g{a1,a2;i1,u1} - 2 g{a2,a1;i1,u1}"},
      {L"R{a1,u1;i1,i2} = 4 g{a1,u1;i1,i2} - 2 g{a1,u1;i2,i1}"},
      {L"R{a1,u1;u2,i1} = 4 g{u1,a1;i1,u2} - 2 g{u1,a1;u2,i1}",
       L"R{a1,u1;i1,u2} = 4 g{a1,u1;i1,u2} - 2 g{u1,a1;i1,u2}"},
      {
          L"R{a1,u1;u2,i1} = -2 f{a1;i1} γ{u1;u2} "
          "- 2 g{a1,u1;i1,u3} γ{u3;u2} "
          "+ 4 g{a1,u1;u3,i1} γ{u3;u2}",
          L"R{a1,u1;i1,u2} = 4 f{a1;i1} γ{u1;u2} "
          "+ 4 g{a1,u1;i1,u3} γ{u3;u2} "
          "- 2 g{a1,u1;u3,i1} γ{u3;u2}",
      },
  };

  REQUIRE(inputs.size() == expected_outputs.size());

  for (std::size_t i = 0; i < inputs.size(); ++i) {
    CAPTURE(inputs.at(i));
    const ResultExpr input = parse_result_expr(inputs.at(i));

    container::svector<ResultExpr> expected;
    for (std::size_t k = 0; k < expected_outputs.at(i).size(); ++k) {
      expected.push_back(parse_result_expr(expected_outputs.at(i).at(k)));
    }

    SECTION("closed_shell" + std::to_string(i)) {
      container::svector<ResultExpr> actual =
          closed_shell_spintrace(input.clone());

      REQUIRE(actual.size() == expected.size());
      for (std::size_t k = 0; k < expected.size(); ++k) {
        CAPTURE(k);
        REQUIRE_THAT(actual.at(k), EquivalentTo(expected.at(k)));
      }
    }
    SECTION("rigorous" + std::to_string(i)) {
      container::svector<ResultExpr> actual = spintrace(input.clone());

      REQUIRE(actual.size() == expected.size());
      for (std::size_t k = 0; k < expected.size(); ++k) {
        CAPTURE(k);
        REQUIRE_THAT(actual.at(k), EquivalentTo(expected.at(k)));
      }
    }
  }
}
}
