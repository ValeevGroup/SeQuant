//
// Created by Eduard Valeyev on 2019-02-19.
//

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/latex.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/parse.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/core/utility/timer.hpp>
#include <SeQuant/domain/mbpt/context.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/op.hpp>
#include <SeQuant/domain/mbpt/rules/df.hpp>
#include <SeQuant/domain/mbpt/rules/thc.hpp>
#include <SeQuant/domain/mbpt/utils.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include "catch2_sequant.hpp"

#include <iostream>
#include <memory>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

#include "SeQuant/core/utility/debug.hpp"

TEST_CASE("mbpt", "[mbpt]") {
  SECTION("nbody_operators") {
    using namespace sequant;

    SECTION("constructor") {
      // tests 1-space quantum number case
      {
        using namespace sequant::mbpt;

        op_t f1([]() -> std::wstring_view { return L"f"; },
                []() -> ExprPtr {
                  return ex<Tensor>(L"f", bra{L"p_1"}, ket{L"p_2"}) *
                         ex<FNOperator>(cre({L"p_1"}), ann({L"p_2"}));
                },
                [](qns_t& qns) { qns += general_type_qns(1); });

        REQUIRE(f1.label() == L"f");

        {
          // exact compare of intervals
          using namespace boost::numeric::interval_lib::compare::possible;
          REQUIRE(operator==(
              f1()[0], general_type_qns(1)[0]));  // produces single replacement
          REQUIRE(operator!=(
              f1()[0],
              general_type_qns(2)[0]));  // cannot produce double replacement
          /// TODO clearly this test does not make sense for context implicit
          /// size of qns. Need help to reimagine this test.
          // REQUIRE(operator==(f1(qns_t{5, 0}), qns_t{{5, 6}, {0, 1}})); //
        }
      }

      // tests 2-space quantum number case
      {
        using namespace sequant::mbpt;

        // this is fock operator in terms of general spaces
        op_t f_gg([]() -> std::wstring_view { return L"f"; },
                  []() -> ExprPtr {
                    return ex<Tensor>(L"f", bra{L"p_1"}, ket{L"p_2"}) *
                           ex<FNOperator>(cre({L"p_1"}), ann({L"p_2"}));
                  },
                  [](qns_t& qns) { qns += mbpt::general_type_qns(1); });
        // excitation part of the Fock operator
        op_t f_uo([]() -> std::wstring_view { return L"f"; },
                  []() -> ExprPtr {
                    return ex<Tensor>(L"f", bra{L"a_2"}, ket{L"i_2"}) *
                           ex<FNOperator>(cre({L"a_1"}), ann({L"i_2"}));
                  },
                  [](qns_t& qns) { qns += mbpt::excitation_type_qns(1); });

        REQUIRE(f_gg.label() == L"f");
        REQUIRE(f_uo.label() == L"f");

        {
          // comparison

          // exact
          REQUIRE((f_uo() == excitation_type_qns(
                                 1)));  // f_uo produces single excitations
          REQUIRE((f_gg() !=
                   excitation_type_qns(
                       1)));  // f_gg does not produce just single excitations
          /* REQUIRE(f_gg().in(excitation_type_qns(1)));  // f_gg can produce
           single excitations REQUIRE(f_gg().in(deexcitation_type_qns(1)));  //
           f_gg can also produce single de-excitations REQUIRE(f_gg().in( {1, 1,
           0, 0}));  // f_gg can produce replacements within occupieds
           REQUIRE(f_gg().in(
               {0, 0, 1, 1}));  // f_gg can produce replacements within virtuals
           REQUIRE(f_gg().in(
               {1, 1, 1, 1}));  // f_gg cannot produce this double replacements,
                                // but this returns true TODO introduce
           constraints
                                // on the total number of creators/annihilators,
                                // the interval logic does not constrain it
           REQUIRE(f_gg().in(
               {0, 0, 0, 0}));  // f_gg cannot produce a null replacement, but
           this
                                // returns true TODO introduce constraints on
           the
                                // total number of creators/annihilators, the
                                // interval logic does not constrain it
                                */
          // most of these seem like artifacts of fixed interval logic. we can
          // add them back if needed

          /*REQUIRE(
              f_uo().in(excitation_type_qns(1)));  // f_uo can produce single
          excitations REQUIRE(!f_uo().in( deexcitation_type_qns(1)));  // f_uo
          cannot produce single de-excitations REQUIRE(!f_uo().in( {1, 1, 0,
          0}));
          // f_uo can produce replacements withing occupieds REQUIRE(!f_uo().in(
              {0, 0, 1, 1}));  // f_uo can produce replacements withing virtuals
          REQUIRE(!f_uo().in(
              {1, 1, 1, 1}));  // f_uo cannot produce double replacements
          REQUIRE(
              !f_uo().in({0, 0, 0, 0}));  // f_uo cannot produce null
          replacements

          REQUIRE(f_gg({0, 1, 1, 0})
                      .in({0, 0, 0, 0}));  // f_gg can produce reference when
                                           // acting on singly-excited
          determinant REQUIRE(f_gg({0, 1, 1, 0}) .in({0, 1, 1, 0}));  // f_gg
          can produce singly-excited determinant
                                  // when acting on singly-excited determinant
          REQUIRE(
              !f_uo({0, 1, 1, 0})
                   .in({0, 0, 0, 0}));  // f_uo can't produce reference when
                                        // acting on singl-y-excited determinant
          REQUIRE(f_uo({0, 1, 1, 0})
                      .in({0, 2, 2,
                           0}));  // f_uo can produce doubly-excited determinant
                                  // when acting on singl-y-excited determinant

          //        REQUIRE(!f1(qns_t{2, 2}).in(0));  // can't produce reference
          //        when
          //                                          // acting on
          doubly-excited
           */
        }
        {
          // equal compare
          // using namespace
          // boost::numeric::interval_lib::compare::lexicographic;
          // REQUIRE(f1(qns_t{0, 0}) == qns_t{-1, 1}); // not same as below due
          // to interaction with Catch could do REQUIRE(operator==(f1(qns_t{0,
          // 0}), qns_t{-1, 1})); but equal is shorter
          //        REQUIRE(equal(f1(qns_t{0, 0}), qns_t{-1, 1}));
          //        REQUIRE(equal(f1(qns_t{-1, 1}), qns_t{-2, 2}));
        }
      }
    }  // SECTION("constructor")

    SECTION("to_latex") {
      using qns_t [[maybe_unused]] = mbpt::qns_t;
      using namespace sequant::mbpt;

      auto f = F();
      auto t1 = T(1);
      auto t2 = T_(2);
      auto lambda1 = Λ_(1);
      auto lambda2 = Λ_(2);
      auto r_2_1 = R_(nₚ(1), nₕ(2));
      auto r_1_2 = R_(nₚ(2), nₕ(1));
      auto theta2 = θ(2);
      REQUIRE(to_latex(theta2) == L"{\\hat{\\theta}_{2}}");
      REQUIRE(to_latex(f) == L"{\\hat{f}}");
      REQUIRE(to_latex(t1) == L"{\\hat{t}_{1}}");
      REQUIRE(to_latex(t2) == L"{\\hat{t}_{2}}");
      REQUIRE(to_latex(lambda1) == L"{\\hat{\\lambda}_{1}}");
      REQUIRE(to_latex(lambda2) == L"{\\hat{\\lambda}_{2}}");
      REQUIRE(to_latex(r_2_1) == L"{\\hat{R}_{2,1}}");
      REQUIRE(to_latex(r_1_2) == L"{\\hat{R}_{1,2}}");

      // projectors
      auto P2 = P(2);
      auto P_neg2 = P(-2);
      auto P_2_1 = P(nₚ(1), nₕ(2));
      auto P_neg_2_1 = P(nₚ(-1), nₕ(-2));
      REQUIRE(to_latex(P2) == L"{\\hat{A}_{-2}}");
      REQUIRE(to_latex(P_neg2) == L"{\\hat{A}_{2}}");
      REQUIRE(to_latex(P_2_1) == L"{\\hat{A}_{-2,-1}}");
      REQUIRE(to_latex(P_neg_2_1) == L"{\\hat{A}_{2,1}}");
    }  // SECTION("to_latex")

    SECTION("canonicalize") {
      using qns_t [[maybe_unused]] = mbpt::qns_t;
      using namespace sequant::mbpt;
      auto f = F();
      auto t1 = T_(1);
      auto l1 = Λ_(1);
      auto t2 = T_(2);
      auto l2 = Λ_(2);
      auto h_pt = Hʼ(1, {.order = 1});
      REQUIRE(to_latex(f * t1 * t2) == to_latex(canonicalize(f * t2 * t1)));
      REQUIRE(to_latex(canonicalize(f * t1 * t2)) ==
              to_latex(canonicalize(f * t2 * t1)));
      REQUIRE(to_latex(t1 * t2 * f * t1 * t2) ==
              to_latex(canonicalize(t2 * t1 * f * t2 * t1)));

      REQUIRE(to_latex(ex<Constant>(3) * f * t1 * t2) ==
              to_latex(simplify(ex<Constant>(2) * f * t2 * t1 + f * t1 * t2)));

      REQUIRE(to_latex(simplify(t1 * l1)) != to_latex(simplify(l1 * t1)));
      REQUIRE(to_latex(simplify(t1 * l2)) != to_latex(simplify(l2 * t1)));

      REQUIRE(to_latex(simplify(l2 * t1)) ==
              L"{{\\hat{\\lambda}_{2}}{\\hat{t}_{1}}}");
      REQUIRE(to_latex(simplify(t1 * l2)) ==
              L"{{\\hat{t}_{1}}{\\hat{\\lambda}_{2}}}");

      REQUIRE(to_latex(simplify(t1 + t1)) ==
              to_latex(simplify(ex<Constant>(2) * t1)));

      REQUIRE(to_latex(simplify(t1 + t1 + t2)) ==
              to_latex(simplify(ex<Constant>(2) * t1 + t2)));

      REQUIRE(to_latex(simplify(f + f)) ==
              to_latex(simplify(ex<Constant>(2) * f)));

      REQUIRE(to_latex(simplify(h_pt + h_pt + f)) ==
              to_latex(simplify(ex<Constant>(2) * h_pt + f)));

      auto t = t1 + t2;

      {
        //      std::wcout << "to_latex(simplify(f * t * t)): "
        //                 << to_latex(simplify(f * t * t)) << std::endl;
        REQUIRE_THAT(simplify(f * t * t),
                     EquivalentTo(f * t2 * t2 + ex<Constant>(2) * f * t1 * t2 +
                                  f * t1 * t1));
      }

      {
        //      std::wcout << "to_latex(simplify(f * t * t * t): "
        //                 << to_latex(simplify(f * t * t * t)) << std::endl;
        REQUIRE_THAT(simplify(f * t * t * t),
                     EquivalentTo(f * t1 * t1 * t1 + f * t2 * t2 * t2 +
                                  ex<Constant>(3) * f * t1 * t2 * t2 +
                                  ex<Constant>(3) * f * t1 * t1 * t2));
      }
    }  // SECTION("canonicalize")

    SECTION("adjoint") {
      using qns_t = mbpt::qns_t;
      using op_t = mbpt::Operator<qns_t>;
      using namespace mbpt;
      op_t f = F()->as<op_t>();
      op_t t1 = T_(1)->as<op_t>();
      op_t lambda2 = Λ_(2)->as<op_t>();
      op_t r_1_2 = R_(nₚ(2), nₕ(1))->as<op_t>();

      REQUIRE_NOTHROW(adjoint(f));
      REQUIRE_NOTHROW(adjoint(t1));
      REQUIRE_NOTHROW(adjoint(lambda2));
      REQUIRE_NOTHROW(adjoint(r_1_2));

      REQUIRE(adjoint(f)() == mbpt::general_type_qns(1));
      REQUIRE(adjoint(t1)() == mbpt::deexcitation_type_qns(1));
      REQUIRE(adjoint(lambda2)() == mbpt::excitation_type_qns(2));
      REQUIRE(adjoint(r_1_2)() == L_(nₚ(2), nₕ(1))->as<op_t>()());

      // adjoint(adjoint(Op)) = Op
      REQUIRE(adjoint(adjoint(t1))() == t1());
      REQUIRE(adjoint(adjoint(r_1_2))() == r_1_2());

      // tensor_form()
      REQUIRE((simplify(adjoint(t1).tensor_form())) ==
              (simplify(adjoint(t1.tensor_form()))));

      REQUIRE(simplify(adjoint(lambda2).tensor_form()) ==
              simplify(adjoint(lambda2.tensor_form())));
      REQUIRE(simplify(adjoint(r_1_2).tensor_form()) ==
              simplify(adjoint(r_1_2.tensor_form())));

      // to_latex()
      REQUIRE(to_latex(adjoint(f).as<Expr>()) == L"{\\hat{f⁺}}");
      REQUIRE(to_latex(adjoint(t1).as<Expr>()) == L"{\\hat{t⁺}^{1}}");
      REQUIRE(to_latex(adjoint(lambda2).as<Expr>()) ==
              L"{\\hat{\\lambda⁺}^{2}}");
      REQUIRE(to_latex(adjoint(r_1_2).as<Expr>()) == L"{\\hat{R⁺}^{1,2}}");

      // adjoint(adjoint(op)) == op
      auto t1_adj = adjoint(t1);
      auto r_1_2_adj = adjoint(r_1_2);
      auto lambda2_adj = adjoint(lambda2);
      REQUIRE(to_latex(adjoint(t1_adj).as<Expr>()) == L"{\\hat{t}_{1}}");
      REQUIRE(to_latex(adjoint(r_1_2_adj).as<Expr>()) == L"{\\hat{R}_{1,2}}");
      REQUIRE(to_latex(adjoint(lambda2_adj).as<Expr>()) ==
              L"{\\hat{\\lambda}_{2}}");
    }  // SECTION("adjoint")

    SECTION("screen") {
      using namespace sequant::mbpt;
      auto g_t2_t2 = H_(2) * T_(2) * T_(2);
      REQUIRE(raises_vacuum_to_rank(g_t2_t2, 2));
      REQUIRE(raises_vacuum_up_to_rank(g_t2_t2, 2));

      auto g_t2 = H_(2) * T_(2);
      REQUIRE(raises_vacuum_to_rank(g_t2, 3));

      auto lambda2_f = Λ_(2) * H_(1);
      REQUIRE(lowers_rank_to_vacuum(lambda2_f, 2));

      auto expr1 = P(nₚ(0), nₕ(1)) * H() * R(nₚ(0), nₕ(1));
      auto expr1_tnsr = lower_to_tensor_form(expr1);
      auto vev1_op = op::vac_av(expr1);
      auto vev1_t = tensor::vac_av(expr1_tnsr);  // no operator level screening
      REQUIRE(to_latex(vev1_op) == to_latex(vev1_t));

      auto expr2 = P(nₚ(2), nₕ(1)) * H() * R(nₚ(1), nₕ(0));
      auto expr2_tnsr = lower_to_tensor_form(expr2);
      auto vev2_op = op::vac_av(expr2);
      auto vev2_t = tensor::vac_av(expr2_tnsr);  // no operator level screening
      REQUIRE(to_latex(vev2_op) == to_latex(vev2_t));

      // Test screen_vac_av
      auto hbar = mbpt::lst(H(), T_(2), 4);  // CCD Hbar
      auto screened_hbar = screen_vac_av(hbar);
      auto expected = H_(2) * T_(2);
      REQUIRE(simplify(screened_hbar - expected) == ex<Constant>(0));

      auto expr3 = P(2) * hbar * R_(nₚ(2), nₕ(2));
      auto screened_expr3 = screen_vac_av(expr3);
      auto expected3 =
          op::P(2) * (H_(2) * T_(2) + H_(2) + H_(1)) * R_(nₚ(2), nₕ(2));
      REQUIRE(simplify(screened_expr3 - expected3) == ex<Constant>(0));

      auto expr4 = P(nₚ(2), nₕ(1)) * hbar * R(nₚ(1), nₕ(0));
      auto screened_expr4 = screen_vac_av(expr4);
      auto expected4 = P(nₚ(2), nₕ(1)) *
                       (H_(1) + H_(1) * T_(2) + H_(2) * T_(2) + H_(2)) *
                       R(nₚ(1), nₕ(0));
      REQUIRE(simplify(screened_expr4 - expected4) == ex<Constant>(0));

    }  // SECTION("screen")

    SECTION("lst") {
      using namespace sequant::mbpt;

      auto commutator = [](const ExprPtr& A, const ExprPtr& B) {
        return A * B - B * A;
      };

      // non-unitary, rank 3
      auto expr1 =
          lst(H(), T_(2), 3, {.unitary = false, .use_commutators = false});
      auto expected1 =
          H() * (ex<Constant>(1) + T_(2) +
                 ex<Constant>(rational{1, 2}) * T_(2) * T_(2) +
                 ex<Constant>(rational{1, 6}) * T_(2) * T_(2) * T_(2));
      REQUIRE(simplify(expr1 - expected1) == ex<Constant>(0));

      auto expr2 =
          lst(H(), T_(2), 3, {.unitary = false, .use_commutators = true});
      auto expected2 =
          H() + commutator(H(), T_(2)) +
          ex<Constant>(rational{1, 2}) *
              commutator(commutator(H(), T_(2)), T_(2)) +
          ex<Constant>(rational{1, 6}) *
              commutator(commutator(commutator(H(), T_(2)), T_(2)), T_(2));
      REQUIRE(simplify(expr2 - expected2) == ex<Constant>(0));

      // unitary, rank 2
      using sequant::adjoint;
      auto expr3 =
          lst(H(), T_(2), 2, {.unitary = true, .use_commutators = false});
      auto expected3 =
          H() + H() * T_(2) + adjoint(T_(2)) * H() +
          adjoint(T_(2)) * H() * T_(2) +
          H() * ex<Constant>(rational{1, 2}) * T_(2) * T_(2) +
          ex<Constant>(rational{1, 2}) * adjoint(T_(2)) * adjoint(T_(2)) * H();
      REQUIRE(simplify(expr3 - expected3) == ex<Constant>(0));

      auto expr4 =
          lst(H(), T_(2), 2, {.unitary = true, .use_commutators = true});
      auto generator = commutator(H(), T_(2)) - commutator(H(), adjoint(T_(2)));
      auto expected4 = H() + generator +
                       ex<Constant>(rational{1, 2}) *
                           (commutator(generator, T_(2)) -
                            commutator(generator, adjoint(T_(2))));
      REQUIRE(simplify(expr4 - expected4) == ex<Constant>(0));
    }  // SECTION("lst")

    SECTION("predefined") {
      // P.S. ref outputs produced with complete canonicalization
      auto ctx = get_default_context();
      ctx.set(CanonicalizeOptions{.method = CanonicalizationMethod::Complete});
      auto _ = set_scoped_default_context(ctx);
      using namespace sequant::mbpt;

      auto theta1 = θ(1)->as<op_t>();
      // std::wcout << "theta1: " << to_latex(simplify(theta1.tensor_form()));
      REQUIRE(to_latex(simplify(theta1.tensor_form())) ==
              L"{{\\theta^{{p_2}}_{{p_1}}}{\\tilde{a}^{{p_1}}_{{p_2}}}}");

      auto R_2 = R_(2)->as<op_t>();
      //    std::wcout << "R_2: " << to_latex(simplify(R_2.tensor_form())) <<
      //    std::endl;
      REQUIRE(
          to_latex(simplify(R_2.tensor_form())) ==
          L"{{{\\frac{1}{4}}}{\\bar{R}^{{i_1}{i_2}}_{{a_1}{a_2}}}{\\tilde{a}^"
          L"{{a_1}{"
          L"a_2}}_{{i_1}{i_2}}}}");

      auto L_3 = L_(3)->as<op_t>();
      //    std::wcout << "L_3: " << to_latex(simplify(L_3.tensor_form())) <<
      //    std::endl;
      REQUIRE(
          to_latex(simplify(L_3.tensor_form())) ==
          L"{{{\\frac{1}{36}}}{\\bar{L}^{{a_1}{a_2}{a_3}}_{{i_1}{i_2}{i_3}}}{"
          L"\\tilde{a}^{{i_1}{i_2}{i_3}}_{{a_1}{a_2}{a_3}}}}");

      auto R_2_3 = R_(nₚ(3), nₕ(2))->as<op_t>();
      //    std::wcout << "R_2_3: " << to_latex(simplify(R_2_3.tensor_form()))
      //    << std::endl;
      REQUIRE(to_latex(simplify(R_2_3.tensor_form())) ==
              L"{{{\\frac{1}{12}}}{\\bar{R}^{{i_1}{i_2}}_{{a_1}{a_2}{a_3}}}{"
              L"\\tilde{a}^{"
              L"{a_1}{a_2}{a_3}}_{\\textvisiblespace\\,{i_1}{i_2}}}}");

      auto L_1_2 = L_(nₚ(1), nₕ(2))->as<op_t>();
      // std::wcout << "L_(1,2): " << to_latex(simplify(L_1_2.tensor_form())) <<
      // std::endl;
      REQUIRE(
          to_latex(simplify(L_1_2.tensor_form())) ==
          L"{{{\\frac{1}{2}}}{\\bar{L}^{{a_1}}_{{i_1}{i_2}}}{\\tilde{a}^{{i_"
          L"1}{i_2}}"
          L"_{\\textvisiblespace\\,{a_1}}}}");

      auto A_2_1 = A(nₚ(2), nₕ(1))->as<op_t>();
      //    std::wcout << "A_2_1: " << to_latex(simplify(A_2_1.tensor_form()))
      //               << std::endl;
      REQUIRE(
          to_latex(simplify(A_2_1.tensor_form())) ==
          L"{{{\\frac{1}{2}}}{\\bar{A}^{{i_1}}_{{a_1}{a_2}}}{\\tilde{a}^{{a_"
          L"1}{a_2}}"
          L"_{\\textvisiblespace\\,{i_1}}}}");

      auto P_0_1 = P(nₚ(0), nₕ(1))->as<op_t>();
      //    std::wcout << "P_0_1: " << to_latex(simplify(P_0_1.tensor_form()))
      //               << std::endl;
      REQUIRE(to_latex(simplify(P_0_1.tensor_form())) ==
              L"{{A^{}_{{i_1}}}{\\tilde{a}^{{i_1}}}}");

      auto P_2_1 = P(nₚ(2), nₕ(1))->as<op_t>();
      //    std::wcout << "P_2_1: " << to_latex(simplify(P_2_1.tensor_form()))
      //               << std::endl;
      REQUIRE(to_latex(simplify(P_2_1.tensor_form())) ==
              L"{{{\\frac{1}{2}}}{\\bar{A}^{{a_1}{a_2}}_{{i_1}}}{\\tilde{a}^{"
              L"\\textvisiblespace\\,{i_1}}_{{a_1}{a_2}}}}");

      auto P_2_3 = P(nₚ(2), nₕ(3))->as<op_t>();
      //    std::wcout << "P_2_3: " << to_latex(simplify(P_3_2.tensor_form()))
      //               << std::endl;
      REQUIRE(to_latex(simplify(P_2_3.tensor_form())) ==
              L"{{{\\frac{1}{12}}}{\\bar{A}^{{a_1}{a_2}}_{{i_1}{i_2}{i_3}}}{"
              L"\\tilde{a}^{"
              L"{i_1}{i_2}{i_3}}_{\\textvisiblespace\\,{a_1}{a_2}}}}");

      auto R33 = R(3);
      lower_to_tensor_form(R33);
      simplify(R33);
      //    std::wcout << "R33: " << to_latex(R33) << std::endl;
      REQUIRE(to_latex(R33) ==
              L"{ "
              L"\\bigl({{R^{{i_1}}_{{a_1}}}{\\tilde{a}^{{a_1}}_{{i_1}}}} + "
              L"{{{\\frac{1}{36}}}{\\bar{R}^{{i_1}{i_2}{i_3}}_{{a_1}{a_2}{a_3}}"
              L"}{\\tilde{a}^{{a_1}{a_2}{a_3}}_{{i_1}{i_2}{i_3}}}} + "
              L"{{{\\frac{1}{4}}}{\\bar{R}^{{i_1}{i_2}}_{{a_1}{a_2}}}{\\tilde{"
              L"a}^{{a_1}{a_2}}_{{i_1}{i_2}}}}\\bigr) }");

      auto R12 = R(nₚ(2), nₕ(1));
      lower_to_tensor_form(R12);
      simplify(R12);
      //    std::wcout << "R12: " << to_latex(R12) << std::endl;
      REQUIRE(to_latex(R12) ==
              L"{ \\bigl({{R^{}_{{a_1}}}{\\tilde{a}^{{a_1}}}} + "
              L"{{{\\frac{1}{2}}}{\\bar{R}^{{i_1}}_{{a_1}{a_2}}}{\\tilde{a}^{{"
              L"a_1}{a_"
              L"2}}_{\\textvisiblespace\\,{i_1}}}}\\bigr) }");

      auto R21 = R(nₚ(1), nₕ(2));
      lower_to_tensor_form(R21);
      simplify(R21);
      //    std::wcout << "R21: " << to_latex(R21) << std::endl;
      REQUIRE(to_latex(R21) ==
              L"{ "
              L"\\bigl({{{\\frac{1}{2}}}{\\bar{R}^{{i_1}{i_2}}_{{a_1}}}{"
              L"\\tilde{a}^{"
              L"\\textvisiblespace\\,{a_1}}_{{i_1}{i_2}}}} + "
              L"{{R^{{i_1}}_{}}{\\tilde{a}_{{i_1}}}}\\bigr) }");

      auto L23 = L(nₚ(2), nₕ(3));
      lower_to_tensor_form(L23);
      simplify(L23);
      // std::wcout << "L23: " << to_latex(L23) << std::endl;
      REQUIRE(to_latex(L23) ==
              L"{ "
              L"\\bigl({{L^{}_{{i_1}}}{\\tilde{a}^{{i_1}}}} + "
              L"{{{\\frac{1}{2}}}{\\bar{L}^{{a_1}}_{{i_1}{i_2}}}{\\tilde{a}^{{"
              L"i_1}{i_2}}_{\\textvisiblespace\\,{a_1}}}} + "
              L"{{{\\frac{1}{12}}}{\\bar{L}^{{a_1}{a_2}}_{{i_1}{i_2}{i_"
              L"3}}}{\\tilde{a}^{{i_1}{i_2}{i_3}}_{\\textvisiblespace\\,{a_1}{"
              L"a_2}}}}\\bigr) }");
    }

    SECTION("batching") {
      // update context to use batching index
      auto isr = sequant::mbpt::make_legacy_spaces();
      mbpt::add_batching_spaces(isr);
      auto ctx_resetter =
          set_scoped_default_context({.index_space_registry_shared_ptr = isr,
                                      .vacuum = Vacuum::SingleProduct});
      REQUIRE_NOTHROW(
          get_default_context().index_space_registry()->retrieve(L"z"));

      using namespace mbpt;
      REQUIRE_NOTHROW(op::Hʼ(1, {.nbatch = 1}));
      REQUIRE_NOTHROW(op::Hʼ(2, {.nbatch = 2}));
      REQUIRE_NOTHROW(op::Λʼ(3, {.batch_ordinals = {3, 4, 5}, .skip1 = true}));
      REQUIRE_NOTHROW(op::Tʼ(1, {.nbatch = 20}));

      // invalid usages
#if SEQUANT_ASSERT_BEHAVIOR == SEQUANT_ASSERT_THROW
      // cannot set both nbatch and batch_ordinals
      REQUIRE_THROWS_AS(op::Hʼ(2, {.nbatch = 2, .batch_ordinals = {1, 2}}),
                        sequant::Exception);
      // all ordinals must be unique
      REQUIRE_THROWS_AS(op::Hʼ(2, {.batch_ordinals = {1, 2, 2}}),
                        sequant::Exception);
      // ordinals must be sorted
      REQUIRE_THROWS_AS(op::Hʼ(1, {.batch_ordinals = {3, 2}}),
                        sequant::Exception);
#endif

      // operations
      auto h0 = op::Hʼ(1);
      REQUIRE(to_latex(h0) == L"{\\hat{h¹}}");

      auto h1 = op::Hʼ(1, {.nbatch = 1});
      auto h1_2 = op::Hʼ(1, {.batch_ordinals = {1, 2}});
      auto pt1 = op::Tʼ(2, {.batch_ordinals = {1}});

      auto sum0 = h0 + h1;
      simplify(sum0);
      REQUIRE(to_latex(sum0) ==
              L"{ \\bigl({\\hat{h¹}}{[{z}_{1}]} + {\\hat{h¹}}\\bigr) }");

      auto sum1 = h1 + h1;
      simplify(sum1);
      REQUIRE(to_latex(sum1) == L"{{{2}}{\\hat{h¹}}{[{z}_{1}]}}");
      auto sum2 = h1 + pt1;
      simplify(sum2);
      // std::wcout << "sum2:  " << to_latex(sum2) << std::endl;
      REQUIRE(to_latex(sum2) ==
              L"{ \\bigl({\\hat{h¹}}{[{z}_{1}]} + {\\hat{t¹}_{2}}{[{z}_{1}]} + "
              L"{\\hat{t¹}_{1}}{[{z}_{1}]}\\bigr) }");

      auto sum3 = h1 + h1_2;
      simplify(sum3);
      // std::wcout << "sum3:  " << to_latex(sum3) << std::endl;
      REQUIRE(to_latex(sum3) ==
              L"{ \\bigl({\\hat{h¹}}{[{z}_{1},{z}_{2}]} + "
              L"{\\hat{h¹}}{[{z}_{1}]}\\bigr) }");

      auto pdt1 = h1 * h1_2;
      simplify(pdt1);
      // std::wcout << "pdt1: " << to_latex(pdt1) << std::endl;
      REQUIRE(to_latex(pdt1) ==
              L"{{\\hat{h¹}}{[{z}_{1}]}{\\hat{h¹}}{[{z}_{1},{z}_{2}]}}");

      auto pdt2 = h1 * pt1;
      simplify(pdt2);
      // std::wcout << "pdt1: " << to_latex(pdt2) << std::endl;
      REQUIRE(to_latex(pdt2) ==
              L"{ \\bigl({{\\hat{h¹}}{[{z}_{1}]}{\\hat{t¹}_{1}}{[{z}_{1}]}} + "
              L"{{\\hat{h¹}}{[{z}_{1}]}{\\hat{t¹}_{2}}{[{z}_{1}]}}\\bigr) }");

      // lowering to tensor form
      auto sum1_t = simplify(lower_to_tensor_form(sum1));
      // std::wcout << "sum1_t: " << to_latex(simplify(sum1_t)) << std::endl;
      REQUIRE_THAT(sum1_t, EquivalentTo(L"2 * h¹{κ1;κ2;z1}:A-C-S * ã{κ2;κ1}"));

      auto sum2_t = simplify(lower_to_tensor_form(sum2));
      // std::wcout << "sum2_t: " << to_latex(sum2_t) << std::endl;
      REQUIRE_THAT(
          sum2_t,
          EquivalentTo(L"h¹{κ2;κ1;z1}:A-C-S * ã{κ1;κ2} + t¹{a1;i1;z1}:A-C-S * "
                       L"ã{i1;a1} + "
                       "(1/4) * t¹{a1,a2;i1,i2;z1}:A-C-S * ã{i1,i2;a1,a2}"));

      auto expr3_t = simplify(lower_to_tensor_form(sum2 * h1_2));
      // std::wcout << "expr3_t: " << to_latex(expr3_t) << std::endl;
      REQUIRE_THAT(
          expr3_t,
          EquivalentTo(
              L"1/4 ã{i1,i2;a1,a2} * t¹{a1,a2;i1,i2;z1}:A-C-S * "
              L"h¹{κ2;κ1;z1,z2}:A-C-S * ã{κ1;κ2} + h¹{κ2;κ1;z1,z2}:A-C-S * "
              L"ã{i1;a1} * ã{κ1;κ2} * t¹{a1;i1;z1}:A-C-S + h¹{κ4;κ3;z1}:A-C-S "
              L"* h¹{κ2;κ1;z1,z2}:A-C-S * ã{κ3;κ4} * ã{κ1;κ2}"));
    }  // SECTION("batching")
  }

  SECTION("wick") {
    using namespace sequant;
    using namespace sequant::mbpt;
    namespace o = sequant::mbpt::op;
    namespace t = sequant::mbpt::tensor;
    TensorCanonicalizer::register_instance(
        std::make_shared<DefaultTensorCanonicalizer>());

    SECTION("SRSO"){
        // H**T12**T12 -> R2
        SECTION("wick(H**T12**T12 -> R2)"){
            auto result = t::vac_av(t::A(nₚ(-2)) * t::H(2) * t::T(2) * t::T(2),
                                    {{1, 2}, {1, 3}});

    //      std::wcout << "H*T12*T12 -> R2 = " << to_latex_align(result, 20)
    //                 << std::endl;
    REQUIRE(result->size() == 15);

    {
      // check against op
      auto result_op = o::vac_av(o::P(nₚ(2)) * o::H() * o ::T(2) * o::T(2));
      REQUIRE(result_op->size() == result->size());  // as compact as result ..
      REQUIRE(simplify(result_op - result) ==
              ex<Constant>(0));  // .. and equivalent to it
    }
  }

  // H2**T3**T3 -> R4
  SECTION("wick(H2**T3**T3 -> R4)") {
    auto result = t::vac_av(t::A(nₚ(-4)) * t::H_(2) * t::T_(3) * t::T_(3),
                            {{1, 2}, {1, 3}});

    // std::wcout << "H2**T3**T3 -> R4 = " << to_latex_align(result, 20)
    //            << std::endl;
    REQUIRE(result->size() == 4);
  }

#ifndef SEQUANT_SKIP_LONG_TESTS
  // the longest term in CCSDTQP
  // H2**T2**T2**T3 -> R5
  {
    ExprPtr ref_result;
    SECTION("wick(H2**T2**T2**T3 -> R5)") {
      ref_result =
          t::vac_av(t::A(-5) * t::H_(2) * t::T_(2) * t::T_(2) * t::T_(3),
                    {{1, 2}, {1, 3}, {1, 4}});
      REQUIRE(ref_result->size() == 7);
    }
  }
#endif  // !defined(SEQUANT_SKIP_LONG_TESTS)
}  // SECTION ("SRSO")

SECTION("SRSO Fock") {
  // <2p1h|H2|1p> ->
  SECTION("wick(<2p1h|H2|1p>)") {
    auto input = t::L_(nₚ(2), nₕ(1)) * t::H_(2) * t::R_(nₚ(1), nₕ(0));
    auto result = t::vac_av(input);

    REQUIRE(result->is<Product>());  // product ...
    REQUIRE(result->size() == 3);    // ... of 3 factors
  }

  // <2p1h|H2|2p1h(c)> ->
  SECTION("wick(<2p1h|H2|2p1h(c)>)") {
    auto input = t::L_(nₚ(2), nₕ(1)) * t::H() * t::R_(nₚ(2), nₕ(1));
    auto result = t::vac_av(input);

    // std::wcout << "<2p1h|H|2p1h(c)> = " << to_latex(result)
    //            << std::endl;
    REQUIRE(result->is<Sum>());    // sub ...
    REQUIRE(result->size() == 4);  // ... of 4 factors
  }
}  // SECTION("SRSO Fock")

SECTION("SRSO-PNO") {
  using sequant::mbpt::Context;
  auto mbpt_ctx =
      sequant::mbpt::set_scoped_default_mbpt_context(Context(mbpt::CSV::Yes));

  // H2**T2**T2 -> R2
  SECTION("wick(H2**T2**T2 -> R2)") {
    auto result = t::vac_av(t::A(nₚ(-2)) * t::H_(2) * t::T_(2) * t::T_(2),
                            {{1, 2}, {1, 3}});

    REQUIRE(result->size() == 4);
  }
}  // SECTION("SRSO-PNO")

SECTION("SRSF") {
  auto ctx = get_default_context();
  ctx.set(SPBasis::Spinfree);
  auto ctx_resetter = set_scoped_default_context(ctx);

  // H2 -> R2
  SECTION("wick(H2 -> R2)") {
    auto result = t::vac_av(t::S(-2) * t::H_(2));

    {
      // std::wcout << "H2 -> R2 = " << to_latex_align(result, 0, 1)
      //            << std::endl;
      REQUIRE(result->is<Sum>());
      REQUIRE(result->size() == 2);
    }
  }

  // H2**T2 -> R2
  SECTION("wick(H2**T2 -> R2)") {
    auto result = t::vac_av(t::S(-2) * t::H_(2) * t::T_(2), {{1, 2}});

    {
      // std::wcout << "H2**T2 -> R2 = " << to_latex_align(result, 0, 1)
      //            << std::endl;
      REQUIRE(result->is<Sum>());
      REQUIRE(result->size() == 12);
    }
  }
}  // SECTION("SRSF")

SECTION("MRSO") {
  auto ctx = get_default_context();
  ctx.set(mbpt::make_mr_spaces());
  auto ctx_resetter = set_scoped_default_context(ctx);

  SECTION("wick(H2**T2 -> 0)") {
    {
      auto result = t::ref_av(t::H_(2) * t::T_(2), {{0, 1}});

      auto result_wo_top = t::ref_av(t::H_(2) * t::T_(2), {{0, 1}},
                                     /* use_topology = */ false);
      REQUIRE(simplify(result - result_wo_top) == ex<Constant>(0));
    }

    // now compute using physical vacuum
    {
      auto ctx = get_default_context();
      ctx.set(mbpt::make_mr_spaces());
      ctx.set(Vacuum::Physical);
      auto ctx_resetter = set_scoped_default_context(ctx);
      auto result_phys = t::ref_av(t::H_(2) * t::T_(2), {{0, 1}});
    }
  }

  // H2 ** T2 ** T2 -> 0
  SECTION("wick(H2**T2**T2 -> 0)") {
    // first without use of topology
    auto result = t::ref_av(t::H_(2) * t::T_(2) * t::T_(2), {{0, 1}},
                            /* use_topology = */ false);
    // now with topology use
    auto result_top = t::ref_av(t::H_(2) * t::T_(2) * t ::T_(2), {{0, 1}},
                                /* use_topology = */ true);

    REQUIRE(simplify(result - result_top) == ex<Constant>(0));
  }

#if 0
    // H**T12 -> R2
    SECTION("wick(H**T2 -> R2)") {
      auto result = t::ref_av(t::A(-2) * t::H() * t::T_(2), {{1, 2}});

      {
        std::wcout << "H*T2 -> R2 = " << to_latex_align(result, 0, 1)
                   << std::endl;
      }
    }
#endif
}  // SECTION("MRSO")

SECTION("MRSF") {
  // now compute using (closed) Fermi vacuum + spinfree basis
  auto ctx = get_default_context();
  ctx.set(mbpt::make_mr_spaces());
  ctx.set(SPBasis::Spinfree);
  auto ctx_resetter = set_scoped_default_context(ctx);

  SECTION("wick(H2**T2 -> 0)") {
    auto result = t::ref_av(t::H_(2) * t::T_(2), {{0, 1}});

    {
      // make sure get same result without use of topology
      auto result_wo_top = t::ref_av(t::H_(2) * t::T_(2), {{0, 1}},
                                     /* use_topology = */ false);

      REQUIRE(simplify(result - result_wo_top) == ex<Constant>(0));
    }

    {
      // make sure get same result using operators
      auto result_op = o::ref_av(o::H_(2) * o::T_(2));

      REQUIRE(result_op->size() == result->size());
      REQUIRE(simplify(result - result_op) == ex<Constant>(0));
    }
  }
}  // SECTION("MRSF")
}

SECTION("rules") {
  using namespace sequant;

  SECTION("density-fit") {
    const std::vector<std::wstring> inputs = {
        L"t{a1,a2;i1,i2} t{a3;i3}",
        L"t{a1,a2;i1,i2} g{i1,i2;a1,a2}",
        L"t{a1,a2;i1,i2} g{i1,i2;a1,a2}:A",
    };
    const std::vector<std::wstring> expected = {
        L"t{a1,a2;i1,i2} t{a3;i3}",
        L"t{a1,a2;i1,i2} B{i1;a1;x_1} B{i2;a2;x_1}",
        L"t{a1,a2;i1,i2} (B{i1;a1;x_1} B{i2;a2;x_1} "
        "- B{i2;a1;x_1} B{i1;a2;x_1})",
    };

    REQUIRE(inputs.size() == expected.size());

    for (std::size_t i = 0; i < inputs.size(); ++i) {
      CAPTURE(inputs.at(i));

      ExprPtr input_expr = parse_expr(inputs.at(i));

      const IndexSpace aux_space =
          get_default_context().index_space_registry()->retrieve(L"x");

      ExprPtr actual = mbpt::density_fit(input_expr, aux_space, L"g", L"B");

      REQUIRE_THAT(actual, EquivalentTo(expected.at(i)));
    }
  }

  SECTION("tensor-hypercontract") {
    const std::vector<std::wstring> inputs = {
        L"t{a1,a2;i1,i2} t{a3;i3}",
        L"t{a1,a2;i1,i2} g{i1,i2;a1,a2}",
        L"t{a1,a2;i1,i2} g{i1,i2;a1,a2}:A",
    };
    const std::vector<std::wstring> expected = {
        L"t{a1,a2;i1,i2} t{a3;i3}",
        L"t{a1,a2;i1,i2} B{i1;;x_1} B{;a1;x_1} C{;;x_1,x_2} B{i2;;x_2} "
        L"B{;a2;x_2}",
        L"t{a1,a2;i1,i2} (B{i1;;x_1} B{;a1;x_1} C{;;x_1,x_2} B{i2;;x_2} "
        L"B{;a2;x_2}"
        " - B{i2;;x_1} B{;a1;x_1} C{;;x_1,x_2} B{i1;;x_2} B{;a2;x_2})",
    };

    REQUIRE(inputs.size() == expected.size());

    for (std::size_t i = 0; i < inputs.size(); ++i) {
      CAPTURE(inputs.at(i));

      ExprPtr input_expr = parse_expr(inputs.at(i));

      const IndexSpace aux_space =
          get_default_context().index_space_registry()->retrieve(L"x");

      ExprPtr actual =
          mbpt::tensor_hypercontract(input_expr, aux_space, L"g", L"B", L"C");

      REQUIRE_THAT(actual, EquivalentTo(expected.at(i)));
    }
  }
}  // SECTION("rules")

SECTION("manuscript-examples") {
  using namespace sequant::mbpt;

  /// When order == 0, returns sim. transformed of zeroth order Hamiltonian.
  /// When order > 0, returns sim. transformed perturbation operator or given
  /// order.
  auto H̅ = [](size_t order = 0) {
    auto hbar0 = lst(op::H(), T(2), 4);
    if (order == 0) return hbar0;
    // only one-body perturbation operator
    auto hbar_pt = lst(op::Hʼ(/*rank*/ 1, {.order = order}), T(2), 2);
    return hbar_pt;
  };

  SECTION("CCD Term") {
    auto expr = ref_av(P(2) * H() * T_(2) * T_(2));
    REQUIRE(expr.size() == 4);
  }

  SECTION("CC LR Function") {
    const int N = 2;  // CC rank

    auto θ̅ = lst(θ(1), T(N), 2);
    auto expr = (1 + Λ(N)) * θ̅ * Tʼ(N) + Λʼ(N) * θ̅;
    auto result = ref_av(expr, {{L"θ", L"t"}, {L"θ", L"t¹"}});
    REQUIRE(result.size() == 21);
  }

  SECTION("EOM-CC Equations") {
    // connectivity info for right and left amplitude equations
    const auto r_connect =
        concat(default_op_connections(),
               OpConnections<OpType>{{OpType::h, OpType::R},
                                     {OpType::f, OpType::R},
                                     {OpType::g, OpType::R}});
    const auto l_connect =
        concat(default_op_connections(),
               OpConnections<OpType>{{OpType::h, OpType::A},
                                     {OpType::f, OpType::A},
                                     {OpType::g, OpType::A}});

    // Note: the correctness of the equations are verified in test_mbpt_cc.cpp
    // here, we just make sure the examples run without any issues

    // EE
    REQUIRE_NOTHROW(ref_av(P(2) * H̅() * R(2), r_connect));
    REQUIRE_NOTHROW(ref_av(L(2) * H̅() * P(-2), l_connect));
    // EA
    REQUIRE_NOTHROW(P(nₚ(2), nₕ(1)) * H̅() * R(nₚ(2), nₕ(1)), r_connect);
    REQUIRE_NOTHROW(L(nₚ(2), nₕ(1)) * H̅() * P(nₚ(-2), nₕ(-1)), l_connect);
    // IP
    REQUIRE_NOTHROW(P(nₚ(1), nₕ(2)) * H̅() * R(nₚ(1), nₕ(2)), r_connect);
    REQUIRE_NOTHROW(L(nₚ(1), nₕ(2)) * H̅() * P(nₚ(-1), nₕ(-2)), l_connect);
  }

  SECTION("CC Perturbed Amplitudes") {
    // connectivity info for perturbed t and λ amplitude equations
    const auto t_connect =
        concat(default_op_connections(),
               OpConnections<OpType>{{OpType::h, OpType::t_1},
                                     {OpType::f, OpType::t_1},
                                     {OpType::g, OpType::t_1},
                                     {OpType::h_1, OpType::t}});

    const auto l_connect =
        concat(default_op_connections(),
               OpConnections<OpType>{{OpType::h, OpType::t_1},
                                     {OpType::f, OpType::t_1},
                                     {OpType::g, OpType::t_1},
                                     {OpType::h_1, OpType::t},
                                     {OpType::h, OpType::A},
                                     {OpType::f, OpType::A},
                                     {OpType::g, OpType::A},
                                     {OpType::h_1, OpType::A}});

    // perturbed t amplitudes (Eq 18 in SQ Manuscript #2)
    auto t = ref_av(P(2) * (H̅(1) + H̅() * Tʼ(2) - "ω" * Tʼ(2)), t_connect);
    REQUIRE(t.size() == 58);

    // perturbed λ amplitudes (Eq 19 in SQ Manuscript #2)
    auto λ = ref_av(
        ((1 + Λ(2)) * (H̅(1) + H̅() * Tʼ(2)) + Λʼ(2) * H̅() + "ω" * Λʼ(2)) * P(-2),
        l_connect);
    REQUIRE(λ.size() == 63);
  }
}  // SECTION("manuscript examples")
}
