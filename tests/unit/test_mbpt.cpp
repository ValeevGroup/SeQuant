//
// Created by Eduard Valeyev on 2019-02-19.
//

#include "SeQuant/core/op.hpp"
#include "SeQuant/core/tensor.hpp"
#include "SeQuant/core/timer.hpp"
#include "SeQuant/domain/mbpt/context.hpp"
#include "SeQuant/domain/mbpt/mr.hpp"
#include "SeQuant/domain/mbpt/op.hpp"
#include "SeQuant/domain/mbpt/spin.hpp"
#include "SeQuant/domain/mbpt/sr.hpp"

#include "catch.hpp"
#include "test_config.hpp"

TEST_CASE("NBodyOp", "[mbpt]") {
  using namespace sequant;

  SECTION("constructor") {
    // tests 1-space quantum number case
    {
      using namespace sequant::mbpt;

      op_t f1([]() -> std::wstring_view { return L"f"; },
              []() -> ExprPtr {
                return ex<Tensor>(L"f", WstrList{L"p_1"}, WstrList{L"p_2"}) *
                       ex<FNOperator>(WstrList{L"p_1"}, WstrList{L"p_2"});
              },
              [](qns_t& qns) {
                qns += qns_t{1, 1};
              });

      REQUIRE(f1.label() == L"f");

      {  // exact compare
        using namespace boost::numeric::interval_lib::compare::possible;
        REQUIRE(operator==(f1(), qns_t{1, 1}));  // produces single replacement
        REQUIRE(operator!=(f1(),
                           qns_t{2, 2}));  // cannot produce double replacement
        REQUIRE(operator==(f1(qns_t{5, 0}), qns_t{{5, 6}, {0, 1}}));
      }
    }

    // tests 2-space quantum number case
    {
      using namespace sequant::mbpt::sr;

      // this is fock operator in terms of general spaces
      op_t f_gg([]() -> std::wstring_view { return L"f"; },
                []() -> ExprPtr {
                  return ex<Tensor>(L"f", WstrList{L"p_1"}, WstrList{L"p_2"}) *
                         ex<FNOperator>(WstrList{L"p_1"}, WstrList{L"p_2"});
                },
                [](qns_t& qns) {
                  qns += qns_t{{0, 1}, {0, 1}, {0, 1}, {0, 1}};
                });
      // excitation part of the Fock operator
      op_t f_uo([]() -> std::wstring_view { return L"f"; },
                []() -> ExprPtr {
                  return ex<Tensor>(L"f", WstrList{L"a_2"}, WstrList{L"i_2"}) *
                         ex<FNOperator>(WstrList{L"a_1"}, WstrList{L"i_2"});
                },
                [](qns_t& qns) {
                  qns += qns_t{0, 1, 1, 0};
                });

      REQUIRE(f_gg.label() == L"f");
      REQUIRE(f_uo.label() == L"f");

      {  // comparison

        // exact
        REQUIRE(
            (f_uo() == qns_t{0, 1, 1, 0}));  // f_uo produces single excitations
        REQUIRE((f_gg() !=
                 qns_t{0, 1, 1,
                       0}));  // f_gg does not produce just single excitations
        REQUIRE((f_gg() !=
                 qns_t{0, 1, 1, 0}));  // f_gg cannot produce double excitations
        REQUIRE(
            f_gg().in({0, 1, 1, 0}));  // f_gg can produce single excitations
        REQUIRE(f_gg().in(
            {1, 0, 0, 1}));  // f_gg can also produce single de-excitations
        REQUIRE(f_gg().in(
            {1, 1, 0, 0}));  // f_gg can produce replacements withing occupieds
        REQUIRE(f_gg().in(
            {0, 0, 1, 1}));  // f_gg can produce replacements withing virtuals
        REQUIRE(f_gg().in(
            {1, 1, 1, 1}));  // f_gg cannot produce this double replacements,
                             // but this returns true TODO introduce constraints
                             // on the total number of creators/annihilators,
                             // the interval logic does not constrain it
        REQUIRE(f_gg().in(
            {0, 0, 0, 0}));  // f_gg cannot produce a null replacement, but this
                             // returns true TODO introduce constraints on the
                             // total number of creators/annihilators, the
                             // interval logic does not constrain it

        REQUIRE(
            f_uo().in({0, 1, 1, 0}));  // f_uo can produce single excitations
        REQUIRE(!f_uo().in(
            {1, 0, 0, 1}));  // f_uo cannot produce single de-excitations
        REQUIRE(!f_uo().in(
            {1, 1, 0, 0}));  // f_uo can produce replacements withing occupieds
        REQUIRE(!f_uo().in(
            {0, 0, 1, 1}));  // f_uo can produce replacements withing virtuals
        REQUIRE(!f_uo().in(
            {1, 1, 1, 1}));  // f_uo cannot produce double replacements
        REQUIRE(
            !f_uo().in({0, 0, 0, 0}));  // f_uo cannot produce null replacements

        REQUIRE(f_gg({0, 1, 1, 0})
                    .in({0, 0, 0, 0}));  // f_gg can produce reference when
                                         // acting on singly-excited determinant
        REQUIRE(f_gg({0, 1, 1, 0})
                    .in({0, 1, 1,
                         0}));  // f_gg can produce singly-excited determinant
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
        //                                          // acting on doubly-excited
      }
      {  // equal compare
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
    using qns_t = mbpt::sr::qns_t;
    using op_t = mbpt::Operator<qns_t>;
    auto f = ex<op_t>([]() -> std::wstring_view { return L"f"; },
                      []() -> ExprPtr {
                        using namespace sequant::mbpt::sr;
                        return F();
                      },
                      [](qns_t& qns) {
                        qns += qns_t{{0, 1}, {0, 1}, {0, 1}, {0, 1}};
                      });
    auto t1 = ex<op_t>([]() -> std::wstring_view { return L"t"; },
                       []() -> ExprPtr {
                         using namespace sequant::mbpt::sr;
                         return T_(1);
                       },
                       [](qns_t& qns) {
                         qns += qns_t{0, 1, 1, 0};
                       });
    auto t2 = ex<op_t>([]() -> std::wstring_view { return L"t"; },
                       []() -> ExprPtr {
                         using namespace sequant::mbpt::sr;
                         return T_(2);
                       },
                       [](qns_t& qns) {
                         qns += qns_t{0, 2, 2, 0};
                       });
    auto lambda1 = ex<op_t>([]() -> std::wstring_view { return L"λ"; },
                            []() -> ExprPtr {
                              using namespace sequant::mbpt::sr;
                              return Λ_(1);
                            },
                            [](qns_t& qns) {
                              qns += qns_t{1, 0, 0, 1};
                            });
    auto lambda2 = ex<op_t>([]() -> std::wstring_view { return L"λ"; },
                            []() -> ExprPtr {
                              using namespace sequant::mbpt::sr;
                              return Λ_(2);
                            },
                            [](qns_t& qns) {
                              qns += qns_t{2, 0, 0, 2};
                            });
    auto r_2_1 = ex<op_t>([]() -> std::wstring_view { return L"R"; },
                          []() -> ExprPtr {
                            using namespace sequant::mbpt::sr;
                            return R_(1, 2);
                          },
                          [](qns_t& qns) {
                            qns += qns_t{0, 2, 1, 0};
                          });
    auto r_1_2 = ex<op_t>([]() -> std::wstring_view { return L"R"; },
                          []() -> ExprPtr {
                            using namespace sequant::mbpt::sr;
                            return R_(2, 1);
                          },
                          [](qns_t& qns) {
                            qns += qns_t{0, 1, 2, 0};
                          });

    REQUIRE(to_latex(f) == L"{\\hat{f}}");
    REQUIRE(to_latex(t1) == L"{\\hat{t}_{1}}");
    REQUIRE(to_latex(t2) == L"{\\hat{t}_{2}}");
    REQUIRE(to_latex(lambda1) == L"{\\hat{\\lambda}_{1}}");
    REQUIRE(to_latex(lambda2) == L"{\\hat{\\lambda}_{2}}");
    REQUIRE(to_latex(r_2_1) == L"{\\hat{R}_{2,1}}");
    REQUIRE(to_latex(r_1_2) == L"{\\hat{R}_{1,2}}");

  }  // SECTION("to_latex")

  SECTION("canonicalize") {
    using qns_t = mbpt::sr::qns_t;
    using op_t = mbpt::Operator<qns_t>;
    auto f = ex<op_t>([]() -> std::wstring_view { return L"f"; },
                      []() -> ExprPtr {
                        using namespace sequant::mbpt::sr;
                        return F();
                      },
                      [](qns_t& qns) {
                        qns += qns_t{{0, 1}, {0, 1}, {0, 1}, {0, 1}};
                      });
    auto t1 = ex<op_t>([]() -> std::wstring_view { return L"t"; },
                       []() -> ExprPtr {
                         using namespace sequant::mbpt::sr;
                         return T_(1);
                       },
                       [](qns_t& qns) {
                         qns += qns_t{0, 1, 1, 0};
                       });
    auto l1 = ex<op_t>([]() -> std::wstring_view { return L"λ"; },
                       []() -> ExprPtr {
                         using namespace sequant::mbpt::sr;
                         return Λ_(1);
                       },
                       [](qns_t& qns) {
                         qns += qns_t{1, 0, 0, 1};
                       });
    auto t2 = ex<op_t>([]() -> std::wstring_view { return L"t"; },
                       []() -> ExprPtr {
                         using namespace sequant::mbpt::sr;
                         return T_(2);
                       },
                       [](qns_t& qns) {
                         qns += qns_t{0, 2, 2, 0};
                       });
    auto l2 = ex<op_t>([]() -> std::wstring_view { return L"λ"; },
                       []() -> ExprPtr {
                         using namespace sequant::mbpt::sr;
                         return Λ_(2);
                       },
                       [](qns_t& qns) {
                         qns += qns_t{2, 0, 0, 2};
                       });

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

    auto t = t1 + t2;

    if constexpr (hash_version() == hash::Impl::BoostPre181) {
      REQUIRE(
          to_latex(simplify(f * t * t)) ==
          to_latex(ex<Constant>(2) * f * t1 * t2 + f * t1 * t1 + f * t2 * t2));
    } else {
      //      std::wcout << "to_latex(simplify(f * t * t)): "
      //                 << to_latex(simplify(f * t * t)) << std::endl;
      REQUIRE(
          to_latex(simplify(f * t * t)) ==
          to_latex(ex<Constant>(2) * f * t1 * t2 + f * t2 * t2 + f * t1 * t1));
    }

    if constexpr (hash_version() == hash::Impl::BoostPre181) {
      REQUIRE(to_latex(simplify(f * t * t * t)) ==
              to_latex(f * t1 * t1 * t1 + ex<Constant>(3) * f * t1 * t2 * t2 +
                       f * t2 * t2 * t2 + ex<Constant>(3) * f * t1 * t1 * t2));
    } else {
      //      std::wcout << "to_latex(simplify(f * t * t * t): "
      //                 << to_latex(simplify(f * t * t * t)) << std::endl;
      REQUIRE(to_latex(simplify(f * t * t * t)) ==
              to_latex(f * t2 * t2 * t2 + f * t1 * t1 * t1 +
                       ex<Constant>(3) * f * t1 * t1 * t2 +
                       ex<Constant>(3) * f * t1 * t2 * t2));
    }

  }  // SECTION("canonicalize")

  SECTION("adjoint") {
    using qns_t = mbpt::sr::qns_t;
    using op_t = mbpt::Operator<qns_t>;
    op_t f([]() -> std::wstring_view { return L"f"; },
           []() -> ExprPtr {
             using namespace sequant::mbpt::sr;
             return F();
           },
           [](qns_t& qns) {
             qns += qns_t{{0, 1}, {0, 1}, {0, 1}, {0, 1}};
           });
    op_t t1([]() -> std::wstring_view { return L"t"; },
            []() -> ExprPtr {
              using namespace sequant::mbpt::sr;
              return T_(1);
            },
            [](qns_t& qns) {
              qns += qns_t{0, 1, 1, 0};
            });
    op_t lambda2([]() -> std::wstring_view { return L"λ"; },
                 []() -> ExprPtr {
                   using namespace sequant::mbpt::sr;
                   return Λ_(2);
                 },
                 [](qns_t& qns) {
                   qns += qns_t{2, 0, 0, 2};
                 });
    op_t r_1_2([]() -> std::wstring_view { return L"R"; },
               []() -> ExprPtr {
                 using namespace sequant::mbpt::sr;
                 return R_(2, 1);
               },
               [](qns_t& qns) {
                 qns += qns_t{0, 1, 2, 0};
               });

    REQUIRE_NOTHROW(adjoint(f));
    REQUIRE_NOTHROW(adjoint(t1));
    REQUIRE_NOTHROW(adjoint(lambda2));
    REQUIRE_NOTHROW(adjoint(r_1_2));

    REQUIRE(adjoint(f)() == qns_t{{0, 1}, {0, 1}, {0, 1}, {0, 1}});
    REQUIRE(adjoint(t1)() == qns_t{{1, 1}, {0, 0}, {0, 0}, {1, 1}});
    REQUIRE(adjoint(lambda2)() == qns_t{{0, 0}, {2, 2}, {2, 2}, {0, 0}});
    REQUIRE(adjoint(r_1_2)() == qns_t{{1, 1}, {0, 0}, {0, 0}, {2, 2}});

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
    REQUIRE(to_latex(adjoint(lambda2).as<Expr>()) == L"{\\hat{\\lambda⁺}^{2}}");
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
    using namespace sequant::mbpt::sr::op;

    auto g_t2_t2 = H_(2) * T_(2) * T_(2);
    REQUIRE(raises_vacuum_to_rank(g_t2_t2, 2));
    REQUIRE(raises_vacuum_up_to_rank(g_t2_t2, 2));

    auto g_t2 = H_(2) * T_(2);
    REQUIRE(raises_vacuum_to_rank(g_t2, 3));

    auto lambda2_f = Λ_(2) * H_(1);
    REQUIRE(lowers_rank_to_vacuum(lambda2_f, 2));

  }  // SECTION("screen")

  SECTION("operators") {
    using namespace sequant::mbpt;

    auto R_2 = sr::op::R_(2, 2)->as<sr::op_t>();
    //    std::wcout << "R_2: " << to_latex(simplify(R_2.tensor_form())) <<
    //    std::endl;
    REQUIRE(to_latex(simplify(R_2.tensor_form())) ==
            L"{{{\\frac{1}{4}}}{R^{{i_1}{i_2}}_{{a_1}{a_2}}}{\\tilde{a}^{{a_1}{"
            L"a_2}}_{{i_1}{i_2}}}}");

    auto L_3 = sr::op::L_(3, 3)->as<sr::op_t>();
    //    std::wcout << "L_3: " << to_latex(simplify(L_3.tensor_form())) <<
    //    std::endl;
    REQUIRE(to_latex(simplify(L_3.tensor_form())) ==
            L"{{{\\frac{1}{36}}}{L^{{a_1}{a_2}{a_3}}_{{i_1}{i_2}{i_3}}}{"
            L"\\tilde{a}^{{i_1}{i_2}{i_3}}_{{a_1}{a_2}{a_3}}}}");

    auto R_2_3 = sr::op::R_(2, 3)->as<sr::op_t>();
    //    std::wcout << "R_2_3: " << to_latex(simplify(R_2_3.tensor_form())) <<
    //    std::endl;
    REQUIRE(to_latex(simplify(R_2_3.tensor_form())) ==
            L"{{{\\frac{1}{12}}}{R^{{i_1}{i_2}}_{{a_1}{a_2}{a_3}}}{\\tilde{a}^{"
            L"{a_1}{a_2}{a_3}}_{\\textvisiblespace\\,{i_1}{i_2}}}}");

    auto L_1_2 = sr::op::L_(1, 2)->as<sr::op_t>();
    //    std::wcout << "L_1_2: " << to_latex(simplify(L_1_2.tensor_form())) <<
    //    std::endl;
    REQUIRE(to_latex(simplify(L_1_2.tensor_form())) ==
            L"{{{\\frac{1}{2}}}{L^{{a_1}{a_2}}_{{i_1}}}{\\tilde{a}^{"
            L"\\textvisiblespace\\,{i_1}}_{{a_1}{a_2}}}}");

    auto A_1_2 = sr::op::A(1, 2)->as<sr::op_t>();
    //    std::wcout << "A_1_2: " << to_latex(simplify(A_1_2.tensor_form()))
    //               << std::endl;
    REQUIRE(to_latex(simplify(A_1_2.tensor_form())) ==
            L"{{{\\frac{1}{2}}}{A^{{i_1}}_{{a_1}{a_2}}}{\\tilde{a}^{{a_1}{a_2}}"
            L"_{\\textvisiblespace\\,{i_1}}}}");

    auto P_1_0 = sr::op::P(1, 0)->as<sr::op_t>();
    //    std::wcout << "P_1_0: " << to_latex(simplify(P_1_0.tensor_form()))
    //               << std::endl;
    REQUIRE(to_latex(simplify(P_1_0.tensor_form())) ==
            L"{{A^{}_{{i_1}}}{\\tilde{a}^{{i_1}}}}");

    auto P_1_2 = sr::op::P(1, 2)->as<sr::op_t>();
    //    std::wcout << "P_1_2: " << to_latex(simplify(P_1_2.tensor_form()))
    //               << std::endl;
    REQUIRE(to_latex(simplify(P_1_2.tensor_form())) ==
            L"{{{\\frac{1}{2}}}{A^{{a_1}{a_2}}_{{i_1}}}{\\tilde{a}^{"
            L"\\textvisiblespace\\,{i_1}}_{{a_1}{a_2}}}}");

    auto P_3_2 = sr::op::P(3, 2)->as<sr::op_t>();
    //    std::wcout << "P_3_2: " << to_latex(simplify(P_3_2.tensor_form()))
    //               << std::endl;
    REQUIRE(to_latex(simplify(P_3_2.tensor_form())) ==
            L"{{{\\frac{1}{12}}}{A^{{a_1}{a_2}}_{{i_1}{i_2}{i_3}}}{\\tilde{a}^{"
            L"{i_1}{i_2}{i_3}}_{\\textvisiblespace\\,{a_1}{a_2}}}}");

    auto R33 = mbpt::sr::op::R(3, 3);
    sr::op::lower_to_tensor_form(R33);
    simplify(R33);
    //    std::wcout << "R33: " << to_latex(R33) << std::endl;
    if constexpr (hash_version() == hash::Impl::Boost181OrLater) {
      REQUIRE(
          to_latex(R33) ==
          L"{ "
          L"\\bigl({{{\\frac{1}{36}}}{R^{{i_1}{i_2}{i_3}}_{{a_1}{a_2}{a_3}}}{"
          L"\\tilde{a}^{{a_1}{a_2}{a_3}}_{{i_1}{i_2}{i_3}}}} + "
          L"{{R^{{i_1}}_{{a_1}}}{\\tilde{a}^{{a_1}}_{{i_1}}}} + "
          L"{{{\\frac{1}{4}}}{R^{{i_1}{i_2}}_{{a_1}{a_2}}}{\\tilde{a}^{{a_1}{"
          L"a_2}}_{{i_1}{i_2}}}}\\bigr) }");
    } else {
      REQUIRE(to_latex(R33) ==
              L"{ \\bigl({{R^{{i_1}}_{{a_1}}}{\\tilde{a}^{{a_1}}_{{i_1}}}} + "
              L"{{{\\frac{1}{4}}}{R^{{i_1}{i_2}}_{{a_1}{a_2}}}{\\tilde{a}^{{a_"
              L"1}{a_2}}_{{i_1}{i_2}}}} + "
              L"{{{\\frac{1}{36}}}{R^{{i_1}{i_2}{i_3}}_{{a_1}{a_2}{a_3}}}{"
              L"\\tilde{a}^{{a_1}{a_2}{a_3}}_{{i_1}{i_2}{i_3}}}}\\bigr) }");
    }

    auto R12 = mbpt::sr::op::R(1, 2);
    sr::op::lower_to_tensor_form(R12);
    simplify(R12);
    //    std::wcout << "R12: " << to_latex(R12) << std::endl;
    if constexpr (hash_version() == hash::Impl::Boost181OrLater) {
      REQUIRE(to_latex(R12) ==
              L"{ \\bigl({{R^{}_{{a_1}}}{\\tilde{a}^{{a_1}}}} + "
              L"{{{\\frac{1}{2}}}{R^{{i_1}}_{{a_1}{a_2}}}{\\tilde{a}^{{a_1}{a_"
              L"2}}_{\\textvisiblespace\\,{i_1}}}}\\bigr) }");
    } else {
      REQUIRE(
          to_latex(R12) ==
          L"{ \\bigl({{R^{}_{{a_1}}}{\\tilde{a}^{{a_1}}}} + "
          L"{{{\\frac{1}{2}}}{R^{{i_1}}_{{a_1}{a_2}}}{\\tilde{a}^{{a_1}{a_2}}"
          L"_{\\textvisiblespace\\,{i_1}}}}\\bigr) }");
    }

    auto R21 = mbpt::sr::op::R(2, 1);
    sr::op::lower_to_tensor_form(R21);
    simplify(R21);
    //    std::wcout << "R21: " << to_latex(R21) << std::endl;

    if constexpr (hash_version() == hash::Impl::Boost181OrLater) {
      REQUIRE(to_latex(R21) ==
              L"{ "
              L"\\bigl({{{\\frac{1}{2}}}{R^{{i_1}{i_2}}_{{a_1}}}{\\tilde{a}^{"
              L"\\textvisiblespace\\,{a_1}}_{{i_1}{i_2}}}} + "
              L"{{R^{{i_1}}_{}}{\\tilde{a}_{{i_1}}}}\\bigr) }");
    } else {
      REQUIRE(to_latex(R21) ==
              L"{ "
              L"\\bigl({{{\\frac{1}{2}}}{R^{{i_1}{i_2}}_{{a_1}}}{\\tilde{a}^{"
              L"\\textvisiblespace\\,{a_1}}_{{i_1}{i_2}}}} + "
              L"{{R^{{i_1}}_{}}{\\tilde{a}_{{i_1}}}}\\bigr) }");
    }

    auto L23 = mbpt::sr::op::L(2, 3);
    sr::op::lower_to_tensor_form(L23);
    simplify(L23);
    //        std::wcout << "L23: " << to_latex(L23) << std::endl;
    if constexpr (hash_version() == hash::Impl::Boost181OrLater) {
      REQUIRE(
          to_latex(L23) ==
          L"{ \\bigl({{{\\frac{1}{12}}}{L^{{a_1}{a_2}{a_3}}_{{i_1}{i_2}}}{"
          L"\\tilde{a}^{\\textvisiblespace\\,{i_1}{i_2}}_{{a_1}{a_2}{a_3}}}} + "
          L"{{{\\frac{1}{2}}}{L^{{a_1}{a_2}}_{{i_1}}}{\\tilde{a}^{"
          L"\\textvisiblespace\\,{i_1}}_{{a_1}{a_2}}}} + "
          L"{{L^{{a_1}}_{}}{\\tilde{a}_{{a_1}}}}\\bigr) }");
    } else {
      REQUIRE(
          to_latex(L23) ==
          L"{ "
          L"\\bigl({{{\\frac{1}{2}}}{L^{{a_1}{a_2}}_{{i_1}}}{\\tilde{a}^{"
          L"\\textvisiblespace\\,{i_1}}_{{a_1}{a_2}}}} + "
          L"{{L^{{a_1}}_{}}{\\tilde{a}_{{a_1}}}} + "
          L"{{{\\frac{1}{12}}}{L^{{a_1}{a_2}{a_3}}_{{i_1}{i_2}}}{\\tilde{a}^{"
          L"\\textvisiblespace\\,{i_1}{i_2}}_{{a_1}{a_2}{a_3}}}}\\bigr) }");
    }

  }  // SECTION("operators")

}  // TEST_CASE("NBodyOp")

TEST_CASE("MBPT", "[mbpt]") {
  using namespace sequant;
  TensorCanonicalizer::register_instance(
      std::make_shared<DefaultTensorCanonicalizer>());

  SECTION("SRSO") {
    using namespace sequant::mbpt::sr;

    // H**T12**T12 -> R2
    SEQUANT_PROFILE_SINGLE("wick(H**T12**T12 -> R2)", {
      auto result = vac_av(A(-2) * H() * T(2) * T(2), {{1, 2}, {1, 3}});

      //      std::wcout << "H*T12*T12 -> R2 = " << to_latex_align(result, 20)
      //                 << std::endl;
      REQUIRE(result->size() == 15);

      {  // check against op
        auto result_op = op::vac_av(op::P(2) * op::H() * op::T(2) * op::T(2));
        REQUIRE(result_op->size() ==
                result->size());  // as compact as result ..
        REQUIRE(simplify(result_op - result) ==
                ex<Constant>(0));  // .. and equivalent to it
      }
    });

    // H2**T3**T3 -> R4
    SEQUANT_PROFILE_SINGLE("wick(H2**T3**T3 -> R4)", {
      auto result = vac_av(A(-4) * H_(2) * T_(3) * T_(3), {{1, 2}, {1, 3}});

      std::wcout << "H2**T3**T3 -> R4 = " << to_latex_align(result, 20)
                 << std::endl;
      REQUIRE(result->size() == 4);
    });

    const std::vector<std::pair<std::wstring, std::wstring>> new_op_connect = {
        {L"h", L"t"}, {L"f", L"t"}, {L"g", L"t"}};

#ifndef SEQUANT_SKIP_LONG_TESTS
    // the longest term in CCSDTQP
    // H2**T2**T2**T3 -> R5
    {
      ExprPtr ref_result;
      SEQUANT_PROFILE_SINGLE("wick(H2**T2**T2**T3 -> R5)", {
        ref_result = op::vac_av(
            op::A(-5) * op::H_(2) * op::T_(2) * op::T_(2) * op::T_(3),
            new_op_connect);
        REQUIRE(ref_result->size() == 7);
      });
      // now computed using specific component of H2
      SEQUANT_PROFILE_SINGLE("wick(H2(oo;vv)**T2**T2**T3 -> R5)", {
        auto result = op::vac_av(
            op::A(-5) * op::H2_oo_vv() * op::T_(2) * op::T_(2) * op::T_(3),
            new_op_connect);
        REQUIRE(result->size() == ref_result->size());
      });
    }
#endif  // !defined(SEQUANT_SKIP_LONG_TESTS)
  }     // SECTION ("SRSO")

  SECTION("SRSO Fock") {
    using namespace sequant::mbpt::sr;

    // <2p1h|H2|1p> ->
    SEQUANT_PROFILE_SINGLE("wick(<2p1h|H2|1p>)", ({
                             auto input = L_(1, 2) * H_(2) * R_(1, 0);
                             auto result = vac_av(input);

                             std::wcout << "<2p1h|H2|1p> = " << to_latex(result)
                                        << std::endl;
                             REQUIRE(result->is<Product>());  // product ...
                             REQUIRE(result->size() == 3);  // ... of 3 factors
                           }));

    // <2p1h|H2|2p1h(c)> ->
    SEQUANT_PROFILE_SINGLE(
        "wick(<2p1h|H2|2p1h(c)>)", ({
          auto input = L_(1, 2) * H() * R_(2, 1);
          auto result = vac_av(input);

          std::wcout << "<2p1h|H|2p1h(c)> = " << to_latex(result) << std::endl;
          REQUIRE(result->is<Sum>());    // sub ...
          REQUIRE(result->size() == 4);  // ... of 4 factors
        }));
  }  // SECTION("SRSO Fock")

  SECTION("SRSO-PNO") {
    using namespace sequant::mbpt;
    using namespace sequant::mbpt::sr;
    using sequant::mbpt::Context;
    auto resetter = set_scoped_default_formalism(Context(CSV::Yes));

    // H2**T2**T2 -> R2
    SEQUANT_PROFILE_SINGLE("wick(H2**T2**T2 -> R2)", {
      auto result = vac_av(A(-2) * H_(2) * T_(2) * T_(2), {{1, 2}, {1, 3}});

      std::wcout << "H2**T2**T2 -> R2 = " << to_latex_align(result, 20)
                 << std::endl;
      REQUIRE(result->size() == 4);
    });
  }  // SECTION("SRSO-PNO")

  SECTION("SRSF") {
    using namespace sequant::mbpt::sr;

    auto ctx_resetter = set_scoped_default_context(
        Context(Vacuum::SingleProduct, IndexSpaceMetric::Unit,
                BraKetSymmetry::conjugate, SPBasis::spinfree));

    // H2 -> R2
    SEQUANT_PROFILE_SINGLE("wick(H2 -> R2)", {
      auto result = vac_av(S(-2) * H_(2));

      {
        std::wcout << "H2 -> R2 = " << to_latex_align(result, 0, 1)
                   << std::endl;
      }
    });

    // H2**T2 -> R2
    SEQUANT_PROFILE_SINGLE("wick(H2**T2 -> R2)", {
      auto result = vac_av(S(-2) * H_(2) * T_(2), {{1, 2}});

      {
        std::wcout << "H2**T2 -> R2 = " << to_latex_align(result, 0, 1)
                   << std::endl;
      }
    });
  }  // SECTION("SRSF")

  SECTION("MRSO") {
    using namespace sequant::mbpt::mr;

    // H2**T2 -> 0
    // std::wcout << "H_(2) * T_(2) = " << to_latex(H_(2) * T_(2)) << std::endl;
    SEQUANT_PROFILE_SINGLE("wick(H2**T2 -> 0)", {
      auto result = vac_av(H_(2) * T_(2), {{0, 1}});

      {
        std::wcout << "H2*T2 -> 0 = " << to_latex_align(result, 0, 1)
                   << std::endl;
      }

      auto result_wo_top =
          vac_av(H_(2) * T_(2), {{0, 1}}, /* use_topology = */ false);

      REQUIRE(simplify(result - result_wo_top) == ex<Constant>(0));

      // now compute using physical vacuum
      {
        auto ctx_resetter = set_scoped_default_context(
            Context(Vacuum::Physical, IndexSpaceMetric::Unit,
                    BraKetSymmetry::conjugate, SPBasis::spinorbital));
        auto result_phys = vac_av(H_(2) * T_(2), {{0, 1}});

        {
          std::wcout << "H2*T2 -> 0 using phys vacuum = "
                     << to_latex_align(result_phys, 0, 1) << std::endl;
        }
      }
    });

    // H2 ** T2 ** T2 -> 0
    SEQUANT_PROFILE_SINGLE("wick(H2**T2**T2 -> 0)", {
      // first without use of topology
      auto result = vac_av(H_(2) * T_(2) * T_(2), {{0, 1}, {0, 2}},
                           /* use_topology = */ false);
      // now with topology use
      auto result_top = vac_av(H_(2) * T_(2) * T_(2), {{0, 1}, {0, 2}},
                               /* use_topology = */ true);

      REQUIRE(simplify(result - result_top) == ex<Constant>(0));
    });

#if 0
    // H**T12 -> R2
    SEQUANT_PROFILE_SINGLE("wick(H**T2 -> R2)", {
      auto result = vac_av(A(-2) * H() * T_(2), {{1, 2}});

      {
        std::wcout << "H*T2 -> R2 = " << to_latex_align(result, 0, 1)
                   << std::endl;
      }
    });
#endif

  }  // SECTION("MRSO")

  SECTION("MRSF") {
    using namespace sequant::mbpt::mr;

    // now compute using (closed) Fermi vacuum + spinfree basis
    auto ctx_resetter = set_scoped_default_context(
        Context(Vacuum::SingleProduct, IndexSpaceMetric::Unit,
                BraKetSymmetry::conjugate, SPBasis::spinfree));

    // H2**T2 -> 0
    std::wcout << "H_(2) * T_(2) = " << to_latex(H_(2) * T_(2)) << std::endl;
    SEQUANT_PROFILE_SINGLE("wick(H2**T2 -> 0)", {
      auto result = vac_av(H_(2) * T_(2), {{0, 1}});

      //      {
      //        std::wcout << "H2*T2 -> 0 = " << to_latex_align(result, 0, 1)
      //                   << std::endl;
      //      }

      {  // make sure get same result without use of topology
        auto result_wo_top =
            vac_av(H_(2) * T_(2), {{0, 1}}, /* use_topology = */ false);

        REQUIRE(simplify(result - result_wo_top) == ex<Constant>(0));
      }

      {  // make sure get same result using operators
        auto result_op = op::vac_av(op::H_(2) * op::T_(2));

        REQUIRE(result_op->size() == result->size());
        REQUIRE(simplify(result - result_op) == ex<Constant>(0));
      }
    });

  }  // SECTION("MRSF")

}  // TEST_CASE("MBPT")
