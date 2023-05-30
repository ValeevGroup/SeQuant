//
// Created by Eduard Valeyev on 2019-02-19.
//

#include "SeQuant/core/op.hpp"
#include "SeQuant/core/tensor.hpp"
#include "SeQuant/core/timer.hpp"
#include "SeQuant/domain/mbpt/formalism.hpp"
#include "SeQuant/domain/mbpt/op.hpp"
#include "SeQuant/domain/mbpt/sr/sr.hpp"

#include "catch.hpp"

TEST_CASE("NBodyOp", "[mbpt]") {
  using namespace sequant;

  SECTION("constructor") {
    auto test1 = [](const auto& null_qns) {
      using qns_t = std::decay_t<decltype(null_qns)>;
      using op_t = mbpt::Operator<qns_t>;
      op_t f1([]() -> std::wstring_view { return L"f"; },
              []() -> ExprPtr {
                return ex<Tensor>(L"f", WstrList{L"p_1"}, WstrList{L"p_2"}) *
                       ex<FNOperator>(WstrList{L"p_1"}, WstrList{L"p_2"});
              },
              [](qns_t& qns) {
                qns += qns_t{-1, 1};
              });
      REQUIRE(f1.label() == L"f");
      {  // possible compare
        using namespace boost::numeric::interval_lib::compare::possible;
        // REQUIRE(f1(qns_t{0, 0}) == 1);   // this is not same as below, due to
        // Catch interference
        REQUIRE(operator==(f1(qns_t{0, 0}),
                           1));  // can produce single excitation
        REQUIRE(operator==(f1(qns_t{0, 0}),
                           -1));  // can produce single de-excitation
        REQUIRE(operator==(
            f1(qns_t{0, 0}),
            0));  // F is normal, so must excite/deexcite, but there is no way
                  // to exclude 0 from boost::interval
        REQUIRE(operator==(
            f1(qns_t{1, 1}),
            0));  // can produce reference when acting on singly-excited
        REQUIRE(operator!=(
            f1(qns_t{2, 2}),
            0));  // cannot produce reference when acting on doubly-excited

        if constexpr (std::is_same_v<qns_t, mbpt::ParticleNumberChange<1>>) {
          REQUIRE(f1(qns_t{0, 0}).in(1));   // can produce single excitation
          REQUIRE(f1(qns_t{0, 0}).in(-1));  // can produce single deexcitation
          REQUIRE(
              f1(qns_t{0, 0})
                  .in(0));  // F is normal, so must excite/deexcite, but there
                            // is no way to exclude 0 from boost::interval
          REQUIRE(f1(qns_t{1, 1}).in(0));   // can produce reference when acting
                                            // on singly-excited
          REQUIRE(!f1(qns_t{2, 2}).in(0));  // can't produce reference when
                                            // acting on doubly-excited
        }
      }
      {  // equal compare
        // using namespace boost::numeric::interval_lib::compare::lexicographic;
        // REQUIRE(f1(qns_t{0, 0}) == qns_t{-1, 1}); // not same as below due to
        // interaction with Catch
        // could do REQUIRE(operator==(f1(qns_t{0, 0}), qns_t{-1, 1})); but
        // equal is shorter
        REQUIRE(equal(f1(qns_t{0, 0}), qns_t{-1, 1}));
        REQUIRE(equal(f1(qns_t{-1, 1}), qns_t{-2, 2}));
      }
    };
    test1(boost::numeric::interval<int>{});
    test1(mbpt::ParticleNumberChange<1>{});

    auto test2 = [](const auto& null_qns) {
      using qns_t = std::decay_t<decltype(null_qns)>;
      using op_t = mbpt::Operator<qns_t>;
      op_t f1([]() -> std::wstring_view { return L"f"; },
              []() -> ExprPtr {
                return ex<Tensor>(L"f", WstrList{L"p_1"}, WstrList{L"p_2"}) *
                       ex<FNOperator>(WstrList{L"p_1"}, WstrList{L"p_2"});
              },
              [](qns_t& qns) {
                qns += qns_t{{0, 0}, {-1, 1}};
              });
      REQUIRE(f1.label() == L"f");
      REQUIRE(
          f1(qns_t{}) ==
          qns_t{
              {0, 0},
              {-1, 1}});  // can produce single de/excitation when applied once
      REQUIRE(f1(f1(qns_t{})) ==
              qns_t{{0, 0}, {-2, 2}});  // can produce up to double
                                        // de/excitation when applied twice
      REQUIRE(f1(qns_t{}).in(std::array{0, 1}));
      REQUIRE(f1(qns_t{}).in(std::array{0, -1}));
      REQUIRE(f1(qns_t{}).in(std::array{
          0, 0}));  // f1 is normal, so must excite/deexcite, but there is no
                    // way to exclude 0 from boost::interbal    }
      REQUIRE(
          f1(qns_t{0, 1}).in(std::array{0, 0}));  // can produce reference when
                                                  // acting on singly-excited
      REQUIRE(!f1(qns_t{0, 2})
                   .in(std::array{0, 0}));  // cannot produce reference when
                                            // acting on doubly-excited
    };
    test2(mbpt::ParticleNumberChange<2>{});
  }  // SECTION("constructor")

  SECTION("to_latex") {
    using qns_t = mbpt::ParticleNumberChange<2>;
    using op_t = mbpt::Operator<qns_t>;
    auto f = ex<op_t>([]() -> std::wstring_view { return L"f"; },
                      []() -> ExprPtr {
                        using namespace sequant::mbpt::sr;
                        return F();
                      },
                      [](qns_t& qns) {
                        qns += qns_t{{0, 0}, {-1, 1}};
                      });
    auto t1 = ex<op_t>([]() -> std::wstring_view { return L"t"; },
                       []() -> ExprPtr {
                         using namespace sequant::mbpt::sr;
                         return T_(1);
                       },
                       [](qns_t& qns) {
                         qns += qns_t{{0, 0}, {1, 1}};
                       });
    auto t2 = ex<op_t>([]() -> std::wstring_view { return L"t"; },
                       []() -> ExprPtr {
                         using namespace sequant::mbpt::sr;
                         return T_(2);
                       },
                       [](qns_t& qns) {
                         qns += qns_t{{0, 0}, {2, 2}};
                       });
    auto lambda1 = ex<op_t>([]() -> std::wstring_view { return L"λ"; },
                            []() -> ExprPtr {
                              using namespace sequant::mbpt::sr;
                              return Lambda_(1);
                            },
                            [](qns_t& qns) {
                              qns += qns_t{{0, 0}, {-1, -1}};
                            });
    auto lambda2 = ex<op_t>([]() -> std::wstring_view { return L"λ"; },
                            []() -> ExprPtr {
                              using namespace sequant::mbpt::sr;
                              return Lambda_(2);
                            },
                            [](qns_t& qns) {
                              qns += qns_t{{0, 0}, {-2, -2}};
                            });
    auto r_2_1 = ex<op_t>([]() -> std::wstring_view { return L"R"; },
                          []() -> ExprPtr {
                            using namespace sequant::mbpt::sr;
                            return R_(1, 2);
                          },
                          [](qns_t& qns) {
                            qns += qns_t{{-1, -1}, {1, 1}};
                          });
    auto r_1_2 = ex<op_t>([]() -> std::wstring_view { return L"R"; },
                          []() -> ExprPtr {
                            using namespace sequant::mbpt::sr;
                            return R_(2, 1);
                          },
                          [](qns_t& qns) {
                            qns += qns_t{{+1, +1}, {2, 2}};
                          });

    REQUIRE(to_latex(f) == L"{\\hat{f}}");
    REQUIRE(to_latex(t1) == L"{\\hat{t}_{1}}");
    REQUIRE(to_latex(t2) == L"{\\hat{t}_{2}}");
    REQUIRE(to_latex(lambda1) == L"{\\hat{\\lambda}_{1}}");
    REQUIRE(to_latex(lambda2) == L"{\\hat{\\lambda}_{2}}");
    //    std::wcout << "to_latex(r_2_1) = " << to_latex(r_2_1) << std::endl;
    //    std::wcout << "to_latex(r_2_1->tensor_form()) = "
    //               << to_latex(r_2_1->as<op_t>().tensor_form()) << std::endl;
    REQUIRE(to_latex(r_2_1) == L"{\\hat{R}_{-2}^{1}}");
    REQUIRE(to_latex(r_1_2) == L"{\\hat{R}_{-1}^{2}}");
  }

  SECTION("canonicalize") {
    using qns_t = mbpt::ParticleNumberChange<2>;
    using op_t = mbpt::Operator<qns_t>;
    auto f = ex<op_t>([]() -> std::wstring_view { return L"f"; },
                      []() -> ExprPtr {
                        using namespace sequant::mbpt::sr;
                        return F();
                      },
                      [](qns_t& qns) {
                        qns += qns_t{{0, 0}, {-1, 1}};
                      });
    auto t1 = ex<op_t>([]() -> std::wstring_view { return L"t"; },
                       []() -> ExprPtr {
                         using namespace sequant::mbpt::sr;
                         return T_(1);
                       },
                       [](qns_t& qns) {
                         qns += qns_t{{0, 0}, {1, 1}};
                       });
    auto t2 = ex<op_t>([]() -> std::wstring_view { return L"t"; },
                       []() -> ExprPtr {
                         using namespace sequant::mbpt::sr;
                         return T_(2);
                       },
                       [](qns_t& qns) {
                         qns += qns_t{{0, 0}, {2, 2}};
                       });

    //    std::wcout << "f1 * t2 * t1=" << to_latex(f1 * t2 * t1) << std::endl;
    //    std::wcout << "canonicalize(f1 * t2 * t1)=" <<
    //    to_latex(canonicalize(f1 * t2 * t1)) << std::endl;
    REQUIRE(to_latex(f * t1 * t2) == to_latex(canonicalize(f * t2 * t1)));

    REQUIRE(to_latex(ex<Constant>(3) * f * t1 * t2) ==
            to_latex(simplify(ex<Constant>(2) * f * t2 * t1 + f * t1 * t2)));
    auto t = t1 + t2;
    // std::wcout << "\\hat{f} \\hat{t} \\hat{t} = " << to_latex(simplify(f1 * t
    // * t)) << std::endl;
    if constexpr (hash_version() == hash::Impl::BoostPre181) {
      REQUIRE(
          to_latex(simplify(f * t * t)) ==
          to_latex(f * t1 * t1 + ex<Constant>(2) * f * t1 * t2 + f * t2 * t2));
    } else {
      REQUIRE(
          to_latex(simplify(f * t * t)) ==
          to_latex(f * t1 * t1 + f * t2 * t2 + ex<Constant>(2) * f * t1 * t2));
    }
    // std::wcout << "\\hat{f} \\hat{t} \\hat{t} \\hat{t} = " <<
    // to_latex(simplify(f * t * t * t)) << std::endl;
    if constexpr (hash_version() == hash::Impl::BoostPre181) {
      REQUIRE(to_latex(simplify(f * t * t * t)) ==
              to_latex(ex<Constant>(3) * f * t1 * t2 * t2 + f * t1 * t1 * t1 +
                       ex<Constant>(3) * f * t1 * t1 * t2 + f * t2 * t2 * t2));
    } else {
      REQUIRE(to_latex(simplify(f * t * t * t)) ==
              to_latex(f * t1 * t1 * t1 + ex<Constant>(3) * f * t1 * t2 * t2 +
                       f * t2 * t2 * t2 + ex<Constant>(3) * f * t1 * t1 * t2));
    }
  }  // SECTION("canonicalize")

  SECTION("adjoint") {
    using qns_t = mbpt::ParticleNumberChange<2>;
    using op_t = mbpt::Operator<qns_t>;
    op_t f([]() -> std::wstring_view { return L"f"; },
           []() -> ExprPtr {
             using namespace sequant::mbpt::sr;
             return F();
           },
           [](qns_t& qns) {
             qns += qns_t{{0, 0}, {-1, 1}};
           });
    op_t t1([]() -> std::wstring_view { return L"t"; },
            []() -> ExprPtr {
              using namespace sequant::mbpt::sr;
              return T_(1);
            },
            [](qns_t& qns) {
              qns += qns_t{{0, 0}, {1, 1}};
            });
    op_t lambda2([]() -> std::wstring_view { return L"λ"; },
                 []() -> ExprPtr {
                   using namespace sequant::mbpt::sr;
                   return Lambda_(2);
                 },
                 [](qns_t& qns) {
                   qns += qns_t{{0, 0}, {-2, -2}};
                 });
    op_t r_1_2([]() -> std::wstring_view { return L"R"; },
               []() -> ExprPtr {
                 using namespace sequant::mbpt::sr;
                 return R_(2, 1);
               },
               [](qns_t& qns) {
                 qns += qns_t{{+1, +1}, {2, 2}};
               });

    REQUIRE_NOTHROW(adjoint(f));
    REQUIRE_NOTHROW(adjoint(t1));
    REQUIRE_NOTHROW(adjoint(lambda2));
    REQUIRE_NOTHROW(adjoint(r_1_2));

    REQUIRE(adjoint(f)(qns_t{}) == qns_t{{0, 0}, {-1, 1}});
    REQUIRE(adjoint(t1)(qns_t{}) == qns_t{{0, 0}, {-1, -1}});
    REQUIRE(adjoint(lambda2)(qns_t{}) == qns_t{{0, 0}, {2, 2}});
    REQUIRE(adjoint(r_1_2)(qns_t{}) == qns_t{{-1, -1}, {-2, -2}});

    // adjoint(adjoint(op)) = op
    REQUIRE(adjoint(adjoint(t1))(qns_t{}) == t1(qns_t{}));
    REQUIRE(adjoint(adjoint(r_1_2))(qns_t{}) == r_1_2(qns_t{}));

    // tensor_form()
    //    std::wcout << to_latex(simplify(r_1_2.tensor_form())) << std::endl;
    //    std::wcout << to_latex(simplify(adjoint(r_1_2).tensor_form())) <<
    //    std::endl; std::wcout << to_latex(simplify(mbpt::sr::T_(1))) <<
    //    std::endl; std::wcout <<
    //    to_latex(simplify(adjoint(r_1_2.tensor_form()))) << std::endl;

    REQUIRE(to_latex(simplify(adjoint(t1).tensor_form())) ==
            L"{{t^{{a_1}}_{{i_1}}}{\\tilde{a}^{{i_1}}_{{a_1}}}}");
    REQUIRE((simplify(adjoint(t1).tensor_form())) ==
            (simplify(adjoint(t1.tensor_form()))));
    REQUIRE(to_latex(simplify(adjoint(lambda2).tensor_form())) ==
            L"{{{\\frac{1}{4}}}{\\bar{λ}^{{i_1}{i_2}}_{{a_1}{a_2}}}{\\tilde{a}^"
            L"{{a_2}{a_1}}_{{i_2}{i_1}}}}");
    REQUIRE(simplify(adjoint(lambda2).tensor_form()) ==
            simplify(adjoint(lambda2.tensor_form())));
    REQUIRE(to_latex(simplify(adjoint(r_1_2).tensor_form())) ==
            L"{{{-\\frac{1}{2}}}{R^{{a_1}{a_2}}_{{i_1}}}{\\tilde{a}^{"
            L"\\textvisiblespace\\,{i_1}}_{{a_2}{a_1}}}}");
    REQUIRE(simplify(adjoint(r_1_2).tensor_form()) ==
            simplify(adjoint(r_1_2.tensor_form())));

    // to_latex()
    REQUIRE(to_latex(f.as<Expr>()) == L"{\\hat{f}}");
    REQUIRE(to_latex(t1.as<Expr>()) == L"{\\hat{t}_{1}}");
    REQUIRE(to_latex(lambda2.as<Expr>()) == L"{\\hat{\\lambda}_{2}}");
    REQUIRE(to_latex(r_1_2.as<Expr>()) == L"{\\hat{R}_{-1}^{2}}");
    //    std::wcout << "t1: " << to_latex(t1.as<Expr>()) << std::endl;
    //    std::wcout << "lambda2: " << to_latex(lambda2.as<Expr>()) <<
    //    std::endl; std::wcout << "r_1_2: " << to_latex(r_1_2.as<Expr>()) <<
    //    std::endl; std::wcout << "f: " << to_latex(f.as<Expr>()) << std::endl;

  }  // SECTION("adjoint")
}

TEST_CASE("MBPT", "[mbpt]") {
  using namespace sequant;
  TensorCanonicalizer::register_instance(
      std::make_shared<DefaultTensorCanonicalizer>());

  SECTION("SRSO") {
    using namespace sequant::mbpt::sr;

    // H**T12**T12 -> R2
    SEQUANT_PROFILE_SINGLE("wick(H**T12**T12 -> R2)", {
      auto result = vac_av(A(2) * H() * T(2) * T(2), {{1, 2}, {1, 3}});

      std::wcout << "H*T12*T12 -> R2 = " << to_latex_align(result, 20)
                 << std::endl;
      REQUIRE(result->size() == 15);
    });

    // H2**T3**T3 -> R4
    SEQUANT_PROFILE_SINGLE("wick(H2**T3**T3 -> R4)", {
      auto result = vac_av(A(4) * H2() * T_(3) * T_(3), {{1, 2}, {1, 3}});

      std::wcout << "H2**T3**T3 -> R4 = " << to_latex_align(result, 20)
                 << std::endl;
      REQUIRE(result->size() == 4);
    });
  }

  SECTION("SRSO Fock") {
    using namespace sequant::mbpt::sr;

    // <2p1h|H2|1p> ->
    SEQUANT_PROFILE_SINGLE("wick(<2p1h|H2|1p>)", ({
                             auto input = L_(1, 2) * H2() * R_(1, 0);
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
  }

  SECTION("SRSO-PNO") {
    using namespace sequant::mbpt::sr;
    using namespace sequant::mbpt;
    auto resetter = set_scoped_default_formalism(
        Formalism::make_default().set(CSVFormalism::CSV));

    // H2**T2**T2 -> R2
    SEQUANT_PROFILE_SINGLE("wick(H2**T2**T2 -> R2)", {
      auto result = vac_av(A(2) * H2() * T_(2) * T_(2), {{1, 2}, {1, 3}});

      std::wcout << "H2**T2**T2 -> R2 = " << to_latex_align(result, 20)
                 << std::endl;
      REQUIRE(result->size() == 4);
    });
  }

}  // TEST_CASE("MBPT")
