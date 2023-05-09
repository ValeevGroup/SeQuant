//
// Created by Eduard Valeyev on 2019-02-19.
//

#include "SeQuant/core/op.hpp"
#include "SeQuant/core/tensor.hpp"
#include "SeQuant/core/timer.hpp"
#include "SeQuant/domain/mbpt/formalism.hpp"
#include "SeQuant/domain/mbpt/op.hpp"
#include "SeQuant/domain/mbpt/sr/sr.hpp"
#include "SeQuant/external/boost/interval.hpp"

#include "catch.hpp"

TEST_CASE("NBodyOp", "[mbpt]") {
  using namespace sequant;

  SECTION("constructor") {
    using qns_t = boost::numeric::interval<int>;
    using op_t = mbpt::Operator<qns_t>;
    op_t f1([]() { return L"F"; },
            []() -> ExprPtr {
              return ex<Tensor>(L"F", WstrList{L"p_1"}, WstrList{L"p_2"}) *
                     ex<FNOperator>(WstrList{L"p_1"}, WstrList{L"p_2"});
            },
            [](qns_t& qns) {
              qns += qns_t{-1, 1};
            });
    REQUIRE(f1.label() == L"F");
    {  // possible compare
      using namespace boost::numeric::interval_lib::compare::possible;
      // REQUIRE(f1(qns_t{0, 0}) == 1);   // this is not same as below, due to
      // Catch interference
      REQUIRE(operator==(f1(qns_t{0, 0}), 1));  // can produce single excitation
      REQUIRE(operator==(f1(qns_t{0, 0}),
                         -1));  // can produce single de-excitation
      REQUIRE(operator==(
          f1(qns_t{0, 0}),
          0));  // F is normal, so must excite/deexcite, but there is no way to
                // exclude 0 from boost::interbal    }
      REQUIRE(operator==(
          f1(qns_t{1, 1}),
          0));  // can produce reference when acting on singly-excited
      REQUIRE(operator!=(
          f1(qns_t{2, 2}),
          0));  // cannot produce reference when acting on doubly-excited
    }
    {  // equal compare
      // using namespace boost::numeric::interval_lib::compare::lexicographic;
      // REQUIRE(f1(qns_t{0, 0}) == qns_t{-1, 1}); // not same as below due to
      // interaction with Catch
      // could do REQUIRE(operator==(f1(qns_t{0, 0}), qns_t{-1, 1})); but equal
      // is shorter
      REQUIRE(equal(f1(qns_t{0, 0}), qns_t{-1, 1}));
      REQUIRE(equal(f1(qns_t{-1, 1}), qns_t{-2, 2}));
    }
  }
}

TEST_CASE("MBPT", "[mbpt]") {
  using namespace sequant;
  TensorCanonicalizer::register_instance(
      std::make_shared<DefaultTensorCanonicalizer>());

  SECTION("SRSO") {
    using namespace sequant::mbpt::sr::so;

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
    using namespace sequant::mbpt::sr::so;

    // <2p1h|H2|1p> ->
    SEQUANT_PROFILE_SINGLE("wick(<2p1h|H2|1p>)", ({
                             auto input = L(1, 2) * H2() * R(1, 0);
                             auto result = vac_av(input);

                             std::wcout << "<2p1h|H2|1p> = " << to_latex(result)
                                        << std::endl;
                             REQUIRE(result->is<Product>());  // product ...
                             REQUIRE(result->size() == 3);  // ... of 3 factors
                           }));

    // <2p1h|H2|2p1h(c)> ->
    SEQUANT_PROFILE_SINGLE(
        "wick(<2p1h|H2|2p1h(c)>)", ({
          auto input = L(1, 2) * H() * R(2, 1);
          auto result = vac_av(input);

          std::wcout << "<2p1h|H|2p1h(c)> = " << to_latex(result) << std::endl;
          REQUIRE(result->is<Sum>());    // sub ...
          REQUIRE(result->size() == 4);  // ... of 4 factors
        }));
  }

  SECTION("SRSO-PNO") {
    using namespace sequant::mbpt::sr::so;
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
