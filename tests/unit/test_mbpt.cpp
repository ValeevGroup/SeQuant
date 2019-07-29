//
// Created by Eduard Valeyev on 2019-02-19.
//

#include "../../src/SeQuant/core/timer.hpp"
#include "../../src/SeQuant/domain/mbpt/sr/sr.hpp"
#include "../../src/SeQuant/core/tensor.hpp"
#include "catch.hpp"

TEST_CASE("MBPT", "[mbpt]") {
  using namespace sequant;
  TensorCanonicalizer::register_instance(
      std::make_shared<DefaultTensorCanonicalizer>());

  SECTION("SRSO") {
    using namespace sequant::mbpt::sr::so;

    // H**T12**T12 -> R2
    SEQUANT_PROFILE_SINGLE("wick(H**T12**T12 -> R2)", {
      auto result = vac_av( A(2) * H() * T(2) * T(2), {{1, 2}, {1, 3}} );

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
      auto result = vac_av( input );

      std::wcout << "<2p1h|H2|1p> = " << to_latex(result)
                 << std::endl;
      REQUIRE(result->is<Product>()); // product ...
      REQUIRE(result->size() == 3); // ... of 3 factors
      }));

    // <2p1h|H2|2p1h(c)> ->
    SEQUANT_PROFILE_SINGLE("wick(<2p1h|H2|2p1h(c)>)", ({
      auto input = L(1, 2) * H() * R(2, 1, true);
      auto result = vac_av( input );

      std::wcout << "<2p1h|H|2p1h(c)> = " << to_latex(result)
                 << std::endl;
      REQUIRE(result->is<Sum>()); // sub ...
      REQUIRE(result->size() == 4); // ... of 4 factors
    }));
  }

  SECTION("SRSO-PNO") {
    using namespace sequant::mbpt::sr::so::csv;

    // H2**T2**T2 -> R2
    SEQUANT_PROFILE_SINGLE("wick(H2**T2**T2 -> R2)", {
      auto result = vac_av(A(2) * H2() * T_(2) * T_(2), {{1, 2}, {1, 3}});

      std::wcout << "H2**T2**T2 -> R2 = " << to_latex_align(result, 20)
                 << std::endl;
      REQUIRE(result->size() == 4);
    });
  }

}  // TEST_CASE("MBPT")