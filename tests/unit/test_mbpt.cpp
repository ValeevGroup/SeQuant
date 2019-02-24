//
// Created by Eduard Valeyev on 2019-02-19.
//

#include "../../src/SeQuant2/mbpt/sr/sr.hpp"
#include "catch.hpp"
#include "timer.hpp"

TEST_CASE("MBPT", "[mbpt]") {
  using namespace sequant2;
  TensorCanonicalizer::set_cardinal_tensor_labels({L"A", L"f", L"g", L"t"});
  TensorCanonicalizer::register_instance(
      std::make_shared<DefaultTensorCanonicalizer>());

  SECTION("SRSO-PNO") {
    using namespace sequant2::mbpt::sr::so::pno;

    // H2**T2**T2 -> R2
    SEQUANT2_PROFILE_SINGLE("wick(H2**T2**T2 -> R2)", {
      auto result = vac_av(A<2>() * H2() * T_<2>() * T_<2>(), {{1, 2}, {1, 3}});

      std::wcout << "H2**T2**T2 -> R2 = " << to_latex_align(result, 20)
                 << std::endl;
      REQUIRE(result->size() == 5);  // but only 4 unique
    });
  }

  SECTION("SRSO") {
    using namespace sequant2::mbpt::sr::so;

#if 0
    // H**T12**T12 -> R2
    SEQUANT2_PROFILE_SINGLE("wick(H**T12**T12 -> R2)", {
      auto result = vac_av( A<2>() * H() * T<2>() * T<2>(), {{1, 2}, {1, 3}} );

      std::wcout << "H*T12*T12 -> R2 = " << to_latex_align(result, 20)
                 << std::endl;
      REQUIRE(result->size() == 0);
    });
#endif

    // H2**T3**T3 -> R4
    SEQUANT2_PROFILE_SINGLE("wick(H2**T3**T3 -> R3)", {
      auto result = vac_av(A<4>() * H2() * T_<3>() * T_<3>(), {{1, 2}, {1, 3}});

      std::wcout << "H2**T3**T3 -> R4 = " << to_latex_align(result, 20)
                 << std::endl;
      REQUIRE(result->size() == 7);  // 2 disconnected + 7 connected terms
    });
  }

  SECTION("SRSO-PNO") {
    using namespace sequant2::mbpt::sr::so::pno;

    // H2**T2**T2 -> R2
    SEQUANT2_PROFILE_SINGLE("wick(H2**T2**T2 -> R2)", {
      auto result = vac_av(A<2>() * H2() * T_<2>() * T_<2>(), {{1, 2}, {1, 3}});

      std::wcout << "H2**T2**T2 -> R2 = " << to_latex_align(result, 20)
                 << std::endl;
      REQUIRE(result->size() == 5);  // but only 4 unique
    });
  }

}  // TEST_CASE("MBPT")