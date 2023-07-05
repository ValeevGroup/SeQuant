//
// Created by Eduard Valeyev on 5/18/23.
//

#include "catch.hpp"

#include <iostream>

#include "SeQuant/core/rational.hpp"
#include "SeQuant/core/wstring.hpp"

TEST_CASE("Rational", "[elements]") {
  using namespace sequant;

  auto print = [](rational r) {
    return sequant::to_wstring(numerator(r)) + L"/" +
           sequant::to_wstring(denominator(r));
  };
  SECTION("to_rational") {
    REQUIRE(to_rational(1. / 3) == rational{1, 3});
    REQUIRE(to_rational(1. / 3, 0.) ==
            rational{6004799503160661ull, 18014398509481984ull});
    REQUIRE(to_rational(1. / 7) == rational{1, 7});
    REQUIRE(to_rational(M_PI) == rational{99023, 31520});
    REQUIRE(to_rational(M_E) == rational{23225, 8544});
    REQUIRE_THROWS_AS(to_rational(std::nan("NaN")), std::invalid_argument);
  }

}  // TEST_CASE("Rational")
