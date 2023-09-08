//
// Created by Eduard Valeyev on 5/18/23.
//

#include "catch.hpp"

#include "SeQuant/core/math.hpp"
#include "SeQuant/core/rational.hpp"
#include "SeQuant/core/runtime.hpp"
#include "SeQuant/core/wstring.hpp"
#include "SeQuant/core/meta.hpp"

#include <iostream>
#include <cmath>
#include <new>
#include <stdexcept>
#include <string>
#include <utility>

#include <range/v3/all.hpp>

TEST_CASE("Rational", "[elements]") {
  using namespace sequant;

  [[maybe_unused]] auto print = [](rational r) {
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

TEST_CASE("Factorial", "[math]") {
  using namespace sequant;
  SECTION("factorial") {
    REQUIRE(sequant::to_string(sequant::factorial(30)) ==
            "265252859812191058636308480000000");
    // 21! has been memoized by now
    REQUIRE(sequant::to_string(sequant::factorial(21)) ==
            "51090942171709440000");

    // try to stress-test reentrancy of memoization
    auto rng = ranges::views::iota(31, 100);
    sequant::for_each(rng, [](const auto& i) {
      REQUIRE_NOTHROW(sequant::to_string(sequant::factorial(i)));
      REQUIRE(sequant::to_string(sequant::factorial(30)) ==
              "265252859812191058636308480000000");
    });

    // 100! has been memoized by now
    REQUIRE(sequant::to_string(sequant::factorial(100)) ==
            "933262154439441526816992388562667004907159682643816214685929638952"
            "175999932299156089414639761565182862536979208272237582511852109168"
            "64000000000000000000000000");
  }
}  // TEST_CASE("Factorial")
