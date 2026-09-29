//
// Created by Eduard Valeyev on 5/18/23.
//

#include <catch2/catch_test_macros.hpp>

#include "catch2_sequant.hpp"

#include <SeQuant/core/math.hpp>
#include <SeQuant/core/meta.hpp>
#include <SeQuant/core/rational.hpp>
#include <SeQuant/core/runtime.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <cmath>
#include <new>
#include <numbers>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <range/v3/view/iota.hpp>

TEST_CASE("math", "[elements]") {
  using namespace sequant;

  SECTION("rational") {
    [[maybe_unused]] auto print = [](rational r) {
      return sequant::to_wstring(numerator(r)) + L"/" +
             sequant::to_wstring(denominator(r));
    };
    SECTION("to_rational") {
      REQUIRE(to_rational(1. / 3) == rational{1, 3});
      REQUIRE(to_rational(1. / 3, 0.) ==
              rational{6004799503160661ull, 18014398509481984ull});
      REQUIRE(to_rational(1. / 7) == rational{1, 7});
      REQUIRE(to_rational(std::numbers::pi) == rational{99023, 31520});
      REQUIRE(to_rational(std::numbers::e) == rational{23225, 8544});
      REQUIRE_THROWS_AS(to_rational(std::nan("NaN")), Exception);
    }
  }

  SECTION("factorial") {
    REQUIRE(sequant::to_string(sequant::factorial(30)) ==
            "265252859812191058636308480000000");
    // 21! has been memoized by now
    REQUIRE(sequant::to_string(sequant::factorial(21)) ==
            "51090942171709440000");

    // try to stress-test reentrancy of memoization; sequant::for_each is
    // parallel and Catch2's assertion macros are not thread-safe, so the loop
    // only records per-position results and the assertions run below
    constexpr int first = 31;
    constexpr int last = 100;
    struct Recorded {
      std::string value;
      std::string memoized;
      std::string error;
    };
    std::vector<Recorded> recorded(static_cast<std::size_t>(last - first));
    auto rng = ranges::views::iota(first, last);
    sequant::for_each(rng, [&recorded](const auto& i) {
      // recorded holds one slot per value of rng, so operator[] is in range and
      // nothing in the worker can throw past this lambda.
      auto& rec = recorded[static_cast<std::size_t>(i - first)];
      try {
        rec.value = sequant::to_string(sequant::factorial(i));
        rec.memoized = sequant::to_string(sequant::factorial(30));
      } catch (const std::exception& e) {
        rec.error = e.what();
      }
    });
    for (std::size_t pos = 0; pos != recorded.size(); ++pos) {
      CAPTURE(first + static_cast<int>(pos));
      REQUIRE(recorded[pos].error.empty());
      REQUIRE(!recorded[pos].value.empty());
      REQUIRE(recorded[pos].memoized == "265252859812191058636308480000000");
    }

    // 100! has been memoized by now
    REQUIRE(sequant::to_string(sequant::factorial(100)) ==
            "933262154439441526816992388562667004907159682643816214685929638952"
            "175999932299156089414639761565182862536979208272237582511852109168"
            "64000000000000000000000000");
  }
}
