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

#include <algorithm>
#include <cmath>
#include <new>
#include <numbers>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <range/v3/view/iota.hpp>

TEST_CASE("math", "[elements]") {
  using namespace sequant;

  SECTION("rational") {
    [[maybe_unused]] auto print = [](rational r) {
      return sequant::to_string(numerator(r)) + "/" +
             sequant::to_string(denominator(r));
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
  SECTION("iroot") {
    using sequant::iroot;
    REQUIRE(iroot(1, 5) == 1);
    REQUIRE(iroot(8, 1) == 8);
    REQUIRE(iroot(8, 3) == 2);
    REQUIRE(iroot(9, 2) == 3);
    REQUIRE(iroot(10, 2) == 3);
    REQUIRE(iroot(26, 3) == 2);
    REQUIRE(iroot(27, 3) == 3);
    REQUIRE(iroot(sequant::intmax_t{1} << 100, 10) == sequant::intmax_t{1}
                                                          << 10);
    REQUIRE(iroot((sequant::intmax_t{1} << 100) - 1, 10) ==
            (sequant::intmax_t{1} << 10) - 1);
    REQUIRE(iroot(5, 1000) == 1);
    if constexpr (sequant::assert_behavior() ==
                  sequant::AssertBehavior::Throw) {
      REQUIRE_THROWS_AS(iroot(8, 0), Exception);
      REQUIRE_THROWS_AS(iroot(-8, 3), Exception);
    }
    const sequant::intmax_t p3_500 = pow(sequant::intmax_t{3}, 500);
    REQUIRE(iroot(p3_500 - 1, 500) == 2);
    REQUIRE(iroot(p3_500, 500) == 3);
    REQUIRE(iroot(p3_500 + 1, 500) == 3);
    const sequant::intmax_t p1009_7 = pow(sequant::intmax_t{1009}, 7);
    REQUIRE(iroot(p1009_7 - 1, 7) == 1008);
    REQUIRE(iroot(p1009_7, 7) == 1009);
    const sequant::intmax_t big = (sequant::intmax_t{1} << 8000) - 1;
    REQUIRE(iroot(big, 2) == (sequant::intmax_t{1} << 4000) - 1);
    REQUIRE(iroot(big, 4000) == 3);
    REQUIRE(iroot(big + 2, 4000) == 4);
  }

  SECTION("coprime_base") {
    using sequant::coprime_base;
    auto sorted = [](container::svector<sequant::intmax_t> base) {
      std::sort(base.begin(), base.end());
      return base;
    };
    REQUIRE(sorted(coprime_base({12, 18})) ==
            container::svector<sequant::intmax_t>{2, 3});
    REQUIRE(sorted(coprime_base({6, 10, 15})) ==
            container::svector<sequant::intmax_t>{2, 3, 5});
    // elements need not be prime, only pairwise coprime
    REQUIRE(sorted(coprime_base({6, 35})) ==
            container::svector<sequant::intmax_t>{6, 35});
    REQUIRE(sorted(coprime_base({1, 1})).empty());
    const sequant::intmax_t p{1000003}, q{1000033};
    REQUIRE(sorted(coprime_base({p * q, p * p, q})) ==
            container::svector<sequant::intmax_t>{p, q});
  }

  SECTION("simplest_radical_form") {
    using sequant::simplest_radical_form;
    const auto half = rational{1, 2};
    const auto third = rational{1, 3};
    auto form = [](auto&&... powers) {
      return simplest_radical_form({powers...});
    };
    auto is = [](const std::optional<RadicalForm>& form, rational rational_part,
                 sequant::intmax_t radicand, sequant::intmax_t root_order) {
      return form && form->rational_part == rational_part &&
             form->radicand == radicand && form->root_order == root_order;
    };
    // only positive bases are admitted, whatever the assert behavior
    REQUIRE_THROWS_AS(form(std::pair{rational{-8}, third}), Exception);
    REQUIRE_THROWS_AS(form(std::pair{rational{0}, third}), Exception);
    // rational values have no root
    REQUIRE(is(form(std::pair{rational{8}, third}), 2, 1, 1));
    REQUIRE(is(form(std::pair{rational{9, 16}, rational{3, 2}}),
               rational{27, 64}, 1, 1));
    REQUIRE(is(form(std::pair{half, half}, std::pair{half, half}), half, 1, 1));
    // the rational part comes out of the radicand, including from a
    // denominator: 48^(1/2) = 4 3^(1/2), (1/2)^(1/2) = (1/2) 2^(1/2)
    REQUIRE(is(form(std::pair{rational{48}, half}), 4, 3, 2));
    REQUIRE(is(form(std::pair{half, half}), half, 2, 2));
    REQUIRE(is(form(std::pair{rational{2}, rational{3, 2}}), 2, 2, 2));
    REQUIRE(is(form(std::pair{rational{2}, -half}), half, 2, 2));
    // roots of different orders combine: 2^(1/2) 3^(1/3) = 72^(1/6)
    REQUIRE(
        is(form(std::pair{rational{2}, half}, std::pair{rational{3}, third}), 1,
           72, 6));
    // a power of a prime beyond trial division is still extracted, and a
    // product of large primes split by gcds: (p^2 3)^(1/2) = p 3^(1/2),
    // (p q)^(2/3) p^(1/3) q^(1/3) = p q
    const rational p{1000003}, q{1000033};
    REQUIRE(is(form(std::pair{p * p * 3, half}), p, 3, 2));
    // ... whatever the order of the root: (p^4)^(1/3) = p p^(1/3),
    // (p^3)^(1/2) = p p^(1/2)
    REQUIRE(is(form(std::pair{p * p * p * p, third}), p, numerator(p), 3));
    REQUIRE(is(form(std::pair{p * p * p, half}), p, numerator(p), 2));
    // the root search is bounded by the size of the cofactor relative to the
    // trial bound, which still admits every root a cofactor can have:
    // 1009 is the smallest prime above the bound, and the 7th root is found
    const sequant::intmax_t q1009{1009};
    REQUIRE(is(form(std::pair{rational{pow(q1009, 7)}, half}), pow(q1009, 3),
               1009, 2));
    // ... and is cheap for a large base with no small factor
    sequant::intmax_t big = (sequant::intmax_t{1} << 2000) + 1;
    for (bool coprime = false; !coprime; big += coprime ? 0 : 2) {
      coprime = true;
      for (unsigned d = 2; d <= 1000 && coprime; ++d) coprime = big % d != 0;
    }
    REQUIRE(is(form(std::pair{rational{big}, half}), 1, big, 2));
    // an integer exponent needs no factoring at all, so a base an order of
    // magnitude larger still folds at once
    const sequant::intmax_t huge = (sequant::intmax_t{1} << 30000) + 1;
    REQUIRE(
        is(form(std::pair{rational{huge}, rational{2}}), huge * huge, 1, 1));
    // ... and a base whose integer part of the exponent alone would exceed
    // the size cap is refused before any factoring
    REQUIRE(!form(std::pair{rational{(sequant::intmax_t{1} << 70000) + 1},
                            rational{3, 2}}));
    REQUIRE(is(form(std::pair{p * q, rational{2, 3}}, std::pair{p, third},
                    std::pair{q, third}),
               p * q, 1, 1));
    // a result too large to build is refused: exponents with coprime
    // denominators multiply the radicand's bits by their lcm, a large
    // integer exponent those of the rational part
    REQUIRE(!form(std::pair{rational{2}, rational{1, 65537}},
                  std::pair{rational{3}, rational{1, 65539}}));
    REQUIRE(is(form(std::pair{rational{2}, rational{1, 65537}}), 1, 2, 65537));
    REQUIRE(!form(std::pair{rational{1000}, rational{10000}}));
    REQUIRE(!form(std::pair{rational{1, 1000}, rational{10000}}));
    REQUIRE(!form(std::pair{rational{1000}, rational{100000, 3}}));
  }
}
