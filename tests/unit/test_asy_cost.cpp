#include "catch.hpp"

#include <SeQuant/core/asy_cost.hpp>
#include <SeQuant/core/rational.hpp>

#include <cstddef>
#include <cmath>
#include <string>

struct MatFlops {
  double occ_range_size;
  double virt_range_size;
  double operator()(unsigned short nocc, unsigned short nvirt) const {
    double ops = 1;
    if (nocc > 0) ops *= std::pow(occ_range_size, nocc);
    if (nvirt > 0) ops *= std::pow(virt_range_size, nvirt);
    return ops > 1 ? ops : 0;
  }
};

TEST_CASE("TEST ASY_COST", "[AsyCost]") {
  using sequant::AsyCost;
  using sequant::rational;

  SECTION("to_text") {
    REQUIRE(AsyCost{0, 0}.text() == "0");

    REQUIRE(AsyCost{}.text() == "0");

    REQUIRE(AsyCost{1, 0}.text() == "O");

    REQUIRE(AsyCost{0, 1}.text() == "V");

    REQUIRE(AsyCost{1, 1}.text() == "OV");

    REQUIRE(AsyCost{2, 1}.text() == "O^2V");

    REQUIRE(AsyCost{1, 2}.text() == "OV^2");

    REQUIRE(AsyCost{2, 2}.text() == "O^2V^2");

    auto c = AsyCost{2, 2} + AsyCost{3, 2} + AsyCost{2, 3} + AsyCost{3, 3};
    REQUIRE(c.text() == "O^3V^3 + O^2V^3 + O^3V^2 + O^2V^2");

    c = AsyCost{1, 1} - AsyCost{2, 3} + AsyCost{2, 2};
    REQUIRE(c.text() == "- O^2V^3 + O^2V^2 + OV");

    REQUIRE(AsyCost{20, 1, 1}.text() == "20*OV");

    REQUIRE(AsyCost{0, 0} == AsyCost::zero());

    REQUIRE(AsyCost{0, 1, 1} == AsyCost::zero());
  }

  SECTION("Comparisons") {
    auto const c1 = AsyCost{0, 0};
    auto const c2 = AsyCost{0, 1};
    auto const c3 = AsyCost{0, 2};
    auto const c4 = AsyCost{0, 2};

    REQUIRE(c1 == AsyCost::zero());

    REQUIRE(c1 < c2);
    REQUIRE(c1 <= c2);

    REQUIRE(c2 > c1);
    REQUIRE(c2 >= c1);

    REQUIRE(c2 < c3);
    REQUIRE(c3 > c2);

    REQUIRE(c3 == c4);

    auto const cc1 = AsyCost{4, 1} + AsyCost{3, 2} + AsyCost{4, 2};
    auto const cc2 = AsyCost{2, 2} + (AsyCost{2, 4} + AsyCost{2, 4});
    REQUIRE(cc1 < cc2);
    REQUIRE_FALSE(cc2 < cc1);
    REQUIRE(cc2 > cc1);
  }

  SECTION("Addition and subtraction") {
    auto const c1 = AsyCost{0, 0};
    auto const c2 = AsyCost{0, 1};
    auto const c3 = AsyCost{1, 2};
    REQUIRE(c1 + c2 == c2);
    REQUIRE(c1 - c2 == -1 * c2);
    REQUIRE((c2 + c3).text() == "OV^2 + V");
    REQUIRE((c2 - c3).text() == "- OV^2 + V");
    REQUIRE(c3 + c3 == 2 * c3);
  }

  SECTION("Assignments") {
    auto c1 = AsyCost{0, 0};
    auto c2 = AsyCost{0, 1};
    auto c3 = AsyCost{1, 2};
    c2 += c1;
    REQUIRE(c2.text() == "V");
    c2 -= c1;
    REQUIRE(c2.text() == "V");
    c2 += c3;
    REQUIRE(c2.text() == "OV^2 + V");
    c2 -= c3;
    REQUIRE(c2.text() == "V");
    c2 += c2;
    REQUIRE(c2.text() == "2*V");
  }

  SECTION("Ops count") {
    const size_t nocc = 2;
    const size_t nvirt = 20;

    auto flops = MatFlops{nocc, nvirt};
    REQUIRE(AsyCost::zero().ops(nocc, nvirt) == 0);

    REQUIRE(AsyCost{2, 3}.ops(nocc, nvirt) == flops(2, 3));

    auto const cost = AsyCost{3, 1} + AsyCost{2, 1};
    REQUIRE(cost.ops(nocc, nvirt) == flops(3, 1) + flops(2, 1));
  }

  SECTION("Fractional costs") {
    auto c0 = AsyCost{{1, 2}, 2, 4};
    REQUIRE(c0.text() == "1/2*O^2V^4");

    auto const c1 = AsyCost{1, 2} * rational{2, 3};
    REQUIRE(c1.text() == "2/3*OV^2");

    auto const c2 = AsyCost{1, 2} / rational{2, 3};
    REQUIRE(c2.text() == "3/2*OV^2");

    auto const c3 = (AsyCost{1, 2} + AsyCost{2, 4}) * 2;
    REQUIRE(c3.text() == "2*O^2V^4 + 2*OV^2");
  }

  SECTION("LaTeX") {
    auto cost = AsyCost{{1, 4}, 2, 3};
    REQUIRE(cost.to_latex() == L"\\frac{1}{4}O^{2}V^{3}");
    cost = AsyCost{2, 3};
    REQUIRE(cost.to_latex() == L"O^{2}V^{3}");
    cost = AsyCost{{1, 1}, 2, 3};
    REQUIRE(cost.to_latex() == L"O^{2}V^{3}");
    cost = AsyCost{{-1, 1}, 2, 3};
    REQUIRE(cost.to_latex() == L"- O^{2}V^{3}");
    cost = AsyCost{{-1, 4}, 2, 3};
    REQUIRE(cost.to_latex() == L"- \\frac{1}{4}O^{2}V^{3}");
  }
}
