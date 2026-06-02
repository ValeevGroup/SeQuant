#include <catch2/catch_test_macros.hpp>

#include <SeQuant/core/asy_cost.hpp>
#include <SeQuant/core/rational.hpp>
#include <SeQuant/core/space.hpp>

#include "catch2_sequant.hpp"

#include <cmath>
#include <cstddef>
#include <string>

namespace {

// Numerical reference for ops() given two-space costs.
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

// Two distinct index spaces used to build occupied/virtual costs. The type
// bits make `virt_space` order above `occ_space`, so the virtual exponent
// dominates cost comparison (and base_key drives how each space prints).
sequant::IndexSpace const occ_space{L"O", 0b01};
sequant::IndexSpace const virt_space{L"V", 0b10};

sequant::AsyCost::ExponentMap ov(std::size_t nocc, std::size_t nvirt) {
  sequant::AsyCost::ExponentMap m;
  if (nocc > 0) m.emplace(occ_space, nocc);
  if (nvirt > 0) m.emplace(virt_space, nvirt);
  return m;
}

}  // namespace

TEST_CASE("asy_cost", "[AsyCost]") {
  using sequant::AsyCost;
  using sequant::rational;

  SECTION("to_text") {
    REQUIRE(AsyCost{ov(0, 0)}.text() == "0");

    REQUIRE(AsyCost{}.text() == "0");

    REQUIRE(AsyCost{ov(1, 0)}.text() == "O");

    REQUIRE(AsyCost{ov(0, 1)}.text() == "V");

    REQUIRE(AsyCost{ov(1, 1)}.text() == "OV");

    REQUIRE(AsyCost{ov(2, 1)}.text() == "O^2V");

    REQUIRE(AsyCost{ov(1, 2)}.text() == "OV^2");

    REQUIRE(AsyCost{ov(2, 2)}.text() == "O^2V^2");

    auto c = AsyCost{ov(2, 2)} + AsyCost{ov(3, 2)} + AsyCost{ov(2, 3)} +
             AsyCost{ov(3, 3)};
    REQUIRE(c.text() == "O^3V^3 + O^2V^3 + O^3V^2 + O^2V^2");

    c = AsyCost{ov(1, 1)} - AsyCost{ov(2, 3)} + AsyCost{ov(2, 2)};
    REQUIRE(c.text() == "- O^2V^3 + O^2V^2 + OV");

    REQUIRE(AsyCost{ov(1, 1), 20}.text() == "20*OV");

    REQUIRE(AsyCost{ov(0, 0)} == AsyCost::zero());

    REQUIRE(AsyCost{ov(1, 1), 0} == AsyCost::zero());
  }

  SECTION("Comparisons") {
    auto const c1 = AsyCost{ov(0, 0)};
    auto const c2 = AsyCost{ov(0, 1)};
    auto const c3 = AsyCost{ov(0, 2)};
    auto const c4 = AsyCost{ov(0, 2)};

    REQUIRE(c1 == AsyCost::zero());

    REQUIRE(c1 < c2);
    REQUIRE_FALSE(c1 == c2);

    REQUIRE(c2 < c3);
    REQUIRE_FALSE(c2 == c3);

    REQUIRE(c1 < c3);
    REQUIRE_FALSE(c1 == c3);

    REQUIRE(c3 == c4);

    auto const cc1 = AsyCost{ov(4, 1)} + AsyCost{ov(3, 2)} + AsyCost{ov(4, 2)};
    auto const cc2 =
        AsyCost{ov(2, 2)} + (AsyCost{ov(2, 4)} + AsyCost{ov(2, 4)});
    REQUIRE(cc1 < cc2);
    REQUIRE_FALSE(cc2 < cc1);
    REQUIRE(cc2 > cc1);
  }

  SECTION("Addition and subtraction") {
    auto const c1 = AsyCost{ov(0, 0)};
    auto const c2 = AsyCost{ov(0, 1)};
    auto const c3 = AsyCost{ov(1, 2)};
    REQUIRE(c1 + c2 == c2);
    REQUIRE(c1 - c2 == -1 * c2);
    REQUIRE((c2 + c3).text() == "OV^2 + V");
    REQUIRE((c2 - c3).text() == "- OV^2 + V");
    REQUIRE(c3 + c3 == 2 * c3);
  }

  SECTION("Assignments") {
    auto c1 = AsyCost{ov(0, 0)};
    auto c2 = AsyCost{ov(0, 1)};
    auto c3 = AsyCost{ov(1, 2)};
    c2 += c1;
    REQUIRE(c2.text() == "V");
    c2 -= c1;
    REQUIRE(c2.text() == "V");
    c2 += c3;
    REQUIRE(c2.text() == "OV^2 + V");
    c2 -= c3;
    REQUIRE(c2.text() == "V");
    c2 -= c3;
    REQUIRE(c2.text() == "- OV^2 + V");
  }

  SECTION("Ops count") {
    const std::size_t nocc = 2;
    const std::size_t nvirt = 20;

    AsyCost::ExtentMap const ext{{occ_space, nocc}, {virt_space, nvirt}};
    auto flops =
        MatFlops{static_cast<double>(nocc), static_cast<double>(nvirt)};
    REQUIRE(AsyCost::zero().ops(ext) == 0);

    REQUIRE(AsyCost{ov(2, 3)}.ops(ext) == flops(2, 3));

    auto const cost = AsyCost{ov(3, 1)} + AsyCost{ov(2, 1)};
    REQUIRE(cost.ops(ext) == flops(3, 1) + flops(2, 1));
  }

  SECTION("Fractional costs") {
    auto c0 = AsyCost{ov(2, 4), rational{1, 2}};
    REQUIRE(c0.text() == "1/2*O^2V^4");

    auto const c1 = AsyCost{ov(1, 2)} * rational{2, 3};
    REQUIRE(c1.text() == "2/3*OV^2");

    auto const c2 = AsyCost{ov(1, 2)} / rational{2, 3};
    REQUIRE(c2.text() == "3/2*OV^2");

    auto const c3 = (AsyCost{ov(1, 2)} + AsyCost{ov(2, 4)}) * 2;
    REQUIRE(c3.text() == "2*O^2V^4 + 2*OV^2");
  }

  SECTION("LaTeX") {
    auto cost = AsyCost{ov(2, 3), rational{1, 4}};
    REQUIRE(cost.to_latex() == L"\\frac{1}{4}O^{2}V^{3}");
    cost = AsyCost{ov(2, 3)};
    REQUIRE(cost.to_latex() == L"O^{2}V^{3}");
    cost = AsyCost{ov(2, 3), rational{1, 1}};
    REQUIRE(cost.to_latex() == L"O^{2}V^{3}");
    cost = AsyCost{ov(2, 3), rational{-1, 1}};
    REQUIRE(cost.to_latex() == L"- O^{2}V^{3}");
    cost = AsyCost{ov(2, 3), rational{-1, 4}};
    REQUIRE(cost.to_latex() == L"- \\frac{1}{4}O^{2}V^{3}");
  }

  SECTION("generalized spaces") {
    using sequant::IndexSpace;
    // Exercise > 2 spaces, including one whose base_key is multi-character.
    // Type bits set the ordering O < U < V < K < Q.
    IndexSpace const O{L"O", 0b00001};
    IndexSpace const U{L"U", 0b00010};
    IndexSpace const V{L"V", 0b00100};
    IndexSpace const K{L"K", 0b01000};
    IndexSpace const Q{L"(qq)", 0b10000};

    AsyCost::ExponentMap m;
    m.emplace(O, 2);
    m.emplace(U, 1);
    m.emplace(V, 3);
    m.emplace(K, 2);
    m.emplace(Q, 1);
    auto const c = AsyCost{m};

    auto const tex = c.to_latex();
    REQUIRE(tex.find(L"O^{2}") != std::wstring::npos);
    REQUIRE(tex.find(L"V^{3}") != std::wstring::npos);
    REQUIRE(tex.find(L"K^{2}") != std::wstring::npos);
    REQUIRE(tex.find(L"(qq)") != std::wstring::npos);
    REQUIRE(tex.find(L"U^{") == std::wstring::npos);     // exp 1: no caret
    REQUIRE(tex.find(L"(qq)^{") == std::wstring::npos);  // exp 1: no caret

    // Equal costs constructed in two different orders compare equal.
    AsyCost::ExponentMap m2;
    m2.emplace(K, 2);
    m2.emplace(Q, 1);
    m2.emplace(V, 3);
    m2.emplace(O, 2);
    m2.emplace(U, 1);
    REQUIRE(AsyCost{m2} == c);

    // Numerical evaluation with explicit extents.
    AsyCost::ExtentMap const ext{{O, 10}, {V, 100}, {U, 5}, {K, 500}, {Q, 7}};
    auto const expected = std::pow(10.0, 2) * std::pow(5.0, 1) *
                          std::pow(100.0, 3) * std::pow(500.0, 2) *
                          std::pow(7.0, 1);
    REQUIRE(c.ops(ext) == expected);

    // Missing extent for K falls back to 1.
    AsyCost::ExtentMap const ext2{{O, 10}, {V, 100}, {U, 5}, {Q, 7}};
    auto const expected2 = std::pow(10.0, 2) * std::pow(5.0, 1) *
                           std::pow(100.0, 3) * 1.0 * std::pow(7.0, 1);
    REQUIRE(c.ops(ext2) == expected2);

    // Ordering follows IndexSpace priority (highest space dominates): with
    // K > V > U > O, a larger K-exponent makes a cost larger regardless of
    // the lower spaces.
    AsyCost const heavy_k{AsyCost::ExponentMap{{O, 4}, {K, 3}}};
    AsyCost const light_k{AsyCost::ExponentMap{{O, 9}, {K, 2}}};
    REQUIRE(light_k < heavy_k);
  }
}
