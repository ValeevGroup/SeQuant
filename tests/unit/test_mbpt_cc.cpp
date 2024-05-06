//
// Created by Eduard Valeyev on 2023-12-06
//

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/timer.hpp>
#include <SeQuant/domain/mbpt/models/cc.hpp>

#include <catch2/catch_test_macros.hpp>
#include "test_config.hpp"

TEST_CASE("SR-TCC", "[mbpt/cc]") {
  using namespace sequant::mbpt::sr;

  SECTION("t") {
    // TCC R1
    SEQUANT_PROFILE_SINGLE("CCSD t", {
      const auto N = 2;
      auto t_eqs = CC{N}.t();
      REQUIRE(t_eqs.size() == N + 1);
      for (auto k = 0; k <= N; ++k) REQUIRE(t_eqs[k]);
      if (N == 2) {
        REQUIRE(size(t_eqs[0]) == 3);
        REQUIRE(size(t_eqs[1]) == 14);
      }
    });

  }  // SECTION("t")

}  // TEST_CASE("SR-TCC")

TEST_CASE("SR-UCC", "[mbpt/cc]") {
  using namespace sequant::mbpt::sr;

  SECTION("t") {
    const auto N = 2;
    const std::size_t C = 3;
    CC::Ansatz ansatz = CC::U;

    // oUCC energy, up to third commutator
    auto t_eqs = CC{N, ansatz}.t(C);
    REQUIRE(t_eqs.size() == N + 1);
    for (auto k = 0; k <= N; ++k) {
      REQUIRE(t_eqs[k]);
    }

    // std::wcout << to_latex_align(t_eqs[0], 5, 3) << std::endl;
    if (C == 3 && ansatz == CC::Ansatz::U) {
      REQUIRE(size(t_eqs[0]) == 56);
    }

  }  // SECTION("t")

}  // TEST_CASE("SR-TCC")
