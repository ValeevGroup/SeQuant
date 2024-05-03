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

TEST_CASE("EOM-CC", "[mbpt/cc]") {
  using namespace sequant::mbpt::sr;

  SECTION("EOM-CCSD") {
    SEQUANT_PROFILE_SINGLE("EE-EOM-CCSD R", {
      const auto N = 2;
      const auto K_occ = 2;
      const auto K_uocc = 2;
      const auto eqs = CC{N}.R(K_occ, K_uocc);
      for (auto k = 1; k < eqs.size(); ++k) REQUIRE(eqs[k]);

      if (N == 2 && K_occ == 2 && K_uocc == 2) {
        REQUIRE(size(eqs[1]) == 21);
        REQUIRE(size(eqs[2]) == 53);
      }
    });

    SEQUANT_PROFILE_SINGLE("EE-EOM-CCSD L", {
      const auto N = 2;
      const auto K_occ = 2;
      const auto K_uocc = 2;
      const auto eqs = CC{N}.L(K_occ, K_uocc);
      for (auto k = 1; k < eqs.size(); ++k) REQUIRE(eqs[k]);

      if (N == 2 && K_occ == 2 && K_uocc == 2) {
        REQUIRE(size(eqs[1]) == 43);
        REQUIRE(size(eqs[2]) == 31);
      }
    });
  }  // SECTION("EOM-CCSD")

  SECTION("EOM-CCSDT") {
    SEQUANT_PROFILE_SINGLE("EE-EOM-CCSDT R", {
      const auto N = 3;
      const auto K_occ = 3;
      const auto K_uocc = 3;
      const auto eqs = CC{N}.R(K_occ, K_uocc);
      for (auto k = 1; k < eqs.size(); ++k) REQUIRE(eqs[k]);

      if (N == 3 && K_occ == 3 && K_uocc == 3) {
        REQUIRE(size(eqs[1]) == 22);
        REQUIRE(size(eqs[2]) == 62);
        REQUIRE(size(eqs[3]) == 99);
      }
    });
  }  // SECTION("EOM-CCSDT")

}  // TEST_CASE("EOM-CC")

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
