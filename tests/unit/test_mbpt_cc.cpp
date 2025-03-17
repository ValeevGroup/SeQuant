//
// Created by Eduard Valeyev on 2023-12-06
//

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/timer.hpp>
#include <SeQuant/domain/mbpt/models/cc.hpp>

#include <catch2/catch_test_macros.hpp>
#include "catch2_sequant.hpp"
#include "test_config.hpp"

TEST_CASE("mbpt_cc", "[mbpt/cc]") {
  using namespace sequant;
  using namespace sequant::mbpt;

  SECTION("sr_tcc") {
    SECTION("t") {
      // TCC R1
      SEQUANT_PROFILE_SINGLE("CCSD t", {
        [[maybe_unused]] auto l = sequant::Logger::instance();
        // l->canonicalize = true;
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
  }

  SECTION("eom_cc") {
    SECTION("EOM-CCSD") {
      const auto N = 2;
      auto cc = CC{N};
      SEQUANT_PROFILE_SINGLE("EE-EOM-CCSD R", {
        const auto np = 2;
        const auto nh = 2;
        const auto eqs = cc.eom_r(nₚ(np), nₕ(nh));
        for (auto k = 1; k < eqs.size(); ++k) REQUIRE(eqs[k]);

        REQUIRE(size(eqs[1]) == 21);
        REQUIRE(size(eqs[2]) == 53);
      });

      SEQUANT_PROFILE_SINGLE("IP-EOM-CCSD R", {
        const auto np = 1;
        const auto nh = 2;
        const auto eqs = cc.eom_r(nₚ(np), nₕ(nh));
        for (auto k = 0; k < eqs.size(); ++k) REQUIRE(eqs[k]);

        REQUIRE(size(eqs[0]) == 9);
        REQUIRE(size(eqs[1]) == 32);
      });

      SEQUANT_PROFILE_SINGLE("EA-EOM-CCSD R", {
        const auto np = 2;
        const auto nh = 1;
        const auto eqs = cc.eom_r(nₚ(np), nₕ(nh));
        for (auto k = 0; k < eqs.size(); ++k) REQUIRE(eqs[k]);

        REQUIRE(size(eqs[0]) == 9);
        REQUIRE(size(eqs[1]) == 32);
      });

      SEQUANT_PROFILE_SINGLE("EE-EOM-CCSD L", {
        const auto np = 2;
        const auto nh = 2;
        const auto eqs = cc.eom_l(nₚ(np), nₕ(nh));
        for (auto k = 1; k < eqs.size(); ++k) REQUIRE(eqs[k]);

        REQUIRE(size(eqs[1]) == 43);
        REQUIRE(size(eqs[2]) == 31);
      });
    }  // SECTION("EOM-CCSD")

    SECTION("EOM-CCSDT") {
      const auto N = 3;
      auto cc = CC{N};
      SEQUANT_PROFILE_SINGLE("EE-EOM-CCSDT R", {
        const auto np = 3;
        const auto nh = 3;
        const auto eqs = cc.eom_r(nₚ(np), nₕ(nh));
        for (auto k = 1; k < eqs.size(); ++k) REQUIRE(eqs[k]);

        REQUIRE(size(eqs[1]) == 22);
        REQUIRE(size(eqs[2]) == 62);
        REQUIRE(size(eqs[3]) == 99);
      });
    }  // SECTION("EOM-CCSDT")
  }

  SECTION("ucc") {
    SECTION("t") {
      const auto N = 2;
      const std::size_t C = 3;
      CC::Ansatz ansatz = CC::Ansatz::U;

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
  }
}
