//
// Created by Eduard Valeyev on 2023-12-06
//

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/logger.hpp>
#include <SeQuant/core/utility/timer.hpp>
#include <SeQuant/domain/mbpt/models/cc.hpp>

#include <catch2/catch_test_macros.hpp>
#include "catch2_sequant.hpp"

TEST_CASE("mbpt_cc", "[mbpt/cc]") {
  using namespace sequant;
  using namespace sequant::mbpt;

  SECTION("sr_tcc") {
    SECTION("t") {
      // TCC R1
      SECTION("CCSD t") {
        [[maybe_unused]] auto& l = sequant::Logger::instance();
        // l.canonicalize = true;
        const auto N = 2;
        auto t_eqs = CC{N}.t();
        REQUIRE(t_eqs.size() == N + 1);
        for (auto k = 0; k <= N; ++k) REQUIRE(t_eqs[k]);
        if (N == 2) {
          REQUIRE(size(t_eqs[0]) == 3);
          REQUIRE(size(t_eqs[1]) == 14);
        }
      }

    }  // SECTION("t")
  }

  SECTION("eom_cc"){SECTION("EOM-CCSD"){const auto N = 2;
  auto cc = CC{N};
  SECTION("EE-EOM-CCSD R") {
    const auto np = 2;
    const auto nh = 2;
    const auto eqs = cc.eom_r(nₚ(np), nₕ(nh));
    for (auto k = 1; k < eqs.size(); ++k) REQUIRE(eqs[k]);

    REQUIRE(size(eqs[1]) == 21);
    REQUIRE(size(eqs[2]) == 53);
  }

  SECTION("IP-EOM-CCSD R") {
    const auto np = 1;
    const auto nh = 2;
    const auto eqs = cc.eom_r(nₚ(np), nₕ(nh));
    for (auto k = 0; k < eqs.size(); ++k) REQUIRE(eqs[k]);

    REQUIRE(size(eqs[0]) == 9);
    REQUIRE(size(eqs[1]) == 32);
  }

  SECTION("EA-EOM-CCSD R") {
    const auto np = 2;
    const auto nh = 1;
    const auto eqs = cc.eom_r(nₚ(np), nₕ(nh));
    for (auto k = 0; k < eqs.size(); ++k) REQUIRE(eqs[k]);

    REQUIRE(size(eqs[0]) == 9);
    REQUIRE(size(eqs[1]) == 32);
  }

  SECTION("EE-EOM-CCSD L") {
    const auto np = 2;
    const auto nh = 2;
    const auto eqs = cc.eom_l(nₚ(np), nₕ(nh));
    for (auto k = 1; k < eqs.size(); ++k) REQUIRE(eqs[k]);

    REQUIRE(size(eqs[1]) == 43);
    REQUIRE(size(eqs[2]) == 31);
  }
}  // SECTION("EOM-CCSD")

#ifndef SEQUANT_SKIP_LONG_TESTS
SECTION("EOM-CCSDT") {
  const auto N = 3;
  auto cc = CC{N};
  SECTION("EE-EOM-CCSDT R") {
    const auto np = 3;
    const auto nh = 3;
    const auto eqs = cc.eom_r(nₚ(np), nₕ(nh));
    for (auto k = 1; k < eqs.size(); ++k) REQUIRE(eqs[k]);

    REQUIRE(size(eqs[1]) == 22);
    REQUIRE(size(eqs[2]) == 62);
    REQUIRE(size(eqs[3]) == 99);
  }
}  // SECTION("EOM-CCSDT")
#endif
}

#ifndef SEQUANT_SKIP_LONG_TESTS
SECTION("ucc") {
  SECTION("t") {
    const auto N = 2;
    const auto C = std::array{2, 3};  // commutator truncation ranks
    CC::Ansatz ansatz = CC::Ansatz::U;

    for (const auto& c : C) {
      // UCC energy, with commutator truncation rank = c
      auto t_eqs = CC(N, {.ansatz = ansatz, .hbar_truncation_rank = c}).t();
      REQUIRE(t_eqs.size() == N + 1);
      for (auto k = 0; k <= N; ++k) {
        REQUIRE(t_eqs[k]);
      }
      // these are numerically verified against http://arxiv.org/abs/2503.00617
      const auto energy_nterms = size(t_eqs[0]);
      if (c == 2) REQUIRE(energy_nterms == 20);
      if (c == 3) REQUIRE(energy_nterms == 74);
    }
  }  // SECTION("t")
}  // SECTION("ucc")
#endif
}
