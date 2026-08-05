//
// Created by Eduard Valeyev on 2023-12-06
//

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/logger.hpp>
#include <SeQuant/core/utility/indices.hpp>
#include <SeQuant/core/utility/timer.hpp>
#include <SeQuant/domain/mbpt/models/cc.hpp>

#include <catch2/catch_test_macros.hpp>
#include "catch2_sequant.hpp"

TEST_CASE("mbpt_cc", "[mbpt/cc][valgrind_skip]") {
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

    SECTION("λ") {
      SECTION("CCSD λ") {
        const auto N = 2;
        auto cc = CC{N};
        auto λ_eqs = cc.λ();
        REQUIRE(λ_eqs.size() == N + 1);
        // element 0 is the λ pseudoenergy: the CC energy with T → Λ⁺, i.e. the
        // CC energy expression with the cluster amplitudes replaced by Λ⁺
        REQUIRE_THAT(λ_eqs[0],
                     EquivalentTo(L"f{i_1;a_1}:A-C-S λ⁺{a_1;i_1}:A-N-S "
                                  L"+ 1/4 g{i_1,i_2;a_1,a_2}:A-C-S "
                                  L"λ⁺{a_1,a_2;i_1,i_2}:A-N-S "
                                  L"+ 1/2 g{i_1,i_2;a_1,a_2}:A-C-S "
                                  L"λ⁺{a_1;i_1}:A-N-S λ⁺{a_2;i_2}:A-N-S"));
      }
    }  // SECTION("λ")
  }

  SECTION("energy") {
    // CC::energy() must equal the p==0 element of CC::t() for both ansätze.
    const auto N = 2;
    SECTION("TCC energy == t()[0]") {
      REQUIRE_THAT(CC{N}.energy(), EquivalentTo(CC{N}.t().at(0)));
    }
    SECTION("UCC energy == t()[0]") {
      const CC::Options opts{.ansatz = CC::Ansatz::U, .hbar_comm_rank = 2};
      REQUIRE_THAT(CC(N, opts).energy(), EquivalentTo(CC(N, opts).t().at(0)));
      // explicit rank override matches an engine built at that rank
      REQUIRE_THAT(
          CC(N, opts).energy(3),
          EquivalentTo(
              CC(N, {.ansatz = CC::Ansatz::U, .hbar_comm_rank = 3}).t().at(0)));
    }
  }  // SECTION("energy")

  SECTION("with_hbar_comm_rank") {
    const auto N = 2;
    const CC::Options opts{.ansatz = CC::Ansatz::U, .hbar_comm_rank = 2};
    auto bch2 = CC(N, opts);
    auto bch3 = bch2.with_hbar_comm_rank(3);
    REQUIRE(bch2.hbar_comm_rank() == 2);
    REQUIRE(bch3.hbar_comm_rank() == 3);
    // rank 0 is a valid truncation (H̄ = H); only CC::λ rejects it, since it
    // derives at rank - 1
    REQUIRE(bch2.with_hbar_comm_rank(0).hbar_comm_rank() == 0);
    if (sequant::assert_behavior() == sequant::AssertBehavior::Throw) {
      // the ctor is the only check
      REQUIRE_THROWS_AS(CC(N, {.ansatz = CC::Ansatz::U}), Exception);
      REQUIRE_THROWS_AS(CC(N, {.hbar_comm_rank = 0}).λ(), Exception);
    }
  }  // SECTION("with_hbar_comm_rank")

  SECTION("rdm") {
    constexpr auto N = 2;
    // Traditional CCSD: γ = <0|(1+Λ) e^{-T} {ã^{p_1..}_{p_{r+1}..}} e^{T}|0>.
    SECTION("TCC RDM") {
      const auto cc = CC{N};
      const auto g = cc.rdm(1);
      const auto G = cc.rdm(2);
      REQUIRE_THAT(
          g,
          EquivalentTo(
              L"t{a_1;i_1}:A-N-S * δ{i_1;p_1}:N-C-S * δ{p_2;a_1}:N-C-S "
              L"+ λ{i_1;a_1}:A-N-S * δ{a_1;p_1}:N-C-S * δ{p_2;i_1}:N-C-S "
              L"- 1/2 t{a_1,a_2;i_1,i_2}:A-N-S * δ{i_2;p_1}:N-C-S * "
              L"δ{p_2;i_3}:N-C-S * λ{i_1,i_3;a_1,a_2}:A-N-S "
              L"- 1 t{a_1;i_1}:A-N-S * t{a_2;i_2}:A-N-S * λ{i_1;a_2}:A-N-S * "
              L"δ{p_2;a_1}:N-C-S * δ{i_2;p_1}:N-C-S "
              L"+ 1/2 t{a_1,a_3;i_1,i_2}:A-N-S * δ{a_2;p_1}:N-C-S * "
              L"δ{p_2;a_3}:N-C-S * λ{i_1,i_2;a_1,a_2}:A-N-S "
              L"- 1 t{a_1;i_1}:A-N-S * λ{i_2;a_1}:A-N-S * δ{i_1;p_1}:N-C-S * "
              L"δ{p_2;i_2}:N-C-S "
              L"- 1/2 t{a_3;i_3}:A-N-S * t{a_1,a_2;i_1,i_2}:A-N-S * "
              L"δ{i_2;p_1}:N-C-S * δ{p_2;a_3}:N-C-S * λ{i_1,i_3;a_1,a_2}:A-N-S "
              L"- 1/2 t{a_2;i_3}:A-N-S * t{a_1,a_3;i_1,i_2}:A-N-S * "
              L"δ{i_3;p_1}:N-C-S * δ{p_2;a_3}:N-C-S * λ{i_1,i_2;a_1,a_2}:A-N-S "
              L"+ t{a_1;i_1}:A-N-S * λ{i_1;a_2}:A-N-S * δ{p_2;a_1}:N-C-S * "
              L"δ{a_2;p_1}:N-C-S "
              L"+ t{a_1,a_2;i_1,i_2}:A-N-S * λ{i_2;a_2}:A-N-S * "
              L"δ{i_1;p_1}:N-C-S * δ{p_2;a_1}:N-C-S"));
      REQUIRE(size(G) == 94);
      // 94 terms is not worth spelling out, so check number of terms and
      // external indices
      const auto ext = get_unique_indices(G);
      REQUIRE(std::ranges::is_permutation(
          ext.ket, container::svector<Index>{L"p_1", L"p_2"}));
      REQUIRE(std::ranges::is_permutation(
          ext.bra, container::svector<Index>{L"p_3", L"p_4"}));
      REQUIRE(ext.aux.empty());
      // an explicit comm_rank truncates early: cutting Γ at the 2nd nested
      // commutator drops the terms the default (4th, exact) picks up
      REQUIRE(size(cc.rdm(2, 2)) == 71);  // 23 terms short of the exact 94
    }
    // raising comm_rank past the default must not add terms; the defaults
    // here are 2 and 4, so these probe three and two orders beyond
    SECTION("TCC RDM default comm_rank") {
      const auto cc = CC{N};
      REQUIRE_THAT(cc.rdm(1), EquivalentTo(cc.rdm(1, 5)));
      REQUIRE_THAT(cc.rdm(2), EquivalentTo(cc.rdm(2, 6)));
    }
    // Unitary CCSD (σ = T − T⁺, no Λ): γ = <0| e^{-σ} {ã^{p_1}_{p_2}} e^{σ}|0>.
    SECTION("UCC RDM") {
      const auto cc = CC(N, {.ansatz = CC::Ansatz::U, .hbar_comm_rank = 2});
      const auto g = cc.rdm(1);
      const auto G = cc.rdm(2);
      REQUIRE_THAT(
          g,
          EquivalentTo(
              L"t{a_1;i_1}:A-N-S * δ{i_1;p_1}:N-C-S * δ{p_2;a_1}:N-C-S "
              L"+ δ{a_1;p_1}:N-C-S * δ{p_2;i_1}:N-C-S * t⁺{i_1;a_1}:A-N-S "
              L"+ t{a_2;i_1}:A-N-S * δ{a_1;p_1}:N-C-S * δ{p_2;a_2}:N-C-S * "
              L"t⁺{i_1;a_1}:A-N-S "
              L"- 1/2 t{a_1,a_2;i_1,i_2}:A-N-S * δ{i_1;p_1}:N-C-S * "
              L"δ{p_2;a_2}:N-C-S * t⁺{i_2;a_1}:A-N-S "
              L"- 1 t{a_1;i_1}:A-N-S * δ{i_1;p_1}:N-C-S * δ{p_2;i_2}:N-C-S * "
              L"t⁺{i_2;a_1}:A-N-S "
              L"+ 1/2 t{a_1,a_3;i_1,i_2}:A-N-S * δ{a_2;p_1}:N-C-S * "
              L"δ{p_2;a_3}:N-C-S * t⁺{i_1,i_2;a_1,a_2}:A-N-S "
              L"- 1/2 t{a_1,a_2;i_1,i_2}:A-N-S * δ{i_2;p_1}:N-C-S * "
              L"δ{p_2;i_3}:N-C-S * t⁺{i_1,i_3;a_1,a_2}:A-N-S "
              L"- 1/2 t{a_2;i_1}:A-N-S * δ{a_1;p_1}:N-C-S * δ{p_2;i_2}:N-C-S * "
              L"t⁺{i_1,i_2;a_1,a_2}:A-N-S"));
      REQUIRE(size(G) == 24);
    }
  }  // SECTION("rdm")

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
      auto t_eqs = CC(N, {.ansatz = ansatz, .hbar_comm_rank = c}).t();
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
