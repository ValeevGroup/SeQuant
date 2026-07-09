//
// Created by Eduard Valeyev on 2023-12-06
//

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/logger.hpp>
#include <SeQuant/core/rational.hpp>
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

  SECTION("hbar_comm_rank") {
    const auto N = 2;
    REQUIRE(CC(N, {.ansatz = CC::Ansatz::U, .hbar_comm_rank = 2})
                .hbar_comm_rank() == 2);
    // rank 0 is a valid truncation (H̄ = H); only CC::λ rejects it, since it
    // derives at rank - 1
    if (sequant::assert_behavior() == sequant::AssertBehavior::Throw) {
      // the ctor is the only check
      REQUIRE_THROWS_AS(CC(N, {.ansatz = CC::Ansatz::U}), Exception);
      REQUIRE_THROWS_AS(CC(N, {.hbar_comm_rank = 0}).λ(), Exception);
    }
  }  // SECTION("hbar_comm_rank")

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
    // regression: op::Λ used to add to a null ExprPtr when singles are skipped
    SECTION("TCC RDM skip_singles") {
      const auto ccd = CC(N, {.skip_singles = true});
      // the doubles-only subset of the CCSD γ pinned above
      REQUIRE_THAT(
          ccd.rdm(1),
          EquivalentTo(L"- 1/2 t{a_1,a_2;i_1,i_2}:A-N-S * δ{i_2;p_1}:N-C-S * "
                       L"δ{p_2;i_3}:N-C-S * λ{i_1,i_3;a_1,a_2}:A-N-S "
                       L"+ 1/2 t{a_1,a_3;i_1,i_2}:A-N-S * δ{a_2;p_1}:N-C-S * "
                       L"δ{p_2;a_3}:N-C-S * λ{i_1,i_2;a_1,a_2}:A-N-S"));
      // Ansatz::oT forces skip_singles, so it reaches the same path
      REQUIRE_THAT(
          ccd.rdm(1),
          EquivalentTo(
              CC(N, {.ansatz = CC::Ansatz::oT, .skip_singles = true}).rdm(1)));
    }
  }  // SECTION("rdm")

  SECTION("extra singles commutators") {
    // hbar_singles_comm_rank wraps H̄_R with a factored t1 similarity transform
    // to order K:  H̄ = Σ_a (1/a!) [H̄_R, t1−t1†]_a.
    const auto N = 2;
    const CC::Options base{.ansatz = CC::Ansatz::U, .hbar_comm_rank = 2};
    const auto hbar_R = CC(N, base).hbar();  // K unset == the R=2 baseline

    // anti-Hermitian singles generator σ1 = t1 − t1†
    const auto sigma1 = op::t(1) - adjoint(op::t(1));
    // commutator built independently of mbpt::lst (plain expression arithmetic)
    const auto comm = [](const ExprPtr& A, const ExprPtr& B) {
      auto r = A * B - B * A;
      simplify(r);
      return r;
    };

    SECTION("K=0/unset reproduces plain hbar") {
      CC::Options k0 = base;
      k0.hbar_singles_comm_rank = 0;
      REQUIRE_THAT(CC(N, k0).hbar(), EquivalentTo(hbar_R));
    }

    SECTION("K=1: extra term is exactly [H̄_R, σ1]") {
      CC::Options k1 = base;
      k1.hbar_singles_comm_rank = 1;
      const auto hbar_k1 = CC(N, k1).hbar();

      // the wrap genuinely adds terms not present in the baseline
      const auto extra = comm(hbar_R, sigma1);
      REQUIRE(size(extra) > 0);
      REQUIRE(size(hbar_k1) > size(hbar_R));

      REQUIRE_THAT(hbar_k1, EquivalentTo(hbar_R + extra));
    }

    SECTION("K=2: increment is (1/2!)[[H̄_R, σ1], σ1]") {
      CC::Options k2 = base;
      k2.hbar_singles_comm_rank = 2;
      const auto expected =
          hbar_R + comm(hbar_R, sigma1) +
          ex<Constant>(rational{1, 2}) * comm(comm(hbar_R, sigma1), sigma1);
      REQUIRE_THAT(CC(N, k2).hbar(), EquivalentTo(expected));
    }
  }  // SECTION("extra singles commutators")

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
