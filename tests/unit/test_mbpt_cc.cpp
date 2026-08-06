//
// Created by Eduard Valeyev on 2023-12-06
//

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/logger.hpp>
#include <SeQuant/core/utility/timer.hpp>
#include <SeQuant/domain/mbpt/bernoulli.hpp>
#include <SeQuant/domain/mbpt/models/cc.hpp>
#include <SeQuant/domain/mbpt/op.hpp>

#include <catch2/catch_test_macros.hpp>
#include "catch2_sequant.hpp"

TEST_CASE("mbpt_cc", "[mbpt/cc][valgrind_skip]") {
  using namespace sequant;
  using namespace sequant::mbpt;

  auto has_tensor = [](const ExprPtr& e, const std::wstring& label) {
    bool found = false;
    e->visit(
        [&](const ExprPtr& n) {
          if (n.is<Tensor>() && n.as<Tensor>().label() == label) found = true;
        },
        /*atoms_only=*/true);
    return found;
  };

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

  SECTION("bernoulli_wick") {
    using namespace sequant;
    using namespace sequant::mbpt;
    // [V, T2] is antisymmetric: [A,B] == -[B,A] after Wick reduction
    const auto V = op::tensor::h(2);
    const auto T2 = op::tensor::t(2);  // rank-2 excitation, tensor form
    const auto ab = bernoulli::detail::wick_commutator(V, T2);
    const auto ba = bernoulli::detail::wick_commutator(T2, V);
    REQUIRE_THAT(ab, EquivalentTo(simplify(ex<Constant>(-1) * ba)));
    // wick_reduce of a bare (already normal-ordered) operator is itself
    REQUIRE_THAT(bernoulli::detail::wick_reduce(V), EquivalentTo(V));
    // Wick reduction adds contractions beyond the naive V*T2 - T2*V
    REQUIRE(bernoulli::detail::wick_commutator(V, T2) != ex<Constant>(0));
    REQUIRE_THAT(bernoulli::detail::wick_commutator(V, T2),
                 !EquivalentTo(simplify(V * T2 - T2 * V)));
  }

  SECTION("bernoulli_expand_to_blocks") {
    using namespace sequant;
    using namespace sequant::mbpt;
    const auto V = op::tensor::h(2);  // general g
    const auto Vx = bernoulli::detail::expand_to_blocks(V);
    // identity on each manifold: the expansion changes no physical content
    for (const auto n : {1, 2})
      REQUIRE_THAT(op::tensor::ref_av(op::tensor::P(nₚ(n)) * Vx),
                   EquivalentTo(op::tensor::ref_av(op::tensor::P(nₚ(n)) * V)));
    REQUIRE(Vx.is<Sum>());
    REQUIRE(Vx.as<Sum>().size() > 1);
    REQUIRE_THAT(bernoulli::detail::expand_to_blocks(Vx),
                 EquivalentTo(Vx));  // idempotent
    // no general index survives: every residual index is occ or uocc
    auto isr = get_default_context().index_space_registry();
    Vx->visit(
        [&](const ExprPtr& n) {
          if (!n.is<NormalOperator<Statistics::FermiDirac>>()) return;
          for (const auto& o :
               n.as<NormalOperator<Statistics::FermiDirac>>().creann()) {
            const auto& sp = o.index().space();
            REQUIRE((isr->is_pure_occupied(sp) || isr->is_pure_unoccupied(sp)));
          }
        },
        /*atoms_only=*/true);
  }

  SECTION("bernoulli_N_R_split") {
    using namespace sequant;
    using namespace sequant::mbpt;
    const auto V = op::tensor::h(2);  // fluctuation potential g (general)
    const auto Vn = bernoulli::detail::N_part(V, 2);
    const auto Vr = bernoulli::detail::R_part(V, 2);
    // N ⊎ R reconstructs V. R stays in compact general-index form, so check the
    // identity on the manifolds rather than symbolically.
    const auto NR = simplify(Vn + Vr);
    for (const auto n : {1, 2})
      REQUIRE_THAT(op::tensor::ref_av(op::tensor::P(nₚ(n)) * NR),
                   EquivalentTo(op::tensor::ref_av(op::tensor::P(nₚ(n)) * V)));
    REQUIRE(Vn != ex<Constant>(0));
    REQUIRE(Vr != ex<Constant>(0));
    // N is idempotent; R has no pure-exc/deexc content
    REQUIRE_THAT(bernoulli::detail::N_part(Vn, 2), EquivalentTo(Vn));
    REQUIRE_THAT(bernoulli::detail::N_part(Vr, 2),
                 EquivalentTo(ex<Constant>(0)));
  }

  SECTION("bernoulli_hbar_structure") {
    using namespace sequant;
    using namespace sequant::mbpt;
    // Equation references are to 10.1063/1.5030344, Sec. III B.
    // The F-cancellation: F appears only in H̄¹, so rank r − rank r−1 is
    // F-free for r ≥ 2.
    auto h0 = bernoulli::hbar(2, 0, false);
    auto h1 = bernoulli::hbar(2, 1, false);
    auto h2 = bernoulli::hbar(2, 2, false);
    auto has_f = [&](const ExprPtr& e) { return has_tensor(e, L"f"); };
    REQUIRE(has_f(simplify(h1 - h0)));  // [F,σ]
    REQUIRE_FALSE(has_f(simplify(h2 - h1)));
    REQUIRE_THAT(h0,  // H̄⁰ = F + V, Eq. (46)
                 EquivalentTo(simplify(op::tensor::F() + op::tensor::h(2))));

    // Reference expectation values of H̄¹ and H̄², Eqs. (47) and (48), taken as
    // successive-rank differences.
    const auto E0 = op::tensor::ref_av(h0);
    const auto E1 = op::tensor::ref_av(h1);
    const auto E2 = op::tensor::ref_av(h2);
    const auto E1_contrib = simplify(E1 - E0);
    const auto E2_contrib = simplify(E2 - E1);

    // <0|H̄¹|0>: g-content is exactly 1/8 <ij||ab> σ_ij^ab + h.c.; the remainder
    // is the [F,σ] Brillouin terms, which vanish at RHF.
    const auto E1_g_closed = deserialize(
        L"1/8 t{a_1,a_2;i_1,i_2}:A-N-S * g{i_1,i_2;a_1,a_2}:A-C-S "
        L"+ 1/8 t⁺{i_1,i_2;a_1,a_2}:A-N-S * g{a_1,a_2;i_1,i_2}:A-C-S");
    const auto E1_brillouin = simplify(E1_contrib - E1_g_closed);
    REQUIRE_FALSE(has_tensor(E1_brillouin, L"g"));
    REQUIRE(has_tensor(E1_brillouin, L"f"));

    // <0|H̄²|0> = 1/12 <ij||ab> σ_i^a σ_j^b + h.c.
    REQUIRE_THAT(E2_contrib,
                 EquivalentTo(L"1/12 t{a_1;i_1}:A-N-S * t{a_2;i_2}:A-N-S "
                              L"* g{i_1,i_2;a_1,a_2}:A-C-S "
                              L"+ 1/12 t⁺{i_1;a_1}:A-N-S * t⁺{i_2;a_2}:A-N-S "
                              L"* g{a_1,a_2;i_1,i_2}:A-C-S"));
  }

  SECTION("bernoulli_config_validation") {
    using namespace sequant;
    using namespace sequant::mbpt;
    // only ranks 0..4 are implemented.
    REQUIRE_THROWS_AS(bernoulli::hbar(2, 5, false), Exception);
  }

  SECTION("bernoulli_quccsd") {
    using namespace sequant;
    using namespace sequant::mbpt;
    const CC::Options opts{.ansatz = CC::Ansatz::U,
                           .hbar_comm_rank = 2,
                           .hbar_expansion = CC::HbarExpansion::Bernoulli};
    CC cc(2, opts);

    // amplitudes through H̄² (hbar_comm_rank)
    const auto amps = cc.t();
    REQUIRE(amps.size() == 3);
    REQUIRE_THAT(amps[1], !EquivalentTo(ex<Constant>(0)));
    REQUIRE_THAT(amps[2], !EquivalentTo(ex<Constant>(0)));

    REQUIRE(size(amps[1]) == 32);
    REQUIRE(size(amps[2]) == 38);

    // one projected equation pinned in full: term counts are blind to the
    // coefficients, which is where a mis-weighted Wick reduction shows up.
    // Doubles at H̄¹, the smallest such equation.
    const auto R2_h1 = CC(2, {.ansatz = CC::Ansatz::U,
                              .hbar_comm_rank = 1,
                              .hbar_expansion = CC::HbarExpansion::Bernoulli})
                           .t()
                           .at(2);
    REQUIRE_THAT(
        R2_h1, EquivalentTo(
                   L"1/2 Â{i_1,i_2;a_1,a_2}:A-C-S * g{a_1,a_2;a_3,a_4}:A-C-S "
                   L"* t{a_3,a_4;i_1,i_2}:A-N-S "
                   L"+ Â{i_1,i_2;a_1,a_2}:A-C-S * g{a_1,a_2;i_1,i_2}:A-C-S "
                   L"+ 1/2 Â{i_1,i_2;a_1,a_2}:A-C-S * g{i_3,i_4;i_1,i_2}:A-C-S "
                   L"* t{a_1,a_2;i_3,i_4}:A-N-S "
                   L"+ 2 Â{i_1,i_2;a_1,a_2}:A-C-S * f{i_3;i_1}:A-C-S "
                   L"* t{a_1,a_2;i_2,i_3}:A-N-S "
                   L"- 2 Â{i_1,i_2;a_1,a_2}:A-C-S * f{a_1;a_3}:A-C-S "
                   L"* t{a_2,a_3;i_1,i_2}:A-N-S "
                   L"+ 2 Â{i_1,i_2;a_1,a_2}:A-C-S * g{a_1,a_2;i_1,a_3}:A-C-S "
                   L"* t{a_3;i_2}:A-N-S "
                   L"+ 2 Â{i_1,i_2;a_1,a_2}:A-C-S * g{i_3,a_1;i_1,i_2}:A-C-S "
                   L"* t{a_2;i_3}:A-N-S "
                   L"- 4 Â{i_1,i_2;a_1,a_2}:A-C-S * g{i_3,a_1;i_1,a_3}:A-C-S "
                   L"* t{a_2,a_3;i_2,i_3}:A-N-S"));

#ifndef SEQUANT_SKIP_LONG_TESTS
    const auto E = cc.energy(3);
    REQUIRE_THAT(E, !EquivalentTo(amps.at(0)));
    REQUIRE(size(E) == 46);
#endif  // !defined(SEQUANT_SKIP_LONG_TESTS)
  }

  SECTION("bernoulli_quccsd_eom") {
    using namespace sequant;
    using namespace sequant::mbpt;
    const CC cc(2, {.ansatz = CC::Ansatz::U,
                    .hbar_comm_rank = 2,
                    .hbar_expansion = CC::HbarExpansion::Bernoulli});
    // qUCCSD block ranks, 10.1063/5.0062090 Sec. II C: SS at the double
    // commutator (Eq. 29), SD/DS at the single (Eqs. 41, 44), DD bare (Eq. 48).
    const std::vector<std::size_t> quccsd = {2, 1, 1, 0};

    const auto ee = cc.eom_r(nₚ(2), nₕ(2), quccsd);
    REQUIRE(ee.size() == 3);
    REQUIRE(!ee[0]);
    REQUIRE(size(ee[1]) == 121);
    REQUIRE(size(ee[2]) == 21);

    // the same ranks drive IP: manifolds are indexed by ascending rank, so
    // {1h, 2h1p} takes the place of {S, D} (10.1021/acs.jctc.5c01991 Table 1)
    const auto ip = cc.eom_r(nₚ(1), nₕ(2), quccsd);
    REQUIRE(ip.size() == 2);
    REQUIRE(size(ip[0]) == 32);
    REQUIRE(size(ip[1]) == 11);

    // one manifold means one block, whose sandwich minus shift is the
    // commutator the uniform path builds; that path shares no code with
    // eom_r_blocked, so this pins the block construction against it
    const CC bch(2, {.ansatz = CC::Ansatz::U, .hbar_comm_rank = 2});
    REQUIRE_THAT(bch.eom_r(nₚ(1), nₕ(1), {2}).at(1),
                 EquivalentTo(bch.eom_r(nₚ(1), nₕ(1)).at(1)));

    if (sequant::assert_behavior() == sequant::AssertBehavior::Throw) {
      // block_ranks must be a K x K matrix over the manifolds ...
      REQUIRE_THROWS_AS(cc.eom_r(nₚ(2), nₕ(2), {2, 1, 0}), Exception);
      // ... and the ansatz must be unitary
      REQUIRE_THROWS_AS(CC(2).eom_r(nₚ(2), nₕ(2), quccsd), Exception);
    }

    // no matrix means uniform truncation at hbar_comm_rank
    REQUIRE_THAT(cc.eom_r(nₚ(2), nₕ(2)).at(1),
                 EquivalentTo(cc.eom_r(nₚ(2), nₕ(2), {2, 2, 2, 2}).at(1)));
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
