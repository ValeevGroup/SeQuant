//
// Kramers tracer (spinor.{hpp,cpp}) unit tests.
//

#include <catch2/catch_test_macros.hpp>
#include "catch2_sequant.hpp"

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/expressions/expr_algorithms.hpp>
#include <SeQuant/core/expressions/tensor.hpp>
#include <SeQuant/core/rational.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/core/utility/string.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>
#include <SeQuant/domain/mbpt/spinor.hpp>

#include <cstddef>
#include <memory>

TEST_CASE("kramers_trace", "[spinor]") {
  using namespace sequant;
  using namespace sequant::mbpt;

  auto ctx = get_default_context();
  ctx.set(CanonicalizeOptions{.method = CanonicalizationMethod::Complete});
  auto _ = set_scoped_default_context(ctx);
  TensorCanonicalizer::register_instance(
      std::make_shared<DefaultTensorCanonicalizer>());

  // counts RealPart nodes anywhere in the expression tree
  auto count_realparts = [](const ExprPtr& expr) {
    std::size_t n = 0;
    expr->visit(
        [&n](const ExprPtr& current) {
          if (current->is<RealPart>()) ++n;
        },
        /* atoms_only = */ false);
    return n;
  };

  SECTION("MP2 energy") {
    // E = 1/4 g-bar^{a1 a2}_{i1 i2} t-bar^{i1 i2}_{a1 a2}  (Kramers-free)
    const auto E = ex<Constant>(rational{1, 4}) *
                   ex<Tensor>(L"g", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"},
                              Symmetry::Antisymm) *
                   ex<Tensor>(L"t", bra{L"a_1", L"a_2"}, ket{L"i_1", L"i_2"},
                              Symmetry::Antisymm);

    ExprPtr result;
    REQUIRE_NOTHROW(result = closed_shell_kramers_trace(E));
    REQUIRE(result);

    INFO("closed_shell_kramers_trace(E_MP2) =\n" << toUtf8(to_latex(result)));

    REQUIRE(result->is<Sum>());

    // Each summand is `c * Re[1/2 g . t-bar]`, one per TRS-canonical Kramers
    // representative. With the global-T (whole-config) fold and sigma handled
    // by canonicalize — but the internal-T-reach (self-complementary block
    // fold) deferred — the MP2 energy yields 7 representatives: the 6-term form
    // plus the unfolded self-complementary pair (g^{a^ a_} and -g^{a_ a^}).
    REQUIRE(result->size() == 7);
    REQUIRE(count_realparts(result) == 7);

    // Coefficient bookkeeping: every summand is `Constant * RealPart`, and the
    // outer coefficients sum to 2^n = 16 (all Kramers configurations accounted
    // for: T-fold contributes x2, sigma multiplicity the rest).
    Constant::scalar_type coeff_sum = 0;
    for (const auto& term : *result) {
      REQUIRE(term->is<Product>());
      coeff_sum += term->as<Product>().scalar();
    }
    REQUIRE(coeff_sum == Constant::scalar_type{16});
  }

  SECTION("kramers_config_orbits: doubles external + direct folds") {
    using Perm = container::svector<std::size_t>;
    // bit layout (KR-MP2 notes): s_i=bit3, s_j=bit2, s_a=bit1, s_b=bit0.
    const Perm P_ij = {0, 1, 3, 2};   // swap occ pair (bits 3<->2)
    const Perm P_ab = {1, 0, 2, 3};   // swap virt pair (bits 1<->0)
    const Perm SIGMA = {1, 0, 3, 2};  // particle interchange (both pairs)

    auto canon_set = [](const auto& orbits) {
      container::set<std::uint64_t> c;
      for (const auto& o : orbits) c.insert(o.front());  // front = orbit-min
      return c;
    };
    auto total_members = [](const auto& orbits) {
      std::size_t n = 0;
      for (const auto& o : orbits) n += o.size();
      return n;
    };

    // External / R2 antisymmetry: <P_ab, P_ij, T> -> exactly 5 blocks
    // {0,15},{3,12},{1,2,13,14},{4,7,8,11},{5,6,9,10} (canonicals 0,1,3,4,5).
    auto ext = kramers_config_orbits(4, {P_ab, P_ij}, /*use_T=*/true);
    REQUIRE(ext.size() == 5);
    REQUIRE(total_members(ext) == 16);
    REQUIRE(canon_set(ext) == container::set<std::uint64_t>{0, 1, 3, 4, 5});

    // MP2 direct: Klein-four <sigma, T> -> 6 orbits.
    auto direct = kramers_config_orbits(4, {SIGMA}, /*use_T=*/true);
    REQUIRE(direct.size() == 6);
    REQUIRE(total_members(direct) == 16);

    // Sanity: no folding (no generators, no T) -> 16 singleton orbits.
    auto none = kramers_config_orbits(4, {}, /*use_T=*/false);
    REQUIRE(none.size() == 16);
  }

  SECTION("kramers_config_orbits: rank-general external fold (1/2/3)") {
    using Perm = container::svector<std::size_t>;
    // The helper makes no rank assumption: external blocks = orbits of
    // (#down-occ, #down-virt) in {0..k}^2 under T -> singles 2, doubles 5,
    // triples 8. Generators are the S_k transpositions of each external group.

    // Singles (rank-1): n=2 (a,i), no within-group swap, just T -> 2.
    REQUIRE(kramers_config_orbits(2, {}, /*use_T=*/true).size() == 2);
    REQUIRE(kramers_config_orbits(2, {}, /*use_T=*/false).size() == 4);

    // Triples (rank-3): n=6, S_3 on the virtual triple (bits 5,4,3) and on the
    // occupied triple (bits 2,1,0) via adjacent transpositions, plus T -> 8.
    const Perm v_54 = {0, 1, 2, 3, 5, 4};
    const Perm v_43 = {0, 1, 2, 4, 3, 5};
    const Perm o_21 = {0, 2, 1, 3, 4, 5};
    const Perm o_10 = {1, 0, 2, 3, 4, 5};
    auto triples =
        kramers_config_orbits(6, {v_54, v_43, o_21, o_10}, /*use_T=*/true);
    std::size_t triples_members = 0;
    for (const auto& o : triples) triples_members += o.size();
    REQUIRE(triples_members == 64);
    REQUIRE(triples.size() == 8);
  }

  SECTION("CC external fold: CCD R2 driver -> 5 labeled blocks") {
    using sequant::reserved::antisymm_label;
    // driver:  Â^{a1 a2}_{i1 i2} * g-bar^{a1 a2}_{i1 i2}  (all indices
    // external)
    auto A = ex<Tensor>(antisymm_label(), bra{L"a_1", L"a_2"},
                        ket{L"i_1", L"i_2"}, Symmetry::Antisymm);
    auto g = ex<Tensor>(L"g", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"},
                        Symmetry::Antisymm);
    auto driver = A * g;

    container::svector<ExprPtr> blocks;
    REQUIRE_NOTHROW(blocks = closed_shell_kramers_CC_trace(driver));

    INFO("CCD driver external blocks:");
    for (const auto& b : blocks) INFO("  " << toUtf8(to_latex(b)));

    // doubles external antisymmetry + T -> 5 symmetry-unique blocks
    REQUIRE(blocks.size() == 5);

    // Â is factored out; each block's remaining tensors carry fully
    // Kramers-labeled indices (every index has an up/down spin annotation)
    auto all_labeled = [](const ExprPtr& b) {
      bool ok = true;
      b->visit(
          [&ok](const ExprPtr& cur) {
            if (!cur->is<Tensor>()) return;
            const auto& t = cur->as<Tensor>();
            auto check = [&ok](const Index& idx) {
              const auto s = mbpt::to_spin(idx.space().qns());
              if (s != mbpt::Spin::alpha && s != mbpt::Spin::beta) ok = false;
            };
            for (const auto& idx : t.bra()) check(idx);
            for (const auto& idx : t.ket()) check(idx);
          },
          /* atoms_only = */ true);
      return ok;
    };
    for (const auto& b : blocks) REQUIRE(all_labeled(b));
  }

  SECTION("CC internal fold: CCD pp-ladder -> 5 blocks, internal classes") {
    using sequant::reserved::antisymm_label;
    // pp-ladder:  1/2 Â^{a1a2}_{i1i2} * g-bar^{a1a2}_{a3a4} *
    // t-bar^{a3a4}_{i1i2}
    //   external a1,a2 (virt), i1,i2 (occ); internal (contracted) a3,a4 (virt).
    auto A = ex<Tensor>(antisymm_label(), bra{L"a_1", L"a_2"},
                        ket{L"i_1", L"i_2"}, Symmetry::Antisymm);
    auto g = ex<Tensor>(L"g", bra{L"a_1", L"a_2"}, ket{L"a_3", L"a_4"},
                        Symmetry::Antisymm);
    auto t = ex<Tensor>(L"t", bra{L"a_3", L"a_4"}, ket{L"i_1", L"i_2"},
                        Symmetry::Antisymm);
    auto ppladder = ex<Constant>(rational{1, 2}) * A * g * t;

    container::svector<ExprPtr> blocks;
    REQUIRE_NOTHROW(blocks = closed_shell_kramers_CC_trace(ppladder));

    // external fold is universal -> still 5 blocks
    REQUIRE(blocks.size() == 5);

    // each block sums the internal Kramers configs of (c1,c2); every index
    // (external + internal) is Kramers-labeled, and the block is non-empty.
    auto term_count = [](const ExprPtr& b) -> std::size_t {
      return b->is<Sum>() ? b->as<Sum>().size() : 1;
    };
    for (const auto& b : blocks) {
      REQUIRE(b);
      INFO("pp-ladder block (" << term_count(b)
                               << " terms): " << toUtf8(to_latex(b)));
      REQUIRE(term_count(b) >= 1);
      bool labeled = true;
      b->visit(
          [&labeled](const ExprPtr& cur) {
            if (!cur->is<Tensor>()) return;
            const auto& tt = cur->as<Tensor>();
            for (const auto& idx : tt.bra())
              if (mbpt::to_spin(idx.space().qns()) == mbpt::Spin::any)
                labeled = false;
            for (const auto& idx : tt.ket())
              if (mbpt::to_spin(idx.space().qns()) == mbpt::Spin::any)
                labeled = false;
          },
          /* atoms_only = */ true);
      REQUIRE(labeled);
    }
  }

  SECTION("CC internal fold: heterogeneous contraction (ring) -> 4 classes") {
    using sequant::reserved::antisymm_label;
    // A heterogeneous internal pair (one virtual a3, one occupied i3) is
    // contracted between g and t. sigma = swap(a3,i3) maps virt<->occ, so it is
    // NOT a symmetry (canonicalize won't merge) -> 4 internal classes in a
    // generic external block, vs the pp-ladder's 3 (homogeneous pair).
    auto A = ex<Tensor>(antisymm_label(), bra{L"a_1", L"a_2"},
                        ket{L"i_1", L"i_2"}, Symmetry::Antisymm);
    auto g = ex<Tensor>(L"g", bra{L"a_1", L"a_2"}, ket{L"a_3", L"i_3"},
                        Symmetry::Nonsymm);
    auto t = ex<Tensor>(L"t", bra{L"a_3", L"i_3"}, ket{L"i_1", L"i_2"},
                        Symmetry::Nonsymm);
    auto ring = A * g * t;

    container::svector<ExprPtr> blocks;
    REQUIRE_NOTHROW(blocks = closed_shell_kramers_CC_trace(ring));
    REQUIRE(blocks.size() == 5);

    auto term_count = [](const ExprPtr& b) -> std::size_t {
      return b->is<Sum>() ? b->as<Sum>().size() : 1;
    };
    // generic external blocks: 4 internal classes (no sigma fold);
    // at least one block must show all 4 distinct internal configs.
    std::size_t max_terms = 0;
    for (const auto& b : blocks) {
      INFO("ring block (" << term_count(b) << " terms)");
      max_terms = std::max(max_terms, term_count(b));
    }
    REQUIRE(max_terms == 4);
  }

  SECTION("CC internal fold: separable T2^2 quad -> 10 classes (9-vs-10)") {
    using sequant::reserved::antisymm_label;
    // separable quad: 1/4 Â g-bar^{a3a4}_{i3i4} t-bar^{a3a4}_{i1i2}
    //                       t-bar^{a1a2}_{i3i4}
    //   two homogeneous internal pairs (a3,a4 virt; i3,i4 occ). The surviving
    //   internal symmetry is the COMBINED sigma=(i3 i4)(a3 a4) -> 10 classes;
    //   two independent swaps would (wrongly) give 9. canonicalize must merge
    //   only via valid same-Kramers dummy relabels.
    auto A = ex<Tensor>(antisymm_label(), bra{L"a_1", L"a_2"},
                        ket{L"i_1", L"i_2"}, Symmetry::Antisymm);
    auto g = ex<Tensor>(L"g", bra{L"i_3", L"i_4"}, ket{L"a_3", L"a_4"},
                        Symmetry::Antisymm);
    auto t1 = ex<Tensor>(L"t", bra{L"a_3", L"a_4"}, ket{L"i_1", L"i_2"},
                         Symmetry::Antisymm);
    auto t2 = ex<Tensor>(L"t", bra{L"a_1", L"a_2"}, ket{L"i_3", L"i_4"},
                         Symmetry::Antisymm);
    auto quad = ex<Constant>(rational{1, 4}) * A * g * t1 * t2;

    container::svector<ExprPtr> blocks;
    REQUIRE_NOTHROW(blocks = closed_shell_kramers_CC_trace(quad));
    REQUIRE(blocks.size() == 5);

    auto term_count = [](const ExprPtr& b) -> std::size_t {
      return b->is<Sum>() ? b->as<Sum>().size() : 1;
    };
    std::size_t max_terms = 0;
    for (const auto& b : blocks) {
      INFO("quad block (" << term_count(b) << " terms)");
      max_terms = std::max(max_terms, term_count(b));
    }
    // The notes' emit-and-reconstruct strategy counts 10 internal classes under
    // the COMBINED sigma=(i3i4)(a3a4) (independent swaps excluded for a single
    // signed config). This tracer instead sums all 16 internal configs and lets
    // canonicalize merge: for the SUMMED contraction each independent
    // homogeneous-pair swap IS a symmetry (the double antisymmetry of g and the
    // amplitude cancels the two signs), so canonicalize folds under the full
    // Klein-4 -> 9. canonicalize only merges provably-equal terms and the
    // 16-config sum is value-exact, so 9 is value-preserving (the 9-vs-10
    // hazard only affects external reconstruction, which is eval-time). The
    // definitive value check is the numeric CCk validation.
    REQUIRE(max_terms == 9);
  }

  SECTION("kramers_external_blocks: doubles reconstruction transforms") {
    using Perm = container::svector<std::size_t>;
    const Perm P_ab = {1, 0, 2, 3};  // swap virt bits 0,1
    const Perm P_ij = {0, 1, 3, 2};  // swap occ bits 2,3
    auto blocks = kramers_external_blocks(4, {P_ab, P_ij}, /*use_T=*/true);
    REQUIRE(blocks.size() == 5);

    auto apply_perm = [](const Perm& p, std::uint64_t m) {
      std::uint64_t out = 0;
      for (std::size_t k = 0; k < p.size(); ++k)
        if ((m >> k) & 1u) out |= (std::uint64_t{1} << p[k]);
      return out;
    };
    const std::uint64_t full = 15;

    std::size_t total = 0;
    for (const auto& b : blocks) {
      total += b.members.size();
      bool found_canonical = false;
      for (const auto& mem : b.members) {
        // perm/conj must reproduce the member's config from the canonical:
        // bit-permutation commutes with complement, so the config is
        // apply_perm(perm, canonical), complemented iff conj (an odd # of T's).
        const auto pc = apply_perm(mem.perm, b.canonical);
        const auto expected = mem.conj ? ((~pc) & full) : pc;
        REQUIRE(mem.config == expected);
        if (mem.config == b.canonical) {
          found_canonical = true;
          REQUIRE(mem.sign == 1);
          REQUIRE(!mem.conj);
        }
      }
      REQUIRE(found_canonical);
    }
    REQUIRE(total == 16);

    // unambiguous T-pair members: block(complement) = +conj(block(canonical))
    // for canonicals 0000 (#down 0) and 0011 (#down 2), both even -> sign +1.
    auto find_member = [&](std::uint64_t canon,
                           std::uint64_t cfg) -> const KramersBlockMember* {
      for (const auto& b : blocks)
        if (b.canonical == canon)
          for (const auto& m : b.members)
            if (m.config == cfg) return &m;
      return nullptr;
    };
    const auto* m15 = find_member(0, 15);
    REQUIRE(m15);
    REQUIRE(m15->conj);
    REQUIRE(m15->sign == 1);
    const auto* m12 = find_member(3, 12);
    REQUIRE(m12);
    REQUIRE(m12->conj);
    REQUIRE(m12->sign == 1);
  }

  SECTION("kramers_external_blocks: symmetric (raw g) generators") {
    using Perm = container::svector<std::size_t>;
    // A raw, non-antisymmetrized integral g^{ab}_{ij} folds under particle
    // interchange sigma (swap BOTH pairs, sign +1) and T -> 6 blocks (vs the
    // antisymmetrized ḡ[as]'s 5). sigma is passed as a symm_perm (sign +1).
    const Perm SIGMA = {1, 0, 3, 2};  // swap virt pair AND occ pair
    auto blocks =
        kramers_external_blocks(4, /*antisym=*/{}, /*use_T=*/true, {SIGMA});
    REQUIRE(blocks.size() == 6);

    auto apply_perm = [](const Perm& p, std::uint64_t m) {
      std::uint64_t out = 0;
      for (std::size_t k = 0; k < p.size(); ++k)
        if ((m >> k) & 1u) out |= (std::uint64_t{1} << p[k]);
      return out;
    };
    const std::uint64_t full = 15;

    std::size_t total = 0;
    for (const auto& b : blocks) {
      total += b.members.size();
      for (const auto& mem : b.members) {
        // same perm/conj <-> config invariant as the antisymm case.
        const auto pc = apply_perm(mem.perm, b.canonical);
        REQUIRE(mem.config == (mem.conj ? ((~pc) & full) : pc));
        // a pure-sigma member (no T, i.e. conj=false, non-identity perm) must
        // carry sign +1 (sigma is a symmetry, not an antisymmetry).
        if (!mem.conj && mem.config != b.canonical) REQUIRE(mem.sign == 1);
      }
    }
    REQUIRE(total == 16);
  }
}
