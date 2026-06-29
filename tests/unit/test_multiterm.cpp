#include <catch2/catch_test_macros.hpp>

#include "catch2_sequant.hpp"

#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/expressions/expr_algorithms.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/optimize/optimize.hpp>
#include <SeQuant/core/space.hpp>

#include <algorithm>
#include <cstddef>
#include <string>
#include <vector>

namespace {

sequant::ExprPtr parse_antisymm(std::wstring_view s) {
  using namespace sequant;
  return deserialize(s, {.def_perm_symm = Symmetry::Antisymm});
}

sequant::ExprPtr parse_nonsymm(std::wstring_view s) {
  using namespace sequant;
  return deserialize(s, {.def_perm_symm = Symmetry::Nonsymm});
}

std::wstring latex(sequant::ExprPtr const& e) { return sequant::to_latex(e); }

/// Every factored-group summand of \p sum, in summand order. A factored group
/// is a Product whose factors include a Sum -- i.e. shared * (partner1 +
/// partner2 +
/// ...).
std::vector<sequant::ExprPtr> find_all_factored(sequant::ExprPtr const& sum) {
  using namespace sequant;
  std::vector<ExprPtr> out;
  if (!sum->is<Sum>()) return out;
  for (auto const& s : sum->as<Sum>().summands()) {
    if (!s->is<Product>()) continue;
    for (auto const& f : s->as<Product>().factors())
      if (f->is<Sum>()) {
        out.push_back(s);
        break;
      }
  }
  return out;
}

/// The first factored-group summand of \p sum (see \ref find_all_factored), or
/// null if none.
sequant::ExprPtr find_factored(sequant::ExprPtr const& sum) {
  auto const all = find_all_factored(sum);
  return all.empty() ? sequant::ExprPtr{} : all.front();
}

}  // namespace

TEST_CASE("multiterm factorization", "[multiterm]") {
  using namespace sequant;

  // Sized spaces so the cost model is well-defined. Set once here, in the
  // outermost scope, rather than per section; the cost-driven-winner section,
  // whose expected fold depends on the occ-vs-virt ordering, overrides them
  // locally via set_sizes.
  auto ctx_resetter = set_scoped_default_context(get_default_context().clone());
  auto reg = get_default_context().mutable_index_space_registry();
  auto set_sizes = [&reg](std::size_t occ, std::size_t virt) {
    reg->retrieve_ptr(L"i")->approximate_size(occ);
    reg->retrieve_ptr(L"a")->approximate_size(virt);
  };
  set_sizes(/*occ=*/10, /*virt=*/20);

  // Multi-term factorization enabled; NoReorder keeps the fold structure
  // observable. The DenseSize and reorder sections build their own options.
  OptimizeOptions const mt{.reorder = ReorderSum::NoReorder,
                           .multiterm = MultiTermFactor::Enable};

  SECTION("disabled is the default and matches an explicit Disable") {
    auto const expr = parse_antisymm(
        L"G{i1,i2;a1,a2} T{a1,a2;i1,i2}"
        L" + G{i1,i2;a1,a2} Z{a1,a2;i1,i2}");
    REQUIRE(*optimize(expr, {.multiterm = MultiTermFactor::Disable}) ==
            *optimize(expr));
  }

  SECTION("enabling the flag does not throw and yields a Sum") {
    auto const expr = parse_antisymm(
        L"G{i1,i2;a1,a2} T{a1,a2;i1,i2}"
        L" + G{i1,i2;a1,a2} Z{a1,a2;i1,i2}");
    ExprPtr res;
    REQUIRE_NOTHROW(res = optimize(expr, mt));
    REQUIRE(res);
    REQUIRE(res->is<Sum>());
  }

  SECTION("one-sided fold: A*B + A*C -> A*(B + C) when it lowers cost") {
    // V contracted with T (or B) over a3,a4 leaves a tensor R{a1,a2;i1,i2}; the
    // contraction (O(a^4 i^2)) dwarfs the partner-sum (O(a^2 i^2)), so folding
    // pays off.
    auto const expr = parse_antisymm(
        L"V{a1,a2;a3,a4} T{a3,a4;i1,i2}"
        L" + V{a1,a2;a3,a4} B{a3,a4;i1,i2}");
    auto const res = optimize(expr, mt);
    // EquivalentTo expands both sides, so this is the round-trip guard *and*
    // pins the exact content (the shared V, the partners, their indices); the
    // factored RHS doubles as documentation of the intended output.
    REQUIRE_THAT(res,
                 EquivalentTo(parse_antisymm(
                     L"V{a1,a2;a3,a4} (T{a3,a4;i1,i2} + B{a3,a4;i1,i2})")));
    // Structural: a fold actually happened (EquivalentTo cannot see this -- it
    // expands the factored form away).
    REQUIRE(res->size() == 1);
    REQUIRE(find_factored(res));
  }

  SECTION("N-ary one-sided fold: A*B + A*C + A*D -> A*(B + C + D)") {
    auto const expr = parse_antisymm(
        L"V{a1,a2;a3,a4} T{a3,a4;i1,i2}"
        L" + V{a1,a2;a3,a4} B{a3,a4;i1,i2}"
        L" + V{a1,a2;a3,a4} U{a3,a4;i1,i2}");
    auto const res = optimize(expr, mt);
    REQUIRE_THAT(res, EquivalentTo(parse_antisymm(
                          L"V{a1,a2;a3,a4} (T{a3,a4;i1,i2}"
                          L" + B{a3,a4;i1,i2} + U{a3,a4;i1,i2})")));
    REQUIRE(res->size() == 1);
    REQUIRE(find_factored(res));
  }

  SECTION("transposed partners fold: V*t + V*t^T -> V*(t + t^T)") {
    // After spin-tracing CC equations a shared factor multiplies two
    // index-transposed copies of a (now non-antisymmetric) amplitude. They are
    // two distinct partner vertices on one side of the biclique, so they fold
    // like any other one-sided pair. (Antisymmetric t would instead collapse
    // the two copies into one term up to sign, so this uses non-antisymmetric
    // tensors -- the realistic post-spintracing case.)
    auto const expr = parse_nonsymm(
        L"V{i1,i2;a1,a2} t{a1,a2;i3,i4}"
        L" + V{i1,i2;a1,a2} t{a2,a1;i4,i3}");
    auto const res = optimize(expr, mt);
    REQUIRE_THAT(res,
                 EquivalentTo(parse_nonsymm(
                     L"V{i1,i2;a1,a2} (t{a1,a2;i3,i4} + t{a2,a1;i4,i3})")));
    REQUIRE(res->size() == 1);
    REQUIRE(find_factored(res));
  }

  SECTION("a bare-tensor summand passes through while its neighbors fold") {
    // V*T + V*B fold to V*(T + B); the bare tensor D{a1,a2;i1,i2} is a leaf (no
    // contraction to split), so extract_core() returns nullopt and it is never
    // interned, scored, or consumed. Reassembly emits untouched summands first
    // (original order), then the folds:
    // D + V*(T + B)
    auto const expr = parse_antisymm(
        L"V{a1,a2;a3,a4} T{a3,a4;i1,i2}"
        L" + V{a1,a2;a3,a4} B{a3,a4;i1,i2}"
        L" + D{a1,a2;i1,i2}");
    auto const res = optimize(expr, mt);
    REQUIRE_THAT(res,
                 EquivalentTo(parse_antisymm(
                     L"D{a1,a2;i1,i2}"
                     L" + V{a1,a2;a3,a4} (T{a3,a4;i1,i2} + B{a3,a4;i1,i2})")));
    REQUIRE(res->size() == 2);  // the untouched bare tensor + one fold
    REQUIRE(find_factored(res));
    // D survived verbatim and, with NoReorder, precedes the fold.
    auto const& smands = res->as<Sum>().summands();
    REQUIRE(smands.front()->is<Tensor>());
    REQUIRE(smands.front()->as<Tensor>().label() == L"D");
  }

  SECTION("two-sided biclique: AX + AY + BX + BY -> (A + B)*(X + Y)") {
    // U,V share a shape (left factors); T,B share a shape (right factors). The
    // four contractions form a complete 2x2 bipartite graph and fold to
    // (U + V)*(T + B).
    auto const expr = parse_antisymm(
        L"U{a1,a2;a3,a4} T{a3,a4;i1,i2}"
        L" + U{a1,a2;a3,a4} B{a3,a4;i1,i2}"
        L" + V{a1,a2;a3,a4} T{a3,a4;i1,i2}"
        L" + V{a1,a2;a3,a4} B{a3,a4;i1,i2}");
    auto const res = optimize(expr, mt);
    REQUIRE_THAT(res, EquivalentTo(parse_antisymm(
                          L"(U{a1,a2;a3,a4} + V{a1,a2;a3,a4})"
                          L" (T{a3,a4;i1,i2} + B{a3,a4;i1,i2})")));
    REQUIRE(res->size() == 1);  // all four folded into one product
    auto const factored = find_factored(res);
    REQUIRE(factored);
    // Both factors are 2-term sums: (U + V) and (T + B).
    for (auto const& f : factored->as<Product>().factors()) {
      REQUIRE(f->is<Sum>());
      REQUIRE(f->size() == 2);
    }
  }

  SECTION("incomplete graph: cost-driven winner flips with extents") {
    // Topology: edges (U,T),(U,B),(V,T) -- no (V,B). Two maximal one-sided
    // folds compete with identical avoided contraction cost C:
    //   {U}x{T,B}  -> U*(T + B),   build cost ~ size(T) = a^2 i^2
    //   {U,V}x{T}  -> (U + V)*T,   build cost ~ size(U) = a^4
    // Greedy keeps the higher-saving fold -- the one whose summed side is
    // cheaper to build: it sums the lighter factors and shares the heavier one.
    // Which factor is heavier is set purely by the occ-vs-virt extent ordering,
    // so the same expression under both orderings must flip the shared factor
    // (U <-> T). The flip shows the choice is cost-driven; the round-trip is
    // the cost-model-independent guard.
    std::wstring const expr_text =
        L"U{a1,a2;a3,a4} T{a3,a4;i1,i2}"
        L" + U{a1,a2;a3,a4} B{a3,a4;i1,i2}"
        L" + V{a1,a2;a3,a4} T{a3,a4;i1,i2}";

    // The shared (non-summed) factor's tensor label, with the structural guards
    // that hold regardless of which side wins: exactly one fold over a 2-term
    // partner sum, plus one untouched leftover summand.
    auto winning_shared_label = [](ExprPtr const& res) -> std::wstring {
      REQUIRE(res->is<Sum>());
      REQUIRE(res->size() == 2);  // one fold + one untouched leftover
      auto const folded = find_all_factored(res);
      REQUIRE(folded.size() == 1);
      ExprPtr shared, partner_sum;
      for (auto const& f : folded.front()->as<Product>().factors())
        (f->is<Sum>() ? partner_sum : shared) = f;
      REQUIRE(partner_sum);
      REQUIRE(partner_sum->size() == 2);
      REQUIRE(shared);
      REQUIRE(shared->is<Tensor>());
      return std::wstring{shared->as<Tensor>().label()};
    };

    // virt-heavy (a=20 > i=10): size(U)=a^4 dominates size(T)=a^2 i^2, so U is
    // the heavier factor and is shared -- summing the lighter T, B into
    // U*(T + B).
    set_sizes(/*occ=*/10, /*virt=*/20);
    auto const expr_vh = parse_antisymm(expr_text);
    auto const res_vh = optimize(expr_vh, mt);
    REQUIRE(winning_shared_label(res_vh) == L"U");
    REQUIRE_THAT(res_vh, EquivalentTo(expr_vh));

    // occ-heavy (i=30 > a=10): now size(T)=a^2 i^2 dominates size(U)=a^4, so T
    // is the heavier factor and is shared -- summing the lighter U, V into
    // (U + V)*T.
    set_sizes(/*occ=*/30, /*virt=*/10);
    auto const expr_oh = parse_antisymm(expr_text);
    auto const res_oh = optimize(expr_oh, mt);
    REQUIRE(winning_shared_label(res_oh) == L"T");
    REQUIRE_THAT(res_oh, EquivalentTo(expr_oh));
  }

  SECTION("full contraction saves nothing: shareable but no fold") {
    // G*T + G*Z both fully contract to a scalar -- exactly the saving()==0
    // boundary, and structurally so (not by extent tuning): for a full
    // contraction the avoided cost C equals the product of the whole index set,
    // and the only build cost is the partner sum (T + Z), whose footprint
    // size(T) is that same whole index set. With the one-sided
    // saving = (m*n - 1)*C - (n - 1)*size(T) at m=1, n=2 this is C - size(T) =
    // 0 for any extents, so nothing folds. The cancellation is exact only under
    // the current (m*n - 1)/(n - 1) coefficients (CostModel::saving,
    // multiterm.cpp); revisit this no-fold check if those change. The
    // partial-contraction sibling below brackets this boundary from above.
    auto const expr = parse_antisymm(
        L"G{i1,i2;a1,a2} T{a1,a2;i1,i2}"
        L" + G{i1,i2;a1,a2} Z{a1,a2;i1,i2}");
    auto const res = optimize(expr, mt);
    REQUIRE(res->is<Sum>());
    REQUIRE(res->size() == 2);
    REQUIRE_FALSE(find_factored(res));
  }

  SECTION("partial contraction clears the threshold: V*(T + Z)") {
    // Contrast with the scalar case: V shares with T and Z, but the contraction
    // is partial (over a3,a4 only), so the result keeps free indices and the
    // avoided cost C = O(a^4 i^2) strictly exceeds the partner-build
    // size(T) = O(a^2 i^2). saving = C - size(T) > 0, so V*(T + Z) folds. With
    // the scalar case sitting exactly on saving()==0, this brackets that
    // boundary from the positive side.
    auto const expr = parse_antisymm(
        L"V{a1,a2;a3,a4} T{a3,a4;i1,i2}"
        L" + V{a1,a2;a3,a4} Z{a3,a4;i1,i2}");
    auto const res = optimize(expr, mt);
    REQUIRE_THAT(res,
                 EquivalentTo(parse_antisymm(
                     L"V{a1,a2;a3,a4} (T{a3,a4;i1,i2} + Z{a3,a4;i1,i2})")));
    REQUIRE(res->size() == 1);  // V * (T + Z)
    REQUIRE(find_factored(res));
  }

  SECTION("DenseSize objective: the fold still pays") {
    // Unlike the other sections (default DenseFLOPs metric), this one sets
    // objective_function to DenseSize, routing CostModel::contraction_cost
    // through memsize_counter instead of flops_counter. V*T + V*B still folds:
    // under DenseSize the avoided contraction's element footprint (O(a^4), the
    // V leg dominating) far exceeds the cost of building the (T + B) partner
    // sum.
    auto const expr = parse_antisymm(
        L"V{a1,a2;a3,a4} T{a3,a4;i1,i2}"
        L" + V{a1,a2;a3,a4} B{a3,a4;i1,i2}");
    auto const res =
        optimize(expr, {.objective_function = ObjectiveFunction::DenseSize,
                        .reorder = ReorderSum::NoReorder,
                        .multiterm = MultiTermFactor::Enable});
    REQUIRE_THAT(res,
                 EquivalentTo(parse_antisymm(
                     L"V{a1,a2;a3,a4} (T{a3,a4;i1,i2} + B{a3,a4;i1,i2})")));
    REQUIRE(res->size() == 1);  // V * (T + B)
    REQUIRE(find_factored(res));
  }

  SECTION("honors reorder independently of multiterm") {
    // reorder and multiterm are independent options, so enabling multiterm must
    // not disable reorder. The input keeps the two effects separable: its two
    // g*t summands fully contract to scalars, so folding saves nothing (the
    // saving()==0 boundary) and the factorizer leaves all three summands alone.
    // reorder, on the other hand, does move them. So latex -- which (unlike
    // EquivalentTo, which canonicalizes summand order away) is order-sensitive
    // -- is the right comparison here.
    auto const expr = parse_antisymm(
        L"g{i3,i4;a1,a2} t{a1,a2;i5,i6} A{i5,i6;i3,i4}"
        L" + p{i1;i2} q{i2;i1}"
        L" + g{i3,i4;a1,a2} t{a1,a2;i5,i6} B{i5,i6;i3,i4}");

    auto const mt_noreorder =
        optimize(expr, {.reorder = ReorderSum::NoReorder,
                        .multiterm = MultiTermFactor::Enable});
    auto const mt_reorder = optimize(
        expr,
        {.reorder = ReorderSum::Reorder, .multiterm = MultiTermFactor::Enable});

    // Nothing folds (full contraction, saving 0), so all three summands survive
    // and reorder is the only thing that can change the output.
    REQUIRE(mt_noreorder->is<Sum>());
    REQUIRE(mt_noreorder->size() == 3);
    REQUIRE_FALSE(find_factored(mt_noreorder));

    // With the bug, reorder got skipped under multiterm and these two came out
    // identical; fixed, they differ.
    REQUIRE(latex(mt_noreorder) != latex(mt_reorder));

    // And the order matches plain reorder: with nothing to fold, multiterm
    // leaves reorder's result untouched.
    auto const reorder_only =
        optimize(expr, {.reorder = ReorderSum::Reorder,
                        .multiterm = MultiTermFactor::Disable});
    REQUIRE(latex(mt_reorder) == latex(reorder_only));
  }

  SECTION("distinct contraction signatures fold independently, not merged") {
    // All four summands share the external indices {a1,a2,i1,i2}. The first two
    // contract over a virtual pair (a3,a4); the last two over an occupied pair
    // (i3,i4). Distinct signatures => two independent buckets => two separate
    // folds V*(T + B) and W*(X + Y), never one cross-bucket merge.
    auto const expr = parse_antisymm(
        L"V{a1,a2;a3,a4} T{a3,a4;i1,i2}"
        L" + V{a1,a2;a3,a4} B{a3,a4;i1,i2}"
        L" + W{a1,a2;i3,i4} X{i3,i4;i1,i2}"
        L" + W{a1,a2;i3,i4} Y{i3,i4;i1,i2}");
    auto const res = optimize(expr, mt);
    REQUIRE_THAT(res,
                 EquivalentTo(parse_antisymm(
                     L"V{a1,a2;a3,a4} (T{a3,a4;i1,i2} + B{a3,a4;i1,i2})"
                     L" + W{a1,a2;i3,i4} (X{i3,i4;i1,i2} + Y{i3,i4;i1,i2})")));
    REQUIRE(res->size() == 2);  // one fold per bucket
    REQUIRE(find_all_factored(res).size() == 2);
  }

  SECTION("sign fold: A*B - A*C -> A*(B - C)") {
    // The factorizer peels the (-1) prefactor off the summand and folds the
    // relative sign into the partner.
    auto const expr = parse_antisymm(
        L"V{a1,a2;a3,a4} T{a3,a4;i1,i2}"
        L" - V{a1,a2;a3,a4} B{a3,a4;i1,i2}");
    auto const res = optimize(expr, mt);
    REQUIRE_THAT(res,
                 EquivalentTo(parse_antisymm(
                     L"V{a1,a2;a3,a4} (T{a3,a4;i1,i2} - B{a3,a4;i1,i2})")));
    REQUIRE(res->size() == 1);  // V * (T - B)
    REQUIRE(find_factored(res));
  }

  SECTION("rational prefactors: 2*A*B + 3*A*C -> A*(2B + 3C)") {
    auto const expr = parse_antisymm(
        L"2 V{a1,a2;a3,a4} T{a3,a4;i1,i2}"
        L" + 3 V{a1,a2;a3,a4} B{a3,a4;i1,i2}");
    auto const res = optimize(expr, mt);
    REQUIRE_THAT(res,
                 EquivalentTo(parse_antisymm(
                     L"V{a1,a2;a3,a4} (2 T{a3,a4;i1,i2} + 3 B{a3,a4;i1,i2})")));
    REQUIRE(res->size() == 1);
    REQUIRE(find_factored(res));
  }

  SECTION("rank-1 two-sided sign matrix folds as (U - V)*(T + B)") {
    // U*T + U*B - V*T - V*B. The sign matrix [[+,+],[-,-]] is multiplicative
    // rank 1, so the whole 2x2 still folds to a single product, the sign riding
    // the left factor-sum.
    auto const expr = parse_antisymm(
        L"U{a1,a2;a3,a4} T{a3,a4;i1,i2}"
        L" + U{a1,a2;a3,a4} B{a3,a4;i1,i2}"
        L" - V{a1,a2;a3,a4} T{a3,a4;i1,i2}"
        L" - V{a1,a2;a3,a4} B{a3,a4;i1,i2}");
    auto const res = optimize(expr, mt);
    REQUIRE_THAT(res, EquivalentTo(parse_antisymm(
                          L"(U{a1,a2;a3,a4} - V{a1,a2;a3,a4})"
                          L" (T{a3,a4;i1,i2} + B{a3,a4;i1,i2})")));
    REQUIRE(res->size() == 1);  // single (U - V)*(T + B)
    auto const factored = find_factored(res);
    REQUIRE(factored);
    for (auto const& f : factored->as<Product>().factors()) {
      REQUIRE(f->is<Sum>());
      REQUIRE(f->size() == 2);
    }
  }

  SECTION("non-rank-1 sign matrix reduces to one-sided folds") {
    // U*T + U*B + V*T - V*B. The sign matrix [[+,+],[+,-]] does not factor, so
    // a single (U +/- V)*(T +/- B) would be wrong. The factorizer must decline
    // the 2x2 and reduce to two one-sided folds. The reduction is not unique
    // (it may share either side), so EquivalentTo on the original sum is the
    // correctness guard while the structural checks pin "two one-sided folds,
    // no bogus 2x2".
    auto const expr = parse_antisymm(
        L"U{a1,a2;a3,a4} T{a3,a4;i1,i2}"
        L" + U{a1,a2;a3,a4} B{a3,a4;i1,i2}"
        L" + V{a1,a2;a3,a4} T{a3,a4;i1,i2}"
        L" - V{a1,a2;a3,a4} B{a3,a4;i1,i2}");
    auto const res = optimize(expr, mt);
    REQUIRE_THAT(res, EquivalentTo(expr));
    REQUIRE(res->size() == 2);  // two one-sided folds
    REQUIRE(find_all_factored(res).size() == 2);
    for (auto const& f : find_all_factored(res)) {
      ExprPtr partner_sum;
      for (auto const& g : f->as<Product>().factors())
        if (g->is<Sum>()) partner_sum = g;
      REQUIRE(partner_sum);
      REQUIRE(partner_sum->size() == 2);
    }
  }

  SECTION("canonicalization-phase fold lands a relative sign on the partner") {
    // Distinct from the scalar-prefactor folds above: there the relative sign
    // was an explicit input coefficient that extract_core peels off. Here
    // neither summand carries any sign or scalar -- every relative sign comes
    // purely from canonicalization. Both summands share the intermediate
    // M = G*W, contracted over the virtual pair a3,a4 (so M is a child factor;
    // only an intermediate can carry a non-trivial canon_phase). The second
    // summand writes G's contracted pair swapped (a4,a3); antisymmetric G makes
    // M canonicalize to the same indices but opposite phase. Phase-relaxed
    // matching interns both Ms into one vertex; the relative -1 rides the edge
    // and lands on the partner: M*(X - Y). size()==1 is the proof the
    // relaxation engaged (phase-strict would leave size()==2).
    auto const expr = parse_antisymm(
        L"G{i1,i2;a3,a4} W{a3,a4;i3,i4} X{i3,i4;a1,a2}"
        L" + G{i1,i2;a4,a3} W{a3,a4;i3,i4} Y{i3,i4;a1,a2}");
    auto const res = optimize(expr, mt);
    REQUIRE(res->is<Sum>());
    REQUIRE(res->size() == 1);  // folded => phase relaxation engaged

    auto const factored = find_factored(res);
    REQUIRE(factored);
    ExprPtr shared, partner_sum;
    for (auto const& f : factored->as<Product>().factors())
      (f->is<Sum>() ? partner_sum : shared) = f;
    // The shared factor is the G*W intermediate (a Product), not a bare leaf.
    REQUIRE(shared);
    REQUIRE(shared->is<Product>());
    REQUIRE(partner_sum);
    REQUIRE(partner_sum->size() == 2);

    // The relative canonicalization sign landed on exactly one partner, even
    // though neither input summand carried a sign or scalar -- the (X - Y)
    // shape. side_expr() emits a +1 partner as a bare Tensor and a non-unit one
    // as a scalar Product, so exactly one partner is a negative Product.
    std::size_t negative_partners = 0;
    for (auto const& g : partner_sum->as<Sum>().summands())
      if (g->is<Product>() && g->as<Product>().scalar().real() < 0)
        ++negative_partners;
    REQUIRE(negative_partners == 1);

    // Round-trip is the cost-model-independent correctness guard.
    REQUIRE_THAT(res, EquivalentTo(expr));
  }

  SECTION(
      "three-factor input yields a single 2-sided fold, not nested factors") {
    // (A + B)(C + D)(E + F) fully expanded into 8 terms. The factorizer only
    // emits 2-sided bicliques and never re-factors the partner sums, so it
    // produces a single fold (AC + AD + BC + BD)*(E + F) -- NOT the nested
    // 3-factor form. This locks the documented scope limit (see multiterm.hpp).
    // Non-antisymmetric tensors so the eight terms stay independent.
    auto const expr = parse_nonsymm(
        L"A{i1,i2;a1,a2} C{a1,a2;a3,a4} E{a3,a4;i3,i4}"
        L" + A{i1,i2;a1,a2} D{a1,a2;a3,a4} E{a3,a4;i3,i4}"
        L" + B{i1,i2;a1,a2} C{a1,a2;a3,a4} E{a3,a4;i3,i4}"
        L" + B{i1,i2;a1,a2} D{a1,a2;a3,a4} E{a3,a4;i3,i4}"
        L" + A{i1,i2;a1,a2} C{a1,a2;a3,a4} F{a3,a4;i3,i4}"
        L" + A{i1,i2;a1,a2} D{a1,a2;a3,a4} F{a3,a4;i3,i4}"
        L" + B{i1,i2;a1,a2} C{a1,a2;a3,a4} F{a3,a4;i3,i4}"
        L" + B{i1,i2;a1,a2} D{a1,a2;a3,a4} F{a3,a4;i3,i4}");
    auto const res = optimize(expr, mt);
    REQUIRE_THAT(res, EquivalentTo(expr));  // correctness guard
    REQUIRE(res->size() == 1);
    auto const factored = find_factored(res);
    REQUIRE(factored);

    // Exactly the two sides of one biclique: a 4-term left sum (still
    // unfactored) and the 2-term (E + F). Both factors are Sums.
    std::vector<std::size_t> sum_sizes;
    for (auto const& f : factored->as<Product>().factors()) {
      REQUIRE(f->is<Sum>());
      sum_sizes.push_back(f->size());
      // No partner is itself factored: no nested 3-factor form.
      for (auto const& g : f->as<Sum>().summands())
        if (g->is<Product>())
          for (auto const& h : g->as<Product>().factors())
            REQUIRE_FALSE(h->is<Sum>());
    }
    std::sort(sum_sizes.begin(), sum_sizes.end());
    REQUIRE(sum_sizes == std::vector<std::size_t>{2, 4});
  }
}
