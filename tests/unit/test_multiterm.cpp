#include <catch2/catch_test_macros.hpp>

#include "catch2_sequant.hpp"

#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/expressions/expr_algorithms.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/optimize/optimize.hpp>
#include <SeQuant/core/space.hpp>

#include <cstddef>
#include <string>
#include <vector>

namespace {

/// Set the occupied (i) and virtual (a) extents so the cost model is
/// well-defined, then run \p body. Tests whose expected fold depends on the
/// occ-vs-virt ordering pass explicit extents rather than coupling to the
/// defaults (see the cost-driven-winner test).
template <typename F>
void with_sized_spaces(std::size_t occ, std::size_t virt, F&& body) {
  using namespace sequant;
  auto ctx_resetter = set_scoped_default_context(get_default_context().clone());
  auto reg = get_default_context().mutable_index_space_registry();
  reg->retrieve_ptr(L"i")->approximate_size(occ);
  reg->retrieve_ptr(L"a")->approximate_size(virt);
  std::forward<F>(body)();
}

template <typename F>
void with_sized_spaces(F&& body) {
  with_sized_spaces(/*occ=*/10, /*virt=*/20, std::forward<F>(body));
}

sequant::ExprPtr parse_antisymm(std::wstring_view s) {
  using namespace sequant;
  return deserialize(s, {.def_perm_symm = Symmetry::Antisymm});
}

std::wstring latex(sequant::ExprPtr const& e) { return sequant::to_latex(e); }

/// A summand of \p sum that is a Product whose factors include a Sum -- i.e. a
/// factored group  shared * (partner1 + partner2 + ...). Returns the matching
/// Product, or null if none.
sequant::ExprPtr find_factored(sequant::ExprPtr const& sum) {
  using namespace sequant;
  if (!sum->is<Sum>()) return {};
  for (auto const& s : sum->as<Sum>().summands()) {
    if (!s->is<Product>()) continue;
    for (auto const& f : s->as<Product>().factors())
      if (f->is<Sum>()) return s;
  }
  return {};
}

/// Every factored-group summand of \p sum (a Product carrying a Sum factor), in
/// summand order. Generalizes \ref find_factored, which returns only the first.
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

}  // namespace

TEST_CASE("multiterm scaffolding", "[multiterm]") {
  using namespace sequant;

  with_sized_spaces([] {
    auto const expr = parse_antisymm(
        L"G{i1,i2;a1,a2} T{a1,a2;i1,i2}"
        L" + G{i1,i2;a1,a2} Z{a1,a2;i1,i2}");

    SECTION("disabled is the default and a no-op vs explicit Disable") {
      auto const def = optimize(expr);
      auto const dis = optimize(expr, {.multiterm = MultiTermFactor::Disable});
      REQUIRE(latex(def) == latex(dis));
    }

    SECTION("enabling the flag does not throw and yields a Sum") {
      ExprPtr res;
      REQUIRE_NOTHROW(
          res = optimize(expr, {.multiterm = MultiTermFactor::Enable}));
      REQUIRE(res);
      REQUIRE(res->is<Sum>());
    }
  });
}

TEST_CASE("multiterm one-sided fold", "[multiterm]") {
  using namespace sequant;

  with_sized_spaces([] {
    SECTION("A*B + A*C -> A*(B + C) when it lowers cost") {
      // V contracted with T (or B) over a3,a4 leaves a tensor R{a1,a2;i1,i2};
      // the contraction (O(a^4 i^2)) dwarfs the partner-sum (O(a^2 i^2)), so
      // folding pays off.
      auto const expr = parse_antisymm(
          L"V{a1,a2;a3,a4} T{a3,a4;i1,i2}"
          L" + V{a1,a2;a3,a4} B{a3,a4;i1,i2}");

      auto const res = optimize(expr, {.reorder = ReorderSum::NoReorder,
                                       .multiterm = MultiTermFactor::Enable});
      REQUIRE(res->is<Sum>());
      // Everything folded into a single product shared * (T + B).
      REQUIRE(res->size() == 1);
      auto const factored = find_factored(res);
      REQUIRE(factored);
      auto const& prod = factored->as<Product>();
      REQUIRE(prod.factors().size() == 2);
      // one factor is the shared tensor V, the other the partner-sum (2 terms)
      ExprPtr shared, partner_sum;
      for (auto const& f : prod.factors())
        (f->is<Sum>() ? partner_sum : shared) = f;
      REQUIRE(shared);
      REQUIRE(shared->is<Tensor>());
      REQUIRE(shared->as<Tensor>().label() == L"V");
      REQUIRE(partner_sum);
      REQUIRE(partner_sum->size() == 2);

      // Round-trip: expanding the factored form reproduces the original sum.
      REQUIRE_THAT(expand(res->clone()), EquivalentTo(expr));
    }

    SECTION("N-ary fold A*B + A*C + A*D -> A*(B + C + D)") {
      auto const expr = parse_antisymm(
          L"V{a1,a2;a3,a4} T{a3,a4;i1,i2}"
          L" + V{a1,a2;a3,a4} B{a3,a4;i1,i2}"
          L" + V{a1,a2;a3,a4} U{a3,a4;i1,i2}");

      auto const res = optimize(expr, {.reorder = ReorderSum::NoReorder,
                                       .multiterm = MultiTermFactor::Enable});
      REQUIRE(res->size() == 1);
      auto const factored = find_factored(res);
      REQUIRE(factored);
      ExprPtr partner_sum;
      for (auto const& f : factored->as<Product>().factors())
        if (f->is<Sum>()) partner_sum = f;
      REQUIRE(partner_sum);
      REQUIRE(partner_sum->size() == 3);
    }
  });
}

TEST_CASE("multiterm pass-through of non-splittable summands", "[multiterm]") {
  using namespace sequant;

  with_sized_spaces([] {
    SECTION(
        "a bare-tensor summand is left untouched while its neighbors fold") {
      // V*T + V*B fold to V*(T + B); the bare tensor D{a1,a2;i1,i2} is a leaf
      // (no contraction to split), so extract_core() returns nullopt and it is
      // never interned, scored, or consumed. Reassembly emits untouched
      // summands first (original order), then the folds, so the result is
      // [D, V*(T + B)].
      auto const expr = parse_antisymm(
          L"V{a1,a2;a3,a4} T{a3,a4;i1,i2}"
          L" + V{a1,a2;a3,a4} B{a3,a4;i1,i2}"
          L" + D{a1,a2;i1,i2}");

      auto const res = optimize(expr, {.reorder = ReorderSum::NoReorder,
                                       .multiterm = MultiTermFactor::Enable});
      REQUIRE(res->is<Sum>());
      REQUIRE(res->size() == 2);  // the untouched bare tensor + one fold

      // The fold is present, with a 2-term partner sum (T + B).
      auto const factored = find_factored(res);
      REQUIRE(factored);
      ExprPtr partner_sum;
      for (auto const& f : factored->as<Product>().factors())
        if (f->is<Sum>()) partner_sum = f;
      REQUIRE(partner_sum);
      REQUIRE(partner_sum->size() == 2);

      // The bare tensor survived verbatim as a standalone summand (not folded
      // in), and -- with NoReorder -- precedes the fold in the output.
      auto const& smands = res->as<Sum>().summands();
      REQUIRE(smands.front()->is<Tensor>());
      REQUIRE(smands.front()->as<Tensor>().label() == L"D");

      REQUIRE_THAT(expand(res->clone()), EquivalentTo(expr));
    }
  });
}

TEST_CASE("multiterm two-sided biclique", "[multiterm]") {
  using namespace sequant;

  with_sized_spaces([] {
    SECTION("(A+B)*(X+Y) recovered from AX + AY + BX + BY") {
      // U,V share a shape (left factors); T,B share a shape (right factors).
      // The four contractions form a complete 2x2 bipartite graph and fold to
      // (U + V) * (T + B).
      auto const expr = parse_antisymm(
          L"U{a1,a2;a3,a4} T{a3,a4;i1,i2}"
          L" + U{a1,a2;a3,a4} B{a3,a4;i1,i2}"
          L" + V{a1,a2;a3,a4} T{a3,a4;i1,i2}"
          L" + V{a1,a2;a3,a4} B{a3,a4;i1,i2}");

      auto const res = optimize(expr, {.reorder = ReorderSum::NoReorder,
                                       .multiterm = MultiTermFactor::Enable});
      REQUIRE(res->is<Sum>());
      REQUIRE(res->size() == 1);  // all four folded into one product

      auto const factored = find_factored(res);
      REQUIRE(factored);
      auto const& prod = factored->as<Product>();
      REQUIRE(prod.factors().size() == 2);
      // Both factors are 2-term sums: (U + V) and (T + B).
      for (auto const& f : prod.factors()) {
        REQUIRE(f->is<Sum>());
        REQUIRE(f->size() == 2);
      }

      // Round-trip: expanding the 2x2 fold reproduces the original four terms.
      REQUIRE_THAT(expand(res->clone()), EquivalentTo(expr));
    }
  });
}

TEST_CASE("multiterm incomplete graph: cost-driven winner flips with extents",
          "[multiterm]") {
  using namespace sequant;

  // Topology: edges (U,T),(U,B),(V,T) -- no (V,B). Two maximal one-sided folds
  // compete, with identical avoided contraction cost C:
  //   {U}x{T,B}  -> U*(T + B),   build cost ~ size(T) = a^2 i^2  (sum T, B)
  //   {U,V}x{T}  -> (U + V)*T,   build cost ~ size(U) = a^4      (sum U, V)
  // Greedy keeps the higher-saving fold -- the one whose summed side is cheaper
  // to build: it sums the lighter factors and shares the heavier one. Which
  // factor is heavier is set purely by the occ-vs-virt extent ordering, so the
  // same expression under both orderings must flip the shared factor (U <-> T).
  // The flip shows the choice is cost-driven, not hardcoded; the round-trip is
  // the cost-model-independent guard and holds in both.
  std::wstring const expr_text =
      L"U{a1,a2;a3,a4} T{a3,a4;i1,i2}"
      L" + U{a1,a2;a3,a4} B{a3,a4;i1,i2}"
      L" + V{a1,a2;a3,a4} T{a3,a4;i1,i2}";

  // The shared (non-summed) factor's tensor label, with the structural guards
  // that hold regardless of which side wins: exactly one fold over a 2-term
  // partner sum, plus one untouched leftover summand (first, under NoReorder).
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

  std::wstring shared_virt_heavy, shared_occ_heavy;

  // virt-heavy (a=20 > i=10): size(U)=a^4 dominates size(T)=a^2 i^2, so U is
  // the heavier factor and is shared -- summing the lighter T, B into U*(T +
  // B).
  with_sized_spaces(/*occ=*/10, /*virt=*/20, [&] {
    auto const expr = parse_antisymm(expr_text);
    auto const res = optimize(expr, {.reorder = ReorderSum::NoReorder,
                                     .multiterm = MultiTermFactor::Enable});
    shared_virt_heavy = winning_shared_label(res);
    REQUIRE_THAT(expand(res->clone()), EquivalentTo(expr));
  });

  // occ-heavy (i=30 > a=10): now size(T)=a^2 i^2 dominates size(U)=a^4, so T is
  // the heavier factor and is shared -- summing the lighter U, V into (U +
  // V)*T.
  with_sized_spaces(/*occ=*/30, /*virt=*/10, [&] {
    auto const expr = parse_antisymm(expr_text);
    auto const res = optimize(expr, {.reorder = ReorderSum::NoReorder,
                                     .multiterm = MultiTermFactor::Enable});
    shared_occ_heavy = winning_shared_label(res);
    REQUIRE_THAT(expand(res->clone()), EquivalentTo(expr));
  });

  REQUIRE(shared_virt_heavy == L"U");
  REQUIRE(shared_occ_heavy == L"T");
  REQUIRE(shared_virt_heavy != shared_occ_heavy);  // the cost-driven flip
}

TEST_CASE("multiterm cost decline", "[multiterm]") {
  using namespace sequant;

  with_sized_spaces([] {
    SECTION("scalar energy: shareable but folding saves nothing") {
      // G*T + G*Z both fully contract to a scalar -- exactly the saving()==0
      // boundary, and structurally so (not by extent tuning): for a full
      // contraction the avoided cost C equals the product of the whole index
      // set, and the only build cost is the partner sum (T + Z), whose
      // footprint size(T) is that same whole index set. With the one-sided
      // saving = (m*n - 1)*C - (n - 1)*size(T) at m=1, n=2 this is C - size(T)
      // = 0 for any extents, so nothing folds. The cancellation is exact only
      // under the current (m*n - 1)/(n - 1) coefficients (CostModel::saving,
      // multiterm.cpp); revisit the no-fold check if those change. The
      // partial-contraction sibling below brackets this boundary from above.
      auto const expr = parse_antisymm(
          L"G{i1,i2;a1,a2} T{a1,a2;i1,i2}"
          L" + G{i1,i2;a1,a2} Z{a1,a2;i1,i2}");

      auto const res = optimize(expr, {.reorder = ReorderSum::NoReorder,
                                       .multiterm = MultiTermFactor::Enable});
      REQUIRE(res->is<Sum>());
      REQUIRE(res->size() == 2);
      REQUIRE_FALSE(find_factored(res));
    }

    SECTION("partial contraction: the same sharing now clears the threshold") {
      // Contrast with the scalar case: V shares with T and Z, but the
      // contraction is partial (over a3,a4 only), so the result keeps free
      // indices and the avoided cost C = O(a^4 i^2) strictly exceeds the
      // partner-build size(T) = O(a^2 i^2). saving = C - size(T) > 0, so
      // V*(T + Z) folds. With the scalar case sitting exactly on saving()==0,
      // this brackets that boundary from the positive side.
      auto const expr = parse_antisymm(
          L"V{a1,a2;a3,a4} T{a3,a4;i1,i2}"
          L" + V{a1,a2;a3,a4} Z{a3,a4;i1,i2}");

      auto const res = optimize(expr, {.reorder = ReorderSum::NoReorder,
                                       .multiterm = MultiTermFactor::Enable});
      REQUIRE(res->is<Sum>());
      REQUIRE(res->size() == 1);  // V * (T + Z)
      auto const factored = find_factored(res);
      REQUIRE(factored);
      ExprPtr partner_sum;
      for (auto const& f : factored->as<Product>().factors())
        if (f->is<Sum>()) partner_sum = f;
      REQUIRE(partner_sum);
      REQUIRE(partner_sum->size() == 2);
      REQUIRE_THAT(expand(res->clone()), EquivalentTo(expr));
    }
  });
}

TEST_CASE("multiterm DenseSize objective", "[multiterm]") {
  using namespace sequant;

  with_sized_spaces([] {
    SECTION("a fold that pays under DenseSize too") {
      // Unlike the other sections (default DenseFLOPs metric), this one sets
      // OptimizeOptions::objective_function to DenseSize, routing
      // CostModel::contraction_cost through memsize_counter instead of
      // flops_counter (multiterm.cpp). V*T + V*B still folds to V*(T + B):
      // under DenseSize the avoided contraction's element footprint (O(a^4),
      // the V leg dominating) far exceeds the cost of building the (T + B)
      // partner sum, so the saving stays positive.
      auto const expr = parse_antisymm(
          L"V{a1,a2;a3,a4} T{a3,a4;i1,i2}"
          L" + V{a1,a2;a3,a4} B{a3,a4;i1,i2}");

      auto const res =
          optimize(expr, {.objective_function = ObjectiveFunction::DenseSize,
                          .reorder = ReorderSum::NoReorder,
                          .multiterm = MultiTermFactor::Enable});
      REQUIRE(res->is<Sum>());
      REQUIRE(res->size() == 1);  // V * (T + B)
      auto const factored = find_factored(res);
      REQUIRE(factored);
      ExprPtr shared, partner_sum;
      for (auto const& f : factored->as<Product>().factors())
        (f->is<Sum>() ? partner_sum : shared) = f;
      REQUIRE(shared);
      REQUIRE(shared->is<Tensor>());
      REQUIRE(shared->as<Tensor>().label() == L"V");
      REQUIRE(partner_sum);
      REQUIRE(partner_sum->size() == 2);
      REQUIRE_THAT(expand(res->clone()), EquivalentTo(expr));
    }
  });
}

TEST_CASE("multiterm honors reorder", "[multiterm]") {
  using namespace sequant;

  // reorder and multiterm are independent options, so enabling multiterm must
  // not disable reorder. The input below keeps the two effects separable: its
  // two g*t summands fully contract to scalars, so folding saves nothing (the
  // saving()==0 boundary from "multiterm cost decline") and the factorizer
  // leaves all three summands alone. reorder, on the other hand, does move
  // them: it groups the two g*t summands and shifts the unrelated p*q term. So
  // "reorder ran" and "reorder skipped" give different output here.
  with_sized_spaces([] {
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

    // Compare serialized forms (order-sensitive, like elsewhere in this file).
    // With the bug, reorder got skipped under multiterm and these two came out
    // identical; fixed, they differ.
    REQUIRE(latex(mt_noreorder) != latex(mt_reorder));

    // And the order matches plain reorder: with nothing to fold, multiterm
    // leaves reorder's result untouched.
    auto const reorder_only =
        optimize(expr, {.reorder = ReorderSum::Reorder,
                        .multiterm = MultiTermFactor::Disable});
    REQUIRE(latex(mt_reorder) == latex(reorder_only));
  });
}

TEST_CASE("multiterm bucketing soundness", "[multiterm]") {
  using namespace sequant;

  with_sized_spaces([] {
    SECTION("different contraction signatures fold independently, not merged") {
      // All four summands share the external indices {a1,a2,i1,i2}. The first
      // two contract over a virtual pair (a3,a4); the last two over an occupied
      // pair (i3,i4). Distinct signatures => two independent buckets => two
      // separate folds (V*(T+B)) and (W*(X+Y)), never one cross-bucket merge.
      auto const expr = parse_antisymm(
          L"V{a1,a2;a3,a4} T{a3,a4;i1,i2}"
          L" + V{a1,a2;a3,a4} B{a3,a4;i1,i2}"
          L" + W{a1,a2;i3,i4} X{i3,i4;i1,i2}"
          L" + W{a1,a2;i3,i4} Y{i3,i4;i1,i2}");

      auto const res = optimize(expr, {.reorder = ReorderSum::NoReorder,
                                       .multiterm = MultiTermFactor::Enable});
      REQUIRE(res->is<Sum>());
      // Two folded products, one per bucket.
      REQUIRE(res->size() == 2);
      auto const folded = find_all_factored(res);
      REQUIRE(folded.size() == 2);  // both summands are independent folds
      for (auto const& f : folded) {
        ExprPtr partner_sum;
        for (auto const& g : f->as<Product>().factors())
          if (g->is<Sum>()) partner_sum = g;
        REQUIRE(partner_sum);
        REQUIRE(partner_sum->size() == 2);  // exactly the two partners
      }
    }
  });
}

TEST_CASE("multiterm scalar and sign folds", "[multiterm]") {
  using namespace sequant;

  // Each section's round-trip (expand the fold, compare to the original) is the
  // correctness guard: signs and prefactors must survive folding intact.
  with_sized_spaces([] {
    SECTION("A*B - A*C -> A*(B - C)") {
      // The factorizer peels the (-1) prefactor off the summand and folds the
      // relative sign into the partner.
      auto const expr = parse_antisymm(
          L"V{a1,a2;a3,a4} T{a3,a4;i1,i2}"
          L" - V{a1,a2;a3,a4} B{a3,a4;i1,i2}");

      auto const res = optimize(expr, {.reorder = ReorderSum::NoReorder,
                                       .multiterm = MultiTermFactor::Enable});
      REQUIRE(res->is<Sum>());
      REQUIRE(res->size() == 1);  // V * (T - B)
      auto const factored = find_factored(res);
      REQUIRE(factored);
      ExprPtr shared, partner_sum;
      for (auto const& f : factored->as<Product>().factors())
        (f->is<Sum>() ? partner_sum : shared) = f;
      REQUIRE(shared);
      REQUIRE(shared->is<Tensor>());
      REQUIRE(shared->as<Tensor>().label() == L"V");
      REQUIRE(partner_sum);
      REQUIRE(partner_sum->size() == 2);  // T and (-B)
      REQUIRE_THAT(expand(res->clone()), EquivalentTo(expr));
    }

    SECTION("rational prefactors: 2*A*B + 3*A*C -> A*(2B + 3C)") {
      auto const expr = parse_antisymm(
          L"2 V{a1,a2;a3,a4} T{a3,a4;i1,i2}"
          L" + 3 V{a1,a2;a3,a4} B{a3,a4;i1,i2}");

      auto const res = optimize(expr, {.reorder = ReorderSum::NoReorder,
                                       .multiterm = MultiTermFactor::Enable});
      REQUIRE(res->size() == 1);
      auto const factored = find_factored(res);
      REQUIRE(factored);
      ExprPtr partner_sum;
      for (auto const& f : factored->as<Product>().factors())
        if (f->is<Sum>()) partner_sum = f;
      REQUIRE(partner_sum);
      REQUIRE(partner_sum->size() == 2);
      REQUIRE_THAT(expand(res->clone()), EquivalentTo(expr));
    }

    SECTION("rank-1 two-sided sign matrix folds as (U - V)*(T + B)") {
      // U*T + U*B - V*T - V*B. The sign matrix [[+,+],[-,-]] is multiplicative
      // rank 1, so the whole 2x2 still folds to a single product, the sign
      // riding the left factor-sum.
      auto const expr = parse_antisymm(
          L"U{a1,a2;a3,a4} T{a3,a4;i1,i2}"
          L" + U{a1,a2;a3,a4} B{a3,a4;i1,i2}"
          L" - V{a1,a2;a3,a4} T{a3,a4;i1,i2}"
          L" - V{a1,a2;a3,a4} B{a3,a4;i1,i2}");

      auto const res = optimize(expr, {.reorder = ReorderSum::NoReorder,
                                       .multiterm = MultiTermFactor::Enable});
      REQUIRE(res->is<Sum>());
      REQUIRE(res->size() == 1);  // single (U - V)*(T + B)
      auto const factored = find_factored(res);
      REQUIRE(factored);
      for (auto const& f : factored->as<Product>().factors()) {
        REQUIRE(f->is<Sum>());
        REQUIRE(f->size() == 2);
      }
      REQUIRE_THAT(expand(res->clone()), EquivalentTo(expr));
    }

    SECTION("non-rank-1 sign matrix reduces to one-sided folds") {
      // U*T + U*B + V*T - V*B. The sign matrix [[+,+],[+,-]] does not factor,
      // so a single (U +/- V)*(T +/- B) would be wrong. The factorizer must
      // decline the 2x2 and reduce to one-sided folds: U*(T + B) and
      // V*(T - B). Round-trip is the correctness guard.
      auto const expr = parse_antisymm(
          L"U{a1,a2;a3,a4} T{a3,a4;i1,i2}"
          L" + U{a1,a2;a3,a4} B{a3,a4;i1,i2}"
          L" + V{a1,a2;a3,a4} T{a3,a4;i1,i2}"
          L" - V{a1,a2;a3,a4} B{a3,a4;i1,i2}");

      auto const res = optimize(expr, {.reorder = ReorderSum::NoReorder,
                                       .multiterm = MultiTermFactor::Enable});
      REQUIRE(res->is<Sum>());
      REQUIRE(res->size() == 2);  // two one-sided folds, no bogus 2x2
      auto const folded = find_all_factored(res);
      REQUIRE(folded.size() == 2);  // both summands are one-sided folds
      for (auto const& f : folded) {
        ExprPtr partner_sum;
        for (auto const& g : f->as<Product>().factors())
          if (g->is<Sum>()) partner_sum = g;
        REQUIRE(partner_sum);
        REQUIRE(partner_sum->size() == 2);
      }
      REQUIRE_THAT(expand(res->clone()), EquivalentTo(expr));
    }
  });
}

TEST_CASE("multiterm canonicalization-phase fold", "[multiterm]") {
  using namespace sequant;

  // Distinct from the scalar-prefactor folds above: there the relative sign was
  // an explicit input coefficient that extract_core peels off (the
  // PrefactoredContraction path). Here neither summand carries any sign or
  // scalar -- every relative sign comes purely from canonicalization.
  // Two factors that are structurally identical but differ by a
  // canonicalization phase must still share one vertex (FactorEq is
  // phase-relaxed), with the relative sign folded onto the edge coefficient
  // (sigma_l/sigma_r in build_buckets, multiterm.cpp). This is the path the
  // scalar cases never reach: canon_phase(factor) != canon_phase(rep).
  with_sized_spaces([] {
    SECTION(
        "antisymmetric index swap inside a shared intermediate folds with a "
        "relative sign") {
      // Both summands share the intermediate M = G*W, contracted over the
      // virtual pair a3,a4 (internal to M) and leaving the cheap occupied
      // result {i1,i2;i3,i4}. With virt > occ this virtual contraction is
      // single_term_opt's first step, so M is a child factor (not a leaf) of
      // the top contraction M*X / M*Y -- only an intermediate can carry a
      // non-trivial canon_phase, so this is what makes the sigma path reachable
      // at all (a leaf factor always canonicalizes to phase +1).
      //
      // The second summand writes G's contracted pair swapped (a4,a3). G is
      // antisymmetric, so M canonicalizes to the same canon_indices but the
      // opposite phase. Phase-relaxed FactorEq interns both Ms into one vertex;
      // the relative -1 rides the edge and lands on the partner: M*(X - Y).
      //
      // size() == 1 is itself the proof that the relaxation engaged: a
      // phase-strict matcher would see two distinct intermediates, fold
      // nothing, and leave size() == 2.
      auto const expr = parse_antisymm(
          L"G{i1,i2;a3,a4} W{a3,a4;i3,i4} X{i3,i4;a1,a2}"
          L" + G{i1,i2;a4,a3} W{a3,a4;i3,i4} Y{i3,i4;a1,a2}");

      auto const res = optimize(expr, {.reorder = ReorderSum::NoReorder,
                                       .multiterm = MultiTermFactor::Enable});
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
      // shape. side_expr() emits a +1 partner as a bare Tensor and a non-unit
      // one as a scalar Product, so exactly one partner is a negative Product.
      std::size_t negative_partners = 0;
      for (auto const& g : partner_sum->as<Sum>().summands())
        if (g->is<Product>() && g->as<Product>().scalar().real() < 0)
          ++negative_partners;
      REQUIRE(negative_partners == 1);

      // Round-trip is the cost-model-independent correctness guard: expanding
      // M*(X - Y) must reproduce the original two summands, the swapped-index
      // (and thus sign-flipped) second term included.
      REQUIRE_THAT(expand(res->clone()), EquivalentTo(expr));
    }
  });
}
