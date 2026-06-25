#include <catch2/catch_test_macros.hpp>

#include "catch2_sequant.hpp"

#include <SeQuant/core/algorithm.hpp>
#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/eval/eval_node.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/optimize/common_subexpression_elimination.hpp>
#include <SeQuant/core/optimize/optimize.hpp>
#include <SeQuant/core/optimize/single_term.hpp>
#include <SeQuant/core/runtime.hpp>
#include <SeQuant/core/space.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>

#include <algorithm>
#include <bit>
#include <cstddef>
#include <functional>
#include <initializer_list>
#include <limits>
#include <memory>
#include <vector>

sequant::ExprPtr extract(sequant::ExprPtr expr,
                         std::initializer_list<size_t> const& idxs) {
  using namespace sequant;
  ExprPtr result = expr;
  for (auto s : idxs) result = result->at(s);
  return result;
}

// number of Tensor leaves in a (binarized) expression tree
size_t count_tensor_leaves(sequant::ExprPtr const& expr) {
  using namespace sequant;
  size_t n = 0;
  expr->visit([&n](auto const& x) { n += x->template is<Tensor>() ? 1 : 0; },
              /*atoms_only=*/true);
  return n;
}

// Minimum peak memory over ALL pairwise contraction sequences of `nt` leaves,
// using S[subset] (subset_footprints) for sizes. Independent oracle for the
// peak DP: enumerates schedules and simulates memory; no recurrence assumed.
static double brute_force_min_peak(std::vector<double> const& S, size_t nt) {
  double const full = static_cast<double>((size_t{1} << nt) - 1);
  // `live` is the set of subset-masks currently materialized (a partition of
  // the full set). Recurse over every pair to merge.
  std::function<double(std::vector<size_t>)> rec =
      [&](std::vector<size_t> live) -> double {
    if (live.size() == 1) return 0.0;  // nothing more to compute
    double live_sum = 0.0;
    for (auto m : live) live_sum += S[m];
    double best = std::numeric_limits<double>::max();
    for (size_t i = 0; i < live.size(); ++i)
      for (size_t j = i + 1; j < live.size(); ++j) {
        size_t merged = live[i] | live[j];
        // instantaneous peak when forming `merged`: all live results plus the
        // new result momentarily co-resident.
        double step_peak = live_sum + S[merged];
        std::vector<size_t> next;
        next.reserve(live.size() - 1);
        for (size_t k = 0; k < live.size(); ++k)
          if (k != i && k != j) next.push_back(live[k]);
        next.push_back(merged);
        best = std::min(best, std::max(step_peak, rec(next)));
      }
    return best;
  };
  std::vector<size_t> leaves;
  for (size_t b = 0; b < nt; ++b) leaves.push_back(size_t{1} << b);
  (void)full;
  return rec(std::move(leaves));
}

// Per-index batch-aware peak oracle.
// T[ctx][subset]: sliced footprints (ctx = bitmask of currently-sliced aux
// indices, subset = bitmask of leaves). open_aux[subset] = bitmask of aux
// indices that are open (external) in that subset. persistent[subset] = 1 iff
// the contraction result at that subset is persistent. nt = number of leaves.
// B = bitmask of aux indices currently being sliced (threading the context).
static double oracle_rec(size_t n, size_t B,
                         std::vector<std::vector<double>> const& T,
                         std::vector<size_t> const& open_aux,
                         std::vector<char> const& persistent, size_t nt) {
  if (std::popcount(n) == 1) return T[B & open_aux[n]][n];
  auto sz = [&](size_t s, size_t ctx) { return T[ctx & open_aux[s]][s]; };
  auto Lof = [&](size_t s, size_t ctx) {
    double r = 0;
    for (size_t b = 0; b < nt; b++)
      if (s & (size_t{1} << b)) r += sz(size_t{1} << b, ctx);
    return r;
  };
  double best = std::numeric_limits<double>::max();
  for (size_t lp = (n - 1) & n; lp; lp = (lp - 1) & n) {
    size_t rp = n ^ lp;
    if (lp > rp) continue;
    // aux contracted at THIS node = open at children but not open at parent
    size_t Acand =
        persistent[n] ? ((open_aux[lp] | open_aux[rp]) & ~open_aux[n]) : 0;
    for (size_t Ap = Acand;; Ap = (Ap - 1) & Acand) {
      size_t C = B | Ap;
      double pl = oracle_rec(lp, C, T, open_aux, persistent, nt);
      double pr = oracle_rec(rp, C, T, open_aux, persistent, nt);
      double both = sz(lp, C) + sz(rp, C) + sz(n, B);
      double lpf = std::max({Lof(rp, C) + pl, sz(lp, C) + pr, both});
      double rpf = std::max({Lof(lp, C) + pr, sz(rp, C) + pl, both});
      best = std::min(best, std::min(lpf, rpf));
      if (Ap == 0) break;
    }
  }
  return best;
}

static double batched_min_peak(std::vector<std::vector<double>> const& T,
                               std::vector<size_t> const& open_aux,
                               std::vector<char> const& persistent, size_t nt) {
  size_t full = (size_t{1} << nt) - 1;
  return oracle_rec(full, 0, T, open_aux, persistent, nt);
}

TEST_CASE("batched_min_peak oracle (per-index)", "[optimize]") {
  // 2 leaves sharing aux F (m=1). open_aux: leaves have F open (bit0=1);
  // the pair has F internal (bit0=0). persistent everywhere.
  // T[ctx][subset]: full (ctx=0) leaves 4, pair result 2; sliced (ctx=1)
  // leaves 2.
  std::vector<std::vector<double>> T(2, std::vector<double>(4, 0.0));
  T[0][0b01] = T[0][0b10] = 4;
  T[0][0b11] = 2;  // full
  T[1][0b01] = T[1][0b10] = 2;
  T[1][0b11] = 2;  // F sliced -> leaves halve
  std::vector<size_t> open_aux(4, 0);
  open_aux[0b01] = 1;
  open_aux[0b10] = 1;
  open_aux[0b11] = 0;  // F open in leaves, internal in pair
  std::vector<char> persistent(4, 1);
  // not batched: 4+4+2=10; batch F at the pair: leaves sliced 2 -> 2+2+2=6.
  // min 6.
  REQUIRE(batched_min_peak(T, open_aux, persistent, /*nt=*/2) == 6.0);
}

TEST_CASE("brute_force_min_peak oracle", "[optimize]") {
  // 3 leaves, hand-checkable sizes by subset mask:
  //   S[001]=2 S[010]=2 S[100]=2 (leaves)
  //   S[011]=1 S[101]=1 S[110]=8 S[111]=1 (pair/full results)
  std::vector<double> S(8, 0.0);
  S[0b001] = S[0b010] = S[0b100] = 2.0;
  S[0b011] = 1.0;
  S[0b101] = 1.0;
  S[0b110] = 8.0;
  S[0b111] = 1.0;
  // Memory model: step_peak = sum(all live subsets) + S[merged result].
  // With 3 leaves all live, live_sum = 2+2+2 = 6.
  // Schedule {001,010}->{011}: step_peak = 6+1 = 7; then {011,100}->{111}:
  //   live_sum=1+2=3, step_peak=3+1=4. Peak = max(7,4) = 7.
  // Schedule {001,100}->{101}: step_peak = 6+1 = 7; then {010,101}->{111}:
  //   live_sum=2+1=3, step_peak=3+1=4. Peak = max(7,4) = 7.
  // Schedule {010,100}->{110}: step_peak = 6+8 = 14. Peak >= 14.
  // min over all schedules = 7.
  REQUIRE(brute_force_min_peak(S, 3) == 7.0);
}

TEST_CASE("optimize", "[optimize]") {
  using namespace sequant;

  // for optimization tests, need to specify index space sizes, so make a clone
  // of the context
  {
    auto ctx_resetter =
        set_scoped_default_context(get_default_context().clone());
    auto reg = get_default_context().mutable_index_space_registry();

    {
      auto occ = reg->retrieve_ptr(L"i");
      auto uocc = reg->retrieve_ptr(L"a");
      auto aux = reg->retrieve_ptr(L"x");
      REQUIRE(occ);
      REQUIRE(uocc);
      REQUIRE(aux);
      occ->approximate_size(10);
      uocc->approximate_size(100);
      aux->approximate_size(4);
      REQUIRE(uocc->approximate_size() == 100);
    }

    auto single_term_opt = [](Product const& prod) {
      return opt::single_term_opt(prod, [](Index const& ix) {
        // null space contributes x1 to the size
        auto sz = ix.nonnull() ? ix.space().approximate_size() : 1;
        return sz;
      });
    };

    auto parse_expr_antisymm = [](auto const& xpr) {
      return deserialize(xpr, {.def_perm_symm = Symmetry::Antisymm});
    };

    SECTION("Single term optimization") {
      const auto prod1 = parse_expr_antisymm(
                             L"g_{i3,i4}^{a3,a4}"     // T1
                             " * t_{a1,a2}^{i3,i4}"   // T2
                             " * t_{a3,a4}^{i1,i2}")  // T3
                             ->as<Product>();
      //
      // Cost of evaluation prod1:
      //
      // ((T1 * T2) * T3)  : 2 * O^2 * V^4  best if nocc > nvirt
      //
      // this is the one we want to find
      // ((T1 * T3) * T2)  : 2 * O^4 * V^2  best if nvirt > nocc
      //
      // (T1 * (T2 * T3))  : 2 * O^4 * V^4  worst sequence of evaluation
      //

      const auto res1 = single_term_opt(prod1);

      REQUIRE(extract(res1, {0, 0}) == prod1.at(0));
      REQUIRE(extract(res1, {0, 1}) == prod1.at(2));
      REQUIRE(extract(res1, {1}) == prod1.at(1));

      const auto prod2 = parse_expr_antisymm(
                             L"   g_{i3,i4}^{a3,a4}"
                             L" * t_{a3,a4}^{i1,i2}"
                             L" * t_{a1}^{i3}"
                             L" * t_{a2}^{i4}")
                             ->as<Product>();

      const auto res2 = single_term_opt(prod2);

      REQUIRE(extract(res2, {0, 0, 0}) == prod2.at(0));
      REQUIRE(extract(res2, {0, 0, 1}) == prod2.at(1));
      REQUIRE(extract(res2, {0, 1}) == prod2.at(2));
      REQUIRE(extract(res2, {1}) == prod2.at(3));

      const auto prod3 = parse_expr_antisymm(
                             L""                   //
                             " g_{i3,i4}^{a3,a4}"  //
                             " t_{a1}^{i3}"        //
                             " t_{a2}^{i4}"        //
                             " t_{a3,a4}^{i1,i2}"  //
                             )
                             ->as<Product>();
      auto res3 = single_term_opt(prod3);

      REQUIRE(extract(res3, {0, 0, 0}) == prod3.at(0));
      REQUIRE(extract(res3, {0, 0, 1}) == prod3.at(3));
      REQUIRE(extract(res3, {0, 1}) == prod3.at(1));
      REQUIRE(extract(res3, {1}) == prod3.at(2));

      //
      // single-term optimization when a dot product occurs in the tensor
      // network
      // ========================

      auto prod4 =
          parse_expr_antisymm(L"1/4 λ{i1;a1} g{i2,i3;a2,a3} t{a2,a3;i2,i3}")
              ->as<Product>();
      auto res4 = single_term_opt(prod4);

      REQUIRE(extract(res4, {0}) == prod4.at(0));
      REQUIRE(extract(res4, {1, 0}) == prod4.at(1));
      REQUIRE(extract(res4, {1, 1}) == prod4.at(2));

      auto prod5 =
          parse_expr_antisymm(L"x{i1,i2;a3,a4} y{a1,a2;i1,i2} z{a3,a4;a1,a2}")
              ->as<Product>();
      auto res5 = single_term_opt(prod5);
      REQUIRE(extract(res5, {0, 0}) == prod5.at(0));
      REQUIRE(extract(res5, {0, 1}) == prod5.at(2));
      REQUIRE(extract(res5, {1}) == prod5.at(1));

      //
      // single-term optimization when sequant::Variables appear in a product
      //
      auto prod6 = deserialize(
                       L"α * β * γ * "
                       "g_{i3,i4}^{a3,a4}"      // T1
                       " * t_{a1,a2}^{i3,i4}"   // T2
                       " * t_{a3,a4}^{i1,i2}")  // T3
                       ->as<Product>();
      auto res6 = single_term_opt(prod6);

      // this is the one we want to find
      // α * β * γ * ((T1 * T3) * T2)  : 2 * O^4 * V^2  best if nvirt > nocc
      REQUIRE(extract(res6, {0}) == prod6.at(0));
      REQUIRE(extract(res6, {1}) == prod6.at(1));
      REQUIRE(extract(res6, {2}) == prod6.at(2));
      REQUIRE(extract(res6, {3, 0}) == prod6.at(3));
      REQUIRE(extract(res6, {3, 1}) == prod6.at(5));
      REQUIRE(extract(res6, {4}) == prod6.at(4));

      //
      // single-term optimization including tensors with auxiliary indices
      //
      auto prod7 = deserialize(
                       L"DF{a_1;a_3;x_1} "  // T1
                       "DF{a_2;i_1;x_1} "   // T2
                       "t{a_3;i_2}"         // T3
                       )
                       ->as<Product>();
      auto res7 = single_term_opt(prod7);

      // this is the one we want to find
      // (T1 T3) T2: V^2 O^1 A^1 + V^2 O^2 A^1 best if nvirt > nocc and nvirt >
      // nact
      REQUIRE(extract(res7, {0, 0}) == prod7.at(0));
      REQUIRE(extract(res7, {0, 1}) == prod7.at(2));
      REQUIRE(extract(res7, {1}) == prod7.at(1));

      auto prod8 =
          deserialize(
              L"T1{i_1;i_2;x_1,x_2,x_3,x_4} T2{i_2;i_1;x_5,x_6,x_7,x_8} "
              L"T3{i_3;;x_1,x_2,x_3,x_4} T4{i_4;;x_5,x_6,x_7,x_8}")
              ->as<Product>();
      auto res8 = single_term_opt(prod8);

      // this is the one we want to find
      // (T1 T3)(T2 T4)
      REQUIRE(extract(res8, {0, 0}) == prod8.at(0));
      REQUIRE(extract(res8, {0, 1}) == prod8.at(2));
      REQUIRE(extract(res8, {1, 0}) == prod8.at(1));
      REQUIRE(extract(res8, {1, 1}) == prod8.at(3));
    }

    SECTION("Single term optimization: n_replay volatility weighting") {
      using namespace sequant;

      // PPL-shaped motif, fully contracted to a scalar:
      //   A = g_{i1,a1}^{x1}   (persistent integral)
      //   B = g_{i2,a2}^{x1}   (persistent integral)
      //   t = t_{a1,a2}^{i1,i2} (VOLATILE amplitude)
      // sizes: i=10 (O), a=100 (V), x=4 (X).
      //
      // (A*B)*t : build I=A*B over x  -> {i1,a1,i2,a2}  cost O^2 V^2 X
      // (persistent)
      //           then I*t            -> scalar         cost O^2 V^2 (volatile)
      // (A*t)*B : build J=A*t over i1,a1 -> {x,i2,a2}   cost O^2 V^2 X
      // (VOLATILE)
      //           then J*B               -> scalar      cost X O V (volatile)
      //
      // n_replay=1  : (A*t)*B wins (O^2 V^2 X + X O V  <  O^2 V^2 X + O^2 V^2)
      //               => t buried in an inner volatile intermediate.
      // n_replay=10 : (A*B)*t wins (persistent build counted once; the only
      //               x10 term is the cheap O^2 V^2 final step)
      //               => t contracted LAST, persistent integral formed first.
      auto idxsz = [](Index const& ix) {
        return ix.nonnull() ? ix.space().approximate_size() : std::size_t{1};
      };

      auto prod = parse_expr_antisymm(
                      L"g_{i1,a1}^{x1}"
                      L" * g_{i2,a2}^{x1}"
                      L" * t_{a1,a2}^{i1,i2}")
                      ->as<Product>();

      auto is_t = [](Tensor const& t) { return t.label() == L"t"; };

      OptimizeOptions base;
      base.idx_to_extent = idxsz;

      // baseline: predicate set but n_replay==1 => weight is 1 everywhere =>
      // reverts to current behavior. (The empty-predicate no-op is checked
      // separately via opts_off below.)
      auto opts1 = base;
      opts1.batch_policy.is_volatile_leaf = is_t;
      opts1.volatile_weight = 1;
      auto res1 = optimize(ex<Product>(prod), opts1);

      auto opts10 = base;
      opts10.batch_policy.is_volatile_leaf = is_t;
      opts10.volatile_weight = 10;
      auto res10 = optimize(ex<Product>(prod), opts10);

      // a bare top-level t leaf means t was contracted last (persistent-first)
      auto top_has_bare_t = [](ExprPtr const& e) {
        if (!e->is<Product>()) return false;
        for (auto const& c : *e)
          if (c->is<Tensor>() && c->as<Tensor>().label() == L"t") return true;
        return false;
      };

      // weighting flips the chosen factorization
      REQUIRE(res1 != res10);
      // volatile_weight=1 reproduces today's behavior: t buried in an inner
      // intermediate
      REQUIRE_FALSE(top_has_bare_t(res1));
      // volatile_weight=10: persistent g*g built first, t contracted last
      REQUIRE(top_has_bare_t(res10));

      // empty predicate => weighting off => identical to volatile_weight=1
      // regardless
      auto opts_off = base;
      opts_off.volatile_weight = 10;  // ignored: predicate empty
      auto res_off = optimize(ex<Product>(prod), opts_off);
      REQUIRE(res_off == res1);
    }

    SECTION("Single term optimization: footprint_weight") {
      using namespace sequant;

      // Network A{i1;x1} * B{x1;i2} * C{i2;a1}: contract x1 (between A,B) and
      // i2 (between B,C); free indices {i1, a1}. With the aux index x LARGE and
      // the virtual index a SMALL, the two viable orders trade FLOPs against
      // the footprint of the single intermediate they materialize:
      //
      //   (A*B)*C : I{i1;i2}  (occ^2 = 100)         FLOPs i*x*i + i*i*a =
      //   100400 A*(B*C) : I{x1;a1}  (aux*virt = 4000)     FLOPs x*i*a + i*x*a
      //   =  80000
      //
      // Pure FLOPs picks A*(B*C): cheaper, but materializes the big
      // aux-carrying intermediate I{x1;a1}. A nonzero footprint_weight
      // penalizes that 4000-element intermediate (vs the 100-element one) and
      // flips the choice to (A*B)*C. Flip threshold here is footprint_weight >
      // ~5.2.
      auto uocc = reg->retrieve_ptr(L"a");
      auto aux = reg->retrieve_ptr(L"x");
      auto const uocc_sz = uocc->approximate_size();
      auto const aux_sz = aux->approximate_size();
      uocc->approximate_size(4);    // virtual: deliberately SMALL
      aux->approximate_size(1000);  // aux: deliberately LARGE

      auto const prod =
          deserialize(L"A{i1;x1} B{x1;i2} C{i2;a1}")->as<Product>();

      // footprint_weight == 0 reproduces the pure-FLOPs choice ...
      auto res0 = optimize(ex<Product>(prod), OptimizeOptions{});
      auto res0_explicit =
          optimize(ex<Product>(prod), OptimizeOptions{.footprint_weight = 0.0});
      REQUIRE(res0 == res0_explicit);  // weight 0 is a no-op

      // ... a large footprint_weight changes the chosen factorization.
      auto resF = optimize(ex<Product>(prod),
                           OptimizeOptions{.footprint_weight = 100.});
      REQUIRE(res0 != resF);

      uocc->approximate_size(uocc_sz);
      aux->approximate_size(aux_sz);
    }

    SECTION("Ensure single-value sums/products are not discarded") {
      auto sum = ex<Sum>();
      sum->as<Sum>().append(
          ex<Product>(ExprPtrList{deserialize(L"f{a_1;i_1}")}));
      REQUIRE(sum->as<Sum>().summand(0).as<Product>().factors().size() == 1);
      auto optimized = optimize(sum);
      REQUIRE(optimized->is<Sum>());
      REQUIRE(optimized->as<Sum>().summands().size() == 1);
      REQUIRE(sum->as<Sum>().summand(0).as<Product>().factors().size() == 1);
    }

    SECTION("Non-covariant indices") {
      auto uocc = reg->retrieve_ptr(L"a");
      auto aux = reg->retrieve_ptr(L"x");
      auto const aux_sz = aux->approximate_size();
      aux->approximate_size(3 * uocc->approximate_size());

      auto const G_abcd_thc =
          deserialize(L"X{a1;;x1} X{;a2;x1} Y{;;x1,x2} X{a3;;x2} X{;a4;x2}")
              ->as<Product>();
      auto const G_abcd_thc_opt =
          deserialize(
              L"((X{a1;;x1} X{;a2;x1}) Y{;;x1,x2})(X{a3;;x2} X{;a4;x2})")
              ->as<Product>();
      REQUIRE(single_term_opt(G_abcd_thc)->as<Product>() == G_abcd_thc_opt);

      auto const GT_abij_thc = deserialize(
                                   L"X{a1;;x1} X{;a2;x1} Y{;;x1,x2} X{a3;;x2} "
                                   L"X{;a4;x2} T{a2,a4;i1,i2}")
                                   ->as<Product>();
      auto const GT_abij_thc_opt = deserialize(
                                       L"(((X{a1;;x1} X{;a2;x1}) Y{;;x1,x2}) ( "
                                       L"X{;a4;x2} T{a2,a4;i1,i2} )) X{a3;;x2}")
                                       ->as<Product>();
      REQUIRE(single_term_opt(GT_abij_thc)->as<Product>() == GT_abij_thc_opt);

      aux->approximate_size(aux_sz);
    }

    SECTION("OptimizeOptions: cost metric and reorder knobs") {
      auto const prod = parse_expr_antisymm(
          L"g_{i3,i4}^{a3,a4} t_{a1,a2}^{i3,i4} t_{a3,a4}^{i1,i2}");

      // both metrics must binarize the 3-tensor product into a binary tree:
      // a 2-factor top product whose leaves are the 3 original tensors
      for (auto objective_function :
           {ObjectiveFunction::DenseFLOPs, ObjectiveFunction::DenseSize}) {
        CAPTURE(static_cast<int>(objective_function));
        auto res = optimize(
            prod, OptimizeOptions{.objective_function = objective_function});
        REQUIRE(res->is<Product>());
        REQUIRE(res->as<Product>().factors().size() == 2);
        REQUIRE(count_tensor_leaves(res) == 3);
      }

      // reorder knob: a two-summand sum is optimized either way, and the
      // optimize() default (reorder) matches an explicit Reorder request
      auto const sum = parse_expr_antisymm(
          L"g_{i3,i4}^{a3,a4} t_{a1,a2}^{i3,i4} t_{a3,a4}^{i1,i2}"
          L" + g_{i3,i4}^{a3,a4} t_{a3,a4}^{i1,i2} t_{a1}^{i3} t_{a2}^{i4}");
      REQUIRE(sum->is<Sum>());

      auto no_reorder =
          optimize(sum, OptimizeOptions{.reorder = ReorderSum::NoReorder});
      auto reorder =
          optimize(sum, OptimizeOptions{.reorder = ReorderSum::Reorder});
      REQUIRE(no_reorder->is<Sum>());
      REQUIRE(reorder->is<Sum>());
      REQUIRE(no_reorder->as<Sum>().size() == sum->as<Sum>().size());
      REQUIRE(reorder->as<Sum>().size() == sum->as<Sum>().size());
      // default options == explicit Reorder
      REQUIRE(*optimize(sum) == *reorder);
    }

    SECTION("Parallel optimization of summands matches sequential") {
      // exercise optimize_impl(..., parallel_outer=true): a multi-summand sum
      // optimized concurrently must yield the same result as single-threaded.
      auto const sum = parse_expr_antisymm(
          L"g_{i3,i4}^{a3,a4} t_{a1,a2}^{i3,i4} t_{a3,a4}^{i1,i2}"
          L" + g_{i3,i4}^{a3,a4} t_{a3,a4}^{i1,i2} t_{a1}^{i3} t_{a2}^{i4}"
          L" + g_{i3,i4}^{a3,a4} t_{a1}^{i3} t_{a2}^{i4} t_{a3,a4}^{i1,i2}");
      REQUIRE(sum->is<Sum>());
      REQUIRE(sum->as<Sum>().size() > 1);

      auto const nthreads_save = num_threads();
      struct ThreadGuard {
        int n;
        ~ThreadGuard() { set_num_threads(n); }
      } guard{nthreads_save};

      set_num_threads(1);
      auto const seq = optimize(sum);
      set_num_threads(4);
      auto const par = optimize(sum);

      REQUIRE(*seq == *par);
    }

    SECTION("subset_footprints") {
      using namespace sequant;
      // i occ (size 2); a virt (size 4). Tensors: g{a1;i1}, g{a2;i2}.
      auto occ_sp = reg->retrieve_ptr(L"i");
      auto virt_sp = reg->retrieve_ptr(L"a");
      auto const occ_sz_save = occ_sp->approximate_size();
      auto const virt_sz_save = virt_sp->approximate_size();
      occ_sp->approximate_size(2);
      virt_sp->approximate_size(4);

      auto idxsz = [](Index const& ix) -> std::size_t {
        return ix.space().approximate_size();
      };
      auto g1 = deserialize(L"g{a1;i1}", {.def_perm_symm = Symmetry::Nonsymm});
      auto g2 = deserialize(L"g{a2;i2}", {.def_perm_symm = Symmetry::Nonsymm});
      TensorNetwork net{std::vector<ExprPtr>{g1, g2}};
      container::svector<Index> targets;  // empty: all indices remain open
      auto S = opt::detail::subset_footprints(net, targets, idxsz);
      REQUIRE(S.size() == 4u);
      REQUIRE(S[0] == 0.0);  // empty subset
      // singleton {T0}: open indices a1,i1 -> 4*2 = 8
      REQUIRE(S[0b01] == 8.0);
      REQUIRE(S[0b10] == 8.0);
      // full {T0,T1}: open a1,i1,a2,i2 -> 4*2*4*2 = 64
      REQUIRE(S[0b11] == 64.0);

      occ_sp->approximate_size(occ_sz_save);
      virt_sp->approximate_size(virt_sz_save);
    }

    SECTION("DensePeakSize DP matches brute-force oracle") {
      using namespace sequant;
      auto idxsz = [](Index const& ix) -> std::size_t {
        return ix.space().approximate_size();
      };
      // A 4-tensor chain whose intermediates differ in size by contraction
      // order.
      auto t0 = deserialize(L"g{a1;i1}", {.def_perm_symm = Symmetry::Nonsymm});
      auto t1 = deserialize(L"g{a1;a2}", {.def_perm_symm = Symmetry::Nonsymm});
      auto t2 = deserialize(L"g{a2;a3}", {.def_perm_symm = Symmetry::Nonsymm});
      auto t3 = deserialize(L"g{a3;i2}", {.def_perm_symm = Symmetry::Nonsymm});
      TensorNetwork net{std::vector<ExprPtr>{t0, t1, t2, t3}};
      container::svector<Index> targets;  // i1,i2 left open
      auto S = opt::detail::subset_footprints(net, targets, idxsz);
      double oracle =
          brute_force_min_peak(std::vector<double>(S.begin(), S.end()), 4);
      double dp = opt::detail::peak_cost(net, targets, idxsz);
      REQUIRE(dp == oracle);
    }

    SECTION("DensePeakSize matches oracle over a battery of small networks") {
      using namespace sequant;
      auto idxsz = [](Index const& ix) -> std::size_t {
        return ix.space().approximate_size();
      };
      std::vector<std::vector<std::wstring>> nets = {
          {L"g{a1;i1}", L"g{a1;a2}", L"g{a2;i2}"},
          {L"g{a1;i1}", L"g{a1;a2}", L"g{a2;a3}", L"g{a3;i2}"},
          {L"g{i1;a1}", L"t{a1,a2;i1,i2}", L"g{a2;i2}"},
          {L"g{a1,a2;i1,i2}", L"t{a1;i1}", L"t{a2;i2}"},
      };
      for (auto const& spec : nets) {
        std::vector<ExprPtr> ts;
        for (auto const& s : spec)
          ts.push_back(deserialize(s, {.def_perm_symm = Symmetry::Nonsymm}));
        TensorNetwork net{ts};
        container::svector<Index> targets;
        auto S = opt::detail::subset_footprints(net, targets, idxsz);
        double oracle = brute_force_min_peak(
            std::vector<double>(S.begin(), S.end()), ts.size());
        double dp = opt::detail::peak_cost(net, targets, idxsz);
        REQUIRE(dp == oracle);
      }
    }

    SECTION("DensePeakSize reconstructed sequence achieves the DP optimum") {
      using namespace sequant;
      auto idxsz = [](Index const& ix) -> std::size_t {
        return ix.space().approximate_size();
      };
      // Simulate the all-co-resident (model A) peak of a reconstructed
      // EvalSequence and confirm it EQUALS peak_cost (the DP's minimum).
      // This proves the emitted contraction order actually realizes the
      // optimum, not just that the DP computed a number.  CRITICAL: all
      // input leaves are resident from the start (model A) -- a naive stack
      // machine that pushes a leaf only when its token appears would compute
      // the Sethi-Ullman (model B) peak and disagree.  A tensor is freed
      // when consumed.
      auto peak_of_sequence = [](EvalSequence const& seq,
                                 container::vector<double> const& S,
                                 size_t nt) -> double {
        container::set<size_t> live;
        for (size_t b = 0; b < nt; ++b) live.insert(size_t{1} << b);
        container::vector<size_t> stack;
        double peak = 0.0;
        for (int tok : seq) {
          if (tok >= 0) {
            stack.push_back(size_t{1} << tok);
          } else {
            size_t rhs = stack.back();
            stack.pop_back();
            size_t lhs = stack.back();
            stack.pop_back();
            size_t merged = lhs | rhs;
            double live_sum = 0.0;
            for (auto m : live) live_sum += S[m];
            peak = std::max(peak, live_sum + S[merged]);
            live.erase(lhs);
            live.erase(rhs);
            live.insert(merged);
            stack.push_back(merged);
          }
        }
        return peak;
      };
      std::vector<std::vector<std::wstring>> nets = {
          {L"g{a1;i1}", L"g{a1;a2}", L"g{a2;i2}"},
          {L"g{a1;i1}", L"g{a1;a2}", L"g{a2;a3}", L"g{a3;i2}"},
          {L"g{a1,a2;i1,i2}", L"t{a1;i1}", L"t{a2;i2}"},
      };
      for (auto const& spec : nets) {
        std::vector<ExprPtr> ts;
        for (auto const& s : spec)
          ts.push_back(deserialize(s, {.def_perm_symm = Symmetry::Nonsymm}));
        TensorNetwork net{ts};
        container::svector<Index> targets;
        auto S = opt::detail::subset_footprints(net, targets, idxsz);
        auto seq = opt::detail::run_single_term_opt(
            opt::detail::PeakModel{idxsz}, net, targets);
        double dp = opt::detail::peak_cost(net, targets, idxsz);
        REQUIRE(peak_of_sequence(seq, S, ts.size()) == dp);
      }
    }

    SECTION("per-index batchability tables") {
      using namespace sequant;
      // Add a dedicated aux/fitting space "F" to the cloned registry so that
      // g{a1;i1;F1} and g{a2;i1;F2} can be deserialized.  Use a fresh type bit
      // (0b10000) that does not overlap with the sr-spaces bits
      // (0b0001..0b1000).
      reg->add(L"F", IndexSpace::Type{0b10000}, 3ul);

      auto idxsz = [](Index const& ix) -> std::size_t {
        return ix.space().approximate_size();
      };
      auto is_batchable = [](Index const& ix) {
        return ix.space().base_key() == L"F";
      };
      auto t0 =
          deserialize(L"g{a1;i1;F1}", {.def_perm_symm = Symmetry::Nonsymm});
      auto t1 =
          deserialize(L"g{a2;i1;F2}", {.def_perm_symm = Symmetry::Nonsymm});
      TensorNetwork net{std::vector<ExprPtr>{t0, t1}};
      container::svector<Index> targets;
      auto aux = opt::detail::batchable_index_list(net, is_batchable);
      REQUIRE(aux.size() == 2u);  // F1, F2 distinct
      auto batch_fn = [](Index const&) -> std::size_t { return 1; };
      auto tables = opt::detail::sliced_footprints(net, targets, idxsz,
                                                   is_batchable, batch_fn, aux);
      REQUIRE(tables.size() == 4u);  // 2^2 sliced-sets
      // B=00 (none sliced) is the full footprint; B=11 (both) the all-sliced.
      REQUIRE(tables[0b00][0b11] > tables[0b11][0b11]);  // full > all-sliced
      // slicing only F1 (bit 0) shrinks the F1-leaf but not the F2-leaf.
      size_t f1bit = 0;  // aux[0]==F1 by appearance order
      REQUIRE(tables[size_t{1} << f1bit][0b01] < tables[0b00][0b01]);
    }

    SECTION(
        "DensePeakSizeBatched all-sliced corner equals Phase-1 batched "
        "peak") {
      using namespace sequant;
      // Dedicated batchable "F" space (fresh type bit, no overlap with
      // sr-spaces bits 0b0001..0b1000); see "per-index batchability tables".
      reg->add(L"F", IndexSpace::Type{0b10000}, 3ul);
      auto idxsz = [](Index const& ix) -> std::size_t {
        return ix.space().approximate_size();
      };
      auto is_batchable = [](Index const& ix) {
        return ix.space().base_key() == L"F";
      };
      std::size_t const batch = 1;
      std::vector<ExprPtr> ts;
      for (auto s : {L"g{a1;i1;F1}", L"g{a2;i1;F1}", L"g{a2;i2;F2}"})
        ts.push_back(deserialize(s, {.def_perm_symm = Symmetry::Nonsymm}));
      TensorNetwork net{ts};
      container::svector<Index> targets;
      auto aux = opt::detail::batchable_index_list(net, is_batchable);
      std::size_t const m = aux.size();
      auto batch_fn = [batch](Index const&) -> std::size_t { return batch; };
      opt::detail::PeakBatchedModel model{idxsz, is_batchable, batch_fn,
                                          /*is_volatile_leaf=*/{}};
      auto ctx = model.build_context(net, targets);
      auto st = opt::detail::solve_single_term(model, net, targets, ctx);
      size_t root = (size_t{1} << ts.size()) - 1;
      size_t allK = (size_t{1} << m) - 1;
      double dp_allsliced = std::numeric_limits<double>::max();
      for (auto const& fp : st[root][allK])
        dp_allsliced = std::min(dp_allsliced, fp.peak);
      // Phase-1 peak with EVERY batchable index sliced: an extent wrapper that
      // slices iff the index is batchable (no batched_extent helper exists).
      auto be = [&](Index const& ix) -> std::size_t {
        std::size_t e = idxsz(ix);
        return is_batchable(ix) ? std::min(e, batch) : e;
      };
      double phase1 = opt::detail::peak_cost(net, targets, be);
      REQUIRE(dp_allsliced == phase1);
    }

    SECTION("DensePeakSizeBatched objective matches per-index oracle") {
      using namespace sequant;
      reg->add(L"F", IndexSpace::Type{0b10000}, 3ul);
      auto idxsz = [](Index const& ix) -> std::size_t {
        return ix.space().approximate_size();
      };
      auto is_batchable = [](Index const& ix) {
        return ix.space().base_key() == L"F";
      };
      std::size_t const batch = 1;
      std::vector<std::vector<std::wstring>> nets = {
          {L"g{a1;i1;F1}", L"g{a2;i1;F1}", L"g{a2;i2;F2}"},  // shared F1
          {L"g{a1;i1;F1}", L"g{a2;i1;F2}", L"g{a2;i2;F2}"},  // two distinct aux
      };
      for (auto const& spec : nets) {
        std::vector<ExprPtr> ts;
        for (auto const& s : spec)
          ts.push_back(deserialize(s, {.def_perm_symm = Symmetry::Nonsymm}));
        TensorNetwork net{ts};
        container::svector<Index> targets;
        auto aux = opt::detail::batchable_index_list(net, is_batchable);
        auto vmask = opt::detail::leaf_volatile_mask(net, {});
        auto batch_fn = [batch](Index const&) -> std::size_t { return batch; };
        auto tables = opt::detail::sliced_footprints(
            net, targets, idxsz, is_batchable, batch_fn, aux);
        // open_aux[s] via the SAME detail helper the DP uses, so DP and oracle
        // index `tables` identically.
        auto open_aux_det = opt::detail::subset_open_aux(net, targets, aux);
        std::vector<size_t> open_aux(open_aux_det.begin(), open_aux_det.end());
        std::vector<std::vector<double>> T(tables.begin(), tables.end());
        std::vector<char> persistent(size_t{1} << ts.size());
        for (size_t n = 0; n < persistent.size(); ++n)
          persistent[n] = ((vmask & n) == 0) ? 1 : 0;
        double oracle = batched_min_peak(T, open_aux, persistent, ts.size());
        double dp = opt::detail::peak_cost_batched(net, targets, idxsz,
                                                   is_batchable, batch_fn, {});
        REQUIRE(dp == oracle);
      }
    }

    SECTION(
        "DensePeakSizeBatched reconstruction achieves the optimum "
        "(numeric)") {
      using namespace sequant;
      // Dedicated batchable "F" space (fresh type bit, no overlap with
      // sr-spaces bits 0b0001..0b1000); see "per-index batchability tables".
      reg->add(L"F", IndexSpace::Type{0b10000}, 3ul);
      auto idxsz = [](Index const& ix) -> std::size_t {
        return ix.space().approximate_size();
      };
      auto is_batchable = [](Index const& ix) {
        return ix.space().base_key() == L"F";
      };
      std::size_t const batch = 1;
      for (auto const& spec : std::vector<std::vector<std::wstring>>{
               {L"g{a1;i1;F1}", L"g{a2;i1;F1}", L"g{a2;i2;F2}"},  // shared F1
               {L"g{a1;i1;F1}", L"g{a2;i1;F2}", L"g{a2;i2;F2}"}}) {  // 2 aux
        std::vector<ExprPtr> ts;
        for (auto const& s : spec)
          ts.push_back(deserialize(s, {.def_perm_symm = Symmetry::Nonsymm}));
        TensorNetwork net{ts};
        container::svector<Index> targets;
        // recompute the chosen tree's peak by simulation over the pr
        // back-pointers (independent of the DP's max/+ recurrence):
        auto batch_fn = [batch](Index const&) -> std::size_t { return batch; };
        double recon = opt::detail::reconstructed_batched_peak(
            net, targets, idxsz, is_batchable, batch_fn, {});
        double dp = opt::detail::peak_cost_batched(net, targets, idxsz,
                                                   is_batchable, batch_fn, {});
        REQUIRE(recon == dp);
      }
    }

    SECTION("optimize() public API dispatches DensePeakSizeBatched") {
      using namespace sequant;
      // Drives Step 4 (the opt_pure_product runtime dispatch for the batched
      // arm). Before Step 4 this hits the exhaustiveness SEQUANT_ASSERT in
      // opt_pure_product and fails; after Step 4 it returns a binarized
      // product over the 3 leaves.
      reg->add(L"F", IndexSpace::Type{0b10000}, 3ul);
      auto expr = deserialize(L"g{a1;i1;F1} * g{a2;i1;F1} * g{a2;i2;F2}",
                              {.def_perm_symm = Symmetry::Nonsymm});
      OptimizeOptions opts;
      opts.objective_function = ObjectiveFunction::DensePeakSizeBatched;
      opts.idx_to_extent = [](Index const& ix) -> std::size_t {
        return ix.space().approximate_size();
      };
      opts.batch_policy.is_batchable_index = [](Index const& ix) {
        return ix.space().base_key() == L"F";
      };
      opts.batch_policy.batch_target_size = [](Index const&) -> std::size_t {
        return 1;
      };
      auto optimized = optimize(expr, opts);
      REQUIRE(optimized);
      REQUIRE(count_tensor_leaves(optimized) == 3u);
    }

    SECTION("optimize() public API dispatches DensePeakSize") {
      using namespace sequant;
      // Drives Step 3 (the opt_pure_product runtime dispatch). Before Step 3
      // this hits the SEQUANT_ASSERT(objective_function == DenseSize) in
      // opt_pure_product and fails; after Step 3 it returns a binarized
      // product.
      auto expr = deserialize(L"g{a1;i1} * g{a1;a2} * g{a2;i2}");
      OptimizeOptions opts;
      opts.objective_function = ObjectiveFunction::DensePeakSize;
      opts.idx_to_extent = [](Index const& ix) -> std::size_t {
        return ix.space().approximate_size();
      };
      auto optimized = optimize(expr, opts);
      REQUIRE(optimized);
      // optimized is a binarized product over the same 3 leaves.
      REQUIRE(count_tensor_leaves(optimized) == 3u);
    }

    SECTION("per-index batch_target_size honored") {
      using namespace sequant;
      // Register F with approximate_size=3; two distinct aux indices F_1, F_2.
      // Index::label() returns base_key + "_" + ordinal, so the F1 ordinal-1
      // index has label L"F_1" (not L"F1").
      reg->add(L"F", IndexSpace::Type{0b10000}, 3ul);
      auto idxsz = [](Index const& ix) -> std::size_t {
        return ix.space().approximate_size();
      };
      auto is_batchable = [](Index const& ix) {
        return ix.space().base_key() == L"F";
      };
      // 3-tensor network where the connector tensor t1 carries BOTH F_1 and
      // F_2: F_1 is contracted between t0 and t1, while F_2 propagates
      // through t1 and is contracted later between t1 and t2.  Because t1
      // carries both aux indices simultaneously, its footprint under a given
      // context depends on the batch sizes of F_1 AND F_2 independently,
      // making c_mixed strictly between c_all1 and c_all2.
      //
      // Network: t0=g{a1;i1;F1}, t1=g{a2;i1;F1,F2}, t2=g{a2;i2;F2}
      //   contracted: i1 (t0-t1), F1 (t0-t1), a2 (t1-t2), F2 (t1-t2)
      //   free (result): a1, i2   -> result size = 100*10 = 1000
      std::vector<ExprPtr> ts;
      for (auto s : {L"g{a1;i1;F1}", L"g{a2;i1;F1,F2}", L"g{a2;i2;F2}"})
        ts.push_back(deserialize(s, {.def_perm_symm = Symmetry::Nonsymm}));
      TensorNetwork net{ts};
      container::svector<Index> targets;
      // Uniform batch sizes for baseline costs.
      auto all1 = [](Index const&) -> std::size_t { return 1; };
      auto all2 = [](Index const&) -> std::size_t { return 2; };
      // mixed: F_1 gets batch=1, F_2 gets batch=2.  The discriminator uses
      // the exact label returned by Index::label(), which is L"F_1" for
      // ordinal-1 (not L"F1").  Because t1 carries both F_1 and F_2, any
      // scalar implementation (same value for all indices) cannot reproduce
      // c_mixed; both per-index values participate in the optimum.
      auto mixed = [](Index const& ix) -> std::size_t {
        return ix.label() == L"F_1" ? std::size_t{1} : std::size_t{2};
      };
      double c_all1 = opt::detail::peak_cost_batched(net, targets, idxsz,
                                                     is_batchable, all1, {});
      double c_all2 = opt::detail::peak_cost_batched(net, targets, idxsz,
                                                     is_batchable, all2, {});
      double c_mixed = opt::detail::peak_cost_batched(net, targets, idxsz,
                                                      is_batchable, mixed, {});
      // c_mixed must differ from both baselines: a scalar batch_target_size
      // (returning the same value for F_1 and F_2) cannot produce c_mixed.
      REQUIRE(c_all1 != c_all2);   // network is sensitive to batch size
      REQUIRE(c_mixed != c_all1);  // mixed is not the same as all-batch-1
      REQUIRE(c_mixed != c_all2);  // mixed is not the same as all-batch-2
    }

    SECTION("CostModel concept conformance + custom model") {
      using namespace sequant;
      // idxsz lambda captures approximate_size() for each index space.
      auto idxsz = [](Index const& ix) {
        return ix.space().approximate_size();
      };

      // --- Static conformance checks ---
      // AdditiveModel (FLOPs variant): two template params (CostFn,
      // FootprintFn).
      static_assert(opt::detail::CostModel<opt::detail::AdditiveModel<
                        decltype(opt::detail::flops_counter(idxsz)),
                        decltype(opt::detail::footprint_counter(idxsz))>>);
      // AdditiveModel (Size variant).
      static_assert(opt::detail::CostModel<opt::detail::AdditiveModel<
                        decltype(opt::detail::memsize_counter(idxsz)),
                        decltype(opt::detail::footprint_counter(idxsz))>>);
      // PeakModel: one template param (IdxToSz).
      static_assert(
          opt::detail::CostModel<opt::detail::PeakModel<decltype(idxsz)>>);
      // PeakBatchedModel: one template param (IdxToSz).
      static_assert(opt::detail::CostModel<
                    opt::detail::PeakBatchedModel<decltype(idxsz)>>);

      // Negative direction: a type lacking State/Context/the six methods must
      // NOT satisfy CostModel (guards against a vacuously-true concept).
      struct NotAModel {};
      static_assert(!opt::detail::CostModel<NotAModel>);

      // --- Custom-model extension-point exercise ---
      // Build an AdditiveModel driven by memsize_counter with a doubled
      // footprint weight, then drive it directly via run_single_term_opt
      // (bypassing the ObjectiveFunction enum).  This proves the public
      // generic entry point is open to user-defined or custom-configured
      // models.
      opt::detail::AdditiveModel custom{opt::detail::memsize_counter(idxsz),
                                        opt::detail::footprint_counter(idxsz),
                                        /*volatile_mask=*/0u,
                                        /*volatile_weight=*/1.0,
                                        /*footprint_weight=*/2.0,
                                        /*subnet_cse=*/false};

      std::vector<ExprPtr> ts;
      for (auto s : {L"g{a1;i1}", L"g{a1;a2}", L"g{a2;i2}"})
        ts.push_back(deserialize(s, {.def_perm_symm = Symmetry::Nonsymm}));
      TensorNetwork net{ts};
      container::svector<Index> targets;
      auto seq = opt::detail::run_single_term_opt(custom, net, targets);

      // A valid binarization of 3 leaves: exactly 3 leaf tokens (>= 0) and
      // 2 merge tokens (-1).
      size_t leaves = 0, merges = 0;
      for (int tok : seq) {
        if (tok >= 0)
          ++leaves;
        else
          ++merges;
      }
      REQUIRE(leaves == 3u);
      REQUIRE(merges == 2u);
    }
  }

  SECTION("CSE") {
    auto ctx_resetter =
        set_scoped_default_context(get_default_context().clone());
    IndexSpaceRegistry registry;
    registry.add("a", 0b001);
    registry.add("i", 0b010);
    registry.add("u", 0b100);
    *get_default_context().mutable_index_space_registry() = registry;

    for (bool force_hash_collisions : {false, true}) {
      CAPTURE(force_hash_collisions);

      for (const auto& [inputs, outputs] : std::vector<
               std::pair<std::vector<std::wstring>, std::vector<std::wstring>>>{
               // Basic example with only scalars
               {{L"R1 = (A B) C", L"R2 = D (A B)"},
                {L"CSE1 = A B", L"R1 = CSE1 C", L"R2 = D CSE1"}},
               // Test case in which the same intermediate is reused but
               // requires different indexing
               {{L"R{a1,a3;i2,i3} = 2 GAM0{a1,a3;a4,a5} T2g{a4,a5;i2,i3} - "
                 L"GAM0{a1,a3;a4,a5} T2g{a4,a5;i3,i2}"},
                {L"CSE1{;;a3,a1,i2,i3} = GAM0{a1,a3;a4,a5} T2g{a4,a5;i2,i3}",
                 L"R{a1,a3;i2,i3} = 2 CSE1{;;a3,a1,i2,i3} - "
                 L"CSE1{;;a3,a1,i3,i2}"}},
               // Scalar CSE with proto-index tensors
               {{L"R1 = (f{i1;a1<i1>} t{a1<i1>;i1}) A",
                 L"R2 = B (f{i1;a1<i1>} t{a1<i1>;i1})"},
                {L"CSE1 = f{i1;a1<i1>} t{a1<i1>;i1}", L"R1 = CSE1 A",
                 L"R2 = B CSE1"}},
               // Tensor CSE with proto-index tensors: reused
               // contraction with different external indexing
               {{L"R{i1;i2} = 2 g{i1;a1<i1>} t{a1<i1>;i2} - "
                 L"g{i2;a1<i2>} t{a1<i2>;i1}"},
                {L"CSE1{;;i1,i2} = g{i1;a1<i1>} t{a1<i1>;i2}",
                 L"R{i1;i2} = 2 CSE1{;;i1,i2} - CSE1{;;i2,i1}"}},
               // ToT CSE: the intermediate itself has proto-indexed
               // indices (tensor-of-tensor)
               {{L"R1{i1;i2} = (g{i1;a1} C{a1;a1<i1>}) h{a1<i1>;i2}",
                 L"R2{i1;i2} = (g{i1;a1} C{a1;a1<i1>}) k{a1<i1>;i2}"},
                {L"CSE1{;;a1<i1>,i1} = g{i1;a1} C{a1;a1<i1>}",
                 L"R1{i1;i2} = CSE1{;;a1<i1>,i1} h{a1<i1>;i2}",
                 L"R2{i1;i2} = CSE1{;;a1<i1>,i1} k{a1<i1>;i2}"}},
               // In this case it is important that the computation of the
               // subexpression isn't simply thrown at the beginning of the
               // expression list as it depends on B, which has to be computed
               // first.
               {{L"B = K J", L"R = (A B) C + (A B) D"},
                {L"B = K J", L"CSE1 = A B", L"R = CSE1 C + CSE1 D"}},
               // CSE in the presence of bra-ket symmetry
               {{L"R2{u2,a1;u1,i1} = -2 f{u3;u4}:N-S Y{u2,u3;u1,u5} "
                 L"t{a1,u5;i1,u4} + f{u3;u4}:N-S Y{u2,u4;u5,u1} t{a1,u5;i1,u3} "
                 L"+ f{u3;u4}:N-S Y{u2,u4;u1,u5} t{a1,u5;u3,i1}"},
                {L"CSE1{;;u4,u2,u1,u5} = f{u3;u4}:N-S Y{u2,u3;u1,u5}",
                 L"R2{u2,a1;u1,i1} = -2 CSE1{;;u4,u2,u1,u5} t{a1,u5;i1,u4}"
                 L" + CSE1{;;u3,u2,u5,u1} t{a1,u5;i1,u3}"
                 L" + CSE1{;;u3,u2,u1,u5} t{a1,u5;u3,i1}"}},
           }) {
        CAPTURE(inputs);

        std::vector<EvalNode<EvalExpr>> expressions;
        std::vector<ResultExpr> expected;

        for (const std::wstring& current : inputs) {
          expressions.push_back(binarize(deserialize<ResultExpr>(
              current, {.def_perm_symm = Symmetry::Nonsymm,
                        .def_braket_symm = BraKetSymmetry::Nonsymm,
                        .def_col_symm = ColumnSymmetry::Nonsymm})));
        }
        for (const std::wstring& current : outputs) {
          expected.push_back(deserialize<ResultExpr>(
              current, {.def_perm_symm = Symmetry::Nonsymm,
                        .def_braket_symm = BraKetSymmetry::Nonsymm,
                        .def_col_symm = ColumnSymmetry::Nonsymm}));
        }

        auto binarizer = [](auto&& expr) {
          // CSE drives binarize() on subexpressions for hash-equivalence
          // detection; positional head is irrelevant here.
          SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
          return binarize(expr);
          SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
        };

        if (force_hash_collisions) {
          // This code path makes all hashes be computed to be zero and hence
          // every pair of objects will yield a hash collision which need to be
          // dealt with by using proper comparison operators.
          static constexpr bool force_collisions = true;
          opt::eliminate_common_subexpressions<
              decltype(expressions), decltype(binarizer),
              opt::cse::AcceptAllPredicate, force_collisions>(expressions,
                                                              binarizer);
        } else {
          opt::eliminate_common_subexpressions(expressions, binarizer);
        }

        std::vector<ResultExpr> actual;
        for (const auto& current : expressions) {
          if (current->is_tensor()) {
            actual.emplace_back(current->expr()->as<Tensor>(),
                                to_expr(current));
          } else {
            REQUIRE(current->is_scalar());
            actual.emplace_back(current->expr()->as<Variable>(),
                                to_expr(current));
          }
        }

        REQUIRE(actual == expected);
      }
    }
  }

  SECTION("Single term optimization with CSE") {
    auto ctx_resetter =
        set_scoped_default_context(get_default_context().clone());
    auto reg = get_default_context().mutable_index_space_registry();
    mbpt::add_df_spaces(reg);
    mbpt::add_pao_spaces(reg);
    mbpt::add_ao_spaces(reg);
    // i 10
    // a 40
    // μ̃ 50
    // Κ 90
    for (auto&& [k, v] :
         std::initializer_list<std::pair<std::wstring_view, size_t>>{
             {L"i", 10}, {L"a", 40}, {L"μ̃", 50}, {L"Κ", 90}}) {
      reg->retrieve_ptr(k)->approximate_size(v);
    }

    auto single_term_opt = [](Product const& prod, bool cse = true) {
      return opt::single_term_opt(
          prod,
          [](Index const& ix) {
            // null space contributes x1 to the size
            auto sz = ix.nonnull() ? ix.space().approximate_size() : 1;
            return sz;
          },
          /*subnet_cse=*/cse);
    };

    auto prod9 =
        deserialize("X{i1;a1} X{i2;a2} Y{a2;i3} Y{a1;i4}")->as<Product>();
    auto res9 = single_term_opt(prod9);
    auto res9_no_cse = single_term_opt(prod9, false);
    // this is the one we want to find
    // (X Y) (X Y)
    REQUIRE(extract(res9, {0, 0}) == prod9.at(0));
    REQUIRE(extract(res9, {0, 1}) == prod9.at(3));
    REQUIRE(extract(res9, {1, 0}) == prod9.at(1));
    REQUIRE(extract(res9, {1, 1}) == prod9.at(2));

    // take a look at res9_no_cse for a result with subnet_cse disabled
    // should give the same result in this case as it's already optimal
    REQUIRE(extract(res9_no_cse, {0, 0}) == prod9.at(0));
    REQUIRE(extract(res9_no_cse, {0, 1}) == prod9.at(3));
    REQUIRE(extract(res9_no_cse, {1, 0}) == prod9.at(1));
    REQUIRE(extract(res9_no_cse, {1, 1}) == prod9.at(2));

    SECTION("CSE effect on optimization result") {
      auto ctx_resetter =
          set_scoped_default_context(get_default_context().clone());
      auto reg = get_default_context().mutable_index_space_registry();
      // Use sizes that make the unbalanced tree better without CSE,
      // but the balanced tree better with CSE.
      // Balanced: ( (X1 Y1) (X2 Y2) )
      // Cost(X1*Y1) = size(i)*size(a)*size(j) = 12*10*12 = 1440.
      // Cost(Inter) = 12^3 = 1728.
      // Total no-CSE: 2*1440 + 1728 = 4608.
      // Total CSE: 1440 + 1728 = 3168.
      // Unbalanced: ( ( (X1 Y1) X2 ) Y2 )
      // Cost(X1*Y1) = 12*10*12 = 1440.
      // Cost((X1*Y1)*X2) = size(i)*size(i)*size(a) = 12*12*10 = 1440.
      // Cost(...) * Y2 = 12*10*12 = 1440.
      // Total Unbalanced: 1440 + 1440 + 1440 = 4320.
      // 3168 < 4320 < 4608.
      reg->retrieve_ptr(L"i")->approximate_size(12);
      reg->retrieve_ptr(L"a")->approximate_size(10);

      auto single_term_opt = [](Product const& prod, bool cse) {
        return opt::single_term_opt(
            prod,
            [](Index const& ix) {
              return ix.nonnull() ? ix.space().approximate_size() : 1;
            },
            cse);
      };

      // X{i1;a1} Y{a1;i2} X{i2;a2} Y{a2;i3}
      auto prod =
          deserialize(L"X{i1;a1} Y{a1;i2} X{i2;a2} Y{a2;i3}")->as<Product>();

      auto res_cse = single_term_opt(prod, true);
      auto res_no_cse = single_term_opt(prod, false);

      // With CSE: Balanced tree
      REQUIRE(res_cse->as<Product>().factors().size() == 2);
      REQUIRE(res_cse->at(0)->is<Product>());
      REQUIRE(res_cse->at(1)->is<Product>());

      // Without CSE: Unbalanced tree
      bool is_unbalanced =
          (res_no_cse->at(0)->is<Tensor>() || res_no_cse->at(1)->is<Tensor>());
      REQUIRE(is_unbalanced);
    }

    SECTION("subnet_cse flows through OptimizeOptions") {
      auto ctx_resetter =
          set_scoped_default_context(get_default_context().clone());
      auto reg = get_default_context().mutable_index_space_registry();
      // Same sizing trick as the section above: CSE prefers balanced,
      // no-CSE prefers unbalanced.
      reg->retrieve_ptr(L"i")->approximate_size(12);
      reg->retrieve_ptr(L"a")->approximate_size(10);

      auto idx_to_extent = [](Index const& ix) -> std::size_t {
        return ix.nonnull() ? ix.space().approximate_size() : 1;
      };

      auto prod =
          deserialize(L"X{i1;a1} Y{a1;i2} X{i2;a2} Y{a2;i3}")->as<Product>();
      auto expr = ex<Product>(prod);

      auto res_cse =
          optimize(expr, OptimizeOptions{.CSE = {.subnet = true},
                                         .idx_to_extent = idx_to_extent});
      auto res_no_cse =
          optimize(expr, OptimizeOptions{.CSE = {.subnet = false},
                                         .idx_to_extent = idx_to_extent});

      // With CSE: balanced tree -- both children are Products.
      REQUIRE(res_cse->is<Product>());
      REQUIRE(res_cse->as<Product>().factors().size() == 2);
      REQUIRE(res_cse->at(0)->is<Product>());
      REQUIRE(res_cse->at(1)->is<Product>());

      // Without CSE: unbalanced tree -- at least one child is a bare Tensor.
      REQUIRE(res_no_cse->is<Product>());
      REQUIRE(res_no_cse->as<Product>().factors().size() == 2);
      bool is_unbalanced =
          res_no_cse->at(0)->is<Tensor>() || res_no_cse->at(1)->is<Tensor>();
      REQUIRE(is_unbalanced);

      // Default OptimizeOptions => subnet_cse Disable => same as no-CSE shape.
      auto res_default =
          optimize(expr, OptimizeOptions{.idx_to_extent = idx_to_extent});
      REQUIRE(res_default->is<Product>());
      REQUIRE(res_default->as<Product>().factors().size() == 2);
      bool default_is_unbalanced =
          res_default->at(0)->is<Tensor>() || res_default->at(1)->is<Tensor>();
      REQUIRE(default_is_unbalanced);
    }
  }

  /// verify that space changes did not leak
  auto reg_check = get_default_context().index_space_registry();
  auto uocc_check = reg_check->retrieve_ptr(L"a");
  REQUIRE(uocc_check);
  REQUIRE(uocc_check->approximate_size() == 10);
}

// ---------------------------------------------------------------------------
// Reproducer: an OSV (proto-indexed) contraction C{a<i>;μ̃} * t{a<i>;i} should
// be done EARLY because it eliminates the OSV index a<i> (rank/size/flops
// drop), producing I{i;μ̃}. Motif distilled from PNO-CCSD residual intermediate
// #1.
// ---------------------------------------------------------------------------
namespace {
std::wstring render_tree(sequant::ExprPtr const& e) {
  using namespace sequant;
  if (e->is<Tensor>()) {
    auto const& t = e->as<Tensor>();
    std::wstring s = std::wstring(t.label()) + L"{";
    bool first = true;
    for (auto const& ix : t.bra()) {
      if (!first) s += L",";
      s += ix.full_label();
      first = false;
    }
    s += L";";
    first = true;
    for (auto const& ix : t.ket()) {
      if (!first) s += L",";
      s += ix.full_label();
      first = false;
    }
    if (t.aux().size()) {
      s += L";";
      first = true;
      for (auto const& ix : t.aux()) {
        if (!first) s += L",";
        s += ix.full_label();
        first = false;
      }
    }
    return s + L"}";
  }
  if (e->is<Product>()) {
    std::wstring s = L"(";
    bool first = true;
    for (auto const& f : e->as<Product>().factors()) {
      if (!first) s += L" * ";
      s += render_tree(f);
      first = false;
    }
    return s + L")";
  }
  return L"?";
}
}  // namespace

TEST_CASE("OSV early-contraction reproducer", "[optimize][osv]") {
  using namespace sequant;
  auto ctx_resetter = set_scoped_default_context(get_default_context().clone());
  auto reg = get_default_context().mutable_index_space_registry();
  mbpt::add_df_spaces(reg);
  mbpt::add_pao_spaces(reg);
  mbpt::add_ao_spaces(reg);
  for (auto&& [k, v] :
       std::initializer_list<std::pair<std::wstring_view, size_t>>{
           {L"i", 10}, {L"a", 40}, {L"μ̃", 50}, {L"Κ", 90}}) {
    reg->retrieve_ptr(k)->approximate_size(v);
  }
  auto idxsz = [](Index const& ix) -> std::size_t {
    return ix.nonnull() ? ix.space().approximate_size() : std::size_t{1};
  };

  // motif from intermediate #1: g(μ̃μ̃Κ) · C(a<i>;μ̃) · C(μ̃;b<i,j>) · t(a<i>;i)
  auto prod =
      deserialize(L"g{μ̃1;μ̃2;Κ1} C{a1<i1>;μ̃1} C{μ̃2;a2<i1,i2>} t{a1<i1>;i1}")
          ->as<Product>();
  std::wcout << L"\n=== INPUT: " << render_tree(ex<Product>(prod)) << L"\n";

  auto show = [&](auto metric, std::wstring name,
                  std::function<double(Index const&, std::size_t)> ip = {}) {
    auto res = opt::single_term_opt<decltype(metric)::value>(
        prod, idxsz, /*subnet_cse=*/false, CostParams{{}, 1.0, 0.0}, {}, {},
        ip);
    std::wcout << name << L":  " << render_tree(res) << L"\n";
  };
  std::wcout
      << L"--- without inner_pow (composite a sized as full uocc=40) ---\n";
  show(std::integral_constant<ObjectiveFunction,
                              ObjectiveFunction::DenseFLOPs>{},
       L"FLOPs");
  show(
      std::integral_constant<ObjectiveFunction, ObjectiveFunction::DenseSize>{},
      L"Size ");
  show(std::integral_constant<ObjectiveFunction,
                              ObjectiveFunction::DensePeakSize>{},
       L"Peak ");

  auto ip = [](Index const&, std::size_t) -> double { return 12.0; };
  std::wcout << L"--- with inner_pow (composite a<i> sized small=12, like a "
                L"PNO/OSV domain) ---\n";
  show(std::integral_constant<ObjectiveFunction,
                              ObjectiveFunction::DenseFLOPs>{},
       L"FLOPs", ip);
  show(
      std::integral_constant<ObjectiveFunction, ObjectiveFunction::DenseSize>{},
      L"Size ", ip);
  show(std::integral_constant<ObjectiveFunction,
                              ObjectiveFunction::DensePeakSize>{},
       L"Peak ", ip);
  std::wcout << L"\n";
}

TEST_CASE("OSV early-contraction reproducer (full term #1)",
          "[optimize][osv]") {
  using namespace sequant;
  auto ctx_resetter = set_scoped_default_context(get_default_context().clone());
  auto reg = get_default_context().mutable_index_space_registry();
  mbpt::add_df_spaces(reg);
  mbpt::add_pao_spaces(reg);
  mbpt::add_ao_spaces(reg);
  for (auto&& [k, v] :
       std::initializer_list<std::pair<std::wstring_view, size_t>>{
           {L"i", 56}, {L"a", 12}, {L"μ̃", 602}, {L"Κ", 1652}}) {
    reg->retrieve_ptr(k)->approximate_size(v);
  }
  auto aux_space = reg->retrieve(L"Κ");
  auto idxsz = [](Index const& ix) -> std::size_t {
    return ix.nonnull() ? ix.space().approximate_size() : std::size_t{1};
  };
  // composite OSV/PNO domain ~ small (like MPQC's ~12-PNO/pair)
  auto ip = [](Index const&, std::size_t) -> double { return 12.0; };
  auto is_batch = [aux_space](Index const& ix) {
    return ix.space() == aux_space;
  };
  std::function<std::size_t(Index const&)> bts = [](Index const&) {
    return std::size_t{236};
  };

  // verbatim term #1 from the water-14 trace
  auto prod = deserialize(
                  L"g{μ̃_1;i_3;Κ_1} C{a_3<i_2>;μ̃_1} g{μ̃_2;μ̃_3;Κ_1} "
                  L"C{a_4<i_1>;μ̃_2} C{μ̃_3;a_1<i_1,i_2>} t{a_4<i_1>;i_1} "
                  L"t{a_5<i_3>;i_3} t{a_3<i_2>;i_2} C{μ̃_5;a_5<i_3>} "
                  L"s{μ̃_4;μ̃_5} C{a_2<i_1,i_2>;μ̃_4}")
                  ->as<Product>();
  std::wcout << L"\n=== FULL TERM #1 (" << prod.factors().size()
             << L" tensors) ===\n";

  {
    auto res = opt::single_term_opt<ObjectiveFunction::DenseFLOPs>(
        prod, idxsz, false, CostParams{{}, 1.0, 0.0}, {}, {}, ip);
    std::wcout << L"FLOPs:        " << render_tree(res) << L"\n";
  }
  {
    auto res = opt::single_term_opt<ObjectiveFunction::DensePeakSize>(
        prod, idxsz, false, CostParams{{}, 1.0, 0.0}, {}, {}, ip);
    std::wcout << L"PeakSize:     " << render_tree(res) << L"\n";
  }
  {
    auto res = opt::single_term_opt<ObjectiveFunction::DensePeakSizeBatched>(
        prod, idxsz, false, CostParams{{}, 1.0, 0.0}, is_batch, bts, ip);
    std::wcout << L"PeakBatched:  " << render_tree(res) << L"\n";
  }
  std::wcout << L"\n  >>> looking for (C{a_4<i_1>;μ̃_1216} * t{a_4<i_1>;i_1}) "
                L"contracted EARLY <<<\n\n";
}

TEST_CASE("OSV deferral reproducer (tetramer term 3)", "[optimize][osv]") {
  using namespace sequant;
  auto ctx_resetter = set_scoped_default_context(get_default_context().clone());
  auto reg = get_default_context().mutable_index_space_registry();
  mbpt::add_df_spaces(reg);
  mbpt::add_pao_spaces(reg);
  mbpt::add_ao_spaces(reg);
  for (auto&& [k, v] :
       std::initializer_list<std::pair<std::wstring_view, size_t>>{
           {L"i", 16}, {L"a", 12}, {L"μ̃", 170}, {L"Κ", 472}}) {
    reg->retrieve_ptr(k)->approximate_size(v);
  }
  auto aux_space = reg->retrieve(L"Κ");
  auto idxsz = [](Index const& ix) -> std::size_t {
    return ix.nonnull() ? ix.space().approximate_size() : std::size_t{1};
  };
  auto ip = [](Index const&, std::size_t) -> double { return 12.0; };
  auto is_batch = [aux_space](Index const& ix) {
    return ix.space() == aux_space;
  };
  std::function<std::size_t(Index const&)> bts = [](Index const&) {
    return std::size_t{236};
  };

  // tetramer trace term 3 (verbatim, indices renumbered to small ordinals)
  auto prod =
      deserialize(
          L"C{a_3<i_1>;μ̃_1} C{μ̃_2;a_1<i_1>} g{μ̃_1;μ̃_2;Κ_1} g{μ̃_3;i_2;Κ_1} "
          L"C{a_2<i_2>;μ̃_3} t{a_2<i_2>;i_2} t{a_3<i_1>;i_1}")
          ->as<Product>();
  std::wcout << L"\n=== TETRAMER TERM 3 (" << prod.factors().size()
             << L" tensors) ===\n";
  std::wcout << L"  trace chose: ...(g{μ̃μ̃Κ} * (...)) * t{a_3<i_1>;i_1}  "
                L"[a_3<i_1> DEFERRED]\n";
  auto fp = opt::detail::footprint_counter(
      idxsz, std::function<double(Index const&, std::size_t)>(ip));
  auto fp_flops_c = opt::detail::flops_counter(
      idxsz, std::function<double(Index const&, std::size_t)>(ip));
  // walk the ExprPtr product tree (no canonicalization); return a subexpr's
  // free indices and update mx with the largest intermediate (result) size
  // seen.
  std::function<std::vector<Index>(ExprPtr const&, double&)> freeix =
      [&](ExprPtr const& e, double& mx) -> std::vector<Index> {
    if (e->is<Tensor>()) {
      std::vector<Index> v;
      for (auto const& ix : e->as<Tensor>().const_braketaux_indices())
        v.push_back(ix);
      return v;
    }
    std::map<std::wstring, std::pair<int, Index>> cnt;
    for (auto const& fct : e->as<Product>().factors()) {
      for (auto const& ix : freeix(fct, mx)) {
        auto k = std::wstring(ix.full_label());
        cnt[k].first++;
        cnt[k].second = ix;
      }
    }
    std::vector<Index>
        result;  // appears in exactly one child => not contracted here
    for (auto const& [k, v] : cnt)
      if (v.first == 1) result.push_back(v.second);
    double here = fp(result);
    if (here > mx) mx = here;
    return result;
  };
  // max LEAF size: the raw DF integral g{μ̃,μ̃,Κ} is a huge leaf that is resident
  // during its own contraction, so it lower-bounds the peak of EVERY schedule.
  double max_leaf = 0.0;
  std::wstring big;
  for (auto const& fct : prod.factors()) {
    std::vector<Index> v;
    for (auto const& ix : fct->as<Tensor>().const_braketaux_indices())
      v.push_back(ix);
    double s = fp(v);
    if (s > max_leaf) {
      max_leaf = s;
      big = render_tree(fct);
    }
  }
  std::wcout << L"max LEAF = " << (long long)max_leaf << L"  (" << big
             << L")\n";
  for (auto const& fct : prod.factors()) {
    std::vector<Index> v;
    for (auto const& ix : fct->as<Tensor>().const_braketaux_indices())
      v.push_back(ix);
    std::wcout << L"   leaf " << render_tree(fct) << L" = " << (long long)fp(v)
               << L"\n";
  }

  // Replicate PeakModel::relax EXACTLY on a fixed binary tree.
  struct TP {
    double peak, leafsum, S;
  };
  std::function<TP(ExprPtr const&)> tree_peak = [&](ExprPtr const& e) -> TP {
    if (e->is<Tensor>()) {
      std::vector<Index> v;
      for (auto const& ix : e->as<Tensor>().const_braketaux_indices())
        v.push_back(ix);
      double s = fp(v);
      return TP{s, s, s};
    }
    auto const& facs = e->as<Product>().factors();
    TP L = tree_peak(facs[0]);
    TP R = tree_peak(facs[1]);
    double dummy = 0.0;
    auto resix = freeix(e, dummy);
    double Snode = fp(resix);
    double both = L.S + R.S + Snode;
    double lfirst = std::max({R.leafsum + L.peak, L.S + R.peak, both});
    double rfirst = std::max({L.leafsum + R.peak, R.S + L.peak, both});
    return TP{std::min(lfirst, rfirst), L.leafsum + R.leafsum, Snode};
  };

  // weighted flops of a fixed binary tree: w=vw when subtree contains a 't'.
  std::function<bool(ExprPtr const&)> has_t = [&](ExprPtr const& e) -> bool {
    if (e->is<Tensor>()) return e->as<Tensor>().label() == L"t";
    for (auto const& fct : e->as<Product>().factors())
      if (has_t(fct)) return true;
    return false;
  };
  auto idxof = [&](ExprPtr const& e) {
    double d = 0.0;
    return freeix(e, d);
  };
  std::function<double(ExprPtr const&, double)> tree_flops =
      [&](ExprPtr const& e, double vw) -> double {
    if (e->is<Tensor>()) return 0.0;
    auto const& facs = e->as<Product>().factors();
    auto a = idxof(facs[0]);
    auto b = idxof(facs[1]);
    auto r = idxof(e);
    double here = fp_flops_c(a, b, r);
    double w = has_t(e) ? vw : 1.0;
    return w * here + tree_flops(facs[0], vw) + tree_flops(facs[1], vw);
  };

  auto report = [&](std::wstring name, sequant::ExprPtr res) -> double {
    double mx = 0.0;
    freeix(res, mx);
    double tpk = tree_peak(res).peak;
    std::wcout << name << L"  recurrence_PEAK=" << (long long)tpk
               << L"  max_imed=" << (long long)mx << L"  wflops(vw1)="
               << (long long)tree_flops(res, 1.0) << L"  wflops(vw100)="
               << (long long)tree_flops(res, 100.0) << L"\n      "
               << render_tree(res) << L"\n";
    return mx;
  };
  double flops_mx =
      report(L"FLOPs:        ",
             opt::single_term_opt<ObjectiveFunction::DenseFLOPs>(
                 prod, idxsz, false, CostParams{{}, 1.0, 0.0}, {}, {}, ip));
  double peak_mx =
      report(L"PeakSize:     ",
             opt::single_term_opt<ObjectiveFunction::DensePeakSize>(
                 prod, idxsz, false, CostParams{{}, 1.0, 0.0}, {}, {}, ip));
  double pbat_mx = report(
      L"PeakBatched:  ",
      opt::single_term_opt<ObjectiveFunction::DensePeakSizeBatched>(
          prod, idxsz, false, CostParams{{}, 1.0, 0.0}, is_batch, bts, ip));
  // Real MPQC config: t is volatile, volatile_weight=100. The tie-break weights
  // volatile (replayed) flops, so it must STILL eliminate the OSV early.
  auto is_t = [](Tensor const& t) { return t.label() == L"t"; };
  // isolate volatile_weight: hold is_volatile_leaf=is_t FIXED, vary only vw.
  double peak_v =
      report(L"PeakSize/is_t,vw1:   ",
             opt::single_term_opt<ObjectiveFunction::DensePeakSize>(
                 prod, idxsz, false, CostParams{is_t, 1.0, 0.0}, {}, {}, ip));
  double peak_v100 =
      report(L"PeakSize/is_t,vw100: ",
             opt::single_term_opt<ObjectiveFunction::DensePeakSize>(
                 prod, idxsz, false, CostParams{is_t, 100.0, 0.0}, {}, {}, ip));
  double pbat_v = report(
      L"PeakBatch/is_t,vw1:  ",
      opt::single_term_opt<ObjectiveFunction::DensePeakSizeBatched>(
          prod, idxsz, false, CostParams{is_t, 1.0, 0.0}, is_batch, bts, ip));
  double pbat_v100 = report(
      L"PeakBatch/is_t,vw100:",
      opt::single_term_opt<ObjectiveFunction::DensePeakSizeBatched>(
          prod, idxsz, false, CostParams{is_t, 100.0, 0.0}, is_batch, bts, ip));
  (void)peak_v100;
  (void)pbat_v100;
  std::wcout << L"\n";
  // The defect this guards against is the OSV-deferred *outer product*: folding
  // the volatile t-amplitude in before contracting the shared subtree forces a
  // ~5.5M-element intermediate (`osv_outer_product`, below). The peak
  // objectives must avoid it; the persistent-only gate must reproduce it.
  // Asserting on that gross structural threshold -- not on exact max_imed
  // equality -- keeps the test robust to the (peak, flops) Pareto frontier's
  // epsilon-tolerant selection, under which a vw-weighted flop reduction may
  // legitimately pick a schedule with a larger single intermediate but a
  // within-tolerance DP peak.
  double const osv_outer_product = 5.0e6;  // actual deferred imed is ~5.5M
  // (0) The order-independent FLOPs objective never forms the OSV outer product
  //     (it is a flop blow-up, not just a memory one) -- the reference
  //     baseline.
  CHECK(flops_mx < osv_outer_product);
  // (1) Default (across-the-board batching, epsilon-tolerant Pareto): every
  // peak
  //     objective avoids the OSV outer product, at every volatile_weight, with
  //     or without a volatility gate (is_t).
  CHECK(peak_mx < osv_outer_product);    // DensePeakSize, no volatile predicate
  CHECK(pbat_mx < osv_outer_product);    // DensePeakSizeBatched, no predicate
  CHECK(peak_v < osv_outer_product);     // DensePeakSize, is_t, vw=1
  CHECK(peak_v100 < osv_outer_product);  // DensePeakSize, is_t, vw=100
  // (2) Batching is applied across the board (not gated on persistence), so the
  //     BATCHED model with a volatility gate (is_t) ALSO avoids the OSV outer
  //     product: slicing reduces footprint regardless of volatility.
  CHECK(pbat_v < osv_outer_product);
  CHECK(pbat_v100 < osv_outer_product);
  // (3) The persistent-only gate is still available as an opt-in: setting
  //     batch_persistent_only restores the old behavior (volatile subtrees not
  //     sliced -> the batched model reverts to deferring the OSV outer
  //     product).
  double pbat_po = report(
      L"PeakBatch/persistent_only:",
      opt::single_term_opt<ObjectiveFunction::DensePeakSizeBatched>(
          prod, idxsz, false, CostParams{is_t, 100.0, 0.0}, is_batch, bts, ip,
          /*batch_persistent_only=*/true));
  CHECK(pbat_po >= osv_outer_product);  // gate restored -> OSV deferred again
}

TEST_CASE("PPL: form 4-PNO W vs fold-t (peak-neutral, flop tie-break)",
          "[optimize][osv]") {
  using namespace sequant;
  auto ctx_resetter = set_scoped_default_context(get_default_context().clone());
  auto reg = get_default_context().mutable_index_space_registry();
  mbpt::add_df_spaces(reg);
  mbpt::add_pao_spaces(reg);
  mbpt::add_ao_spaces(reg);
  for (auto&& [k, v] :
       std::initializer_list<std::pair<std::wstring_view, size_t>>{
           {L"i", 16}, {L"a", 12}, {L"μ̃", 170}, {L"Κ", 472}})
    reg->retrieve_ptr(k)->approximate_size(v);
  auto aux = reg->retrieve(L"Κ");
  auto idxsz = [](Index const& ix) -> std::size_t {
    return ix.nonnull() ? ix.space().approximate_size() : std::size_t{1};
  };
  auto ip = [](Index const&, std::size_t) -> double { return 12.0; };
  auto is_batch = [aux](Index const& ix) { return ix.space() == aux; };
  std::function<std::size_t(Index const&)> bts = [](Index const&) {
    return std::size_t{236};
  };
  auto is_t = [](Tensor const& t) { return t.label() == L"t"; };

  // PPL: R_ij^{a1 a2} = (a1 a3 | a2 a4) t_ij^{a3 a4}, DF: (a1 a3|K)=gCC, (a2
  // a4|K)=gCC
  auto prod = deserialize(
                  L"C{a_1<i_1,i_2>;μ̃_1} g{μ̃_1;μ̃_2;Κ_1} C{μ̃_2;a_3<i_1,i_2>} "
                  L"C{a_2<i_1,i_2>;μ̃_3} g{μ̃_3;μ̃_4;Κ_1} C{μ̃_4;a_4<i_1,i_2>} "
                  L"t{a_3<i_1,i_2>,a_4<i_1,i_2>;i_1,i_2}")
                  ->as<Product>();

  auto fp = opt::detail::footprint_counter(
      idxsz, std::function<double(Index const&, std::size_t)>(ip));
  auto fc = opt::detail::flops_counter(
      idxsz, std::function<double(Index const&, std::size_t)>(ip));
  std::function<std::vector<Index>(ExprPtr const&, double&)> freeix =
      [&](ExprPtr const& e, double& mx) -> std::vector<Index> {
    if (e->is<Tensor>()) {
      std::vector<Index> v;
      for (auto const& ix : e->as<Tensor>().const_braketaux_indices())
        v.push_back(ix);
      return v;
    }
    std::map<std::wstring, std::pair<int, Index>> c;
    for (auto const& f : e->as<Product>().factors())
      for (auto const& ix : freeix(f, mx)) {
        auto k = std::wstring(ix.full_label());
        c[k].first++;
        c[k].second = ix;
      }
    std::vector<Index> r;
    for (auto const& [k, v] : c)
      if (v.first == 1) r.push_back(v.second);
    double here = fp(r);
    if (here > mx) mx = here;
    return r;
  };
  std::function<bool(ExprPtr const&)> has_t = [&](ExprPtr const& e) -> bool {
    if (e->is<Tensor>()) return e->as<Tensor>().label() == L"t";
    for (auto const& f : e->as<Product>().factors())
      if (has_t(f)) return true;
    return false;
  };
  std::function<double(ExprPtr const&, double)> wflops =
      [&](ExprPtr const& e, double vw) -> double {
    if (e->is<Tensor>()) return 0.0;
    auto const& f = e->as<Product>().factors();
    double d = 0;
    auto a = freeix(f[0], d), b = freeix(f[1], d), r = freeix(e, d);
    double w = has_t(e) ? vw : 1.0;
    return w * fc(a, b, r) + wflops(f[0], vw) + wflops(f[1], vw);
  };
  auto rep = [&](std::wstring name, ExprPtr res) -> double {
    double mx = 0;
    freeix(res, mx);
    double wf = wflops(res, 100.0);
    std::wcout << name << L"  max_imed=" << (long long)mx << L"  wflops(vw100)="
               << (long long)wf << L"\n      " << render_tree(res) << L"\n";
    return wf;
  };
  std::wcout << L"\n=== PPL term: form-W vs fold-t ===\n";
  // The flop-optimal schedule forms the persistent 4-PNO integral W=(ac|bd)=gCC
  // *gCC once, then contracts the volatile amplitude t into it -- the volatile
  // (replayed) flops are the cheap step. The "fold-t" alternative folds t into
  // a gCC half-transform first; that recomputes the ladder on every replay, a
  // ~7x larger volatile-weighted flop count. The 4-PNO W carries a free PAO leg
  // so it is slightly LARGER in peak than fold-t; only the peak_flops_tolerance
  // (default 0.10) lets the peak objectives accept that within-tolerance peak
  // bump in exchange for the large volatile-flop win, i.e. form W like FLOPs.
  double flops_w =
      rep(L"FLOPs/vw100:  ",
          opt::single_term_opt<ObjectiveFunction::DenseFLOPs>(
              prod, idxsz, false, CostParams{is_t, 100.0, 0.0}, {}, {}, ip));
  double peak_w =
      rep(L"PeakSize/vw100:",
          opt::single_term_opt<ObjectiveFunction::DensePeakSize>(
              prod, idxsz, false, CostParams{is_t, 100.0, 0.0}, {}, {}, ip));
  double pbat_w = rep(
      L"PeakB/vw100:  ",
      opt::single_term_opt<ObjectiveFunction::DensePeakSizeBatched>(
          prod, idxsz, false, CostParams{is_t, 100.0, 0.0}, is_batch, bts, ip));
  std::wcout << L"\n";
  // The epsilon-tolerant Pareto selection makes both peak objectives form W,
  // matching the flop-optimal volatile-weighted flop count.
  CHECK(peak_w == flops_w);
  CHECK(pbat_w == flops_w);
  // Strict peak-min (peak_flops_tolerance=0) instead defers to fold-t: it
  // refuses the W peak bump, paying the much larger volatile-flop ladder.
  double peak_w_strict =
      rep(L"PeakSize/strict:",
          opt::single_term_opt<ObjectiveFunction::DensePeakSize>(
              prod, idxsz, false,
              CostParams{is_t, 100.0, 0.0, /*peak_flops_tolerance=*/0.0}, {},
              {}, ip));
  CHECK(peak_w_strict > flops_w);
  // At volatile_weight=1 (the caching-off regime, where persistent
  // intermediates cannot be amortized across replays) the persistent 4-PNO W is
  // strictly dominated by fold-t on BOTH peak and flops: the batched root
  // frontier collapses to the single fold-t point, so every objective folds the
  // amplitude t into a Kappa-batched half-transform ladder instead of
  // forming/caching W. (mpqc's SeQuantEngine pins volatile_weight to 1 when
  // eval:cache is off for exactly this reason.) A form-W tree would reproduce
  // flops_w; fold-t does not.
  double pbat_w_vw1 = rep(
      L"PeakB/vw1:    ",
      opt::single_term_opt<ObjectiveFunction::DensePeakSizeBatched>(
          prod, idxsz, false, CostParams{is_t, 1.0, 0.0}, is_batch, bts, ip));
  CHECK(pbat_w_vw1 > flops_w);
}

// Quadratic-bubble (g·t2·t2) exchange term in PNO/CSV basis: the two competing
// factorizations of one residual contribution.
//
//   early-K: contract the shared DF aux Κ between the two dressed integrals
//            g{i_4;a_3<i_1,i_3>;Κ}·g{i_3;a_4<i_2,i_4>;Κ} FIRST, forming the
//            held-whole 4-occ/2-PNO integral I{i_1..i_4; a_3,a_4} (Κ-free,
//            peak ≈ o⁴·p²), then bring in the amplitudes.
//   late-K : build each half M_x = t·(gC) FIRST (Κ retained, sliced to K_b:
//            M_x{i_1..i_4,Κ; a_1<i_1,i_2>}, peak ≈ o⁴·p·K_b), contract Κ LAST;
//            the two halves co-reside at that final node (peak ≈ 2·o⁴·p·K_b).
//
// Crossover (pure peak): early-K wins iff K_b > p/2. With accumulation_factor λ
// charging the held-whole accumulator's co-resident batch contribution, early-K
// is priced (1+λ)·o⁴·p², so late-K wins iff K_b < (1+λ)·p/2 — i.e. raising λ
// favors the (batchable, memory-bounded) late-K route.
TEST_CASE("quadratic bubble: early-K integral vs late-K t·(gC)",
          "[optimize][bubble]") {
  using namespace sequant;
  auto ctx_resetter = set_scoped_default_context(get_default_context().clone());
  auto reg = get_default_context().mutable_index_space_registry();
  mbpt::add_df_spaces(reg);
  mbpt::add_pao_spaces(reg);
  mbpt::add_ao_spaces(reg);
  // water-20-scale extents (≈ water-14 OSV extents scaled by 20/14).
  for (auto&& [k, v] :
       std::initializer_list<std::pair<std::wstring_view, size_t>>{
           {L"i", 80}, {L"a", 12}, {L"μ̃", 860}, {L"Κ", 2360}}) {
    reg->retrieve_ptr(k)->approximate_size(v);
  }
  auto aux_space = reg->retrieve(L"Κ");
  auto idxsz = [](Index const& ix) -> std::size_t {
    return ix.nonnull() ? ix.space().approximate_size() : std::size_t{1};
  };
  // PNO domain per pair (composite inner extent). Crossover is K_b = p/2 = 6.
  double const p = 12.0;
  auto ip = [p](Index const&, std::size_t) -> double { return p; };
  auto is_batch = [aux_space](Index const& ix) {
    return ix.space() == aux_space;
  };

  // Full EXCHANGE quadratic bubble, parenthesized per half then flattened to
  // 12 leaves. S(PNO-PNO overlap) is written CsC (s = PAO-PAO overlap), the
  // form MPQC exposes. Externals: i_1,i_2,a_1<i_1,i_2>,a_2<i_1,i_2>; the halves
  // share/contract Κ_1, i_3, i_4.
  auto nested = deserialize(
      L"( g{i_4;μ̃_1;Κ_1} * C{μ̃_1;a_3<i_1,i_3>}"
      L"  * t{a_1<i_1,i_3>,a_3<i_1,i_3>;i_1,i_3}"
      L"  * C{a_1<i_1,i_3>;μ̃_2} * s{μ̃_2;μ̃_3} * C{μ̃_3;a_1<i_1,i_2>} )"
      L"* ( g{i_3;μ̃_4;Κ_1} * C{μ̃_4;a_4<i_2,i_4>}"
      L"  * t{a_2<i_2,i_4>,a_4<i_2,i_4>;i_2,i_4}"
      L"  * C{a_2<i_2,i_4>;μ̃_5} * s{μ̃_5;μ̃_6} * C{μ̃_6;a_2<i_1,i_2>} )",
      {.def_perm_symm = Symmetry::Nonsymm});
  Product flatp{};
  for (auto const& half : nested->as<Product>().factors())
    flatp.append(1, half, Product::Flatten::Yes);
  REQUIRE(flatp.factors().size() == 12);
  Product const& prod = flatp;

  // Subtree predicates (robust to scalar/constant factors).
  std::function<bool(ExprPtr const&)> has_g = [&](ExprPtr const& e) -> bool {
    if (e->is<Tensor>()) return e->as<Tensor>().label() == L"g";
    if (e->is<Product>())
      for (auto const& f : e->as<Product>().factors())
        if (has_g(f)) return true;
    return false;
  };
  std::function<bool(ExprPtr const&)> has_t = [&](ExprPtr const& e) -> bool {
    if (e->is<Tensor>()) return e->as<Tensor>().label() == L"t";
    if (e->is<Product>())
      for (auto const& f : e->as<Product>().factors())
        if (has_t(f)) return true;
    return false;
  };
  // True iff the chosen tree builds the held-whole (g·C)(g·C) integral: some
  // node joins the two DF g's (each in a different factor) with NO amplitude t
  // anywhere below it. In late-K each g is fused with its t first, so the g's
  // never share a t-free subtree.
  std::function<bool(ExprPtr const&)> forms_integral =
      [&](ExprPtr const& e) -> bool {
    if (!e->is<Product>()) return false;
    auto const& facs = e->as<Product>().factors();
    int gcount = 0;
    for (auto const& f : facs)
      if (has_g(f)) ++gcount;
    if (gcount == 2 && !has_t(e)) return true;
    for (auto const& f : facs)
      if (forms_integral(f)) return true;
    return false;
  };

  // inner_pow-aware, Kappa-sliced peak of a chosen binary tree (mirrors the
  // DP's co-resident recurrence). Composites sized by ip; batchable Kappa
  // capped at K_b. Returns element count; *8 bytes for memory.
  auto batched_peak = [&](ExprPtr const& root, std::size_t Kb) -> double {
    auto ext = [&, Kb](Index const& ix) -> std::size_t {
      std::size_t e = idxsz(ix);
      return is_batch(ix) ? std::min(e, Kb) : e;
    };
    auto fp = opt::detail::footprint_counter(
        ext, std::function<double(Index const&, std::size_t)>(ip));
    // result (free) indices of a subtree: those appearing in exactly one child.
    std::function<std::vector<Index>(ExprPtr const&)> freeix =
        [&](ExprPtr const& e) -> std::vector<Index> {
      if (e->is<Tensor>()) {
        std::vector<Index> v;
        for (auto const& ix : e->as<Tensor>().const_braketaux_indices())
          v.push_back(ix);
        return v;
      }
      if (!e->is<Product>()) return {};
      std::map<std::wstring, std::pair<int, Index>> cnt;
      for (auto const& fct : e->as<Product>().factors())
        for (auto const& ix : freeix(fct)) {
          auto k = std::wstring(ix.full_label());
          cnt[k].first++;
          cnt[k].second = ix;
        }
      std::vector<Index> result;
      for (auto const& [k, v] : cnt)
        if (v.first == 1) result.push_back(v.second);
      return result;
    };
    struct TP {
      double peak, leafsum, S;
    };
    std::function<TP(ExprPtr const&)> tp = [&](ExprPtr const& e) -> TP {
      if (e->is<Tensor>()) {
        double s = fp(freeix(e));
        return TP{s, s, s};
      }
      std::vector<ExprPtr> kids;
      for (auto const& f : e->as<Product>().factors())
        if (f->is<Tensor>() || f->is<Product>()) kids.push_back(f);
      if (kids.size() == 1) return tp(kids[0]);
      TP L = tp(kids[0]), R = tp(kids[1]);
      double Snode = fp(freeix(e));
      double both = L.S + R.S + Snode;
      double lf = std::max({R.leafsum + L.peak, L.S + R.peak, both});
      double rf = std::max({L.leafsum + R.peak, R.S + L.peak, both});
      return TP{std::min(lf, rf), L.leafsum + R.leafsum, Snode};
    };
    return tp(root).peak;
  };

  auto choose = [&](std::size_t Kb, double lambda) -> bool {
    std::function<std::size_t(Index const&)> bts = [Kb](Index const&) {
      return Kb;
    };
    // CostParams: {is_volatile_leaf, volatile_weight, footprint_weight,
    //              peak_flops_tolerance, roofline, accumulation_factor}.
    // peak_flops_tolerance = 0 => strict peak-min (no tolerance band).
    auto res = opt::single_term_opt<ObjectiveFunction::DensePeakSizeBatched>(
        prod, idxsz, false, CostParams{{}, 1.0, 0.0, 0.0, {}, lambda}, is_batch,
        bts, ip);
    bool integral = forms_integral(res);
    double gb = batched_peak(res, Kb) * 8.0 / 1e9;
    std::wcout << L"  K_b=" << Kb << L"\tlambda=" << lambda
               << (integral ? L"\tEARLY-K (integral)" : L"\tLATE-K  (t.(gC))")
               << L"\tpeak=" << gb << L" GB\n";
    return integral;
  };

  std::wcout << L"\n=== quadratic bubble exchange: i=80 mu=860 K=2360 p=" << p
             << L" (crossover K_b=p/2=" << p / 2 << L") ===\n";
  for (double lam : {0.0, 1.0, 10.0, 40.0}) {
    std::wcout << L"-- lambda=" << lam << L" --\n";
    for (std::size_t Kb : {std::size_t{2}, std::size_t{6}, std::size_t{12},
                           std::size_t{72}, std::size_t{236}})
      choose(Kb, lam);
  }
  // Real MPQC config: t is volatile (replayed), default peak tolerance 0.10
  // lets the (replay-weighted) flop tie-break trade peak for building the
  // persistent, t-free, Kappa-free integral ONCE. This is the suspected real
  // driver of the held-whole 4-occ/2-PNO object (the C60 OOM), independent of
  // accumulation.
  auto is_t = [](Tensor const& t) { return t.label() == L"t"; };
  std::wcout
      << L"-- volatile t, vw=100, tolerance=0.10 (real config), lambda=0 "
         L"--\n";
  for (std::size_t Kb :
       {std::size_t{2}, std::size_t{12}, std::size_t{72}, std::size_t{236}}) {
    std::function<std::size_t(Index const&)> bts = [Kb](Index const&) {
      return Kb;
    };
    auto res = opt::single_term_opt<ObjectiveFunction::DensePeakSizeBatched>(
        prod, idxsz, false, CostParams{is_t, 100.0, 0.0, 0.10, {}, 0.0},
        is_batch, bts, ip);
    // peak of the actually-chosen tree, plus the rejected alternative's peak
    // (force pure-peak to recover the late-K tree) to see the traded premium.
    double gb = batched_peak(res, Kb) * 8.0 / 1e9;
    auto alt = opt::single_term_opt<ObjectiveFunction::DensePeakSizeBatched>(
        prod, idxsz, false, CostParams{{}, 1.0, 0.0, 0.0, {}, 0.0}, is_batch,
        bts, ip);
    double gb_alt = batched_peak(alt, Kb) * 8.0 / 1e9;
    std::wcout << L"  K_b=" << Kb
               << (forms_integral(res) ? L"\tEARLY-K (integral)"
                                       : L"\tLATE-K  (t.(gC))")
               << L"\tpeak=" << gb << L" GB   (pure-peak alt="
               << (forms_integral(alt) ? L"early" : L"late") << L" " << gb_alt
               << L" GB)\n";
  }
  std::wcout
      << L"  chosen tree (K_b=236, lambda=0):\n      "
      << render_tree(
             opt::single_term_opt<ObjectiveFunction::DensePeakSizeBatched>(
                 prod, idxsz, false, CostParams{{}, 1.0, 0.0, 0.0, {}, 0.0},
                 is_batch, [](Index const&) { return std::size_t{236}; }, ip))
      << L"\n";

  // Robust sanity: with Kappa sliced to a tiny batch the held-whole integral
  // cannot win (o^4*p^2 > 2*o^4*p*K_b for K_b < p/2), so the optimizer must
  // take the batchable late-K route.
  CHECK_FALSE(choose(/*K_b=*/2, /*lambda=*/0.0));
}
