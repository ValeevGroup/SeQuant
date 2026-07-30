#include <catch2/catch_approx.hpp>
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
#include <cstdlib>
#include <functional>
#include <initializer_list>
#include <limits>
#include <memory>
#include <unordered_map>
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
// indices, subset = bitmask of leaves). open_modes[subset] = bitmask of aux
// indices that are open (external) in that subset. persistent[subset] = 1 iff
// the contraction result at that subset is persistent. nt = number of leaves.
// B = bitmask of aux indices currently being sliced (threading the context).
static double oracle_rec(size_t n, size_t B,
                         std::vector<std::vector<double>> const& T,
                         std::vector<size_t> const& open_modes,
                         std::vector<char> const& persistent, size_t nt) {
  if (std::popcount(n) == 1) return T[B & open_modes[n]][n];
  auto sz = [&](size_t s, size_t ctx) { return T[ctx & open_modes[s]][s]; };
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
    size_t contracted_here =
        persistent[n] ? ((open_modes[lp] | open_modes[rp]) & ~open_modes[n])
                      : 0;
    for (size_t Ap = contracted_here;; Ap = (Ap - 1) & contracted_here) {
      size_t C = B | Ap;
      double pl = oracle_rec(lp, C, T, open_modes, persistent, nt);
      double pr = oracle_rec(rp, C, T, open_modes, persistent, nt);
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
                               std::vector<size_t> const& open_modes,
                               std::vector<char> const& persistent, size_t nt) {
  size_t full = (size_t{1} << nt) - 1;
  return oracle_rec(full, 0, T, open_modes, persistent, nt);
}

TEST_CASE("batched_min_peak oracle (per-index)", "[optimize]") {
  // 2 leaves sharing aux F (m=1). open_modes: leaves have F open (bit0=1);
  // the pair has F internal (bit0=0). persistent everywhere.
  // T[ctx][subset]: full (ctx=0) leaves 4, pair result 2; sliced (ctx=1)
  // leaves 2.
  std::vector<std::vector<double>> T(2, std::vector<double>(4, 0.0));
  T[0][0b01] = T[0][0b10] = 4;
  T[0][0b11] = 2;  // full
  T[1][0b01] = T[1][0b10] = 2;
  T[1][0b11] = 2;  // F sliced -> leaves halve
  std::vector<size_t> open_modes(4, 0);
  open_modes[0b01] = 1;
  open_modes[0b10] = 1;
  open_modes[0b11] = 0;  // F open in leaves, internal in pair
  std::vector<char> persistent(4, 1);
  // not batched: 4+4+2=10; batch F at the pair: leaves sliced 2 -> 2+2+2=6.
  // min 6.
  REQUIRE(batched_min_peak(T, open_modes, persistent, /*nt=*/2) == 6.0);
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

    SECTION("ObjectiveFunction aliases share enum values") {
      using sequant::ObjectiveFunction;
      // The deprecated names are aliases (same underlying value), so all
      // existing `== DensePeakSize` guards keep catching the renamed constant.
      STATIC_REQUIRE(ObjectiveFunction::DensePeakSize ==
                     ObjectiveFunction::DenseSpaceTime);
      STATIC_REQUIRE(ObjectiveFunction::DensePeakSizeBatched ==
                     ObjectiveFunction::DenseSpaceTimeBatched);
      // The new perf-first values are distinct from the peak-first ones.
      STATIC_REQUIRE(ObjectiveFunction::DenseTimeSpace !=
                     ObjectiveFunction::DenseSpaceTime);
      STATIC_REQUIRE(ObjectiveFunction::DenseTimeSpaceBatched !=
                     ObjectiveFunction::DenseSpaceTimeBatched);
      // ABI: DenseSpaceTime keeps DensePeakSize's old numeric value (2).
      STATIC_REQUIRE(static_cast<int>(ObjectiveFunction::DenseSpaceTime) == 2);
    }

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
          optimize(ex<Product>(prod),
                   OptimizeOptions{.inner_pow = {}, .footprint_weight = 0.0});
      REQUIRE(res0 == res0_explicit);  // weight 0 is a no-op

      // ... a large footprint_weight changes the chosen factorization.
      auto resF =
          optimize(ex<Product>(prod),
                   OptimizeOptions{.inner_pow = {}, .footprint_weight = 100.});
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
            prod, OptimizeOptions{.objective_function = objective_function,
                                  .inner_pow = {}});
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

      auto no_reorder = optimize(
          sum,
          OptimizeOptions{.reorder = ReorderSum::NoReorder, .inner_pow = {}});
      auto reorder = optimize(
          sum,
          OptimizeOptions{.reorder = ReorderSum::Reorder, .inner_pow = {}});
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
      auto aux = opt::detail::batchable_mode_list(net, is_batchable);
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
      auto aux = opt::detail::batchable_mode_list(net, is_batchable);
      std::size_t const m = aux.size();
      auto batch_fn = [batch](Index const&) -> std::size_t { return batch; };
      opt::detail::PeakBatchedModel model{idxsz, batch_fn,
                                          /*is_volatile_leaf=*/{}};
      model.is_batchable_contracted_index = is_batchable;
      model.is_batchable_external_index =
          is_batchable;  // external role (Task-4)
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
        auto aux = opt::detail::batchable_mode_list(net, is_batchable);
        auto vmask = opt::detail::leaf_volatile_mask(net, {});
        auto batch_fn = [batch](Index const&) -> std::size_t { return batch; };
        auto tables = opt::detail::sliced_footprints(
            net, targets, idxsz, is_batchable, batch_fn, aux);
        // open_modes[s] via the SAME detail helper the DP uses, so DP and
        // oracle index `tables` identically.
        auto open_aux_det = opt::detail::subset_open_aux(net, targets, aux);
        std::vector<size_t> open_modes(open_aux_det.begin(),
                                       open_aux_det.end());
        std::vector<std::vector<double>> T(tables.begin(), tables.end());
        std::vector<char> persistent(size_t{1} << ts.size());
        for (size_t n = 0; n < persistent.size(); ++n)
          persistent[n] = ((vmask & n) == 0) ? 1 : 0;
        double oracle = batched_min_peak(T, open_modes, persistent, ts.size());
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
      opts.batch_policy.is_batchable_contracted_index = [](Index const& ix) {
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

    auto binarizer = [](auto&& expr) {
      // CSE drives binarize() on subexpressions for hash-equivalence
      // detection; positional head is irrelevant here.
      SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
      return binarize(expr);
      SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    };

    auto collect_as_expr = [](auto&& expressions) {
      std::vector<ResultExpr> actual;
      for (const auto& current : expressions) {
        if (current->is_tensor()) {
          actual.emplace_back(current->expr()->template as<Tensor>(),
                              to_expr(current));
        } else {
          REQUIRE(current->is_scalar());
          actual.emplace_back(current->expr()->template as<Variable>(),
                              to_expr(current));
        }
      }

      return actual;
    };

    auto parse_inputs = [](auto&& inputs) {
      std::vector<EvalNode<EvalExpr>> expressions;
      for (const std::wstring& current : inputs) {
        expressions.push_back(binarize(deserialize<ResultExpr>(
            current, {.def_perm_symm = Symmetry::Nonsymm,
                      .def_braket_symm = BraKetSymmetry::Nonsymm,
                      .def_col_symm = ColumnSymmetry::Nonsymm})));
      }
      return expressions;
    };

    auto parse_expected = [](auto&& outputs) {
      std::vector<ResultExpr> expected;
      for (const std::wstring& current : outputs) {
        expected.push_back(deserialize<ResultExpr>(
            current, {.def_perm_symm = Symmetry::Nonsymm,
                      .def_braket_symm = BraKetSymmetry::Nonsymm,
                      .def_col_symm = ColumnSymmetry::Nonsymm}));
      }
      return expected;
    };

    SECTION("standard") {
      for (bool force_hash_collisions : {false, true}) {
        CAPTURE(force_hash_collisions);

        for (const auto& [inputs, outputs] :
             std::vector<std::pair<std::vector<std::wstring>,
                                   std::vector<std::wstring>>>{
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
                   L"t{a1,u5;i1,u4} + f{u3;u4}:N-S Y{u2,u4;u5,u1} "
                   L"t{a1,u5;i1,u3} "
                   L"+ f{u3;u4}:N-S Y{u2,u4;u1,u5} t{a1,u5;u3,i1}"},
                  {L"CSE1{;;u4,u2,u1,u5} = f{u3;u4}:N-S Y{u2,u3;u1,u5}",
                   L"R2{u2,a1;u1,i1} = -2 CSE1{;;u4,u2,u1,u5} t{a1,u5;i1,u4}"
                   L" + CSE1{;;u3,u2,u5,u1} t{a1,u5;i1,u3}"
                   L" + CSE1{;;u3,u2,u1,u5} t{a1,u5;u3,i1}"}},
             }) {
          CAPTURE(inputs);

          std::vector<EvalNode<EvalExpr>> expressions = parse_inputs(inputs);
          const std::vector<ResultExpr> expected = parse_expected(outputs);

          if (force_hash_collisions) {
            // This code path makes all hashes be computed to be zero and hence
            // every pair of objects will yield a hash collision which need to
            // be dealt with by using proper comparison operators.
            static constexpr bool force_collisions = true;
            opt::eliminate_common_subexpressions<
                decltype(expressions), decltype(binarizer), force_collisions>(
                expressions, binarizer);
          } else {
            opt::eliminate_common_subexpressions(expressions, binarizer);
          }

          REQUIRE(collect_as_expr(expressions) == expected);
        }
      }
    }
    SECTION("batch indices") {
      const opt::CSEOptions<EvalNode<EvalExpr>> opts = {.batch_indices = {
                                                            "i5",
                                                            "i6",
                                                            "a5",
                                                            "a6",
                                                        }};

      for (const auto& [inputs, outputs] : std::vector<
               std::pair<std::vector<std::wstring>, std::vector<std::wstring>>>{
               // Can't eliminate any CSE due to differences in batching indices
               {{L"R = (A{i1;i5} B{i5;i2}) C{i2;i1} + "
                 L"(A{i1;i6} B{i6;i2}) D{i2;i1}"},
                {L"R = (A{i1;i5} B{i5;i2}) C{i2;i1} + "
                 L"(A{i1;i6} B{i6;i2}) D{i2;i1}"}},
               {{L"R = (A{i1;i5} B{i5;i2}) C{i2;i1} + "
                 L"(A{i1;i3} B{i3;i2}) D{i2;i1}"},
                {L"R = (A{i1;i5} B{i5;i2}) C{i2;i1} + "
                 L"(A{i1;i3} B{i3;i2}) D{i2;i1}"}},
               // Can eliminate if batched index is same
               {{L"R = (A{i1;i5} B{i5;i2}) C{i2;i1} + "
                 L"(A{i3;i5} B{i5;i4}) D{i4;i3}"},
                {L"CSE1{;;i1,i2} = A{i1;i5} B{i5;i2}",
                 L"R = CSE1{;;i1,i2} C{i2;i1} + "
                 L"CSE1{;;i3,i4} D{i4;i3}"}},
           }) {
        CAPTURE(inputs);

        std::vector<EvalNode<EvalExpr>> expressions = parse_inputs(inputs);
        const std::vector<ResultExpr> expected = parse_expected(outputs);

        opt::eliminate_common_subexpressions(expressions, binarizer, opts);

        REQUIRE(collect_as_expr(expressions) == expected);
      }
    }
  }

  SECTION("Single term optimization with CSE") {
    auto ctx_resetter =
        set_scoped_default_context(get_default_context().clone());
    auto reg = get_default_context().mutable_index_space_registry();
    mbpt::add_df_spaces(reg);
    mbpt::add_pao_spaces(reg, mbpt::Spin::any);
    mbpt::add_ao_spaces(reg, mbpt::Spin::any);
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
                                         .idx_to_extent = idx_to_extent,
                                         .inner_pow = {}});
      auto res_no_cse =
          optimize(expr, OptimizeOptions{.CSE = {.subnet = false},
                                         .idx_to_extent = idx_to_extent,
                                         .inner_pow = {}});

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
      auto res_default = optimize(
          expr,
          OptimizeOptions{.idx_to_extent = idx_to_extent, .inner_pow = {}});
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

// Role-based batchability (domain-generic): the caller supplies two space sets
// -- is_batchable admits spaces batchable in the CONTRACTED role,
// is_batchable_external admits spaces batchable in the EXTERNAL role. In
//   R{i_1;a_1} = f{i_1;i_2} t{i_2;a_1}
// i_1 is external (open on the network root) while i_2 is contracted between f
// and t. With the occupied space admitted in the EXTERNAL role only,
// batchable_mode_list admits both as CANDIDATES (a mode's role is not known
// there) and build_context's role filter keeps i_1 but drops i_2 -- so the
// contracted occurrences of an external-only space never enter the 2^m search.
TEST_CASE(
    "batchable roles: external-only space keeps external, drops contracted",
    "[optimize][batch]") {
  using namespace sequant;
  auto ctx_resetter = set_scoped_default_context(get_default_context().clone());
  auto reg = get_default_context().mutable_index_space_registry();
  {
    auto occ = reg->retrieve_ptr(L"i");
    auto uocc = reg->retrieve_ptr(L"a");
    REQUIRE(occ);
    REQUIRE(uocc);
    occ->approximate_size(10);
    uocc->approximate_size(100);
  }

  auto is_batchable = [](Index const& ix) {
    return ix.space().base_key() == L"a";
  };
  std::function<bool(Index const&)> is_batchable_external =
      [](Index const& ix) { return ix.space().base_key() == L"i"; };

  auto f = deserialize(L"f{i_1;i_2}", {.def_perm_symm = Symmetry::Nonsymm});
  auto t = deserialize(L"t{i_2;a_1}", {.def_perm_symm = Symmetry::Nonsymm});
  TensorNetwork net{std::vector<ExprPtr>{f, t}};

  auto i1 = Index{L"i_1"};
  auto i2 = Index{L"i_2"};

  // Candidate list: roles unresolved, so BOTH occupied modes are admitted.
  auto aux = opt::detail::batchable_mode_list(net, is_batchable,
                                              is_batchable_external);
  REQUIRE(ranges::find(aux, i1) != ranges::end(aux));
  REQUIRE(ranges::find(aux, i2) != ranges::end(aux));

  // build_context resolves each mode's role and applies the role filter.
  auto idxsz = [](Index const& ix) -> std::size_t {
    return ix.nonnull() ? ix.space().approximate_size() : std::size_t{1};
  };
  auto batch_fn = [](Index const&) -> std::size_t { return 4; };
  opt::detail::PeakBatchedModel<decltype(idxsz)> model{idxsz, batch_fn, {}};
  model.is_batchable_contracted_index = is_batchable;
  model.is_batchable_external_index = is_batchable_external;
  container::svector<Index> const targets{};
  auto ctx = model.build_context(net, targets);
  // i_1 occurs external => kept by the external-role set.
  CHECK(ranges::find(ctx.batchable_modes, i1) !=
        ranges::end(ctx.batchable_modes));
  // i_2 occurs contracted, and the occupied space is not contracted-batchable.
  CHECK(ranges::find(ctx.batchable_modes, i2) ==
        ranges::end(ctx.batchable_modes));
}

// Ground-truth for the batchability ROLE SPLIT (the behavioral delta callers
// rely on). A batchable mode that occurs in the CONTRACTED role is DP-sliced
// only when its space is admitted in the CONTRACTED role -- NOT when its space
// is admitted only in the EXTERNAL role. This is exactly what lets an
// external-batching caller (occ-like spectator) declare its mode external
// WITHOUT the DP slicing that mode in a contracted cell the runtime never
// realizes.
//
// Network g{a1;F1} h{F1;F2} t{F2;a2}: F1 (shared by g,h) and F2 (shared by h,t)
// are CONTRACTED batchable modes of space "F"; a1,a2 are external virtual (not
// batchable). Two arms, ONE test:
//   (control) F admitted in the CONTRACTED role => F1/F2 are DP modes and some
//             DP cell slices them: the DP CAN and DOES slice a contracted F.
//   (fix)     F admitted in the EXTERNAL role ONLY (contracted predicate empty)
//             => the role filter drops F1/F2, so NO DP cell slices them.
// The control arm is what makes the fix assertion NON-VACUOUS: it proves the
// modes are genuinely sliceable, so their absence under the external-only role
// is a real role-gating effect, not a "nothing was ever sliceable" artifact.
// Equivalently: leaving a contracted mode in the contracted predicate (the
// pre-migration routing) DOES slice it (control, > 0); moving it to the
// external role does NOT (fix, == 0).
TEST_CASE("role filter: contracted mode sliced only in the contracted role",
          "[optimize][role-filter]") {
  using namespace sequant;
  namespace o = sequant::opt::detail;
  auto ctx_clone = get_default_context().clone();
  ctx_clone.mutable_index_space_registry()->add(L"F", IndexSpace::Type{0b10000},
                                                8ul);
  auto ctx_resetter = set_scoped_default_context(std::move(ctx_clone));

  auto idxsz = [](Index const& ix) { return ix.space().approximate_size(); };
  auto batch_fn = [](Index const&) -> std::size_t { return 2; };
  std::function<bool(Index const&)> is_F = [](Index const& ix) {
    return ix.space().base_key() == L"F";
  };

  std::vector<ExprPtr> ts;
  for (auto s : {L"g{a1;F1}", L"h{F1;F2}", L"t{F2;a2}"})
    ts.push_back(deserialize(s, {.def_perm_symm = Symmetry::Nonsymm}));
  TensorNetwork net{ts};
  container::svector<Index> const
      targets{};  // a1,a2 external; F1,F2 contracted

  auto count_F_modes = [](auto const& ctx) -> std::size_t {
    std::size_t n = 0;
    for (auto const& ix : ctx.batchable_modes)
      if (ix.space().base_key() == L"F") ++n;
    return n;
  };
  // Number of DP cells whose sliced-set carries at least one "F" mode.
  auto n_F_sliced_cells = [](auto const& ctx) -> std::size_t {
    std::size_t mask = 0;
    for (std::size_t k = 0; k < ctx.m; ++k)
      if (ctx.batchable_modes[k].space().base_key() == L"F")
        mask |= (std::size_t{1} << k);
    std::size_t count = 0;
    for (std::size_t id = 0; id < ctx.nCells; ++id)
      if (ctx.cell_union(id) & mask) ++count;
    return count;
  };

  // control: F batchable in the CONTRACTED role.
  o::PeakBatchedModel<decltype(idxsz)> control{idxsz, batch_fn, {}};
  control.is_batchable_contracted_index = is_F;
  control.order_aware_recompute = true;  // ordered cells (the production path)
  auto cctx = control.build_context(net, targets);
  CHECK(count_F_modes(cctx) == 2);    // F1,F2 are contracted DP modes
  CHECK(n_F_sliced_cells(cctx) > 0);  // ... and the DP slices them

  // fix: F batchable in the EXTERNAL role ONLY (contracted predicate empty).
  o::PeakBatchedModel<decltype(idxsz)> fixed{idxsz, batch_fn, {}};
  fixed.is_batchable_external_index = is_F;
  fixed.order_aware_recompute = true;
  auto fctx = fixed.build_context(net, targets);
  CHECK(count_F_modes(fctx) == 0);     // dropped by the role filter
  CHECK(n_F_sliced_cells(fctx) == 0);  // ... so NO cell slices an F mode
}

// Task 4 [role-filter]: the historical external->contracted fallback is GONE.
// A space admitted ONLY in the contracted role, whose mode occurs EXTERNALLY
// (open on the term root), must NOT be batchable in the external role: with the
// external building block at its default (decline) there is no fallback to the
// contracted predicate, so build_context's role filter drops that external
// occurrence.
//
// RED verification (fallback present, i.e. BASE): the empty-default external
// predicate fell back to the contracted predicate, so the external F mode WAS
// admitted and this CHECK saw n_F == 1. After the fallback removal +
// default-decline the mode is dropped (n_F == 0). The RED run also proves the
// assertion is non-vacuous: the mode is a genuine batchable-space mode that the
// fallback DID admit.
TEST_CASE("role filter: external mode needs the external role, not a fallback",
          "[optimize][role-filter]") {
  using namespace sequant;
  namespace o = sequant::opt::detail;
  auto ctx_clone = get_default_context().clone();
  ctx_clone.mutable_index_space_registry()->add(L"F", IndexSpace::Type{0b10000},
                                                8ul);
  auto ctx_resetter = set_scoped_default_context(std::move(ctx_clone));

  auto idxsz = [](Index const& ix) { return ix.space().approximate_size(); };
  auto batch_fn = [](Index const&) -> std::size_t { return 2; };
  std::function<bool(Index const&)> is_F = [](Index const& ix) {
    return ix.space().base_key() == L"F";
  };

  // F1 occurs once (open on the root => EXTERNAL); a1 is shared (contracted).
  std::vector<ExprPtr> ts;
  for (auto s : {L"g{a1;F1}", L"h{a1;a2}"})
    ts.push_back(deserialize(s, {.def_perm_symm = Symmetry::Nonsymm}));
  TensorNetwork net{ts};
  container::svector<Index> const targets{};  // a2, F1 external; a1 contracted

  // F admitted ONLY in the contracted role; the external role stays at its
  // default (decline). The single external F1 occurrence must be dropped -- the
  // external role no longer falls back to the contracted predicate.
  o::PeakBatchedModel<decltype(idxsz)> model{idxsz, batch_fn, {}};
  model.is_batchable_contracted_index = is_F;
  auto ctx = model.build_context(net, targets);
  std::size_t n_F = 0;
  for (auto const& ix : ctx.batchable_modes)
    if (ix.space().base_key() == L"F") ++n_F;
  CHECK(n_F == 0);  // external F1 NOT admitted via a (now-removed) fallback
}

TEST_CASE("OSV early-contraction reproducer", "[optimize][osv]") {
  using namespace sequant;
  auto ctx_resetter = set_scoped_default_context(get_default_context().clone());
  auto reg = get_default_context().mutable_index_space_registry();
  mbpt::add_df_spaces(reg);
  mbpt::add_pao_spaces(reg, mbpt::Spin::any);
  mbpt::add_ao_spaces(reg, mbpt::Spin::any);
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
        prod, idxsz, /*subnet_cse=*/false,
        CostParams{.volatile_weight = 1.0, .inner_pow = ip});
    std::wcout << name << L":  " << render_tree(res) << L"\n";
  };
  // Omitting inner_pow on a term with composite (CSV/PNO) indices used to
  // silently size the composite `a` at its full uocc extent (40), inverting
  // the factorization (this is the 4-PAO-integral class of bug). That path is
  // now forbidden: the sizing code throws instead of guessing.
  std::wcout
      << L"--- without inner_pow: now REJECTED (composite present) ---\n";
  CHECK_THROWS_AS(show(std::integral_constant<ObjectiveFunction,
                                              ObjectiveFunction::DenseFLOPs>{},
                       L"FLOPs"),
                  std::invalid_argument);

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
  mbpt::add_pao_spaces(reg, mbpt::Spin::any);
  mbpt::add_ao_spaces(reg, mbpt::Spin::any);
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
        prod, idxsz, false,
        CostParams{.volatile_weight = 1.0, .inner_pow = ip});
    std::wcout << L"FLOPs:        " << render_tree(res) << L"\n";
  }
  {
    auto res = opt::single_term_opt<ObjectiveFunction::DensePeakSize>(
        prod, idxsz, false,
        CostParams{.volatile_weight = 1.0, .inner_pow = ip});
    std::wcout << L"PeakSize:     " << render_tree(res) << L"\n";
  }
  {
    auto res = opt::single_term_opt<ObjectiveFunction::DensePeakSizeBatched>(
        prod, idxsz, false,
        CostParams{.volatile_weight = 1.0,
                   .is_batchable_contracted_index = is_batch,
                   .batch_target_size = bts,
                   .inner_pow = ip});
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
  mbpt::add_pao_spaces(reg, mbpt::Spin::any);
  mbpt::add_ao_spaces(reg, mbpt::Spin::any);
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
                 prod, idxsz, false,
                 CostParams{.volatile_weight = 1.0, .inner_pow = ip}));
  double peak_mx =
      report(L"PeakSize:     ",
             opt::single_term_opt<ObjectiveFunction::DensePeakSize>(
                 prod, idxsz, false,
                 CostParams{.volatile_weight = 1.0, .inner_pow = ip}));
  double pbat_mx =
      report(L"PeakBatched:  ",
             opt::single_term_opt<ObjectiveFunction::DensePeakSizeBatched>(
                 prod, idxsz, false,
                 CostParams{.volatile_weight = 1.0,
                            .is_batchable_contracted_index = is_batch,
                            .batch_target_size = bts,
                            .inner_pow = ip}));
  // Real MPQC config: t is volatile, volatile_weight=100. The tie-break weights
  // volatile (replayed) flops, so it must STILL eliminate the OSV early.
  auto is_t = [](Tensor const& t) { return t.label() == L"t"; };
  // isolate volatile_weight: hold is_volatile_leaf=is_t FIXED, vary only vw.
  double peak_v = report(L"PeakSize/is_t,vw1:   ",
                         opt::single_term_opt<ObjectiveFunction::DensePeakSize>(
                             prod, idxsz, false,
                             CostParams{.is_volatile_leaf = is_t,
                                        .volatile_weight = 1.0,
                                        .inner_pow = ip}));
  double peak_v100 =
      report(L"PeakSize/is_t,vw100: ",
             opt::single_term_opt<ObjectiveFunction::DensePeakSize>(
                 prod, idxsz, false,
                 CostParams{.is_volatile_leaf = is_t,
                            .volatile_weight = 100.0,
                            .inner_pow = ip}));
  double pbat_v =
      report(L"PeakBatch/is_t,vw1:  ",
             opt::single_term_opt<ObjectiveFunction::DensePeakSizeBatched>(
                 prod, idxsz, false,
                 CostParams{.is_volatile_leaf = is_t,
                            .volatile_weight = 1.0,
                            .is_batchable_contracted_index = is_batch,
                            .batch_target_size = bts,
                            .inner_pow = ip}));
  double pbat_v100 =
      report(L"PeakBatch/is_t,vw100:",
             opt::single_term_opt<ObjectiveFunction::DensePeakSizeBatched>(
                 prod, idxsz, false,
                 CostParams{.is_volatile_leaf = is_t,
                            .volatile_weight = 100.0,
                            .is_batchable_contracted_index = is_batch,
                            .batch_target_size = bts,
                            .inner_pow = ip}));
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
  //     product). This specifically probes peak-driven (memory-constrained)
  //     selection under the gate, so it needs a finite (near-zero)
  //     peak_threshold: the default +inf would pick the min-flops schedule
  //     regardless of the gate (flops do not depend on which nodes may be
  //     sliced), masking the effect the gate is meant to demonstrate.
  CostParams po_cost{is_t, 100.0, 0.0};
  po_cost.peak_threshold = 1.0;  // near-zero => infeasible => min-peak fallback
  po_cost.is_batchable_contracted_index = is_batch;
  po_cost.batch_target_size = bts;
  po_cost.inner_pow = ip;
  po_cost.batch_persistent_only = true;
  double pbat_po =
      report(L"PeakBatch/persistent_only:",
             opt::single_term_opt<ObjectiveFunction::DensePeakSizeBatched>(
                 prod, idxsz, false, po_cost));
  CHECK(pbat_po >= osv_outer_product);  // gate restored -> OSV deferred again
}

// C60 PNO-CCSD residual "member 2": the batched intermediate that OOMs (~185 GB
// materialized) is the double-proto I(i_1,i_2,i_3,Κ; a_3<i_2,i_3>,
// a_1<i_1,i_2>).
//   flat product:  g{i_3;i_1;Κ} · C{a_1<i_1,i_2>;μ̃} · g{μ̃;μ̃;Κ} ·
//   C{μ̃;a_3<i_2,i_3>}
// Two 2-occ PNO composites sharing i_2, no t amplitude. This probe answers:
//  (1) does STO form the double-proto intermediate, and under which objective;
//  (2) how does the (block-sparse) cost model size it (max_imed);
//  (3) do the rejected brackets avoid it, i.e. does a cheaper schedule exist?
TEST_CASE("C60 member-2 double-proto probe", "[optimize][osv][c60]") {
  using namespace sequant;
  auto ctx_resetter = set_scoped_default_context(get_default_context().clone());
  auto reg = get_default_context().mutable_index_space_registry();
  mbpt::add_df_spaces(reg);
  mbpt::add_pao_spaces(reg, mbpt::Spin::any);
  mbpt::add_ao_spaces(reg, mbpt::Spin::any);
  // C60 cc-pVDZ-F12-ish extents: active occ 120, PNO domain 42, PAO 1800,
  // DF aux 6000 (sliced to batch 30 by the batched objective).
  for (auto&& [k, v] :
       std::initializer_list<std::pair<std::wstring_view, size_t>>{
           {L"i", 120}, {L"a", 42}, {L"μ̃", 1800}, {L"Κ", 6000}}) {
    reg->retrieve_ptr(k)->approximate_size(v);
  }
  auto aux_space = reg->retrieve(L"Κ");
  auto idxsz = [](Index const& ix) -> std::size_t {
    return ix.nonnull() ? ix.space().approximate_size() : std::size_t{1};
  };
  // constant PNO domain 42 => all power means equal, neutralizing the
  // RMS-vs-mean magnitude question; here we test STRUCTURE (which tree, is the
  // double-proto avoidable), not the heavy-tail correction.
  auto ip = [](Index const&, std::size_t) -> double { return 42.0; };
  auto is_batch = [aux_space](Index const& ix) {
    return ix.space() == aux_space;
  };
  std::function<std::size_t(Index const&)> bts = [](Index const&) {
    return std::size_t{30};
  };

  auto prod = deserialize(
                  L"g{i_3;i_1;Κ_1} C{a_1<i_1,i_2>;μ̃_1} "
                  L"g{μ̃_1;μ̃_2;Κ_1} C{μ̃_2;a_3<i_2,i_3>}")
                  ->as<Product>();

  auto fp = opt::detail::footprint_counter(
      idxsz, std::function<double(Index const&, std::size_t)>(ip));

  std::function<std::vector<Index>(ExprPtr const&, double&)> freeix =
      [&](ExprPtr const& e, double& mx) -> std::vector<Index> {
    if (e->is<Tensor>()) {
      std::vector<Index> v;
      for (auto const& ix : e->as<Tensor>().const_braketaux_indices())
        v.push_back(ix);
      return v;
    }
    std::map<std::wstring, std::pair<int, Index>> cnt;
    for (auto const& fct : e->as<Product>().factors())
      for (auto const& ix : freeix(fct, mx)) {
        auto k = std::wstring(ix.full_label());
        cnt[k].first++;
        cnt[k].second = ix;
      }
    std::vector<Index> result;
    for (auto const& [k, v] : cnt)
      if (v.first == 1) result.push_back(v.second);
    double here = fp(result);
    if (here > mx) mx = here;
    return result;
  };
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
    double Snode = fp(freeix(e, dummy));
    double both = L.S + R.S + Snode;
    double lfirst = std::max({R.leafsum + L.peak, L.S + R.peak, both});
    double rfirst = std::max({L.leafsum + R.peak, R.S + L.peak, both});
    return TP{std::min(lfirst, rfirst), L.leafsum + R.leafsum, Snode};
  };
  auto report = [&](std::wstring name, ExprPtr res) {
    double mx = 0.0;
    freeix(res, mx);
    double tpk = tree_peak(res).peak;
    std::wcout << name << L"  max_imed(elems)=" << (long long)mx
               << L"  recurrence_PEAK(elems)=" << (long long)tpk << L"\n      "
               << render_tree(res) << L"\n";
  };

  std::wcout
      << L"\n=== C60 MEMBER-2 (i=120, a=42, μ̃=1800, Κ=6000, batch=30) ===\n";
  report(L"FLOPs:       ",
         opt::single_term_opt<ObjectiveFunction::DenseFLOPs>(
             prod, idxsz, false,
             CostParams{.volatile_weight = 1.0, .inner_pow = ip}));
  report(L"PeakSize:    ",
         opt::single_term_opt<ObjectiveFunction::DensePeakSize>(
             prod, idxsz, false,
             CostParams{.volatile_weight = 1.0, .inner_pow = ip}));
  report(L"PeakBatched: ",
         opt::single_term_opt<ObjectiveFunction::DensePeakSizeBatched>(
             prod, idxsz, false,
             CostParams{.volatile_weight = 1.0,
                        .is_batchable_contracted_index = is_batch,
                        .batch_target_size = bts,
                        .inner_pow = ip}));
  // Experiment: make PAO (μ̃) batchable TOO (batch 100). Prediction: STO
  // switches away from the double-proto to a g·g-first tree whose worst
  // intermediate carries μ̃ and is now sliceable.
  auto mu_space = reg->retrieve(L"μ̃");
  auto is_batch2 = [aux_space, mu_space](Index const& ix) {
    return ix.space() == aux_space || ix.space() == mu_space;
  };
  std::function<std::size_t(Index const&)> bts2 = [aux_space](Index const& ix) {
    return ix.space() == aux_space ? std::size_t{30} : std::size_t{100};
  };
  report(L"PeakBatched(+μ̃):",
         opt::single_term_opt<ObjectiveFunction::DensePeakSizeBatched>(
             prod, idxsz, false,
             CostParams{.volatile_weight = 1.0,
                        .is_batchable_contracted_index = is_batch2,
                        .batch_target_size = bts2,
                        .inner_pow = ip}));
  // Reference: the double-proto intermediate as the model sizes it, dense at
  // d=42 (== the block-sparse upper bound), with Κ at one batch of 30:
  //   i^3 * Κ_batch * a * a  elements.
  std::wcout << L"double-proto imed (i^3*Κbatch*a^2) elems = "
             << (long long)(120.0 * 120 * 120 * 30 * 42 * 42) << L"  (="
             << (long long)(120.0 * 120 * 120 * 30 * 42 * 42 * 8)
             << L" B dense)\n\n";
}

TEST_CASE("PPL: form 4-PNO W vs fold-t (peak-neutral, flop tie-break)",
          "[optimize][osv]") {
  using namespace sequant;
  auto ctx_resetter = set_scoped_default_context(get_default_context().clone());
  auto reg = get_default_context().mutable_index_space_registry();
  mbpt::add_df_spaces(reg);
  mbpt::add_pao_spaces(reg, mbpt::Spin::any);
  mbpt::add_ao_spaces(reg, mbpt::Spin::any);
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
  double flops_w = rep(L"FLOPs/vw100:  ",
                       opt::single_term_opt<ObjectiveFunction::DenseFLOPs>(
                           prod, idxsz, false,
                           CostParams{.is_volatile_leaf = is_t,
                                      .volatile_weight = 100.0,
                                      .inner_pow = ip}));
  double peak_w = rep(L"PeakSize/vw100:",
                      opt::single_term_opt<ObjectiveFunction::DensePeakSize>(
                          prod, idxsz, false,
                          CostParams{.is_volatile_leaf = is_t,
                                     .volatile_weight = 100.0,
                                     .inner_pow = ip}));
  double pbat_w =
      rep(L"PeakB/vw100:  ",
          opt::single_term_opt<ObjectiveFunction::DensePeakSizeBatched>(
              prod, idxsz, false,
              CostParams{.is_volatile_leaf = is_t,
                         .volatile_weight = 100.0,
                         .is_batchable_contracted_index = is_batch,
                         .batch_target_size = bts,
                         .inner_pow = ip}));
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
              CostParams{.is_volatile_leaf = is_t,
                         .volatile_weight = 100.0,
                         .peak_flops_tolerance = 0.0,
                         .inner_pow = ip}));
  CHECK(peak_w_strict > flops_w);
  // At volatile_weight=1 (the caching-off regime, where persistent
  // intermediates cannot be amortized across replays) the persistent 4-PNO W is
  // strictly dominated by fold-t on BOTH peak and flops: the batched root
  // frontier collapses to the single fold-t point, so every objective folds the
  // amplitude t into a Kappa-batched half-transform ladder instead of
  // forming/caching W. (mpqc's SeQuantEngine pins volatile_weight to 1 when
  // eval:cache is off for exactly this reason.) A form-W tree would reproduce
  // flops_w; fold-t does not.
  double pbat_w_vw1 =
      rep(L"PeakB/vw1:    ",
          opt::single_term_opt<ObjectiveFunction::DensePeakSizeBatched>(
              prod, idxsz, false,
              CostParams{.is_volatile_leaf = is_t,
                         .volatile_weight = 1.0,
                         .is_batchable_contracted_index = is_batch,
                         .batch_target_size = bts,
                         .inner_pow = ip}));
  CHECK(pbat_w_vw1 > flops_w);

  // Perf-first (DenseTimeSpaceBatched) is flops-primary, so it forms the
  // flop-optimal W REGARDLESS of peak_flops_tolerance and regardless of a tight
  // peak_threshold: its replay-weighted flop count equals the FLOPs optimum.
  // (Peak-first at strict tolerance did NOT: peak_w_strict > flops_w above.)
  double tbat_w =
      rep(L"TimeSpaceB/vw100:",
          opt::single_term_opt<ObjectiveFunction::DenseTimeSpaceBatched>(
              prod, idxsz, false,
              CostParams{.is_volatile_leaf = is_t,
                         .volatile_weight = 100.0,
                         .is_batchable_contracted_index = is_batch,
                         .batch_target_size = bts,
                         .inner_pow = ip}));
  CHECK(tbat_w == flops_w);
  // Even with a near-zero peak budget, perf-first stays flop-optimal (threshold
  // is not a feasibility gate), unlike peak-first's min-peak fallback.
  CostParams tight{is_t, 100.0, 0.0};
  tight.peak_threshold = 1.0;
  tight.is_batchable_contracted_index = is_batch;
  tight.batch_target_size = bts;
  tight.inner_pow = ip;
  double tbat_w_tight =
      rep(L"TimeSpaceB/tight:",
          opt::single_term_opt<ObjectiveFunction::DenseTimeSpaceBatched>(
              prod, idxsz, false, tight));
  CHECK(tbat_w_tight == flops_w);
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
// Exercises single-term optimization on a 12-leaf network. The DP is ~100x
// slower in Debug (~100 s per single_term_opt call, tens of minutes total) and
// times out CI's Debug/Valgrind/Sanitizer jobs; gate it to optimized builds,
// where the whole test runs in seconds. NDEBUG is defined in all SeQuant build
// types (asserts are governed separately by SEQUANT_ASSERT_BEHAVIOR_), so it
// cannot distinguish Debug from Release; __OPTIMIZE__ (defined by GCC and Clang
// only at -O1 and above, never at -O0) is the reliable discriminator. See
// PR #559.
#ifdef __OPTIMIZE__
TEST_CASE("quadratic bubble: early-K integral vs late-K t·(gC)",
          "[optimize][bubble]") {
  using namespace sequant;
  auto ctx_resetter = set_scoped_default_context(get_default_context().clone());
  auto reg = get_default_context().mutable_index_space_registry();
  mbpt::add_df_spaces(reg);
  mbpt::add_pao_spaces(reg, mbpt::Spin::any);
  mbpt::add_ao_spaces(reg, mbpt::Spin::any);
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
    // peak_flops_tolerance is no longer consulted by DensePeakSizeBatched's
    // final selection (Task 1.2: threshold-gated). This probe is inherently
    // about PEAK (which factorization is smaller in memory, not flops), so it
    // needs a near-zero peak_threshold to force the min-peak fallback path
    // (the default +inf would instead pick purely by flops, masking the
    // crossover this test demonstrates).
    CostParams cost{{}, 1.0, 0.0, 0.0, {}, lambda};
    cost.peak_threshold = 1.0;
    cost.is_batchable_contracted_index = is_batch;
    cost.batch_target_size = bts;
    cost.inner_pow = ip;
    auto res = opt::single_term_opt<ObjectiveFunction::DensePeakSizeBatched>(
        prod, idxsz, false, cost);
    bool integral = forms_integral(res);
    double gb = batched_peak(res, Kb) * 8.0 / 1e9;
    std::wcout << L"  K_b=" << Kb << L"\tlambda=" << lambda
               << (integral ? L"\tEARLY-K (integral)" : L"\tLATE-K  (t.(gC))")
               << L"\tpeak=" << gb << L" GB\n";
    return integral;
  };

  std::wcout << L"\n=== quadratic bubble exchange: i=80 mu=860 K=2360 p=" << p
             << L" (crossover K_b=p/2=" << p / 2 << L") ===\n";
  // Crossover at K_b = p/2 = 6. With Kappa sliced to a tiny batch the
  // held-whole integral cannot win (o^4*p^2 > 2*o^4*p*K_b for K_b < p/2), so
  // the optimizer takes the batchable late-K route; with a peak-trading premium
  // (lambda) and K_b above the crossover the held-whole (g·C)(g·C) integral
  // (early-K) wins.
  CHECK_FALSE(choose(/*K_b=*/2, /*lambda=*/0.0));  // below crossover -> late-K
  CHECK_FALSE(choose(/*K_b=*/6, /*lambda=*/1.0));  // at crossover   -> late-K
  CHECK(choose(/*K_b=*/72, /*lambda=*/10.0));      // above crossover -> early-K

  // Real MPQC config: t is volatile (replayed). This probes the same PEAK
  // crossover under replay-weighted flops (is_t, volatile_weight=100), so it
  // needs the same near-zero peak_threshold as choose() above: the default
  // +inf would instead pick purely by (replay-weighted) flops, which favor
  // forming the persistent, t-free, Kappa-free integral ONCE REGARDLESS of
  // K_b (flops do not depend on the aux batch size) -- exactly the behavior
  // exercised by the "threshold gates batching" test (Task 1.2, below), which
  // reuses this same motif. With peak_threshold forced near-zero (min-peak
  // fallback), this instead reproduces the peak-driven crossover: below the
  // crossover it stays late-K; above it flips to early-K -- the suspected
  // real driver of the held-whole 4-occ/2-PNO object (the C60 OOM),
  // independent of accumulation.
  auto is_t = [](Tensor const& t) { return t.label() == L"t"; };
  auto real_config_integral = [&](std::size_t Kb) -> bool {
    std::function<std::size_t(Index const&)> bts = [Kb](Index const&) {
      return Kb;
    };
    CostParams cost{is_t, 100.0, 0.0, 0.10, {}, 0.0};
    cost.peak_threshold = 1.0;
    cost.is_batchable_contracted_index = is_batch;
    cost.batch_target_size = bts;
    cost.inner_pow = ip;
    auto res = opt::single_term_opt<ObjectiveFunction::DensePeakSizeBatched>(
        prod, idxsz, false, cost);
    return forms_integral(res);
  };
  CHECK_FALSE(real_config_integral(/*K_b=*/2));  // below crossover -> late-K
  CHECK(real_config_integral(/*K_b=*/236));      // above crossover -> early-K
}
#endif  // __OPTIMIZE__

// Task 1.2: threshold-gated root selection in PeakBatchedModel::reconstruct.
// Reuses the "quadratic bubble" motif (early-K integral vs late-K t.(gC),
// above) since it is a proven, already-tuned case where the OLD
// (peak-first, epsilon-tolerant) root selection and the flop-optimal choice
// diverge: at K_b=2 (below the peak crossover), the replay-weighted flops
// (is_t volatile, volatile_weight=100) still favor forming the persistent
// Kappa-free integral ONCE (early-K) -- flops do not depend on K_b at all,
// only peak does -- but the old peak_flops_tolerance=0.10 band rejects it
// there because its peak exceeds the tolerance band around the (much
// smaller) late-K peak (see CHECK_FALSE(real_config_integral(2)) above).
// The new threshold-gated selection no longer consults peak_flops_tolerance:
// with the default peak_threshold = +inf every root-frontier point is
// feasible, so pure min-flops wins, forming the integral even at this small
// K_b. A near-zero peak_threshold instead makes every point infeasible,
// triggering the min-peak fallback, which reproduces the old (peak-driven)
// answer.
TEST_CASE(
    "threshold gates batching: default (+inf) picks min-flops regardless of "
    "K_b; near-zero threshold falls back to min-peak",
    "[optimize][threshold]") {
  using namespace sequant;
  auto ctx_resetter = set_scoped_default_context(get_default_context().clone());
  auto reg = get_default_context().mutable_index_space_registry();
  mbpt::add_df_spaces(reg);
  mbpt::add_pao_spaces(reg, mbpt::Spin::any);
  mbpt::add_ao_spaces(reg, mbpt::Spin::any);
  // water-20-scale extents (matches the "quadratic bubble" test above).
  for (auto&& [k, v] :
       std::initializer_list<std::pair<std::wstring_view, size_t>>{
           {L"i", 80}, {L"a", 12}, {L"μ̃", 860}, {L"Κ", 2360}}) {
    reg->retrieve_ptr(k)->approximate_size(v);
  }
  auto aux_space = reg->retrieve(L"Κ");
  auto idxsz = [](Index const& ix) -> std::size_t {
    return ix.nonnull() ? ix.space().approximate_size() : std::size_t{1};
  };
  double const p = 12.0;
  auto ip = [p](Index const&, std::size_t) -> double { return p; };
  auto is_batch = [aux_space](Index const& ix) {
    return ix.space() == aux_space;
  };

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

  auto is_t = [](Tensor const& t) { return t.label() == L"t"; };
  auto integral_at = [&](std::size_t Kb, double peak_threshold) -> bool {
    std::function<std::size_t(Index const&)> bts = [Kb](Index const&) {
      return Kb;
    };
    CostParams cost{is_t, 100.0, 0.0};
    cost.peak_threshold = peak_threshold;
    cost.is_batchable_contracted_index = is_batch;
    cost.batch_target_size = bts;
    cost.inner_pow = ip;
    auto res = opt::single_term_opt<ObjectiveFunction::DensePeakSizeBatched>(
        prod, idxsz, false, cost);
    return forms_integral(res);
  };

  double const inf = std::numeric_limits<double>::infinity();
  // Below the peak crossover (K_b=2): default peak_threshold (+inf) => every
  // root point is feasible => pure min-(replay-weighted)-flops selection =>
  // forms the persistent integral (early-K), unlike the old peak-tolerance-
  // band selection (CHECK_FALSE(real_config_integral(2)) in the test above).
  CHECK(integral_at(/*Kb=*/2, inf));
  // A near-zero peak_threshold makes every point infeasible => min-peak
  // fallback => reproduces the old (peak-driven) answer: late-K.
  CHECK_FALSE(integral_at(/*Kb=*/2, /*peak_threshold(bytes)=*/1.0));

  // Perf-first (DenseTimeSpaceBatched): peak_threshold is a CEILING, not a
  // min-peak objective. At the two ENDPOINTS the selection is min-flops either
  // way -- +inf makes every point feasible (min-flops among all), and a
  // near-zero threshold makes every point infeasible (perf-first best-effort
  // fallback = GLOBAL min-flops, NOT peak-first's min-peak). So both endpoints
  // form the persistent Kappa-free integral. Contrast the peak-first checks
  // above, whose near-zero fallback is min-peak (CHECK_FALSE). Only an
  // INTERMEDIATE budget (below the min-flops peak, above a batched schedule's
  // peak) makes perf-first batch -- exercised by
  // reconstruct_batched_modes_emits_external_per_node's Assertion 1b.
  auto integral_at_perf = [&](std::size_t Kb, double peak_threshold) -> bool {
    std::function<std::size_t(Index const&)> bts = [Kb](Index const&) {
      return Kb;
    };
    CostParams cost{is_t, 100.0, 0.0};
    cost.peak_threshold = peak_threshold;
    cost.is_batchable_contracted_index = is_batch;
    cost.batch_target_size = bts;
    cost.inner_pow = ip;
    auto res = opt::single_term_opt<ObjectiveFunction::DenseTimeSpaceBatched>(
        prod, idxsz, false, cost);
    return forms_integral(res);
  };
  // Endpoint-insensitive: forms the integral at +inf AND at near-zero (both
  // resolve to min-flops -- see the ceiling note above).
  CHECK(integral_at_perf(/*Kb=*/2, inf));
  CHECK(integral_at_perf(/*Kb=*/2, /*peak_threshold(bytes)=*/1.0));
}

// Task 2.2: validate that the batched DP's accumulation_factor charge is
// priced correctly under NESTED accumulation -- i.e. more than one batchable
// index, contracted at *different* nodes of the binarized tree -- and that
// the ctx.m <= 1 restriction that used to guard accumulation_factor != 0 can
// be lifted. Network: T0=g{i;a;K}, T1=g{mu1;mu2;K}, T2=f{mu1;a2},
// T3=f{mu2;a3}. T1 is a hub carrying K, mu1 AND mu2 simultaneously; K is
// contracted between T0/T1, mu1 between T1/T2, mu2 between T1/T3 -- three
// distinct batchable Index instances (K, mu1, mu2) spanning the aux (Κ) and
// PAO (μ̃) spaces, forcing the DP to slice more than one mode, at more than
// one node, to reach its minimum-peak schedule.
TEST_CASE("batched DP peak matches oracle with two modes and accumulation",
          "[optimize][batched-accum]") {
  using namespace sequant;
  auto ctx_resetter = set_scoped_default_context(get_default_context().clone());
  auto reg = get_default_context().mutable_index_space_registry();
  mbpt::add_df_spaces(reg);
  mbpt::add_pao_spaces(reg, mbpt::Spin::any);
  mbpt::add_ao_spaces(reg, mbpt::Spin::any);
  for (auto&& [k, v] :
       std::initializer_list<std::pair<std::wstring_view, size_t>>{
           {L"i", 20}, {L"a", 20}, {L"μ̃", 200}, {L"Κ", 300}}) {
    reg->retrieve_ptr(k)->approximate_size(v);
  }
  auto aux = reg->retrieve(L"Κ");
  auto pao = reg->retrieve(L"μ̃");
  auto idxsz = [](Index const& ix) -> std::size_t {
    return ix.nonnull() ? ix.space().approximate_size() : std::size_t{1};
  };
  auto is_batch = [aux, pao](Index const& ix) {
    return ix.space() == aux || ix.space() == pao;
  };
  auto is_batch_K = [aux](Index const& ix) { return ix.space() == aux; };
  auto is_batch_mu = [pao](Index const& ix) { return ix.space() == pao; };
  std::function<std::size_t(Index const&)> bts = [](Index const&) {
    return std::size_t{10};
  };
  std::function<bool(Tensor const&)> novol = {};
  // Two batchable contracted indices K and mu-tilde, contracted at different
  // nodes (nested accumulation).
  auto prod = deserialize(
                  L"g{i_1;a_1;Κ_1} g{μ̃_1;μ̃_2;Κ_1} f{μ̃_1;a_2} "
                  L"f{μ̃_2;a_3}")
                  ->as<Product>();
  TensorNetwork net(prod.factors());
  container::svector<Index> tidxs{};
  double const acc = 1.0;
  double const dp = opt::detail::peak_cost_batched(net, tidxs, idxsz, is_batch,
                                                   bts, novol, acc);
  double const oracle = opt::detail::reconstructed_batched_peak(
      net, tidxs, idxsz, is_batch, bts, novol, acc);
  // The correctness gate: the DP's max/+ recurrence for the accumulation
  // charge must agree with the independent memory-simulation oracle even
  // when more than one batchable index is sliced along the chosen tree.
  CHECK(dp == Catch::Approx(oracle));

  // Confirm the identity is not vacuous: the chosen 2-mode schedule must
  // actually be lower-peak than slicing either mode alone, proving the DP
  // engaged (and nested) both K and mu-tilde rather than degenerating to a
  // single-mode schedule.
  double const dp_K_only = opt::detail::peak_cost_batched(
      net, tidxs, idxsz, is_batch_K, bts, novol, acc);
  double const dp_mu_only = opt::detail::peak_cost_batched(
      net, tidxs, idxsz, is_batch_mu, bts, novol, acc);
  CHECK(dp < dp_K_only);
  CHECK(dp < dp_mu_only);

  // Resident-scan (order_aware_recompute): with the flag ON a batching node's
  // accumulator is charged as resident across its loop, so the modeled peak
  // RISES -- and the DP recurrence and the independent memory-simulation oracle
  // must still agree (the res term is mirrored in both). Flag OFF is
  // byte-identical (the CHECK(dp == oracle) above, default flag=false).
  double const dp_oar = opt::detail::peak_cost_batched(
      net, tidxs, idxsz, is_batch, bts, novol, acc, /*order_aware=*/true);
  double const oracle_oar = opt::detail::reconstructed_batched_peak(
      net, tidxs, idxsz, is_batch, bts, novol, acc, /*order_aware=*/true);
  CHECK(dp_oar == Catch::Approx(oracle_oar));  // parity holds with res on
  // res is monotone -- it can only raise the peak (adds to staged terms of a
  // max), never lower it. On THIS network the contraction moment dominates so
  // it is inert; the strict rise is exercised in [resident-scan] below.
  CHECK(dp_oar >= dp);
}

// Resident-scan demonstrator: a network where a batching node's accumulator,
// held across its loop while a child subtree evaluates, is the peak-setting
// co-residency -- so turning order_aware_recompute ON strictly RAISES the
// modeled peak, and the DP and oracle agree on the raised value.
TEST_CASE("resident-scan raises the batched peak by the accumulator",
          "[optimize][resident-scan]") {
  using namespace sequant;
  auto ctx_clone = get_default_context().clone();
  ctx_clone.mutable_index_space_registry()->add(L"F", IndexSpace::Type{0b10000},
                                                4ul);
  auto ctx_resetter = set_scoped_default_context(std::move(ctx_clone));
  auto idxsz = [](Index const& ix) { return ix.space().approximate_size(); };
  auto is_batchable = [](Index const& ix) {
    return ix.space().base_key() == L"F";
  };
  std::function<std::size_t(Index const&)> bts = [](Index const&) {
    return std::size_t{1};
  };
  std::function<bool(Tensor const&)> novol = {};
  std::vector<ExprPtr> ts;
  for (auto str : {L"g{a4;F1}", L"h{F1;a1}", L"s{a1;a2}", L"t{a2;a3}"})
    ts.push_back(deserialize(str, {.def_perm_symm = Symmetry::Nonsymm}));
  TensorNetwork net{ts};
  container::svector<Index> tidxs{};
  double const acc = 0.0;
  double const off = opt::detail::peak_cost_batched(
      net, tidxs, idxsz, is_batchable, bts, novol, acc);
  double const on = opt::detail::peak_cost_batched(
      net, tidxs, idxsz, is_batchable, bts, novol, acc, /*order_aware=*/true);
  double const oracle_on = opt::detail::reconstructed_batched_peak(
      net, tidxs, idxsz, is_batchable, bts, novol, acc, /*order_aware=*/true);
  CHECK(on > off);                        // the resident term bites here
  CHECK(on == Catch::Approx(oracle_on));  // DP == independent oracle, res on
}

// Ordered-key gate: on the Carr != 0 network, {s,t} carries the aux F_2 but not
// F_1. The set-keyed cell (flag off) can only charge F_1's recompute (rf=4 =>
// 1600) because it cannot express "F_1 inner". With order_aware_recompute on,
// the DP also holds the [F_2 outer, F_1 inner] ordered cell, where {s,t} hoists
// above F_1's loop and is charged nothing (rf=1 => 400). Read directly off the
// cells, minimizing over all orderings with the same sliced-SET {F_1,F_2}.
TEST_CASE("ordered key prices the hoistable order (Carr != 0)",
          "[optimize][ordered-key]") {
  using namespace sequant;
  auto ctx_clone = get_default_context().clone();
  ctx_clone.mutable_index_space_registry()->add(L"F", IndexSpace::Type{0b10000},
                                                4ul);
  auto ctx_resetter = set_scoped_default_context(std::move(ctx_clone));
  auto idxsz = [](Index const& ix) { return ix.space().approximate_size(); };
  auto is_batchable = [](Index const& ix) {
    return ix.space().base_key() == L"F";
  };
  auto batch_fn = [](Index const&) -> std::size_t { return 1; };
  std::vector<ExprPtr> ts;
  for (auto str : {L"g{a4;F1}", L"h{F1;F2}", L"s{F2;a2}", L"t{a2;a3}"})
    ts.push_back(deserialize(str, {.def_perm_symm = Symmetry::Nonsymm}));
  TensorNetwork net{ts};
  container::svector<Index> targets;
  std::size_t const st_set = 0b1100;  // {s,t}
  std::size_t const both_F = 0b11;    // union {F_1, F_2}

  auto min_flops_over_union = [&](bool oar) {
    opt::detail::PeakBatchedModel model{idxsz, batch_fn, {}};
    model.is_batchable_contracted_index = is_batchable;
    model.charge_batch_recompute = true;
    model.order_aware_recompute = oar;
    auto ctx = model.build_context(net, targets);
    auto st = opt::detail::solve_single_term(model, net, targets, ctx);
    double m = std::numeric_limits<double>::max();
    for (std::size_t id = 0; id < ctx.nCells; ++id)
      if (ctx.cell_union(id) == both_F)
        for (auto const& fp : st[st_set][id]) m = std::min(m, fp.flops);
    return m;
  };

  CHECK(min_flops_over_union(false) == 1600.0);  // set-keyed: over-charged
  CHECK(min_flops_over_union(true) == 400.0);  // ordered: hoistable order found
}

// S3.2: the ordered-cell enumeration must exclude EXTERNAL batchable modes
// (open on the root, contracted nowhere) -- only contracted-only modes are
// nestable, so an external mode must never appear in any cell's union. Network
// g{a1;F1} h{F1;F2} with empty targets: a1 (virtual space, not batchable, not
// occupied so batchable_mode_list's occ-external pass does not admit it) and
// F2 (space "F", batchable) each occur in exactly one tensor, so both are open
// at the root; F1 is shared by both tensors, so it is contracted. F2 is
// therefore EXTERNAL and F1 is the sole contracted batchable mode. Before the
// fix, build_cells enumerates ordered sequences over ALL batchable modes
// (including F2), so some cell's union carries F2's bit -- this CHECK fails
// today.
TEST_CASE("ordered cells exclude external modes", "[optimize][ext-place]") {
  using namespace sequant;
  auto ctx_clone = get_default_context().clone();
  ctx_clone.mutable_index_space_registry()->add(L"F", IndexSpace::Type{0b10000},
                                                4ul);
  auto ctx_resetter = set_scoped_default_context(std::move(ctx_clone));
  auto idxsz = [](Index const& ix) { return ix.space().approximate_size(); };
  auto is_batchable = [](Index const& ix) {
    return ix.space().base_key() == L"F";
  };
  auto batch_fn = [](Index const&) -> std::size_t { return 1; };
  std::vector<ExprPtr> ts;
  for (auto s : {L"g{a1;F1}", L"h{F1;F2}"})
    ts.push_back(deserialize(s, {.def_perm_symm = Symmetry::Nonsymm}));
  TensorNetwork net{ts};
  container::svector<Index> targets;  // empty: a1 and F2 both single-occurrence

  opt::detail::PeakBatchedModel model{idxsz, batch_fn, {}};

  model.is_batchable_contracted_index = is_batchable;
  model.is_batchable_external_index = is_batchable;  // external role (Task-4)
  model.charge_batch_recompute = true;
  model.order_aware_recompute = true;
  auto ctx = model.build_context(net, targets);

  REQUIRE(ctx.m == 2);
  std::size_t x_bit = ctx.m;
  for (std::size_t k = 0; k < ctx.m; ++k)
    if (model.is_external_mode(ctx, k)) x_bit = k;
  REQUIRE(x_bit < ctx.m);  // F2 is the external mode

  for (std::size_t id = 0; id < ctx.nCells; ++id)
    CHECK((ctx.cell_union(id) & (std::size_t{1} << x_bit)) == 0);
}

// S3.3: phase-2 node-level external-mode placement. Network
// g{F1;a1} h{a1;a2} t{a2;a3}: F1 (space "F", batchable, LARGE=100) sits on g
// only and stays open on the root, so it is a genuine EXTERNAL, is_batchable
// mode; a1/a2 are contracted (shared), a3 is a small external. Min-flops picks
// g*(h*t) (927 vs 1800 flops), whose (h*t) node does NOT carry F. Because that
// node escapes F's loop, the OLD root-level forest seed (seeded_forest_peak) is
// DECLINED as non-work-neutral, so before this task ON == OFF and the CHECKs
// below are RED. The phase-2 pass prices by PEAK only: at the over-budget root
// (which carries F) it slices F, shrinking every F-carrying node by
// block/extent (1/100). Expected root peak OFF = stage_form
// sz(g)+sz(ht)+sz(root) = 300+9+300 = 609 elems => 4872 B; ON with F sliced =
// 3+9+3 = 15 elems for the form stage but the F-free (h*t) subtree (peak 27
// elems) now dominates the staged max at 30 elems => 240 B.
TEST_CASE("phase-2 places an external mode on an over-budget node",
          "[optimize][ext-place]") {
  using namespace sequant;
  namespace o = sequant::opt::detail;
  auto ctx_clone = get_default_context().clone();
  ctx_clone.mutable_index_space_registry()->add(L"F", IndexSpace::Type{0b10000},
                                                100ul);
  ctx_clone.mutable_index_space_registry()
      ->retrieve_ptr(L"a")
      ->approximate_size(3ul);
  auto ctx_resetter = set_scoped_default_context(std::move(ctx_clone));
  auto idxsz = [](Index const& ix) { return ix.space().approximate_size(); };
  auto is_batchable = [](Index const& ix) {
    return ix.space().base_key() == L"F";
  };
  auto batch_fn = [](Index const&) -> std::size_t { return 1; };
  std::vector<ExprPtr> ts;
  for (auto s : {L"g{F1;a1}", L"h{a1;a2}", L"t{a2;a3}"})
    ts.push_back(deserialize(s, {.def_perm_symm = Symmetry::Nonsymm}));
  TensorNetwork net{ts};
  container::svector<Index> targets;  // empty: F1 and a3 single-occurrence

  // Run reconstruct_batched_modes and report (root peak bytes, #External
  // stamps) with batch_spectator_indices ON vs OFF; everything else identical.
  auto run = [&](bool spectator) {
    o::PeakBatchedModel model{idxsz, batch_fn, {}};
    model.is_batchable_contracted_index = is_batchable;
    // F is external-only here; admit it in the external role explicitly so it
    // survives the fallback removal (Task 4). Byte-identical to the fallback.
    model.is_batchable_external_index = is_batchable;
    model.charge_batch_recompute = true;
    model.order_aware_recompute = true;
    model.perf_first = true;
    model.batch_spectator_indices = spectator;
    model.node_level_placement = spectator;  // node-level emit (decoupled)
    // Finite budget strictly between the F-sliced (240 B) and unsliced (4872 B)
    // root peaks => the root is over budget and slicing F brings it under.
    model.peak_threshold = 1000.0;
    auto ctx = model.build_context(net, targets);
    REQUIRE(ctx.m == 1);                      // only F is batchable
    REQUIRE(model.is_external_mode(ctx, 0));  // and it is external
    auto st = o::solve_single_term(model, net, targets, ctx);
    double peak_bytes = 0.0;
    auto emitted =
        model.reconstruct_batched_modes(ctx, st, net, targets, &peak_bytes);
    std::size_t ext_stamps = 0;
    for (auto const& modes : emitted.second)
      for (auto const& [ix, knd] : modes.axes)
        if (knd == BatchModeType::External) ++ext_stamps;
    return std::pair<double, std::size_t>{peak_bytes, ext_stamps};
  };

  auto const [peak_off, ext_off] = run(false);
  auto const [peak_on, ext_on] = run(true);

  // OFF: pass does not fire => unsliced footprint, no External stamps.
  CHECK(peak_off == Catch::Approx(4872.0));
  CHECK(ext_off == 0u);
  // ON: F placed on the over-budget root => reported peak drops; F stamped.
  CHECK(peak_on == Catch::Approx(240.0));
  CHECK(ext_on >= 1u);
  // The point of S3.3: node-level external placement fired and cut the peak.
  CHECK(peak_on < peak_off);
}

// Task 3.3: binarize() must stamp EvalExpr::batched_here() from the optimizer's
// per-node sliced-sets (OptimizeOptions::term_batch_axes ->
// BinarizationOptions::node_batch_axes), and the two post-orders (the
// optimizer's DP reconstruction and binarize's Product recursion) must line
// up exactly -- this round-trips a real optimize() -> binarize() call and
// checks that (a) some node actually got annotated and (b) the annotated
// mode is the aux index the batch policy was forced to slice.
TEST_CASE("binarize stamps per-node batch modes from optimize()",
          "[optimize][annotate]") {
  using namespace sequant;
  auto ctx_resetter = set_scoped_default_context(get_default_context().clone());
  auto reg = get_default_context().mutable_index_space_registry();
  mbpt::add_df_spaces(reg);
  for (auto&& [k, v] :
       std::initializer_list<std::pair<std::wstring_view, size_t>>{
           {L"i", 30}, {L"a", 30}, {L"Κ", 500}}) {
    reg->retrieve_ptr(k)->approximate_size(v);
  }
  auto aux = reg->retrieve(L"Κ");
  auto idxsz = [](Index const& ix) -> std::size_t {
    return ix.nonnull() ? ix.space().approximate_size() : std::size_t{1};
  };
  auto is_batch = [aux](Index const& ix) { return ix.space() == aux; };
  std::function<std::size_t(Index const&)> bts = [](Index const&) {
    return std::size_t{20};
  };

  // 4-tensor network with Κ_1 shared between the two g's (the only batchable
  // contraction) so a forced-batching run must slice Κ at that node.
  auto expr =
      deserialize(L"g{a_1;i_1;Κ_1} g{a_2;i_2;Κ_1} f{i_1;i_3} f{i_2;i_4}");

  auto axes_map = std::make_shared<std::unordered_map<
      Expr const*, container::vector<NodeBatchAnnotation>>>();

  OptimizeOptions opts;
  opts.objective_function = ObjectiveFunction::DensePeakSizeBatched;
  opts.idx_to_extent = idxsz;
  opts.batch_policy.is_batchable_contracted_index = is_batch;
  opts.batch_policy.batch_target_size = bts;
  opts.batch_policy.peak_threshold = 1.0;  // tiny budget => force batching
  opts.term_batch_axes = axes_map;

  auto optimized = optimize(expr, opts);
  REQUIRE(optimized);

  // optimize() on a bare Product returns opt_pure_product's result directly
  // (no enclosing Sum), so `optimized` itself is the key term_batch_axes was
  // recorded under.
  auto it = axes_map->find(optimized.get());
  REQUIRE(it != axes_map->end());
  auto const& node_axes = it->second;
  REQUIRE(!node_axes.empty());

  BinarizationOptions bopts;
  bopts.node_batch_axes = node_axes;

  SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
  // binarize(ExprPtr) is deprecated in favor of binarize(ResultExpr); using it
  // here for ordering only -- the positional head layout it warns about does
  // not matter for this internal-node annotation check.
  auto node = binarize(optimized, {}, bopts);
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_END

  bool any_annotated = false;
  bool aux_found = false;
  node.visit([&](auto const& n) {
    if (n->batched_here().empty()) return;
    any_annotated = true;
    for (auto const& entry : n->batched_here())
      if (entry.first.space() == aux) aux_found = true;
  });
  // The essential assertion: the round-trip actually annotated a node. If the
  // optimizer's and binarize's post-orders had diverged, either this would be
  // false (silently under-annotated) or the SEQUANT_ASSERT inside binarize()
  // would have already fired above.
  CHECK(any_annotated);
  CHECK(aux_found);
}

// Task 3/4: reconstruct_batched_modes must emit BatchModeType::External entries
// for a genuine external (external, never-contracted) mode at every node whose
// subset carries it, gated by BOTH BatchPolicy::batch_spectator_indices AND
// the term-level emit_external gate (perf-first objective selected AND the
// selected root's unseeded byte peak exceeds peak_threshold). i_1 lives only
// on the first tensor of a 4-tensor chain -- a genuine external
// that is never contracted anywhere, so is_external_mode(i_1) holds
// regardless of the chosen factorization -- while Kappa_1 (the DF aux shared
// by the first two tensors) is a genuine CONTRACTED mode, never an external.
// Distinguishes the two kinds even though both get admitted into
// ctx.batchable_modes.
TEST_CASE("reconstruct_batched_modes_emits_external_per_node",
          "[optimize][batch]") {
  using namespace sequant;
  auto ctx_resetter = set_scoped_default_context(get_default_context().clone());
  auto reg = get_default_context().mutable_index_space_registry();
  mbpt::add_df_spaces(reg);
  for (auto&& [k, v] :
       std::initializer_list<std::pair<std::wstring_view, size_t>>{
           {L"i", 8}, {L"a", 8}, {L"Κ", 400}})
    reg->retrieve_ptr(k)->approximate_size(v);
  auto aux_space = reg->retrieve(L"Κ");
  auto occ_space = reg->retrieve(L"i");

  auto idxsz = [](Index const& ix) -> std::size_t {
    return ix.nonnull() ? ix.space().approximate_size() : std::size_t{1};
  };

  // Chain network T1-T2-T3-T4: T1--(a_1,Kappa_1)--T2--(a_2)--T3--(a_4)--T4.
  // i_1 (on T1 only) and a_5 (on T4 only) are the root-external indices.
  auto expr =
      deserialize(L"g{i_1;a_1;Κ_1} g{a_1;a_2;Κ_1} f{a_2;a_4} h{a_4;a_5}");

  auto axes_map = std::make_shared<std::unordered_map<
      Expr const*, container::vector<NodeBatchAnnotation>>>();

  OptimizeOptions opts;
  opts.objective_function = ObjectiveFunction::DenseTimeSpaceBatched;
  opts.reorder = ReorderSum::NoReorder;
  opts.idx_to_extent = idxsz;
  opts.batch_policy.is_batchable_contracted_index =
      [aux_space](Index const& ix) { return ix.space() == aux_space; };
  opts.batch_policy.batch_target_size = [](Index const&) -> std::size_t {
    return 20;
  };
  // The external mode is admitted in the EXTERNAL role: this space is batchable
  // where it is open on the term root, not where it is contracted. (Previously
  // such modes were auto-admitted by an occupancy-specific special case in
  // batchable_mode_list; batchability is now role-based and caller-declared.)
  opts.batch_policy.is_batchable_external_index = [](Index const& ix) {
    return ix.space().base_key() == L"i";
  };
  opts.batch_policy.batch_spectator_indices = true;
  // This case characterizes the LEGACY emit_external regime -- root-level
  // forest seed via the `else if (emit_external)` branch, NOT node-level
  // placement (see the Assertion 1b comment below). Pin order_aware_recompute
  // OFF explicitly (it is also the BatchPolicy default) so this legacy regime
  // is pinned regardless of any future default change.
  opts.batch_policy.order_aware_recompute = false;
  opts.term_batch_axes = axes_map;

  // Determine this term's actual perf-first root peak (bytes) directly, via
  // the same PeakBatchedModel construction single_term.hpp's
  // DenseTimeSpaceBatched arm uses (idxsz, is_batchable_index,
  // batch_target_size, perf_first=true), so peak_threshold values below and
  // above it can be chosen without hardcoding a magic number tied to the
  // extents above. peak_threshold is NOT a factorization gate under
  // perf_first (see PeakBatchedModel::select_root), so this value does not
  // depend on which threshold the model below is eventually given.
  namespace o = sequant::opt::detail;
  container::svector<ExprPtr> tensors;
  for (auto const& f : expr->as<Product>().factors())
    if (f->is<Tensor>()) tensors.push_back(f);
  TensorNetwork const net{tensors};
  o::PeakBatchedModel<decltype(idxsz)> peek_model{
      idxsz, opts.batch_policy.batch_target_size, /*is_volatile_leaf=*/{}};
  peek_model.is_batchable_contracted_index =
      opts.batch_policy.is_batchable_contracted_index;
  peek_model.perf_first = true;
  // Same role split as the optimize() path, so the peeked root peak (and hence
  // the threshold derived from it) is computed over the same mode set.
  peek_model.is_batchable_external_index =
      opts.batch_policy.is_batchable_external_index;
  container::svector<Index> const tidxs{};
  auto pctx = peek_model.build_context(net, tidxs);
  auto pst = o::solve_single_term(peek_model, net, tidxs, pctx);
  std::size_t const proot = (std::size_t{1} << pctx.nt) - 1;
  int const pbest = peek_model.select_root(pctx, pst);
  REQUIRE(pbest >= 0);
  double const root_peak_bytes =
      pst[proot][0][pbest].peak * peek_model.numeric_size;
  REQUIRE(root_peak_bytes > 0.0);

  // Assertion 1: peak_threshold strictly BELOW the term's actual peak =>
  // gate (iv) holds => External entries ARE emitted. (For this small chain no
  // batched schedule fits root_peak/2 either -- the peak shrinks less than 2x
  // under Kappa slicing -- so select_root's perf-first ceiling falls back to
  // the unbatched min-flops schedule here, same as before the ceiling existed;
  // the ceiling ENGAGING is covered by select_root_perf_first_ceiling below.)
  opts.batch_policy.peak_threshold = root_peak_bytes / 2.0;

  auto optimized = optimize(expr, opts);
  REQUIRE(optimized);

  auto it = axes_map->find(optimized.get());
  REQUIRE(it != axes_map->end());
  auto const& node_axes = it->second;
  REQUIRE(!node_axes.empty());

  bool node_with_external = false;
  bool node_without_external = false;
  for (auto const& modes : node_axes) {
    bool has_i1_external = false;
    for (auto const& entry : modes.axes) {
      // Kappa (a genuine contracted DF aux, never an external) must never be
      // tagged External.
      if (entry.first.space() == aux_space)
        CHECK(entry.second == BatchModeType::Contracted);
      if (entry.first.space() == occ_space &&
          entry.second == BatchModeType::External)
        has_i1_external = true;
    }
    if (has_i1_external)
      node_with_external = true;
    else
      node_without_external = true;
  }
  CHECK(node_with_external);
  CHECK(node_without_external);

  // Assertion 1b (review follow-up to D1 fix (ii), commit 601fe3205): this run
  // is the LEGACY `emit_external` regime -- batch_spectator_indices on,
  // order_aware_recompute at its default (off), over budget -- where
  // emit_mask_n comes from the `else if (emit_external)` branch
  // (chosen_seed_mask & ctx.open_modes[n]), NOT node-level placement. The
  // build lambda in cost_model.hpp pushes External entries before Contracted
  // ones regardless of which branch fed emit_mask_n, so the T1-T2 node here --
  // which co-carries the adopted external i_1 AND the contracted Κ_1 -- must
  // realize the External axis before the Contracted axis in `ann.axes`. This
  // characterizes the now-shipping legacy-regime behavior (it would have
  // FAILED before 601fe3205, when Contracted entries were pushed first).
  // Guard with found_co_carrying so the compare cannot pass vacuously.
  bool found_co_carrying = false;
  for (auto const& modes : node_axes) {
    int first_ext = -1, first_con = -1;
    for (int p = 0; p < static_cast<int>(modes.axes.size()); ++p) {
      auto const& entry = modes.axes[static_cast<std::size_t>(p)];
      if (entry.second == BatchModeType::External && first_ext < 0)
        first_ext = p;
      if (entry.second == BatchModeType::Contracted && first_con < 0)
        first_con = p;
    }
    if (first_ext >= 0 && first_con >= 0) {
      found_co_carrying = true;
      CHECK(first_ext < first_con);
    }
  }
  REQUIRE(found_co_carrying);

  // Assertion 2: peak_threshold = +infinity (the term never "needs" batching
  // under any budget) with batch_spectator_indices still true => the
  // unseeded-peak-over-threshold gate fails => ZERO External entries
  // anywhere. Before the Task 4 gate, reconstruct_batched_modes emitted
  // External regardless of peak_threshold, so this assertion fails without the
  // fix and passes with it.
  axes_map->clear();
  opts.batch_policy.peak_threshold = std::numeric_limits<double>::infinity();
  auto optimized_inf = optimize(expr, opts);
  REQUIRE(optimized_inf);
  auto it_inf = axes_map->find(optimized_inf.get());
  REQUIRE(it_inf != axes_map->end());
  for (auto const& modes : it_inf->second)
    for (auto const& entry : modes.axes)
      CHECK(entry.second != BatchModeType::External);

  // Assertion 3: with the external gate off (batch_spectator_indices=false)
  // -- threshold restored below the actual peak so only the flag is under
  // test -- the default optimize() path stays byte-identical: no External
  // entries anywhere (existing [optimize] tests are its witness).
  axes_map->clear();
  opts.batch_policy.peak_threshold = root_peak_bytes / 2.0;
  opts.batch_policy.batch_spectator_indices = false;
  auto optimized_off = optimize(expr, opts);
  REQUIRE(optimized_off);
  auto it_off = axes_map->find(optimized_off.get());
  REQUIRE(it_off != axes_map->end());
  for (auto const& modes : it_off->second)
    for (auto const& entry : modes.axes)
      CHECK(entry.second != BatchModeType::External);
}

// select_root under perf_first (DenseTimeSpaceBatched) treats peak_threshold
// as a CEILING: with a finite budget it picks the fewest-flops schedule whose
// peak fits, so a term is batched only when its unbatched (min-flops) schedule
// would exceed the budget -- flops are never traded for peak below the ceiling.
// At +inf every point is feasible (min-flops), and below even the min-peak
// point nothing fits so it falls back to GLOBAL min-flops (perf-first
// character, NOT peak-first's min-peak). Directly probes select_root on the DP
// frontier.
TEST_CASE("select_root_perf_first_ceiling", "[optimize][batch]") {
  using namespace sequant;
  namespace o = sequant::opt::detail;
  auto ctx_resetter = set_scoped_default_context(get_default_context().clone());
  auto reg = get_default_context().mutable_index_space_registry();
  mbpt::add_df_spaces(reg);
  mbpt::add_pao_spaces(reg, mbpt::Spin::any);
  mbpt::add_ao_spaces(reg, mbpt::Spin::any);
  // water-20-scale extents (matches the "threshold gates batching" test).
  for (auto&& [k, v] :
       std::initializer_list<std::pair<std::wstring_view, size_t>>{
           {L"i", 80}, {L"a", 12}, {L"μ̃", 860}, {L"Κ", 2360}})
    reg->retrieve_ptr(k)->approximate_size(v);
  auto aux_space = reg->retrieve(L"Κ");
  auto idxsz = [](Index const& ix) -> std::size_t {
    return ix.nonnull() ? ix.space().approximate_size() : std::size_t{1};
  };
  auto is_batch = [aux_space](Index const& ix) {
    return ix.space() == aux_space;
  };
  std::function<std::size_t(Index const&)> bts = [](Index const&) {
    return std::size_t{118};
  };

  // Two halves each carry g{.;.;Κ_1} sharing the SAME aux Κ_1, joined only
  // through the amplitude/PAO chains, so contracting each half down while Κ_1
  // stays open forms a large Κ_1-CARRYING intermediate -- the exact water-20
  // g(μ̃,μ̃,Κ) situation where slicing Κ lowers the peak (unlike a chain that
  // contracts Κ immediately). Flatten to a single 12-tensor product.
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
  container::svector<ExprPtr> tensors;
  for (auto const& f : flatp.factors())
    if (f->is<Tensor>()) tensors.push_back(f);
  TensorNetwork const net{tensors};

  // composite (CSV/PNO) domain size ~ small, like an OSV/PNO per-pair domain
  auto ip = [](Index const&, std::size_t) -> double { return 12.0; };
  o::PeakBatchedModel<decltype(idxsz)> model{idxsz, bts,
                                             /*is_volatile_leaf=*/{}, ip};
  model.is_batchable_contracted_index = is_batch;
  model.perf_first = true;
  container::svector<Index> const tidxs{};
  auto ctx = model.build_context(net, tidxs);
  auto st = o::solve_single_term(model, net, tidxs, ctx);
  std::size_t const root = (std::size_t{1} << ctx.nt) - 1;
  auto peak_bytes = [&](double p) { return p * model.numeric_size; };
  double const inf = std::numeric_limits<double>::infinity();

  // +inf budget => every point feasible => min-flops schedule (the unbatched
  // one, since batching only adds recompute flops).
  model.peak_threshold = inf;
  int const unb = model.select_root(ctx, st);
  REQUIRE(unb >= 0);
  double const p_unb = st[root][0][unb].peak;

  // Min achievable peak over the whole root frontier (the fully-sliced point).
  double p_min = std::numeric_limits<double>::max();
  for (auto const& fp : st[root][0]) p_min = std::min(p_min, fp.peak);
  // Batching must actually reduce the peak here, else there is no ceiling
  // window to exercise.
  REQUIRE(p_min < p_unb);

  // A budget strictly between p_min and p_unb: the unbatched min-flops point no
  // longer fits, so the ceiling picks a LOWER-peak batched point that does.
  double const mid_bytes = peak_bytes((p_min + p_unb) / 2.0);
  model.peak_threshold = mid_bytes;
  int const mid = model.select_root(ctx, st);
  REQUIRE(mid >= 0);
  CHECK(peak_bytes(st[root][0][mid].peak) <= mid_bytes);    // fits the ceiling
  CHECK(st[root][0][mid].peak < p_unb);                     // batched down
  CHECK(st[root][0][mid].flops >= st[root][0][unb].flops);  // >= min-flops

  // Below even p_min: nothing fits => perf-first best effort = GLOBAL min flops
  // (NOT peak-first's min-peak), i.e. the same flops as the +inf selection.
  model.peak_threshold = peak_bytes(p_min) / 2.0;
  int const none = model.select_root(ctx, st);
  REQUIRE(none >= 0);
  CHECK(st[root][0][none].flops == st[root][0][unb].flops);
}

TEST_CASE("contractible_adjacency", "[optimize][pruning]") {
  using namespace sequant;
  namespace o = sequant::opt::detail;
  auto ctx_resetter = set_scoped_default_context(get_default_context().clone());
  auto reg = get_default_context().mutable_index_space_registry();
  reg->retrieve_ptr(L"i")->approximate_size(10);
  reg->retrieve_ptr(L"a")->approximate_size(100);
  reg->retrieve_ptr(L"x")->approximate_size(4);
  auto parse = [](auto const& s) {
    return deserialize(s, {.def_perm_symm = Symmetry::Antisymm});
  };
  auto net_of = [](ExprPtr const& prod) {
    container::vector<ExprPtr> v;
    for (auto&& f : prod->as<Product>().factors())
      if (f->is<Tensor>()) v.push_back(f);
    return TensorNetwork{v};
  };
  auto edge_count = [](container::vector<std::size_t> const& adj) {
    std::size_t s = 0;
    for (auto m : adj) s += static_cast<std::size_t>(std::popcount(m));
    return s / 2;  // each undirected edge counted twice
  };

  // Chain: f-g share a1; g-h share a2. Targets {i1,i2}. 2 edges.
  auto chain = net_of(parse(L"f_{i1}^{a1} g_{a1}^{a2} h_{a2}^{i2}"));
  std::vector<Index> tgt2{Index{L"i_1"}, Index{L"i_2"}};
  auto adj_chain = o::contractible_adjacency(chain, tgt2);
  REQUIRE(adj_chain.size() == 3);
  CHECK(edge_count(adj_chain) == 2);

  // Hyperedge: p,q,r all carry summed a5 -> clique (3 edges).
  auto hyper = net_of(parse(L"p_{i1}^{a5} q_{i2}^{a5} r_{i3}^{a5}"));
  std::vector<Index> tgt3{Index{L"i_1"}, Index{L"i_2"}, Index{L"i_3"}};
  auto adj_hyper = o::contractible_adjacency(hyper, tgt3);
  CHECK(edge_count(adj_hyper) == 3);

  // External-only: two tensors share only the target index i1 -> no edges.
  auto spec = net_of(parse(L"u_{i1}^{a1} v_{i1}^{a2}"));
  std::vector<Index> tgt_spec{Index{L"i_1"}, Index{L"a_1"}, Index{L"a_2"}};
  auto adj_spec = o::contractible_adjacency(spec, tgt_spec);
  CHECK(edge_count(adj_spec) == 0);

  // Composite (CSV/PNO) indices: u and v carry DIFFERENT top-level composite
  // indices (distinct base labels a vs g) that happen to share the SAME
  // occupied protoindices {i1,i2}. contractible_adjacency only iterates
  // top-level bra/ket/aux indices, so a1<i1,i2> and g1<i1,i2> are two
  // distinct entries in the carrier map (different base label/space) even
  // though their protoindices coincide; those shared protoindices i1,i2 are
  // never themselves iterated as standalone indices, so no edge must form.
  // Targets are the plain externals i3,i4 only -- neither the composite
  // indices nor their protoindices are targets, so a zero edge count here is
  // not merely an artifact of exclusion-via-target.
  auto composite_diff = net_of(parse(L"u_{i3}^{a1<i1,i2>} v_{i4}^{g1<i1,i2>}"));
  std::vector<Index> tgt_composite_diff{Index{L"i_3"}, Index{L"i_4"}};
  auto adj_composite_diff =
      o::contractible_adjacency(composite_diff, tgt_composite_diff);
  CHECK(edge_count(adj_composite_diff) == 0);

  // Positive control: u and v now genuinely share the SAME top-level
  // composite index a1<i1,i2> (not merely the same protoindices), which
  // must produce an edge. This confirms the zero count above is because the
  // composites differ, not because composite indices can never be
  // adjacency carriers at all.
  auto composite_same = net_of(parse(L"u_{i5}^{a1<i1,i2>} v_{a1<i1,i2>}^{i6}"));
  std::vector<Index> tgt_composite_same{Index{L"i_5"}, Index{L"i_6"}};
  auto adj_composite_same =
      o::contractible_adjacency(composite_same, tgt_composite_same);
  CHECK(edge_count(adj_composite_same) == 1);
}

TEST_CASE("connected_subsets and outer_product_connectivity",
          "[optimize][pruning]") {
  using namespace sequant;
  namespace o = sequant::opt::detail;

  // Path 0-1-2 : adj[0]={1}, adj[1]={0,2}, adj[2]={1}.
  container::vector<std::size_t> path{0b010, 0b101, 0b010};
  auto c = o::connected_subsets(path, 3);
  REQUIRE(c.size() == 8);
  CHECK(c[0b001] == 1);  // singleton
  CHECK(c[0b011] == 1);  // {0,1} share edge
  CHECK(c[0b110] == 1);  // {1,2} share edge
  CHECK(c[0b101] == 0);  // {0,2} NOT adjacent -> disconnected
  CHECK(c[0b111] == 1);  // {0,1,2} connected via 1

  // Two disjoint components {0,1} and {2,3}: full set disconnected.
  container::vector<std::size_t> prod{0b0010, 0b0001, 0b1000, 0b0100};
  auto cp = o::connected_subsets(prod, 4);
  CHECK(cp[0b0011] == 1);  // {0,1}
  CHECK(cp[0b1100] == 1);  // {2,3}
  CHECK(cp[0b1111] == 0);  // full: disconnected product

  // outer_product_connectivity: env-disabled -> all ones.
  auto ctx_resetter = set_scoped_default_context(get_default_context().clone());
  auto reg = get_default_context().mutable_index_space_registry();
  reg->retrieve_ptr(L"i")->approximate_size(10);
  reg->retrieve_ptr(L"a")->approximate_size(100);
  reg->retrieve_ptr(L"x")->approximate_size(4);
  auto prod_expr = deserialize(L"f_{i1}^{a1} g_{a1}^{i2}",
                               {.def_perm_symm = Symmetry::Antisymm});
  container::vector<ExprPtr> v;
  for (auto&& f : prod_expr->as<Product>().factors())
    if (f->is<Tensor>()) v.push_back(f);
  TensorNetwork tn{v};
  std::vector<Index> tgt{Index{L"i_1"}, Index{L"i_2"}};
  setenv("SEQUANT_DISABLE_OUTER_PRODUCT_PRUNING", "1", 1);
  auto m_off = o::outer_product_connectivity(tn, tgt);
  unsetenv("SEQUANT_DISABLE_OUTER_PRODUCT_PRUNING");
  for (auto val : m_off) CHECK(val == 1);
}

TEST_CASE("outer-product pruning parity (pruned == unpruned)",
          "[optimize][pruning]") {
  using namespace sequant;
  auto ctx_resetter = set_scoped_default_context(get_default_context().clone());
  auto reg = get_default_context().mutable_index_space_registry();
  reg->retrieve_ptr(L"i")->approximate_size(10);
  reg->retrieve_ptr(L"a")->approximate_size(100);
  reg->retrieve_ptr(L"x")->approximate_size(4);

  auto opts_for = [&](ObjectiveFunction obj) {
    OptimizeOptions o;
    o.objective_function = obj;
    o.idx_to_extent = [](Index const& ix) -> std::size_t {
      return ix.nonnull() ? ix.space().approximate_size() : 1;
    };
    o.batch_policy.is_batchable_contracted_index = [](Index const&) {
      return false;
    };
    o.batch_policy.batch_target_size = [](Index const&) -> std::size_t {
      return 1;
    };
    return o;
  };

  std::vector<ObjectiveFunction> const objs{
      ObjectiveFunction::DenseFLOPs, ObjectiveFunction::DenseSize,
      ObjectiveFunction::DensePeakSize,
      ObjectiveFunction::DensePeakSizeBatched};

  std::vector<std::wstring> const terms{
      L"1/4 g_{i2,i3}^{a2,a3} t_{a2,a3}^{i2,i3} t_{a1}^{i1}",
      L"x_{i1,i2}^{a3,a4} y_{a1,a2}^{i1,i2} z_{a3,a4}^{a1,a2}",
      L"g_{i1,i2}^{a1,a2} t_{a1}^{i1} t_{a2}^{i2} t_{a3}^{i3}",
      // hyperedge-flavored: three tensors sharing a common summed index a5
      L"p_{i1}^{a5} q_{i2}^{a5} r_{i3}^{a5} s_{a5}^{i4}",
      // "star": g bridges both t's (whole network connected), but {t,t}
      // shares no index, so the {t,t} subset is a genuinely prunable outer
      // product. This is the fixture that actually exercises the pruning
      // skip -- see the non-vacuousness guard below.
      L"g_{i1,i2}^{a1,a2} t_{a1}^{i1} t_{a2}^{i2}",
  };

  // Non-vacuousness guard: the parity checks below only prove pruning is
  // loss-free if the pruned DP actually skips subsets the unpruned DP visits.
  // A fixture whose full network is disconnected falls back to an all-connected
  // mask (pruning is a no-op), so parity would hold even with a broken skip.
  // Assert that the "star" fixture yields a connected full network yet a
  // disconnected proper subset -- i.e. its mask really does prune something.
  {
    namespace o = sequant::opt::detail;
    unsetenv("SEQUANT_DISABLE_OUTER_PRODUCT_PRUNING");
    auto star = deserialize(L"g_{i1,i2}^{a1,a2} t_{a1}^{i1} t_{a2}^{i2}",
                            {.def_perm_symm = Symmetry::Antisymm});
    container::vector<ExprPtr> sv;
    for (auto&& f : star->as<Product>().factors())
      if (f->is<Tensor>()) sv.push_back(f);
    TensorNetwork star_tn{sv};
    std::vector<Index> const star_tgt{};  // fully contracted -> empty target
    auto star_mask = o::outer_product_connectivity(star_tn, star_tgt);
    REQUIRE(star_mask.back() == 1);  // full network connected (not fallback)
    bool has_pruned_subset = false;
    for (auto v : star_mask)
      if (v == 0) has_pruned_subset = true;
    REQUIRE(has_pruned_subset);  // at least one subset is genuinely pruned
  }

  auto run = [&](std::wstring const& term, ObjectiveFunction obj, bool prune) {
    if (prune)
      unsetenv("SEQUANT_DISABLE_OUTER_PRODUCT_PRUNING");
    else
      setenv("SEQUANT_DISABLE_OUTER_PRODUCT_PRUNING", "1", 1);
    auto expr = deserialize(term, {.def_perm_symm = Symmetry::Antisymm});
    auto out = optimize(expr, opts_for(obj));
    unsetenv("SEQUANT_DISABLE_OUTER_PRODUCT_PRUNING");
    REQUIRE(out);
    return to_latex(out);
  };

  for (std::size_t ti = 0; ti < terms.size(); ++ti)
    for (auto obj : objs) {
      auto pruned = run(terms[ti], obj, true);
      auto unpruned = run(terms[ti], obj, false);
      CAPTURE(ti, static_cast<int>(obj));
      CHECK(pruned == unpruned);
    }
}

TEST_CASE("prune_outer_products option controls pruning (default on)",
          "[optimize][pruning]") {
  using namespace sequant;
  // The requirement: pruning is user-controllable and ON by default.
  CHECK(OptimizeOptions{}.prune_outer_products == true);

  auto ctx_resetter = set_scoped_default_context(get_default_context().clone());
  auto reg = get_default_context().mutable_index_space_registry();
  reg->retrieve_ptr(L"i")->approximate_size(10);
  reg->retrieve_ptr(L"a")->approximate_size(100);
  reg->retrieve_ptr(L"x")->approximate_size(4);
  // "star" term: g bridges two mutually-disconnected t's, so the {t,t} subset
  // is a genuinely prunable outer product (the option's code path fires).
  auto expr = deserialize(L"g_{i1,i2}^{a1,a2} t_{a1}^{i1} t_{a2}^{i2}",
                          {.def_perm_symm = Symmetry::Antisymm});
  auto opts_for = [&](bool prune) {
    OptimizeOptions o;
    o.objective_function = ObjectiveFunction::DensePeakSize;
    o.idx_to_extent = [](Index const& ix) -> std::size_t {
      return ix.nonnull() ? ix.space().approximate_size() : 1;
    };
    o.prune_outer_products = prune;
    return o;
  };
  // Pruning is loss-free, so both option values give the same factorization.
  auto with_prune = to_latex(optimize(expr, opts_for(true)));
  auto no_prune = to_latex(optimize(expr, opts_for(false)));
  CHECK(with_prune == no_prune);
  // prune_outer_products == false must reproduce the env force-disable path.
  setenv("SEQUANT_DISABLE_OUTER_PRODUCT_PRUNING", "1", 1);
  auto env_disabled = to_latex(optimize(expr, opts_for(true)));
  unsetenv("SEQUANT_DISABLE_OUTER_PRODUCT_PRUNING");
  CHECK(no_prune == env_disabled);
}

// PeakBatchedModel's relax tie-break reads the precomputed per-subset atom
// tables (Context::fast_flops) instead of rebuilding index sets and calling
// flops_of per bipartition. fast_flops is loss-free by construction (atom IDs
// in FullLabelCompare order => same multiply order as inner_aware_volume). Gate
// that construction directly: for every connected bipartition of a composite
// (PNO/OSV) term and a plain term, with and without inner_pow engaged,
// fast_flops(lp, rp) must equal flops_of(idx[lp], idx[rp], idx[lp|rp]) exactly.
TEST_CASE("fast_flops equals flops_of over all bipartitions (parity)",
          "[optimize][fast-flops-parity]") {
  using namespace sequant;
  namespace o = opt::detail;
  auto ctx_resetter = set_scoped_default_context(get_default_context().clone());
  auto reg = get_default_context().mutable_index_space_registry();
  mbpt::add_df_spaces(reg);
  mbpt::add_pao_spaces(reg, mbpt::Spin::any);
  mbpt::add_ao_spaces(reg, mbpt::Spin::any);
  for (auto&& [k, v] :
       std::initializer_list<std::pair<std::wstring_view, size_t>>{
           {L"i", 10}, {L"a", 40}, {L"μ̃", 50}, {L"Κ", 90}})
    reg->retrieve_ptr(k)->approximate_size(v);
  auto idxsz = [](Index const& ix) -> std::size_t {
    return ix.nonnull() ? ix.space().approximate_size() : std::size_t{1};
  };
  auto is_batchable = [](Index const& ix) {
    return ix.space().base_key() == L"Κ";
  };
  auto batch = [](Index const&) -> std::size_t { return 30; };

  // fast_flops vs flops_of must agree for ANY inner_pow. Check two distinct
  // non-empty moment functions (empty inner_pow with composites is now
  // forbidden -- inner_aware_volume throws -- so it is not a valid parity
  // mode).
  std::function<double(Index const&, std::size_t)> const ip_on =
      [](Index const&, std::size_t) -> double { return 12.0; };
  std::function<double(Index const&, std::size_t)> const ip_off =
      [](Index const&, std::size_t k) -> double {
    return 5.0 * static_cast<double>(k);
  };

  bool composite_inner_checked = false;
  for (std::wstring const term :
       {std::wstring(L"g{μ̃1;μ̃2;Κ1} C{a1<i1>;μ̃1} C{μ̃2;a2<i1,i2>} t{a1<i1>;i1}"),
        std::wstring(L"g{i1;a1;Κ1} g{i2;a2;Κ1} t{a1;i1} t{a2;i2}")}) {
    for (auto const* ip : {&ip_on, &ip_off}) {
      auto prod = deserialize(term)->as<Product>();
      std::vector<ExprPtr> ts;
      for (auto&& f : prod.factors())
        if (f->is<Tensor>()) ts.push_back(f);
      TensorNetwork net{ts};
      container::svector<Index> const targets;
      o::PeakBatchedModel<decltype(idxsz)> model{idxsz, batch,
                                                 /*is_volatile_leaf=*/{}, *ip};
      model.is_batchable_contracted_index = is_batchable;
      auto ctx = model.build_context(net, targets);
      REQUIRE(ctx.use_fast_flops);
      auto const connected = o::outer_product_connectivity(net, targets);
      std::size_t const root = (std::size_t{1} << ts.size()) - 1;
      bool checked_bipart = false;
      for (std::size_t n = 1; n <= root; ++n) {
        if (std::popcount(n) < 2 || !connected[n]) continue;
        // enumerate every proper non-empty submask lp of n (rp = n ^ lp)
        for (std::size_t lp = (n - 1) & n; lp; lp = (lp - 1) & n) {
          std::size_t const rp = n ^ lp;
          if (rp == 0 || !connected[lp] || !connected[rp]) continue;
          checked_bipart = true;
          double const fast = ctx.fast_flops(lp, rp);
          double const slow =
              ctx.flops_of(ctx.idx[lp], ctx.idx[rp], ctx.idx[n]);
          CAPTURE(term, (ip == &ip_on), n, lp, rp, fast, slow);
          CHECK(fast == slow);  // bit-identical by construction
          if (!ctx.f_inner[lp].empty() || !ctx.f_inner[rp].empty())
            composite_inner_checked = true;
        }
      }
      REQUIRE(checked_bipart);  // non-vacuous
    }
  }
  // The composite term must actually exercise the inner (CSV/PNO) atom path.
  REQUIRE(composite_inner_checked);
}

TEST_CASE("outer-product pruning: large connected term optimizes quickly",
          "[optimize][pruning]") {
  using namespace sequant;
  auto ctx_resetter = set_scoped_default_context(get_default_context().clone());
  auto reg = get_default_context().mutable_index_space_registry();
  reg->retrieve_ptr(L"i")->approximate_size(10);
  reg->retrieve_ptr(L"a")->approximate_size(100);
  reg->retrieve_ptr(L"x")->approximate_size(4);
  // A connected chain: each adjacent pair shares one summed index, so the whole
  // term is one connected component. The pruned DP explores only connected
  // subsets and finishes fast; the unpruned 3^n enumeration is far slower.
  auto expr = deserialize(
      L"g_{a0,a1}^{a2,a3} t_{a2}^{a4} t_{a3}^{a5} t_{a4}^{a6} t_{a5}^{a7} "
      L"t_{a6}^{a8} t_{a7}^{a9} v_{a8,a9}^{a0,a1}",
      {.def_perm_symm = Symmetry::Antisymm});
  OptimizeOptions o;
  o.objective_function = ObjectiveFunction::DensePeakSize;
  o.idx_to_extent = [](Index const& ix) -> std::size_t {
    return ix.nonnull() ? ix.space().approximate_size() : 1;
  };
  // Pruning ON (default): must complete (the assertion is that it returns).
  auto out = optimize(expr, o);
  CHECK(out);
}

TEST_CASE("outer-product pruning: multi-component product falls back unpruned",
          "[optimize][pruning]") {
  using namespace sequant;
  namespace o = sequant::opt::detail;
  auto ctx_resetter = set_scoped_default_context(get_default_context().clone());
  auto reg = get_default_context().mutable_index_space_registry();
  reg->retrieve_ptr(L"i")->approximate_size(10);
  reg->retrieve_ptr(L"a")->approximate_size(100);
  reg->retrieve_ptr(L"x")->approximate_size(4);
  auto net_of = [](ExprPtr const& p) {
    container::vector<ExprPtr> v;
    for (auto&& f : p->as<Product>().factors())
      if (f->is<Tensor>()) v.push_back(f);
    return TensorNetwork{v};
  };

  // Two independent contractions sharing NO summed index: {f,g} over a1,
  // {p,q} over a2. The full adjacency graph has two components.
  auto prod = deserialize(L"f_{i1}^{a1} g_{a1}^{i2} p_{i3}^{a2} q_{a2}^{i4}",
                          {.def_perm_symm = Symmetry::Antisymm});
  auto tn = net_of(prod);
  std::vector<Index> tgt{Index{L"i_1"}, Index{L"i_2"}, Index{L"i_3"},
                         Index{L"i_4"}};
  auto mask = o::outer_product_connectivity(tn, tgt);
  for (auto v : mask) CHECK(v == 1);  // disconnected full net -> all-connected

  // And optimize() on the product term must match the env-disabled result.
  OptimizeOptions opt;
  opt.objective_function = ObjectiveFunction::DenseFLOPs;
  opt.idx_to_extent = [](Index const& ix) -> std::size_t {
    return ix.nonnull() ? ix.space().approximate_size() : 1;
  };
  auto with = to_latex(optimize(prod, opt));
  setenv("SEQUANT_DISABLE_OUTER_PRODUCT_PRUNING", "1", 1);
  auto without = to_latex(optimize(prod, opt));
  unsetenv("SEQUANT_DISABLE_OUTER_PRODUCT_PRUNING");
  CHECK(with == without);
}

// A2 PROBE (Phase A, order-aware multilevel batching). Measures what the DP
// charges TODAY for the gC/middle-gap shape, on a small hand-built network, so
// the RED assertion is written against ground truth rather than a predicted
// failure direction.
//
// Network (spec's example): R{i1,a2} = (g{i1;F1} * h{F1;a1}) * (s{a1;F2} *
// t{F2;a2}).  Batched set {F, i}. The right sub-intermediate carries F but NOT
// i, so it is invariant to any i-loop -- the I2/gC class.
//
// NOTE the premise conflict this probe exists to settle: plan Task A2 asserts
// "today's esc prices rf=1 (the phantom)" via the DP DROPPING a mode from B,
// but Task A1 established B is tree-faithful (B at a node is exactly the set
// of batched modes contracted at strict ancestors, by the C = B | aprime
// descent), so that mechanism does not exist. If rf=1 shows up here it is for
// a DIFFERENT reason -- i is EXTERNAL (free on the root, contracted nowhere),
// so it never enters any aprime and hence never enters B on the emit walk,
// which starts at B=0. Print the per-node charge and let the numbers decide.
TEST_CASE("loop-tree recompute charge prices the middle gap",
          "[.][loop-tree]") {
  using namespace sequant;
  // Scoped context clone: these TEST_CASEs are file-scope, so mutating the
  // shared registry would make them collide on the second reg->add(L"F").
  auto ctx_clone = get_default_context().clone();
  ctx_clone.mutable_index_space_registry()->add(L"F", IndexSpace::Type{0b10000},
                                                4ul);
  auto ctx_resetter = set_scoped_default_context(std::move(ctx_clone));

  auto idxsz = [](Index const& ix) -> std::size_t {
    return ix.space().approximate_size();
  };
  // batch BOTH the aux-like F and the occupied i
  auto is_batchable = [](Index const& ix) {
    auto const bk = ix.space().base_key();
    return bk == L"F" || bk == L"i";
  };
  auto batch_fn = [](Index const& ix) -> std::size_t {
    return ix.space().base_key() == L"F" ? std::size_t{1} : std::size_t{1};
  };

  std::vector<ExprPtr> ts;
  for (auto s : {L"g{i1;F1}", L"h{F1;a1}", L"s{a1;F2}", L"t{F2;a2}"})
    ts.push_back(deserialize(s, {.def_perm_symm = Symmetry::Nonsymm}));
  TensorNetwork net{ts};
  container::svector<Index> targets;

  auto aux = opt::detail::batchable_mode_list(net, is_batchable);
  std::wcerr << L"\n[loop-tree-probe] batchable modes (" << aux.size() << L"):";
  for (auto const& ix : aux) std::wcerr << L" " << ix.full_label();
  std::wcerr << L"\n";

  opt::detail::PeakBatchedModel model{idxsz, batch_fn,
                                      /*is_volatile_leaf=*/{}};
  model.is_batchable_contracted_index = is_batchable;
  // The batchable mode occurs EXTERNALLY here (open on the term root); admit it
  // in the external role explicitly. Byte-identical to the old external->
  // contracted fallback that Task 4 removed -- this hidden [.][loop-tree] test
  // was a fallback-reliant caller missed by Task 4 (its migration ran
  // [optimize], not [loop-tree]).
  model.is_batchable_external_index = is_batchable;
  model.charge_batch_recompute = true;
  auto ctx = model.build_context(net, targets);
  for (std::size_t k = 0; k < ctx.m; ++k)
    std::wcerr << L"[loop-tree-probe] mode "
               << ctx.batchable_modes[k].full_label() << L" nbatches="
               << ctx.nbatches[k] << L" external="
               << (model.is_external_mode(ctx, k) ? 1 : 0) << L"\n";

  auto st = opt::detail::solve_single_term(model, net, targets, ctx);
  // Emit walk with SEQUANT_DP_RECOMPUTE_DEBUG=1 set externally prints the
  // per-node carried/inside/escaped/rf triple that answers the question.
  // What the DP PRICES for the emitted schedule (emit walk starts at B=0)
  // vs what it would price with the external mode SEEDED into B at the root
  // (the enclosing loop the runtime actually executes). The seeded frontier is
  // already computed by the B-loop -- seeded_forest_peak reads it as a guard --
  // so both numbers are available today.
  std::size_t const root = (std::size_t{1} << ts.size()) - 1;
  auto min_flops_at = [&](std::size_t B) {
    double m = std::numeric_limits<double>::max();
    for (auto const& fp : st[root][B]) m = std::min(m, fp.flops);
    return m;
  };
  std::size_t i_bit = ctx.m;
  for (std::size_t k = 0; k < ctx.m; ++k)
    if (model.is_external_mode(ctx, k)) i_bit = k;
  REQUIRE(i_bit < ctx.m);
  double const unseeded = min_flops_at(0);
  double const seeded = min_flops_at(std::size_t{1} << i_bit);
  std::wcerr << L"[loop-tree-probe] min flops  B=0 (what emit prices) = "
             << unseeded << L"\n[loop-tree-probe] min flops  B={"
             << ctx.batchable_modes[i_bit].full_label()
             << L"} (runtime truth)   = " << seeded << L"\n"
             << L"[loop-tree-probe] ratio = " << (seeded / unseeded)
             << L"  (nbatches=" << ctx.nbatches[i_bit] << L")\n";

  auto const emitted = model.reconstruct_batched_modes(ctx, st, net, targets);
  std::wcerr << L"[loop-tree-probe] emitted " << emitted.second.size()
             << L" node stamps\n";
  for (std::size_t n = 0; n < emitted.second.size(); ++n) {
    std::wcerr << L"  node " << n << L":";
    for (auto const& [ix, knd] : emitted.second[n].axes)
      std::wcerr << L" " << ix.full_label()
                 << (knd == BatchModeType::External ? L":ext" : L":con");
    std::wcerr << L"\n";
  }
  // ---- RED (Phase A gate, Task A2) -------------------------------------
  // i_1 is EXTERNAL: free on the root, contracted nowhere. At runtime it is an
  // outer block-loop over the whole term, so every node that does NOT carry
  // i_1 is rebuilt nBatch(i_1) times. The DP must price that.
  //
  // It already CAN: st[root][{i_1}] is computed by the B-loop and carries the
  // charge (seeded_forest_peak reads exactly this frontier as a work-neutrality
  // guard). The defect is that the emit walk starts at B=0 -- build(root, 0,
  // best) -- so the schedule is PRICED as if the loop did not exist.
  //
  // Measured here: emit prices 720; the DP's own seeded frontier says 1200.
  // An order-aware cell key cannot close this gap -- ordering B changes how
  // modes ALREADY IN B are charged, and i_1 never enters B at all (B grows
  // only by aprime, and an external mode is in no contracted_here). The fix is
  // to seed external batched modes into B at the root. See
  // .superpowers/sdd/oamb-a0-note.md section 15.
  CHECK(seeded > unseeded);  // the DP's seeded frontier does carry the charge

  // Does the EMIT actually produce an unpriced i_1 loop? Turn on spectator
  // batching with a budget small enough to make the term over-budget, and see
  // whether an External stamp for i_1 is emitted. seeded_forest_peak declines a
  // seed whose slice is not work-neutral, which is exactly this case
  // (1200 != 720), so the expectation is NO stamp -- i.e. the runtime never
  // loops over i_1 and there is nothing unpriced.
  std::size_t n_ext_stamps = 0;
  {
    opt::detail::PeakBatchedModel m2{idxsz, batch_fn, {}};
    m2.is_batchable_contracted_index = is_batchable;
    m2.is_batchable_external_index = is_batchable;  // external role (Task-4)
    m2.charge_batch_recompute = true;
    m2.batch_spectator_indices = true;
    m2.perf_first = true;
    m2.peak_threshold = 1.0;  // force over_budget
    auto c2 = m2.build_context(net, targets);
    auto s2 = opt::detail::solve_single_term(m2, net, targets, c2);
    auto e2 = m2.reconstruct_batched_modes(c2, s2, net, targets);
    for (auto const& nd : e2.second)
      for (auto const& [ix, knd] : nd.axes)
        if (knd == BatchModeType::External) ++n_ext_stamps;
    std::wcerr << L"[loop-tree] spectator-on, over-budget: External stamps = "
               << n_ext_stamps << L"\n";
  }

  // RETRACTED GATE. An earlier revision asserted
  //   CHECK(unseeded >= seeded);   // "the emit under-prices the i_1 loop"
  // on the reasoning that i_1, being external, never enters B, so nodes
  // invariant to it are never charged. The first half is true (see the
  // inside_batch column above) but the conclusion is WRONG, and the check
  // immediately above is why: seeded_forest_peak DECLINES a seed whose slice
  // is not work-neutral, which is exactly this case (1200 != 720). No
  // External stamp is emitted, so the runtime never opens an i_1 loop and
  // there is nothing to under-price. The external path is safe BY
  // CONSTRUCTION -- the guard admits only seeds carried on every node, which
  // are legitimately rf=1.
  //
  // Phase A's real target is therefore the CONTRACTED case, where a mode DOES
  // enter B and the order-blind `esc` charge mis-prices it: per
  // oamb-a1-note.md section 1.3 today's charge is a systematic OVER-charge,
  // billing nBatch for escaped loops the node could hoist above for free. That
  // needs a network exhibiting the free-hoist (I2) shape and an assertion that
  // such a node is NOT charged. See oamb-a0-note.md section 16.
  CHECK(n_ext_stamps == 0u);  // no external loop is opened, so none is unpriced
}

// Phase A RED gate (Task A2, re-aimed per oamb-a0-note.md section 16).
//
// The free-hoist / I2 shape: a subset that carries NO batchable mode, sitting
// inside an enclosing loop over a CONTRACTED batchable mode. Such a node is
// invariant to that loop AND slicing the mode would not shrink it, so hoisting
// it above the loop is FREE -- unchanged footprint, computed once. The correct
// recompute factor is rf = 1.
//
// Today's charge (cost_model.hpp ~909-913) is order-blind: it bills
// nBatch(x) for EVERY x in `esc = B & ~open_modes[n]`, with no notion of where
// the node sits relative to those loops. So it over-charges this node by
// nBatch(F). A1 section 1.3 calls the Carr == 0 case "today's worst
// over-charge".
//
// Asserted directly on the DP cell st[n][B] rather than on an emitted schedule,
// so the test pins the CHARGE and does not depend on which tree the DP happens
// to select.
TEST_CASE("loop-tree charge must not bill a free hoist", "[.][loop-tree]") {
  using namespace sequant;
  // Scoped context clone: these TEST_CASEs are file-scope, so mutating the
  // shared registry would make them collide on the second reg->add(L"F").
  auto ctx_clone = get_default_context().clone();
  ctx_clone.mutable_index_space_registry()->add(L"F", IndexSpace::Type{0b10000},
                                                4ul);
  auto ctx_resetter = set_scoped_default_context(std::move(ctx_clone));

  auto idxsz = [](Index const& ix) -> std::size_t {
    return ix.space().approximate_size();
  };
  auto is_batchable = [](Index const& ix) {
    return ix.space().base_key() == L"F";
  };
  auto batch_fn = [](Index const&) -> std::size_t { return 1; };

  // g,h carry F_1 (contracted between them); s,t carry no F at all, so the
  // subset {s,t} is the free-hoist node.
  std::vector<ExprPtr> ts;
  // NB: no occupied index anywhere -- batchable_mode_list's external-occ pass
  // admits a pure-occupied index that is open on the root REGARDLESS of
  // is_batchable, which would silently make ctx.m == 2.
  for (auto str : {L"g{a4;F1}", L"h{F1;a1}", L"s{a1;a2}", L"t{a2;a3}"})
    ts.push_back(deserialize(str, {.def_perm_symm = Symmetry::Nonsymm}));
  TensorNetwork net{ts};
  container::svector<Index> targets;

  opt::detail::PeakBatchedModel model{idxsz, batch_fn, {}};

  model.is_batchable_contracted_index = is_batchable;
  model.charge_batch_recompute = true;
  model.order_aware_recompute = true;  // the fix under test
  auto ctx = model.build_context(net, targets);
  REQUIRE(ctx.m == 1u);  // exactly F_1
  REQUIRE(ctx.nbatches[0] > 1.0);

  auto st = opt::detail::solve_single_term(model, net, targets, ctx);

  std::size_t const n_st = 0b1100;                    // subset {s, t}
  REQUIRE(((ctx.open_modes[n_st] >> 0) & 1u) == 0u);  // carries no F_1

  auto min_flops = [&](std::size_t n, std::size_t B) {
    double m = std::numeric_limits<double>::max();
    for (auto const& fp : st[n][B]) m = std::min(m, fp.flops);
    return m;
  };
  double const outside = min_flops(n_st, 0);   // no enclosing loop
  double const inside = min_flops(n_st, 0b1);  // inside the F_1 loop
  std::wcerr << L"[loop-tree] free-hoist node {s,t}: flops outside F_1 loop = "
             << outside << L", inside = " << inside << L"  (nBatch(F_1) = "
             << ctx.nbatches[0] << L")\n";

  // {s,t} carries no batchable mode, so slicing F_1 cannot shrink it and it can
  // be built once above the F_1 loop at unchanged footprint. rf must be 1.
  CHECK(inside == outside);  // FAILS today: inside == nBatch(F_1) * outside
}

// Task A4: reconstruct_batched_modes must emit, per contraction node, its
// effective use count (effective_count) alongside the batch axes, on the
// order-aware path. The oracle is the note's formula evaluated straight from
// the Context at each node's (n, B): effective count = the product of
// nbatches over the escaped-outer set. We re-walk the SAME chosen
// back-pointer tree the emit walks and cross-check every emitted node, so the
// assertion is independent of which factorization the DP happens to select.
TEST_CASE("loop-tree emit: per-node effective_count", "[.][loop-tree]") {
  using namespace sequant;
  auto ctx_clone = get_default_context().clone();
  ctx_clone.mutable_index_space_registry()->add(L"F", IndexSpace::Type{0b10000},
                                                4ul);
  auto ctx_resetter = set_scoped_default_context(std::move(ctx_clone));

  auto idxsz = [](Index const& ix) -> std::size_t {
    return ix.space().approximate_size();
  };
  // Both F1 and F2 are batchable and contracted (nestable) => two ordered
  // enclosing modes.
  auto is_batchable = [](Index const& ix) {
    return ix.space().base_key() == L"F";
  };
  auto batch_fn = [](Index const&) -> std::size_t {
    return 1;  // batch tile 1 => nbatches == extent (4) > 1
  };

  // slot index (in ctx.batchable_modes) of a batchable Index, by full label.
  auto slot_of = [](auto const& ctx, Index const& ix) -> int {
    for (std::size_t k = 0; k < ctx.m; ++k)
      if (ctx.batchable_modes[k].full_label() == ix.full_label())
        return static_cast<int>(k);
    return -1;
  };

  // Re-walk the chosen back-pointer tree in the emit's left-first post-order,
  // recording (n, B) per contraction node so node_nb[j] pairs with the j-th
  // emitted annotation.
  auto walk_nb = [](auto const& ctx, auto const& st, int best) {
    container::vector<std::pair<std::size_t, std::size_t>> node_nb;
    std::size_t const root = (std::size_t{1} << ctx.nt) - 1;
    std::function<void(std::size_t, std::size_t, int)> go =
        [&](std::size_t n, std::size_t B, int idx) {
          if (std::popcount(n) == 1) return;
          auto const& r = st[n][B][idx];
          std::size_t const C = ctx.descend(B, r.aprime);
          std::size_t const fs = r.lp_first ? r.lp : r.rp;
          int const fi = r.lp_first ? r.lp_idx : r.rp_idx;
          std::size_t const ss = r.lp_first ? r.rp : r.lp;
          int const si = r.lp_first ? r.rp_idx : r.lp_idx;
          go(fs, C, fi);
          go(ss, C, si);
          node_nb.push_back({n, B});
        };
    go(root, 0, best);
    return node_nb;
  };

  // ---- Assertion 1: nested chain -- cross-check every emitted node ---------
  {
    std::vector<ExprPtr> ts;
    for (auto s : {L"g{a4;F1}", L"h{F1;F2}", L"s{F2;a2}", L"t{a2;a3}"})
      ts.push_back(deserialize(s, {.def_perm_symm = Symmetry::Nonsymm}));
    TensorNetwork net{ts};
    container::svector<Index> targets;

    opt::detail::PeakBatchedModel model{idxsz, batch_fn, {}};

    model.is_batchable_contracted_index = is_batchable;
    model.charge_batch_recompute = true;
    model.order_aware_recompute = true;
    auto ctx = model.build_context(net, targets);
    REQUIRE(ctx.m == 2u);

    auto st = opt::detail::solve_single_term(model, net, targets, ctx);
    int const best = model.select_root(ctx, st);
    REQUIRE(best >= 0);
    auto const node_nb = walk_nb(ctx, st, best);

    auto const emitted = model.reconstruct_batched_modes(ctx, st, net, targets);
    REQUIRE(emitted.second.size() == node_nb.size());
    REQUIRE(!emitted.second.empty());

    for (std::size_t j = 0; j < emitted.second.size(); ++j) {
      auto const& ann = emitted.second[j];
      auto const [n, B] = node_nb[j];

      // Expected effective_count = prod nbatches over the escaped-outer set.
      std::size_t const esc = ctx.escaped_outer(B, ctx.open_modes[n]);
      double rf = 1.0;
      for (std::size_t k = 0; k < ctx.m; ++k)
        if (esc & (std::size_t{1} << k)) rf *= ctx.nbatches[k];
      CHECK(ann.effective_count == static_cast<std::size_t>(std::llround(rf)));

      // Contracted axes emitted in ascending batchable-slot order (descend).
      int last = -1;
      for (auto const& [ix, knd] : ann.axes)
        if (knd == BatchModeType::Contracted) {
          int const slot = slot_of(ctx, ix);
          CHECK(slot > last);
          last = slot;
        }
    }

    // The term root carries no enclosing batched mode (its cell is the empty
    // sequence), so it has unit use count.
    CHECK(emitted.second.back().effective_count == 1u);

    // OFF path (order_aware_recompute=false) leaves effective_count at its
    // default -- byte-identical to before this task.
    opt::detail::PeakBatchedModel m_off{idxsz, batch_fn, {}};
    m_off.is_batchable_contracted_index = is_batchable;
    m_off.charge_batch_recompute = true;
    auto c_off = m_off.build_context(net, targets);
    auto s_off = opt::detail::solve_single_term(m_off, net, targets, c_off);
    auto const e_off =
        m_off.reconstruct_batched_modes(c_off, s_off, net, targets);
    for (auto const& ann : e_off.second) {
      CHECK(ann.effective_count == 1u);
    }
  }

  // ---- Assertion 2: a node contracting two batched modes emits them in ----
  // ---- ascending batched-slot order (matches descend). -------------------
  {
    std::vector<ExprPtr> ts;
    for (auto s : {L"x{a1;F1,F2}", L"y{F1,F2;a2}", L"z{a2;a3}"})
      ts.push_back(deserialize(s, {.def_perm_symm = Symmetry::Nonsymm}));
    TensorNetwork net{ts};
    container::svector<Index> targets;

    opt::detail::PeakBatchedModel model{idxsz, batch_fn, {}};

    model.is_batchable_contracted_index = is_batchable;
    model.charge_batch_recompute = true;
    model.order_aware_recompute = true;
    auto ctx = model.build_context(net, targets);
    REQUIRE(ctx.m == 2u);

    auto st = opt::detail::solve_single_term(model, net, targets, ctx);
    auto const emitted = model.reconstruct_batched_modes(ctx, st, net, targets);

    bool saw_two_mode = false;
    for (auto const& ann : emitted.second) {
      int last = -1, n_con = 0;
      for (auto const& [ix, knd] : ann.axes)
        if (knd == BatchModeType::Contracted) {
          int const slot = slot_of(ctx, ix);
          CHECK(slot > last);  // strictly ascending slot order
          last = slot;
          ++n_con;
        }
      if (n_con >= 2) saw_two_mode = true;
    }
    // x.y contracts BOTH F1 and F2 at one node.
    CHECK(saw_two_mode);
  }
}

// D1 fix (i) gate (external-placement propagation). When an external mode is
// adopted at an ancestor (the outermost over-budget node), the phase-2 place
// walk historically stamped `placed_at_node` ONLY at that ancestor -- the
// (tiny) residual root -- while the giant DESCENDANTS that carry the external
// mode FREE were left `(none)`. The runtime slices a node ONLY from that node's
// OWN `batched_here()` External stamp, so those descendants were never sliced
// and materialized/cached at full extent (the C60 4-occ giants). The fix
// propagates the adopted placement DOWN to every descendant carrying the mode.
//
// Same network idiom as the neighbouring [.][loop-tree] tests: R{i1,i2} chain
// g{i1;a1} * h{a1;a2} * k{a2;i2}. i1,i2 are EXTERNAL (free on the root,
// contracted nowhere); a1,a2 are CONTRACTED but non-batchable. With
// peak_threshold forcing over-budget, an external mode is adopted at the root
// and an interior carrier node (e.g. {g,h} carrying i1 free) is a descendant of
// that root carrying the external mode. RED pre-fix: the descendant's emitted
// `batched_here()` (axes) has NO External entry. GREEN post-fix: it carries an
// External entry for the adopted mode.
TEST_CASE(
    "loop-tree emit: external placement propagates to carrying descendants",
    "[.][loop-tree]") {
  using namespace sequant;
  auto ctx_clone = get_default_context().clone();
  auto reg = ctx_clone.mutable_index_space_registry();
  reg->retrieve_ptr(L"i")->approximate_size(10);  // occupied, external, batched
  reg->retrieve_ptr(L"a")->approximate_size(20);  // virtual, contracted, inert
  auto ctx_resetter = set_scoped_default_context(std::move(ctx_clone));

  auto idxsz = [](Index const& ix) -> std::size_t {
    return ix.space().approximate_size();
  };
  auto is_batchable = [](Index const& ix) {
    return ix.space().base_key() == L"i";
  };
  auto batch_fn = [](Index const&) -> std::size_t { return 1; };

  std::vector<ExprPtr> ts;
  for (auto s : {L"g{i1;a1}", L"h{a1;a2}", L"k{a2;i2}"})
    ts.push_back(deserialize(s, {.def_perm_symm = Symmetry::Nonsymm}));
  TensorNetwork net{ts};
  container::svector<Index> targets;

  opt::detail::PeakBatchedModel model{idxsz, batch_fn, {}};

  model.is_batchable_contracted_index = is_batchable;
  model.is_batchable_external_index = is_batchable;  // external role (Task-4)
  model.charge_batch_recompute = true;
  model.order_aware_recompute = true;    // order-aware cost model (selection)
  model.batch_spectator_indices = true;  // external batching
  model.node_level_placement = true;     // engage node-level placement (emit)
  model.peak_threshold = 1.0;            // force every node over budget

  auto ctx = model.build_context(net, targets);

  std::size_t ext_mask = 0;
  for (std::size_t k = 0; k < ctx.m; ++k)
    if (model.is_external_mode(ctx, k)) ext_mask |= (std::size_t{1} << k);
  REQUIRE(ext_mask != 0);

  auto st = opt::detail::solve_single_term(model, net, targets, ctx);
  int const best = model.select_root(ctx, st);
  REQUIRE(best >= 0);

  // Re-walk the chosen back-pointer tree in the emit's left-first post-order,
  // recording n so node_n[j] pairs with the j-th emitted annotation.
  std::size_t const root = (std::size_t{1} << ctx.nt) - 1;
  container::vector<std::size_t> node_n;
  std::function<void(std::size_t, std::size_t, int)> go =
      [&](std::size_t n, std::size_t B, int idx) {
        if (std::popcount(n) == 1) return;
        auto const& r = st[n][B][idx];
        std::size_t const C = ctx.descend(B, r.aprime);
        std::size_t const fs = r.lp_first ? r.lp : r.rp;
        int const fi = r.lp_first ? r.lp_idx : r.rp_idx;
        std::size_t const ss = r.lp_first ? r.rp : r.lp;
        int const si = r.lp_first ? r.rp_idx : r.lp_idx;
        go(fs, C, fi);
        go(ss, C, si);
        node_n.push_back(n);
      };
  go(root, 0, best);

  auto const emitted = model.reconstruct_batched_modes(ctx, st, net, targets);
  REQUIRE(emitted.second.size() == node_n.size());

  // Locate an interior (non-root) DESCENDANT that carries an external mode
  // FREE in its own result -- this is the giant-class node the runtime must
  // slice. It is a descendant of the root, which is where the external mode is
  // adopted (peak_threshold=1.0 forces the root over budget).
  int desc_j = -1;
  std::size_t desc_mode_mask = 0;
  for (std::size_t j = 0; j < node_n.size(); ++j) {
    std::size_t const n = node_n[j];
    if (n == root) continue;
    std::size_t const carried_ext = ctx.open_modes[n] & ext_mask;
    if (carried_ext) {
      desc_j = static_cast<int>(j);
      desc_mode_mask = carried_ext;
      break;
    }
  }
  REQUIRE(desc_j >= 0);  // the factorization must expose a descendant carrier

  // The adopted external mode must be stamped External on this descendant's
  // own annotation. RED pre-fix: no External entry (only the root was stamped);
  // GREEN post-fix: the placement propagated down to the carrier.
  auto const& desc = emitted.second[static_cast<std::size_t>(desc_j)];
  bool desc_has_ext = false;
  for (auto const& [ix, knd] : desc.axes) {
    if (knd != BatchModeType::External) continue;
    for (std::size_t k = 0; k < ctx.m; ++k)
      if ((desc_mode_mask & (std::size_t{1} << k)) &&
          ctx.batchable_modes[k].full_label() == ix.full_label())
        desc_has_ext = true;
  }
  CHECK(desc_has_ext);
}

// D1 fix (ii) gate (emit External BEFORE Contracted at co-carrying nodes). A
// node that co-carries a stamped External (occ) mode AND a Contracted (aux)
// mode it contracts at itself must realize the External loop OUTER of the
// Contracted one, else the runtime scatter widens the external mode to full
// extent per contracted block (the C60 co-carrying 4-occ giants). The emitted
// annotation's `axes` == the node's own realized-loop order == the runtime
// pick order, so "External outer" means the External entry appears BEFORE the
// Contracted entry in `ann.axes`.
//
// Network: g{i1;a1} * k{i2;a1}. Both `i` (external, free on the root R{i1,i2},
// contracted nowhere) and `a` (contracted at the single {g,k} node, appears on
// no target free) are made BATCHABLE with tile size 1 => nbatches == extent
// (> 1). The lone internal node {g,k} therefore co-carries: it contracts a1
// (a batched Contracted mode) AND carries i1,i2 free external, and with
// peak_threshold forcing over budget the external modes are adopted/stamped
// there. RED pre-fix: `axes` is [Contracted a1, External i...] (contracted
// pushed first). GREEN post-fix: [External i..., Contracted a1].
TEST_CASE(
    "loop-tree emit: External axis precedes Contracted at co-carrying node",
    "[.][loop-tree]") {
  using namespace sequant;
  auto ctx_clone = get_default_context().clone();
  auto reg = ctx_clone.mutable_index_space_registry();
  reg->retrieve_ptr(L"i")->approximate_size(10);  // occupied, external, batched
  reg->retrieve_ptr(L"a")->approximate_size(
      20);  // virtual, contracted, batched
  auto ctx_resetter = set_scoped_default_context(std::move(ctx_clone));

  auto idxsz = [](Index const& ix) -> std::size_t {
    return ix.space().approximate_size();
  };
  // batch BOTH the external occupied i and the contracted virtual a, so the
  // co-carrying node emits an External axis (i) and a Contracted axis (a).
  auto is_batchable = [](Index const& ix) {
    auto const key = ix.space().base_key();
    return key == L"i" || key == L"a";
  };
  auto batch_fn = [](Index const&) -> std::size_t {
    return 1;  // batch tile 1 => nbatches == extent (> 1 for both i and a)
  };

  std::vector<ExprPtr> ts;
  for (auto s : {L"g{i1;a1}", L"k{i2;a1}"})
    ts.push_back(deserialize(s, {.def_perm_symm = Symmetry::Nonsymm}));
  TensorNetwork net{ts};
  container::svector<Index> targets;

  opt::detail::PeakBatchedModel model{idxsz, batch_fn, {}};

  model.is_batchable_contracted_index = is_batchable;
  // i is external, a is contracted; both are admitted in both roles here, so
  // the role filter keeps i external and a contracted. Byte-identical to the
  // fallback but survives its removal (Task 4).
  model.is_batchable_external_index = is_batchable;
  model.charge_batch_recompute = true;
  model.order_aware_recompute = true;    // order-aware cost model (selection)
  model.batch_spectator_indices = true;  // external batching
  model.node_level_placement = true;     // engage node-level placement (emit)
  model.peak_threshold = 1.0;            // force every node over budget

  auto ctx = model.build_context(net, targets);

  auto st = opt::detail::solve_single_term(model, net, targets, ctx);
  int const best = model.select_root(ctx, st);
  REQUIRE(best >= 0);

  auto const emitted = model.reconstruct_batched_modes(ctx, st, net, targets);

  // Find an emitted annotation that co-carries BOTH an External and a
  // Contracted axis, and assert the FIRST External entry precedes the FIRST
  // Contracted entry. REQUIRE such a node exists (guard against a vacuous
  // pass: if the setup failed to produce a co-carrying node the invariant is
  // untested).
  bool found_co_carrying = false;
  for (auto const& ann : emitted.second) {
    int first_ext = -1, first_con = -1;
    for (int p = 0; p < static_cast<int>(ann.axes.size()); ++p) {
      auto const knd = ann.axes[static_cast<std::size_t>(p)].second;
      if (knd == BatchModeType::External && first_ext < 0) first_ext = p;
      if (knd == BatchModeType::Contracted && first_con < 0) first_con = p;
    }
    if (first_ext >= 0 && first_con >= 0) {
      found_co_carrying = true;
      // RED pre-fix: Contracted pushed first => first_con < first_ext.
      // GREEN post-fix: External pushed first => first_ext < first_con.
      CHECK(first_ext < first_con);
    }
  }
  REQUIRE(found_co_carrying);
}

// PROBE (not a gate): the resident-scan peak. When {s,t} is free-hoisted above
// the F_1 loop (order_aware_recompute), it is RESIDENT across that loop, so its
// footprint should appear in the peak of everything evaluated inside -- the
// top-node peak under B={F_1}. A3a adds no peak term, so today the hoist looks
// free on peak. Measure the cells before writing the resident-scan RED gate.
TEST_CASE("loop-tree probe: resident-scan peak of a hoisted node",
          "[.][loop-tree-peak]") {
  using namespace sequant;
  auto ctx_clone = get_default_context().clone();
  ctx_clone.mutable_index_space_registry()->add(L"F", IndexSpace::Type{0b10000},
                                                4ul);
  auto ctx_resetter = set_scoped_default_context(std::move(ctx_clone));

  auto idxsz = [](Index const& ix) -> std::size_t {
    return ix.space().approximate_size();
  };
  auto is_batchable = [](Index const& ix) {
    return ix.space().base_key() == L"F";
  };
  auto batch_fn = [](Index const&) -> std::size_t { return 1; };

  std::vector<ExprPtr> ts;
  for (auto str : {L"g{a4;F1}", L"h{F1;a1}", L"s{a1;a2}", L"t{a2;a3}"})
    ts.push_back(deserialize(str, {.def_perm_symm = Symmetry::Nonsymm}));
  TensorNetwork net{ts};
  container::svector<Index> targets;

  std::size_t const top = 0b1111, st_set = 0b1100;

  auto min_peak = [](auto const& frontier) {
    double m = std::numeric_limits<double>::max();
    for (auto const& fp : frontier) m = std::min(m, fp.peak);
    return m;
  };

  for (bool oar : {false, true}) {
    opt::detail::PeakBatchedModel model{idxsz, batch_fn, {}};
    model.is_batchable_contracted_index = is_batchable;
    model.charge_batch_recompute = true;
    model.order_aware_recompute = oar;
    auto ctx = model.build_context(net, targets);
    auto st = opt::detail::solve_single_term(model, net, targets, ctx);
    std::wcerr << L"[loop-tree-peak] order_aware=" << oar << L"  sz{s,t}="
               << ctx.sz(st_set, 0) << L"  top.peak[B=0]="
               << min_peak(st[top][0]) << L"  top.peak[B={F_1}]="
               << min_peak(st[top][0b1]) << L"\n";
  }
}

// PROBE (not a gate): the Carr != 0 order-dependent case. A node that carries
// one batched contracted mode (F_2) and not another (F_1) sits in a cell keyed
// by the SET {F_1,F_2}. The correct charge depends on the ORDER: if F_1's loop
// is OUTER to F_2's the node is legitimately recomputed per F_1 block, but if
// F_1's loop is INNER to F_2's the node can hoist above it for free. Both
// orders map to the same cell, so one answer must be wrong -- that is the
// representability defect. Measure what the cells actually hold before
// asserting anything.
TEST_CASE("loop-tree probe: order-dependent Carr != 0 cells",
          "[.][loop-tree-order]") {
  using namespace sequant;
  auto ctx_clone = get_default_context().clone();
  ctx_clone.mutable_index_space_registry()->add(L"F", IndexSpace::Type{0b10000},
                                                4ul);
  auto ctx_resetter = set_scoped_default_context(std::move(ctx_clone));

  auto idxsz = [](Index const& ix) -> std::size_t {
    return ix.space().approximate_size();
  };
  auto is_batchable = [](Index const& ix) {
    return ix.space().base_key() == L"F";
  };
  auto batch_fn = [](Index const&) -> std::size_t { return 1; };

  std::vector<ExprPtr> ts;
  for (auto str : {L"g{a4;F1}", L"h{F1;F2}", L"s{F2;a2}", L"t{a2;a3}"})
    ts.push_back(deserialize(str, {.def_perm_symm = Symmetry::Nonsymm}));
  TensorNetwork net{ts};
  container::svector<Index> targets;

  opt::detail::PeakBatchedModel model{idxsz, batch_fn, {}};

  model.is_batchable_contracted_index = is_batchable;
  model.charge_batch_recompute = true;
  auto ctx = model.build_context(net, targets);
  std::wcerr << L"\n[loop-tree-order] modes:";
  for (std::size_t k = 0; k < ctx.m; ++k)
    std::wcerr << L" " << ctx.batchable_modes[k].full_label() << L"(bit" << k
               << L",nB=" << ctx.nbatches[k] << L")";
  std::wcerr << L"\n";

  auto st = opt::detail::solve_single_term(model, net, targets, ctx);
  std::size_t const n_st = 0b1100;  // {s,t}
  std::wcerr << L"[loop-tree-order] subset {s,t} open modes bitmask = "
             << ctx.open_modes[n_st] << L"\n";
  auto min_flops = [&](std::size_t n, std::size_t B) {
    double m = std::numeric_limits<double>::max();
    for (auto const& fp : st[n][B]) m = std::min(m, fp.flops);
    return m;
  };
  for (std::size_t B = 0; B < ctx.nB; ++B)
    std::wcerr << L"[loop-tree-order]   st[{s,t}][B=" << B << L"] min flops = "
               << min_flops(n_st, B) << L"\n";
  SUCCEED("probe");
}

// PROBE (not a gate): resident-scan peak on the nested Carr != 0 network. Reads
// the top-node peak per B to see whether an enclosing accumulator (the F_1
// contraction's result, resident during the F_2 loop and vice versa) is
// counted.
TEST_CASE("loop-tree probe: resident-scan peak, nested",
          "[.][loop-tree-peak-nested]") {
  using namespace sequant;
  auto ctx_clone = get_default_context().clone();
  ctx_clone.mutable_index_space_registry()->add(
      L"F", IndexSpace::Type{0b10000},
      1000ul);  // large: dominates peak
  auto ctx_resetter = set_scoped_default_context(std::move(ctx_clone));

  auto idxsz = [](Index const& ix) { return ix.space().approximate_size(); };
  auto is_batchable = [](Index const& ix) {
    return ix.space().base_key() == L"F";
  };
  auto batch_fn = [](Index const&) -> std::size_t { return 1; };

  std::vector<ExprPtr> ts;
  for (auto str : {L"g{a4;F1}", L"h{F1;F2}", L"s{F2;a2}", L"t{a2;a3}"})
    ts.push_back(deserialize(str, {.def_perm_symm = Symmetry::Nonsymm}));
  TensorNetwork net{ts};
  container::svector<Index> targets;

  auto min_peak = [](auto const& frontier) {
    double m = std::numeric_limits<double>::max();
    for (auto const& fp : frontier) m = std::min(m, fp.peak);
    return m;
  };

  // batch F to tiles of 1 so nbatches == extent; large F so slicing it
  // dominates
  auto batch1 = [](Index const&) -> std::size_t { return 1; };
  for (bool oar : {false, true}) {
    opt::detail::PeakBatchedModel model{idxsz, batch1, {}};
    model.is_batchable_contracted_index = is_batchable;
    model.charge_batch_recompute = true;
    model.order_aware_recompute = oar;
    auto ctx = model.build_context(net, targets);
    auto st = opt::detail::solve_single_term(model, net, targets, ctx);
    std::size_t const top = (std::size_t{1} << ctx.nt) - 1;
    std::wcerr << L"[peak-nested] oar=" << oar << L"  m=" << ctx.m
               << L"  sz{s,t}[B=0]=" << ctx.sz(0b1100, 0) << L"  sz{g,h}[B=0]="
               << ctx.sz(0b0011, 0);
    for (std::size_t B = 0; B < ctx.nB; ++B)
      std::wcerr << L"  top.peak[B=" << B << L"]=" << min_peak(st[top][B]);
    std::wcerr << L"\n";
  }
}

// Task 1 [role-api]: the two batchability building-block predicates + their
// derived "any role" accessors. BatchPolicy / CostParams / PeakBatchedModel
// each expose is_batchable_contracted_index and is_batchable_external_index as
// settable fields. The "batchable in any role" query is DERIVED, never a
// settable field: BatchPolicy::is_batchable_index() returns the std::function
// union; PeakBatchedModel::is_batchable(ix) is the bool method == OR of the two
// building blocks.
TEST_CASE("batchability role-split building-block predicates",
          "[optimize][role-api]") {
  using namespace sequant;
  auto ctx_resetter = set_scoped_default_context(get_default_context().clone());
  auto reg = get_default_context().mutable_index_space_registry();
  {
    auto occ = reg->retrieve_ptr(L"i");
    auto uocc = reg->retrieve_ptr(L"a");
    REQUIRE(occ);
    REQUIRE(uocc);
  }
  Index const a{L"a_1"};  // unoccupied space "a"
  Index const i{L"i_1"};  // occupied space "i"

  auto is_a = [](Index const& ix) { return ix.space().base_key() == L"a"; };
  auto is_i = [](Index const& ix) { return ix.space().base_key() == L"i"; };

  SECTION("BatchPolicy::is_batchable_index() is the union of both roles") {
    BatchPolicy p;
    p.is_batchable_contracted_index = is_a;
    p.is_batchable_external_index = is_i;
    std::function<bool(Index const&)> const u = p.is_batchable_index();
    CHECK(u(a));
    CHECK(u(i));
    // contracted-role-only => the derived union == the contracted predicate.
    BatchPolicy p2;
    p2.is_batchable_contracted_index = is_a;
    std::function<bool(Index const&)> const u2 = p2.is_batchable_index();
    CHECK(u2(a));
    CHECK_FALSE(u2(i));
  }

  SECTION("PeakBatchedModel::is_batchable(ix) == OR of both building blocks") {
    auto idxsz = [](Index const& ix) -> std::size_t {
      return ix.nonnull() ? ix.space().approximate_size() : std::size_t{1};
    };
    auto batch_fn = [](Index const&) -> std::size_t { return 4; };
    opt::detail::PeakBatchedModel<decltype(idxsz)> model{idxsz, batch_fn, {}};
    model.is_batchable_contracted_index = is_a;
    model.is_batchable_external_index = is_i;
    CHECK(model.is_batchable(a));
    CHECK(model.is_batchable(i));
    // external-role-only mode is still batchable via the OR.
    opt::detail::PeakBatchedModel<decltype(idxsz)> m2{idxsz, batch_fn, {}};
    m2.is_batchable_external_index = is_i;
    CHECK(m2.is_batchable(i));
    CHECK_FALSE(m2.is_batchable(a));
  }

  SECTION("CostParams carries both role predicates + batch_target_size") {
    CostParams cost;
    cost.is_batchable_contracted_index = is_a;
    cost.is_batchable_external_index = is_i;
    cost.batch_target_size = [](Index const&) -> std::size_t { return 8; };
    CHECK(cost.is_batchable_contracted_index(a));
    CHECK(cost.is_batchable_external_index(i));
    CHECK(cost.batch_target_size(a) == 8u);
  }
}
