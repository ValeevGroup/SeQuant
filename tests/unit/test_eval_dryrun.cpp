// Task 1 (DryRun eval backend, Phase 1): infrastructure + DP-inspection
// mechanics for the eventual go/no-go probe on the C60 batched-schedule
// factorization.
//
// SCOPE NOTE (two course corrections; see task-1-report.md for the full
// history): the committed fixture `data/csv_ccsd_residual.txt` (mpqc's
// `csv_eqn_Rs.serialized`) is a PRE-transform residual: every two-electron
// integral `g` is still 4-index (`g{i,i;a,a}`-shaped), there is no auxiliary
// DF ("K"/Kappa) index anywhere in it, and the non-proto `a_NNNNN` indices are
// the GENERIC BASE VIRTUAL space, NOT the PAO ("mu~") space -- PAO and DF aux
// are introduced downstream by mpqc's CSV->base-integral transform and
// density-fitting refactorization, which run in mpqc's own pipeline against
// mpqc's extended LCAO space registry, not the standalone SeQuant mbpt
// registry this test has access to. So:
//   - the two K-based memsize anchors from the cluster trace (~1.87 GB /
//     ~186 GB) are NOT reproducible from this fixture, and this file does not
//     attempt them (see the self-consistency case below instead);
//   - a faithful PAO/K batch-axis GO/NO-GO VERDICT is NOT obtainable from this
//     pre-transform fixture either, since the factorizer that produced the
//     C60 trace ran on the POST-transform (mu~ + K) form. That verdict is
//     deferred to a follow-up once a post-transform fixture is available.
//
// What this file DOES deliver and verify, against the real residual:
//   1. the deserialize path works under the standalone mbpt context (mirrors
//      test_optimize.cpp's context setup), and reports the term count;
//   2. the DP-inspection MECHANICS work end-to-end: running optimize() with
//      objective_function = DensePeakSizeBatched and a finite peak_threshold,
//      using a generic (stand-in) batchable predicate over the base virtual
//      space "a", and reading back per-node batch_axes() via
//      OptimizeOptions::term_batch_axes -> BinarizationOptions::
//      node_batch_axes -> binarize() -> EvalExpr::batch_axes(). This is the
//      exact plumbing later tasks (and the eventual real verdict) will reuse;
//      it is explicitly NOT a PAO/mu~ finding (see the log labelling below);
//   3. a harness self-check with a tiny forced peak_threshold, proving the
//      plumbing CAN produce a non-empty batch_axes annotation when the DP is
//      forced to batch (rules out "the harness can't express batching" as an
//      explanation for whatever the mechanics case observes);
//   4. a self-consistency memsize check against SizeRegime's own extents/
//      moments (replaces the unreachable cluster K-anchors).

// Flip Trace::Default -> Trace::On for this ENTIRE translation unit (must
// precede every SeQuant/core/eval/eval.hpp inclusion, directly or
// transitively -- header include guards mean the Trace enum's Default member
// is fixed, for this whole TU, by whether this macro is defined before the
// header's FIRST inclusion). This is load-bearing for the [dryrun-eval]
// replay below: make_batched_custom_evaluator's inner per-batch/per-member
// replay calls `evaluate(*mem, le_g, bs.cache)` WITHOUT an explicit `<Trace::
// ...>` argument, so those nested calls (and their note_working_set() calls,
// which are gated at COMPILE TIME on the EvalTrace template argument, not
// just Logger's runtime level) only fire if Trace::Default resolves to On.
// Without this, working_set_hwmark() would reflect only the outermost
// custom-evaluator interception and stay blind to everything the batched
// replay does inside it -- exactly the visibility Task 6's witness needs.
// mpqc's own C60 trace (614336.log, cited throughout this file) was captured
// the same way (SEQUANT_EVAL_TRACE defined in that build).
#define SEQUANT_EVAL_TRACE 1

#include <SeQuant/core/eval/backends/dryrun/cost_model_object.hpp>
#include <SeQuant/core/eval/backends/dryrun/eval_expr.hpp>
#include <SeQuant/core/eval/backends/dryrun/result.hpp>
#include <SeQuant/core/eval/backends/dryrun/size_regime.hpp>

#include <SeQuant/core/batch_policy.hpp>
#include <SeQuant/core/eval/eval.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/expressions/tensor.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/logger.hpp>
#include <SeQuant/core/optimize/optimize.hpp>
#include <SeQuant/core/optimize/options.hpp>
#include <SeQuant/core/optimize/single_term_detail.hpp>
#include <SeQuant/core/runtime.hpp>
#include <SeQuant/core/space.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <chrono>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <vector>

using namespace sequant;
using sequant::eval::dryrun::SizeRegime;

namespace {

std::string slurp(std::string const& path) {
  std::ifstream in(path);
  std::stringstream ss;
  ss << in.rdbuf();
  return ss.str();
}

// The fixture is mpqc's raw debug dump: three "==== Rs[k] ====" header lines,
// each immediately followed by one line holding the serialized residual
// member k (a Sum of Products). Strip the headers; return the three
// expression lines in order.
std::vector<std::string> split_residual_members(std::string const& text) {
  std::vector<std::string> members;
  std::istringstream in(text);
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line.rfind("====", 0) == 0) continue;
    members.push_back(line);
  }
  return members;
}

// Stand-in regime for the pre-transform residual: extents are plausible
// order-of-magnitude values (the real per-atom/per-basis extents live in
// mpqc's post-transform pipeline, not reachable from a bare SeQuant unit
// test -- see the file header). `a` is dual-purpose exactly as in mpqc's
// ctx.idx_to_extent (cck.ipp): a bare (non-proto) occurrence is the generic
// base virtual (sized by space_extent); a proto-indexed occurrence is a PNO
// (2 protos) or OSV (1 proto) domain composite (sized by the moment tables).
// NOTE: bare "a" here is NOT the PAO ("mu~") space (see file header) -- it is
// used only as a stand-in batchable axis to exercise the DP-inspection
// mechanics below.
SizeRegime c60_regime() {
  SizeRegime r;
  r.space_extent = {
      {L"i", 120},   // active occupied (order-of-magnitude stand-in)
      {L"a", 1800},  // generic base virtual (order-of-magnitude stand-in)
  };
  // CC-average PNO ~58, OSV ~20 (order-of-magnitude stand-ins; see report).
  double const pno = 58.0, osv = 20.0;
  for (std::size_t k = 0; k <= 4; ++k) {
    r.csv_pno_moment[k] = std::pow(pno, double(k));
    r.csv_osv_moment[k] = std::pow(osv, double(k));
  }
  return r;
}

// Stand-in batchable predicate for the mechanics check: a non-proto-indexed
// index in the (generic base virtual) space "a". This is NOT a claim that "a"
// is the PAO space -- see the file header. It only needs to pick SOME
// contracted index so the DP-inspection round-trip (OptimizeOptions ->
// term_batch_axes -> BinarizationOptions -> binarize -> EvalExpr::
// batch_axes()) can be exercised end-to-end against the real residual.
bool is_stand_in_batchable(Index const& ix) {
  return ix.space().base_key() == L"a" && !ix.has_proto_indices();
}

// Free (external) indices of a binarized-tree node's result, read off its
// EvalExpr tensor (scalar-result nodes return an empty vector).
std::vector<Index> node_free_indices(EvalExpr const& n) {
  std::vector<Index> v;
  if (!n.is_tensor()) return v;
  for (auto const& ix : n.as_tensor().const_braketaux_indices())
    v.push_back(ix);
  return v;
}

std::wstring describe_indices(std::vector<Index> const& ixs) {
  std::wstring s;
  for (auto const& ix : ixs) {
    s += std::wstring(ix.full_label());
    s += L" ";
  }
  return s;
}

}  // namespace

TEST_CASE("dryrun size regime basic extents", "[dryrun-probe]") {
  auto r = c60_regime();
  // A bare occ index resolves to its space extent.
  auto i = Index{L"i_1"};
  CHECK(r.extent(i) == 120);
  auto a = Index{L"a_1"};
  CHECK(r.extent(a) == 1800);
  CHECK_THROWS_AS(r.extent(Index{L"x_1"}), std::out_of_range);
}

TEST_CASE("dryrun regime self-consistency (no cluster anchors needed)",
          "[dryrun-probe]") {
  // Replaces the (unreachable) two K-anchor checks: pick one concrete index
  // set, hand-compute its expected element count from SizeRegime's own
  // extents/moments, and confirm memsize_counter agrees. This validates the
  // size MODEL (extent lookup + inner_pow composite dispatch), independent of
  // whether the specific numeric extents match any particular cluster run.
  auto r = c60_regime();
  auto memsize = sequant::opt::detail::memsize_counter(r.idx_to_extent(),
                                                       r.inner_pow_fn());

  // A dense operand: two occ + two bare (base-virtual) indices, e.g. the
  // ladder term's g{i,i;a,a}.
  {
    std::vector<Index> ixs{Index{L"i_1"}, Index{L"i_2"}, Index{L"a_1"},
                           Index{L"a_2"}};
    double const expected = 120.0 * 120.0 * 1800.0 * 1800.0;  // i * i * a * a
    double const got = memsize(ixs, std::vector<Index>{}, std::vector<Index>{});
    CHECK(got == Catch::Approx(expected));
    INFO("dense operand elems=" << got << " (" << (got * 8 / 1e9)
                                << " GB dense)");
  }

  // A CSV/PNO composite operand: one occ-pair-domain composite (2 protos),
  // e.g. an intermediate carrying a_1<i_1,i_2>. tot_indices() pulls the
  // composite's own proto indices (i_1, i_2) into the outer set too (they are
  // genuinely free tensor slots the composite depends on), so the expected
  // volume is i_1 * i_2 * <#PNO^1>, not just one occ factor.
  {
    Index i1{L"i_1"}, i2{L"i_2"};
    Index pno_leg{L"a_1", {i1, i2}};
    std::vector<Index> ixs{pno_leg};
    double const expected = 120.0 * 120.0 * 58.0;  // i * i * <#PNO^1>
    double const got = memsize(ixs, std::vector<Index>{}, std::vector<Index>{});
    CHECK(got == Catch::Approx(expected));
  }
}

// Sanity check for the probe's own plumbing (is_stand_in_batchable /
// batch_policy / term_batch_axes / binarize round-trip), independent of the
// real residual: with a tiny peak_threshold that makes every schedule
// infeasible except the sliced one, the DP MUST annotate a batch axis
// somewhere. If this failed it would mean the harness itself cannot express
// batching (rather than the mechanics case below observing a genuine
// non-annotation for some other reason).
TEST_CASE("dryrun harness sanity: forced threshold does annotate a batch axis",
          "[dryrun-probe]") {
  auto ctx = get_default_context().clone();
  ctx.set_first_dummy_index_ordinal(1000000);
  auto ctx_resetter = set_scoped_default_context(std::move(ctx));

  auto r = c60_regime();

  // 4-tensor network with a_1 shared between the two g's (the only batchable
  // contraction), mirroring test_optimize.cpp's "binarize stamps per-node
  // batch axes from optimize()" case but using the stand-in base-virtual
  // ("a", non-proto) axis instead of a DF aux axis.
  auto expr = deserialize<ExprPtr>(
      "g{a_100;i_1;a_1} g{a_101;i_2;a_1} f{i_1;i_3} f{i_2;i_4}");
  REQUIRE(static_cast<bool>(expr));

  auto axes_map = std::make_shared<std::unordered_map<
      Expr const*, container::vector<container::svector<Index>>>>();

  OptimizeOptions opts;
  opts.objective_function = ObjectiveFunction::DensePeakSizeBatched;
  opts.idx_to_extent = r.idx_to_extent();
  opts.inner_pow = r.inner_pow_fn();
  opts.batch_policy.is_batchable_index = [](Index const& ix) {
    return is_stand_in_batchable(ix);
  };
  opts.batch_policy.batch_target_size = [](Index const&) {
    return std::size_t{20};
  };
  opts.batch_policy.peak_threshold = 1.0;  // tiny budget => force batching
  opts.term_batch_axes = axes_map;

  auto optimized = optimize(expr, opts);
  REQUIRE(static_cast<bool>(optimized));

  auto it = axes_map->find(optimized.get());
  REQUIRE(it != axes_map->end());
  auto const& node_axes = it->second;
  REQUIRE(!node_axes.empty());

  BinarizationOptions bopts;
  bopts.node_batch_axes = node_axes;

  SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
  auto node = binarize(optimized, {}, bopts);
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_END

  bool any_annotated = false;
  bool axis_found = false;
  node.visit([&](auto const& n) {
    if (n->batch_axes().empty()) return;
    any_annotated = true;
    for (auto const& ix : n->batch_axes())
      if (is_stand_in_batchable(ix)) axis_found = true;
  });
  CHECK(any_annotated);
  CHECK(axis_found);
}

TEST_CASE("dryrun residual fixture deserializes", "[dryrun-probe]") {
  // The fixture's dummy-index ordinals (a_21514, ...) come from mpqc's
  // globally-incrementing internal counter and can reach into the tens of
  // thousands; raise first_dummy_index_ordinal well above that so the parser
  // does not mistake them for SeQuant's own reserved temporary-index range
  // (default threshold 100, set by test_main.cpp).
  auto ctx = get_default_context().clone();
  ctx.set_first_dummy_index_ordinal(1000000);
  auto ctx_resetter = set_scoped_default_context(std::move(ctx));

  auto const residual = slurp(std::string(SEQUANT_UNIT_TESTS_SOURCE_DIR) +
                              "/data/csv_ccsd_residual.txt");
  REQUIRE(!residual.empty());
  auto members = split_residual_members(residual);
  REQUIRE(members.size() == 3);

  for (std::size_t k = 0; k < members.size(); ++k) {
    auto expr = deserialize<ExprPtr>(members[k]);
    bool const have_expr = static_cast<bool>(expr);
    REQUIRE(have_expr);
    std::size_t nterms =
        expr->is<Sum>() ? expr->as<Sum>().summands().size() : 1;
    INFO("Rs[" << k << "] term count = " << nterms);
    CHECK(nterms > 0);
    std::wcout << L"Rs[" << k << L"] terms=" << nterms << L"\n";
  }
}

// DP-inspection MECHANICS check on the real (pre-transform) residual: run
// optimize() with objective_function = DensePeakSizeBatched and a finite
// peak_threshold, using the stand-in batchable predicate over the generic
// base-virtual space "a" (NOT the PAO/mu~ space -- see file header), purely
// to exercise the API end-to-end and confirm per-node batch_axes() can be
// read back after optimize(). This is NOT a PAO/mu~ go/no-go verdict; that
// verdict needs a post-transform (mu~ + K) fixture (follow-up).
TEST_CASE(
    "dryrun DP-inspection mechanics on pre-transform residual "
    "(verdict pending post-transform fixture)",
    "[dryrun-probe]") {
  // See the ordinal-range comment in the "residual fixture deserializes"
  // case above.
  auto ctx = get_default_context().clone();
  ctx.set_first_dummy_index_ordinal(1000000);
  auto ctx_resetter = set_scoped_default_context(std::move(ctx));

  auto const residual = slurp(std::string(SEQUANT_UNIT_TESTS_SOURCE_DIR) +
                              "/data/csv_ccsd_residual.txt");
  auto members = split_residual_members(residual);
  REQUIRE(members.size() == 3);

  auto r = c60_regime();

  double biggest_bytes = 0.0;
  std::wstring biggest_desc;
  bool biggest_has_axis = false;
  bool any_axis_anywhere = false;
  std::size_t n_carrying_nodes = 0;
  // Confirms the term_batch_axes read-back mechanism itself (not just that
  // optimize() ran): at least one optimized summand must have a non-empty
  // entry in the map, i.e. the DP actually recorded a per-node sliced-set
  // vector for it (see SeQuantEngine::make_optimize_options /
  // sequant.cpp's to_node for the same read-back shape in mpqc).
  std::size_t n_summands_with_nonempty_axes_entry = 0;

  auto memsize = sequant::opt::detail::memsize_counter(r.idx_to_extent(),
                                                       r.inner_pow_fn());

  for (std::size_t k = 0; k < members.size(); ++k) {
    auto expr = deserialize<ExprPtr>(members[k]);
    bool const have_expr = static_cast<bool>(expr);
    REQUIRE(have_expr);

    // One summand per Sum entry. IMPORTANT: mirror mpqc's own to_node
    // (sequant.cpp, process_for_evaluation) and call optimize() SEPARATELY
    // on each individual summand, rather than once on the whole Sum. Calling
    // optimize() on a multi-summand Sum with the default
    // opts.reorder == ReorderSum::Reorder makes optimize_impl's Sum branch
    // run opt::reorder(new_sum, nodes) at the end, which returns a NEW Sum
    // whose summand ExprPtr identities differ from the ones opt_pure_product
    // recorded term_batch_axes against -- so a post-hoc
    // axes_map->find(optimized_summand.get()) on the whole-Sum result always
    // misses (I hit exactly this: it silently returned 0 non-empty entries
    // for all 86 real summands, no crash, no assert -- see report). Per-term
    // optimize() calls take the Product branch of optimize_impl directly, so
    // the returned ExprPtr's identity IS the term_batch_axes key.
    // deserialize() preserves the literal (deeply nested) parenthesization
    // of the source text, e.g. "2 ((g * (C * t)) * C) * t" stays a Product
    // whose own factors can themselves be Products; it is NOT auto-flattened
    // to a single flat factor list. single_term_opt (and thus
    // TensorNetwork::tensors().size(), which run_single_term_opt_axes's
    // nt==1 shortcut keys off) only sees the network of GENUINE Tensor
    // leaves after flattening -- test_optimize.cpp's own "threshold gates
    // batching" test hits the exact same thing and re-flattens explicitly
    // (Product::Flatten::Yes) before calling single_term_opt. Do the same
    // here; without it every term looks like a trivial 1-or-2-tensor network
    // to the DP (which is exactly what I originally observed: node_axes came
    // back empty for all 86 real summands, even ones with genuine multi-GB
    // PAO-role-carrying intermediates found by the manual node walk below).
    auto flatten_product = [](ExprPtr const& e) -> ExprPtr {
      if (!e->is<Product>()) return e;
      auto const& p = e->as<Product>();
      return ex<Product>(p.scalar(), p.factors(), Product::Flatten::Yes);
    };

    std::vector<ExprPtr> input_terms;
    if (expr->is<Sum>()) {
      for (auto const& s : expr->as<Sum>()) input_terms.push_back(s);
    } else {
      input_terms.push_back(expr);
    }
    for (auto& t : input_terms) t = flatten_product(t);

    for (auto const& input_term : input_terms) {
      auto axes_map = std::make_shared<std::unordered_map<
          Expr const*, container::vector<container::svector<Index>>>>();

      OptimizeOptions opts;
      opts.objective_function = ObjectiveFunction::DensePeakSizeBatched;
      opts.idx_to_extent = r.idx_to_extent();
      opts.inner_pow = r.inner_pow_fn();
      opts.batch_policy.is_batchable_index = [](Index const& ix) {
        return is_stand_in_batchable(ix);
      };
      opts.batch_policy.batch_target_size = [](Index const&) {
        return std::size_t{100};
      };
      opts.batch_policy.peak_threshold = 40e9;  // 40 GB, matches the C60 run
      opts.term_batch_axes = axes_map;

      auto term = optimize(input_term, opts);
      bool const have_optimized = static_cast<bool>(term);
      REQUIRE(have_optimized);

      auto it = axes_map->find(term.get());
      container::vector<container::svector<Index>> node_axes;
      if (it != axes_map->end()) {
        node_axes = it->second;
        if (!node_axes.empty()) ++n_summands_with_nonempty_axes_entry;
      }

      BinarizationOptions bopts;
      bopts.node_batch_axes = node_axes;

      SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
      auto node = binarize(term, {}, bopts);
      SEQUANT_PRAGMA_IGNORE_DEPRECATED_END

      // visit_internal only: leaves (original tensor factors, e.g. "g"/"C")
      // are never DP decision points and can never carry a batch-axis
      // annotation, so including them would understate the DP's actual
      // annotation rate on genuine contraction (intermediate) results.
      node.visit_internal([&](auto const& n) {
        auto free_ixs = node_free_indices(*n);
        if (free_ixs.empty()) return;
        bool has_axis_candidate = false;
        for (auto const& ix : free_ixs)
          if (is_stand_in_batchable(ix)) has_axis_candidate = true;
        if (!has_axis_candidate) return;
        ++n_carrying_nodes;

        double const bytes =
            memsize(free_ixs, std::vector<Index>{}, std::vector<Index>{}) * 8.0;

        bool node_has_axis = false;
        for (auto const& ax : n->batch_axes())
          if (is_stand_in_batchable(ax)) node_has_axis = true;
        if (node_has_axis) any_axis_anywhere = true;

        std::wcout << L"[dryrun-probe] Rs[" << k << L"] node free={"
                   << describe_indices(free_ixs) << L"} bytes=" << bytes
                   << L" (" << (bytes / 1e9) << L" GB) batch_axes={"
                   << describe_indices(container::vector<Index>(
                          n->batch_axes().begin(), n->batch_axes().end()))
                   << L"}\n";

        if (bytes > biggest_bytes) {
          biggest_bytes = bytes;
          biggest_desc = describe_indices(free_ixs);
          biggest_has_axis = node_has_axis;
        }
      });
    }
  }

  INFO("largest stand-in-axis-carrying node: "
       << std::string(biggest_desc.begin(), biggest_desc.end())
       << " bytes=" << biggest_bytes);
  std::wcout
      << L"\n=== DP-INSPECTION MECHANICS (pre-transform residual; NOT a "
         L"PAO/mu~ verdict) ===\n"
      << L"stand-in-axis-carrying contraction nodes seen: " << n_carrying_nodes
      << L"\n"
      << L"summands with a non-empty term_batch_axes entry: "
      << n_summands_with_nonempty_axes_entry << L"\n"
      << L"any node annotated with the stand-in batch axis: "
      << (any_axis_anywhere ? L"YES" : L"NO") << L"\n"
      << L"largest stand-in-axis-carrying node: " << biggest_desc << L" bytes="
      << biggest_bytes << L" (" << (biggest_bytes / 1e9) << L" GB)"
      << L" -- got stand-in batch axis: " << (biggest_has_axis ? L"YES" : L"NO")
      << L"\n"
      << L"NOTE: 'a' here is the generic base virtual space, not PAO/mu~; "
         L"the real PAO/K go/no-go verdict needs a post-transform fixture "
         L"(follow-up).\n";

  // Mechanics check only: confirms the round-trip runs and per-node
  // batch_axes() is readable; does not assert a pass/fail verdict on whether
  // any particular node got annotated (see file header).
  REQUIRE(n_carrying_nodes > 0);
  // Confirms the term_batch_axes read-back mechanism specifically: at least
  // one real summand must have produced a non-empty per-node axes vector.
  REQUIRE(n_summands_with_nonempty_axes_entry > 0);
  SUCCEED();
}

// ===========================================================================
// POST-TRANSFORM VERDICT: the actual PAO(mu~)/DF-aux(K) go/no-go probe.
//
// Fixture data/csv_ccsd_doubles_residual_df.txt is the CSV CCSD DOUBLES
// residual dumped from mpqc (repro/w8-batch-min.json) at the EXACT point it is
// handed to sequant::optimize() -- i.e. AFTER the CSV->PAO base transform and
// DF refactorization. It therefore contains the real PAO index "mu~"
// (μ̃_NNNN) and DF-aux index "K" (Κ_N), the 3-center DF integrals
// g{μ̃;i;Κ}, and the CSV coefficients C{a<i>;μ̃}. Standalone SeQuant's default
// mbpt registry has neither mu~ nor K, so we augment a cloned context with
// mbpt::add_pao_spaces (mu~) + mbpt::add_df_spaces (K).
//
// This case reproduces, offline and in milliseconds, the batched-objective DP
// decision that on the cluster left the free-mu~ giant un-sliced (single-mode
// aux batching -> OOM). It reports, per contraction node carrying a free mu~,
// whether the DP annotated a mu~ (or K) batch axis -- the localizing signal:
//   - if the giant mu~-carrying intermediate never gets a mu~ axis even as its
//     modeled size dwarfs the threshold => the COST MODEL / DP is the gap;
//   - if it does get one here => the gap is downstream (binarize/runtime).
// The DryRun harness's value is that the SAME fixture can be swept across size
// regimes (water-8 vs C60) by changing only the SizeRegime extents.
// ===========================================================================

namespace {

// Regime for the post-transform (mu~ + K) residual. Extents are per-index
// domain sizes fed to the DP's idx_to_extent; proto-indexed "a" legs are
// PNO/OSV composites sized by the moment tables (same dispatch as the
// pre-transform regime). Defaults are water-8-scale (from the w8-batch-min
// run: active occ = 32, DF aux cc-pvdz-ri = 672, avg PNOs/pair ~19); the
// giant-relevant knob is the mu~ (PAO domain) extent, swept below.
SizeRegime df_regime(std::size_t mu_tilde, std::size_t aux, std::size_t i_occ,
                     double pno, double osv) {
  SizeRegime r;
  r.space_extent = {
      {L"i", i_occ},
      {L"μ̃", mu_tilde},
      {L"Κ", aux},
      {L"a", mu_tilde},  // bare "a" should not appear (all a are proto here)
  };
  // csv_{pno,osv}_moment[k] must be the k-th POWER MEAN M_k=(<d^k>)^(1/k) of
  // the per-pair domain, NOT the k-th moment <d^k>. inner_aware_volume sizes a
  // k-group of proto-sharing composites as the PRODUCT of inner_pow(c,k) over
  // its k members, so M_k^k telescopes to <d^k> and (times outer occ^N) equals
  // the true block volume Sum_pairs d^k. With only the scalar mean (no per-pair
  // histogram), approximate the domain as constant -> M_k = mean for ALL k.
  // pow(pno,k) (the k-th MOMENT) is wrong: for a k=2 result like
  // R{a_1<i,i>,a_2<i,i>;i,i} it gives occ^2*pno^4 not occ^2*pno^2 (the 358 GB
  // root artifact). k=1 single-composite giants are unaffected (pno^1==pno).
  for (std::size_t k = 0; k <= 4; ++k) {
    r.csv_pno_moment[k] = pno;
    r.csv_osv_moment[k] = osv;
  }
  return r;
}

// Batchable = the two axes mpqc's runtime batches on the CSV path: PAO (mu~)
// and DF aux (K). Both are non-proto base spaces (mu~ = PAO, K = DFBS aux).
bool is_df_batchable(Index const& ix) {
  auto const k = ix.space().base_key();
  return k == L"μ̃" || k == L"Κ";
}

}  // namespace

TEST_CASE("dryrun POST-transform PAO/K batch-axis verdict", "[dryrun-df]") {
  // Augment the default mbpt registry with PAO (mu~) and DF-aux (K) so the
  // post-transform fixture deserializes; raise the dummy-ordinal ceiling for
  // mpqc's high internal ordinals (mu~_1152, a_21674, ...).
  auto ctx = get_default_context().clone();
  ctx.set_first_dummy_index_ordinal(1000000);
  auto isr = ctx.mutable_index_space_registry();
  REQUIRE(isr != nullptr);
  sequant::mbpt::add_pao_spaces(isr);  // mu~
  sequant::mbpt::add_df_spaces(isr);   // K
  auto ctx_resetter = set_scoped_default_context(std::move(ctx));

  auto const body = slurp(std::string(SEQUANT_UNIT_TESTS_SOURCE_DIR) +
                          "/data/csv_ccsd_doubles_residual_df.txt");
  REQUIRE(!body.empty());

  // The fixture body is a single line: one Sum of Products (the whole doubles
  // residual). Deserialize, then split into summands and re-flatten each
  // (deserialize keeps literal nesting; single_term_opt needs a flat factor
  // list or term_batch_axes silently comes back empty -- see the mechanics
  // case above).
  std::string line = body;
  if (auto nl = line.find('\n'); nl != std::string::npos)
    line = line.substr(0, nl);
  auto expr = deserialize<ExprPtr>(line);
  bool const parsed = static_cast<bool>(expr);
  REQUIRE(parsed);

  auto flatten_product = [](ExprPtr const& e) -> ExprPtr {
    if (!e->is<Product>()) return e;
    auto const& p = e->as<Product>();
    return ex<Product>(p.scalar(), p.factors(), Product::Flatten::Yes);
  };
  std::vector<ExprPtr> terms;
  if (expr->is<Sum>())
    for (auto const& s : expr->as<Sum>()) terms.push_back(flatten_product(s));
  else
    terms.push_back(flatten_product(expr));
  REQUIRE(!terms.empty());
  std::wcout << L"\n=== POST-TRANSFORM DOUBLES RESIDUAL: " << terms.size()
             << L" summands (mu~ + K present) ===\n";

  // One regime (C60-scale mu~ domain) at the C60 peak_threshold. Per-term
  // progress + timing to std::cerr (unbuffered) so a slow/pathological term is
  // visible; each term reports whether its largest free-mu~ intermediate got a
  // mu~/K batch axis.
  double const peak_threshold = 40e9;  // 40 GB, matches the C60 run
  // Env overrides for the root-cause BISECT (default = faithful real C60
  // config). DRYRUN_OCC/PNO/OSV vary extents; DRYRUN_PAO_TS/AUX_TS vary batch
  // target sizes; DRYRUN_VW varies volatile_weight; DRYRUN_ROOFLINE=0 disables
  // the roofline tie-break. Absent env var => real-config default.
  auto env_d = [](char const* k, double dflt) {
    char const* v = std::getenv(k);
    return v ? std::atof(v) : dflt;
  };
  auto env_u = [](char const* k, std::size_t dflt) {
    char const* v = std::getenv(k);
    return v ? static_cast<std::size_t>(std::atoll(v)) : dflt;
  };
  std::size_t const occ_ext = env_u("DRYRUN_OCC", 120u);
  double const pno_mom = env_d("DRYRUN_PNO", 42.0);
  double const osv_mom = env_d("DRYRUN_OSV", 310.0);
  std::size_t const pao_ts = env_u("DRYRUN_PAO_TS", 256u);
  std::size_t const aux_ts = env_u("DRYRUN_AUX_TS", 72u);
  double const vol_weight = env_d("DRYRUN_VW", 20.0);
  bool const use_roofline = env_u("DRYRUN_ROOFLINE", 1u) != 0;
  // DRYRUN_AUX_ONLY=1 makes ONLY the DF aux (K) sliceable, NOT PAO (mu~) --
  // reproduces MPQC's aux-only batching to check the proper (gC)^2 PPL
  // factorization keeps the mu~-full giant (large realized peak = OOM) and
  // does NOT form the fully-sliceable 4-PAO integral.
  bool const aux_only = env_u("DRYRUN_AUX_ONLY", 0u) != 0;
  auto is_batchable = [aux_only](Index const& ix) {
    auto const k = ix.space().base_key();
    return aux_only ? (k == L"Κ") : (k == L"μ̃" || k == L"Κ");
  };
  std::wcerr << L"[dryrun-df] config: occ=" << occ_ext << L" pno=" << pno_mom
             << L" osv=" << osv_mom << L" pao_ts=" << pao_ts << L" aux_ts="
             << aux_ts << L" vw=" << vol_weight << L" roofline="
             << (use_roofline ? 1 : 0) << L" aux_only=" << (aux_only ? 1 : 0)
             << L"\n";
  // C60 pVDZ-F12 scale, REAL run dimensions (from the 614336 Owl job log):
  //   active occupied = 120 (tiles [0,120), elements [60,180))
  //   #PAO = #AO = 1800 (872/pair is a sparsity metric, NOT the dense cost)
  //   DF aux (aug-cc-pVDZ-RI) = 4320
  //   Average PNOs per pair = 41.92 ; Average OSVs per pair = 309.67
  // The earlier water-8 occ/PNO (i=32, PNO=19, OSV=57) were as wrong as the
  // K=672; they under-scaled every intermediate and (via peak_threshold
  // pruning) changed the DP's axis choice.
  auto regime = df_regime(/*mu_tilde=*/1800u, /*aux=*/4320u, /*i_occ=*/occ_ext,
                          /*pno=*/pno_mom, /*osv=*/osv_mom);
  auto memsize = sequant::opt::detail::memsize_counter(regime.idx_to_extent(),
                                                       regime.inner_pow_fn());
  auto has_free_mu_tilde = [](std::vector<Index> const& ixs) {
    for (auto const& ix : ixs)
      if (ix.space().base_key() == L"μ̃") return true;
    return false;
  };

  std::size_t total_mu_nodes = 0, n_terms_with_giant = 0;
  std::size_t total_mu_nodes_with_mu_axis = 0, total_mu_nodes_with_k_only = 0;
  double overall_biggest = 0.0;
  bool overall_biggest_has_mu = false, overall_biggest_has_k = false;
  std::wstring overall_biggest_desc, overall_biggest_axes,
      overall_biggest_sources;

  for (std::size_t ti = 0; ti < terms.size(); ++ti) {
    auto t0 = std::chrono::steady_clock::now();
    std::cerr << "[dryrun-df] optimizing term " << (ti + 1) << "/"
              << terms.size() << " ..." << std::flush;

    auto axes_map = std::make_shared<std::unordered_map<
        Expr const*, container::vector<container::svector<Index>>>>();
    OptimizeOptions opts;
    opts.objective_function = ObjectiveFunction::DensePeakSizeBatched;
    opts.idx_to_extent = regime.idx_to_extent();
    opts.inner_pow = regime.inner_pow_fn();
    opts.batch_policy.is_batchable_index = is_batchable;
    // FAITHFUL replica of the real C60 run's OptimizeOptions
    // (make_optimize_options + cck.ipp batch policy, 614336 job log):
    //   batch:pao_target_size=256, batch:aux_target_size=72
    //   optimize:volatile_weight=20, machine_balance=200, fast_mem_elems=1e6
    //   is_volatile_leaf = (label == "t"), accumulation_factor = 1.0 (MPQC
    //   dflt)
    opts.batch_policy.batch_target_size =
        [pao_ts, aux_ts](Index const& ix) -> std::size_t {
      return ix.space().base_key() == L"μ̃" ? pao_ts   // pao_target_size
                                           : aux_ts;  // aux_target_size
    };
    opts.batch_policy.is_volatile_leaf = [](Tensor const& t) {
      return t.label() == L"t";
    };
    opts.batch_policy.accumulation_factor = 1.0;
    opts.batch_policy.peak_threshold = peak_threshold;
    opts.volatile_weight = vol_weight;
    if (use_roofline) {
      opts.roofline.machine_balance = 200.0;
      opts.roofline.fast_mem_elems = 1000000.0;
    }
    opts.term_batch_axes = axes_map;

    auto term = optimize(terms[ti], opts);
    if (!static_cast<bool>(term)) {
      std::cerr << " (skipped)\n";
      continue;
    }
    auto it = axes_map->find(term.get());
    container::vector<container::svector<Index>> node_axes;
    if (it != axes_map->end()) node_axes = it->second;
    BinarizationOptions bopts;
    bopts.node_batch_axes = node_axes;
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto node = binarize(term, {}, bopts);
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END

    // CORRECT metric: a node can only slice a batchable index it CONTRACTS
    // (relax()'s Acand = open at children & closed at parent). A free index on
    // an intermediate is sliceable only at the ANCESTOR that contracts it, so
    // the giant's REALIZED size is its nominal size with each free index that
    // some ancestor sliced reduced to batch_target_size. Walk top-down carrying
    // active = union of ancestor batch_axes (by FULL label, so mu~_1241 sliced
    // above only reduces the SAME mu~_1241 below, not a different mu~_j).
    auto keyof = [](Index const& ix) { return std::wstring(ix.full_label()); };
    auto ext_of = [](Index const& ix) -> double {
      auto bk = ix.space().base_key();
      if (bk == L"μ̃") return 1800.0;
      if (bk == L"Κ") return 4320.0;
      return 0.0;
    };
    auto tgt_of = [pao_ts, aux_ts](Index const& ix) -> double {
      return ix.space().base_key() == L"μ̃" ? double(pao_ts) : double(aux_ts);
    };
    double term_biggest = 0.0;  // max REALIZED bytes over free-mu~ nodes
    std::wstring term_biggest_desc, term_biggest_axes, term_biggest_sources;
    bool term_biggest_has_mu = false;  // its free mu~ ESCAPED slicing
    bool term_biggest_has_k = false;   // its free K ESCAPED slicing
    std::size_t term_mu_nodes = 0;
    // active maps a sliced index's FULL label -> a descriptor of the ANCESTOR
    // node that slices it (its result free-index signature + its batch_axes).
    // This is the node-dump: it turns "escaped={}" from an inference into a
    // concrete "mu~_X is sliced by ancestor <node>" (or ESCAPED).
    std::function<void(std::remove_cvref_t<decltype(node)> const&,
                       std::map<std::wstring, std::wstring>)>
        walk = [&](auto const& n, std::map<std::wstring, std::wstring> active) {
          auto free_ixs = node_free_indices(*n);
          if (has_free_mu_tilde(free_ixs)) {
            ++term_mu_nodes;
            ++total_mu_nodes;
            double nominal =
                memsize(free_ixs, std::vector<Index>{}, std::vector<Index>{}) *
                8.0;
            double factor = 1.0;
            std::vector<Index> escaped;  // batchable free ixs NOT sliced above
            std::wstring sources;        // per-free-batchable-index provenance
            for (auto const& ix : free_ixs) {
              auto bk = ix.space().base_key();
              if (bk != L"μ̃" && bk != L"Κ") continue;
              auto it = active.find(keyof(ix));
              if (it != active.end()) {
                factor *= tgt_of(ix) / ext_of(ix);
                sources += L"\n      " + keyof(ix) + L": SLICED by ancestor " +
                           it->second;
              } else {
                escaped.push_back(ix);
                sources += L"\n      " + keyof(ix) +
                           L": ESCAPED (no ancestor "
                           L"slices it)";
              }
            }
            double const realized = nominal * factor;
            bool mu_esc = false, k_esc = false;
            for (auto const& ix : escaped)
              (ix.space().base_key() == L"μ̃" ? mu_esc : k_esc) = true;
            if (mu_esc) ++total_mu_nodes_with_mu_axis;  // repurposed: escaped
            if (k_esc && !mu_esc) ++total_mu_nodes_with_k_only;
            if (realized > 5e10) {
              std::wcerr << L"\n  [GIANT realized=" << (realized / 1e9)
                         << L"GB nominal=" << (nominal / 1e9) << L"GB] free={"
                         << describe_indices(free_ixs) << L"}" << sources
                         << L"\n";
            }
            if (realized > term_biggest) {
              term_biggest = realized;
              term_biggest_desc = describe_indices(free_ixs);
              term_biggest_has_mu = mu_esc;
              term_biggest_has_k = k_esc;
              term_biggest_sources = sources;
              term_biggest_axes = describe_indices(
                  container::vector<Index>(escaped.begin(), escaped.end()));
            }
          }
          if (!n.leaf()) {
            std::map<std::wstring, std::wstring> child_active = active;
            std::wstring const self_desc =
                L"[free={" + describe_indices(node_free_indices(*n)) +
                L"} batch_axes={" +
                describe_indices(container::vector<Index>(
                    n->batch_axes().begin(), n->batch_axes().end())) +
                L"}]";
            for (auto const& ax : n->batch_axes())
              child_active[keyof(ax)] = self_desc;
            walk(n.left(), child_active);
            walk(n.right(), child_active);
          }
        };
    walk(node, {});
    if (term_mu_nodes > 0) ++n_terms_with_giant;
    if (term_biggest > overall_biggest) {
      overall_biggest = term_biggest;
      overall_biggest_desc = term_biggest_desc;
      overall_biggest_axes = term_biggest_axes;
      overall_biggest_has_mu = term_biggest_has_mu;
      overall_biggest_has_k = term_biggest_has_k;
      overall_biggest_sources = term_biggest_sources;
    }
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                  std::chrono::steady_clock::now() - t0)
                  .count();
    std::cerr << " " << ms << "ms  free-mu~ nodes=" << term_mu_nodes
              << " biggest_realized=" << (term_biggest / 1e9)
              << "GB mu~escaped=" << (term_biggest_has_mu ? "YES" : "NO")
              << " Kescaped=" << (term_biggest_has_k ? "YES" : "NO") << "\n";
  }

  std::wcerr << L"\n=== POST-TRANSFORM VERDICT (REALIZED peak, mu~=1800, "
                L"thr=40GB) ===\n"
             << L"terms with a free-mu~ intermediate: " << n_terms_with_giant
             << L"/" << terms.size() << L"\n"
             << L"free-mu~ contraction nodes: " << total_mu_nodes << L"\n"
             << L"  ... with a free mu~ that ESCAPED slicing: "
             << total_mu_nodes_with_mu_axis << L"\n"
             << L"  ... with only a free K that escaped:      "
             << total_mu_nodes_with_k_only << L"\n"
             << L"LARGEST REALIZED free-mu~ intermediate: {"
             << overall_biggest_desc << L"} = " << (overall_biggest / 1e9)
             << L" GB (realized, after ancestor slices)\n"
             << L"  escaped (un-sliced) indices = {" << overall_biggest_axes
             << L"}\n"
             << L"  -> free mu~ escaped: "
             << (overall_biggest_has_mu ? L"YES" : L"NO")
             << L" | free K escaped: "
             << (overall_biggest_has_k ? L"YES" : L"NO") << L"\n"
             << L"  NODE-DUMP (per free batchable index, which ancestor slices "
                L"it):"
             << overall_biggest_sources << L"\n"
             << L"INTERPRETATION: realized peak = what the runtime actually "
                L"materializes. If it is >> 40GB with mu~ escaped=YES, the DP "
                L"left the giant's free mu~ un-sliced at its contracting "
                L"ancestor -> reproduces the C60 OOM in predicted space.\n";

  REQUIRE(total_mu_nodes > 0);
  SUCCEED();
}

// ===========================================================================
// Task 2-6: DryRun eval BACKEND (zero-data Result + eval_expr) and the
// end-to-end replay harness that WITNESSES the runtime's batch-axis
// realization on the real post-transform giant term. The POST-TRANSFORM
// VERDICT case above established the DP side of this story (the DP DOES
// annotate a mu~ batch axis on the giant, surviving binarize). These new
// cases do NOT modify anything above; [dryrun-probe]/[dryrun-df] stay exactly
// as committed.
// ===========================================================================

namespace {

using sequant::eval::dryrun::CostModel;
using sequant::eval::dryrun::DryRunLeafEvaluator;
using sequant::eval::dryrun::EvalExprDryRun;
using sequant::eval::dryrun::EvalNodeDryRun;
using sequant::eval::dryrun::ExtentOverrides;
using sequant::eval::dryrun::make_dryrun_result;
using sequant::eval::dryrun::ResultDryRun;
using sequant::eval::dryrun::ResultDryRunNested;

// A small, self-consistent regime for the backend unit tests below (distinct
// from c60_regime()/df_regime() above -- this one just needs a couple of
// named spaces plus non-trivial PNO moments).
SizeRegime backend_test_regime() {
  SizeRegime r;
  r.space_extent = {
      {L"i", 10},
      {L"a", 20},
  };
  double const pno = 4.0;
  for (std::size_t k = 0; k <= 4; ++k)
    r.csv_pno_moment[k] = std::pow(pno, double(k));
  r.csv_osv_moment = r.csv_pno_moment;
  return r;
}

std::array<std::any, 3> annot3(container::svector<Index> l,
                               container::svector<Index> r,
                               container::svector<Index> res) {
  return {std::any{std::move(l)}, std::any{std::move(r)},
          std::any{std::move(res)}};
}

}  // namespace

TEST_CASE("dryrun cost model memsize matches memsize_counter",
          "[dryrun-costmodel]") {
  auto r = backend_test_regime();
  CostModel cm{r};
  container::svector<Index> idx{Index{L"i_1"}, Index{L"i_2"}, Index{L"a_4"}};
  auto direct = sequant::opt::detail::memsize_counter(
      r.idx_to_extent(), r.inner_pow_fn())(idx, container::svector<Index>{},
                                           container::svector<Index>{});
  CHECK(cm.memsize(idx) == static_cast<std::size_t>(direct * 8.0));
}

TEST_CASE("dryrun cost model memsize honors an extent override",
          "[dryrun-costmodel]") {
  auto r = backend_test_regime();
  CostModel cm{r};
  container::svector<Index> idx{Index{L"a_3"}, Index{L"i_1"}};
  auto const full = cm.memsize(idx);
  ExtentOverrides ov;
  ov[Index{L"a_3"}] = 5;  // narrowed from 20 to 5
  auto const sliced = cm.memsize(idx, ov);
  CHECK(sliced < full);
  CHECK(full == sliced * 4);  // linear in a_3's extent
}

TEST_CASE("dryrun cost model flops and exec_cost are finite/positive",
          "[dryrun-costmodel]") {
  auto r = backend_test_regime();
  CostModel cm{r};
  container::svector<Index> out{Index{L"a_3"}, Index{L"i_1"}};
  container::svector<Index> contracted{Index{L"i_2"}};
  auto const f = cm.flops(out, contracted);
  CHECK(f > 0.0);
  CHECK(cm.exec_cost(f, cm.memsize(out), 4096) > 0.0);
}

TEST_CASE("dryrun flat result size delegates to cost model",
          "[dryrun-result]") {
  // Result::prod/sum/permute/slice_mode/mode_batches/size_in_bytes are all
  // overridden PRIVATE in the concrete DryRun classes (mirroring
  // ResultTensorTAPP), since real callers only ever reach a Result through a
  // ResultPtr/Result const& (base-class access); tests that want to call
  // them on a concrete object must do the same -- via a `Result const&`.
  auto r = backend_test_regime();
  auto cm = std::make_shared<CostModel const>(r);
  container::svector<Index> idx{Index{L"i_1"}, Index{L"i_2"}, Index{L"a_4"}};
  ResultDryRun t{idx, cm};
  Result const& rt = t;
  CHECK(rt.size_in_bytes() == cm->memsize(idx));
}

TEST_CASE("dryrun flat result prod yields the result annotation index set",
          "[dryrun-result]") {
  auto r = backend_test_regime();
  auto cm = std::make_shared<CostModel const>(r);
  ResultDryRun l{{Index{L"a_3"}, Index{L"i_2"}}, cm};
  ResultDryRun rr{{Index{L"i_2"}, Index{L"i_1"}}, cm};
  Result const& rl = l;
  container::svector<Index> res{Index{L"a_3"}, Index{L"i_1"}};
  auto out = rl.prod(rr, annot3(l.indices(), rr.indices(), res), DeNest::False);
  REQUIRE(out);
  bool const is_flat = out->is<ResultDryRun>();
  CHECK(is_flat);
  auto const& ot = out->as<ResultDryRun>();
  CHECK(ot.indices() == res);
  CHECK(out->size_in_bytes() == cm->memsize(res));
}

TEST_CASE("dryrun flat result full contraction yields a scalar",
          "[dryrun-result]") {
  auto r = backend_test_regime();
  auto cm = std::make_shared<CostModel const>(r);
  ResultDryRun l{{Index{L"i_1"}}, cm};
  ResultDryRun rr{{Index{L"i_1"}}, cm};
  Result const& rl = l;
  auto out = rl.prod(rr, annot3(l.indices(), rr.indices(), {}), DeNest::False);
  REQUIRE(out);
  bool const is_scalar = out->is<ResultScalar<double>>();
  CHECK(is_scalar);
}

TEST_CASE("dryrun flat result slice_mode shrinks the sliced mode",
          "[dryrun-result]") {
  auto r = backend_test_regime();
  auto cm = std::make_shared<CostModel const>(r);
  ResultDryRun t{{Index{L"a_3"}, Index{L"i_2"}}, cm};  // mu~ extent 20
  Result const& rt = t;
  auto const full = rt.size_in_bytes();
  auto sliced = rt.slice_mode(0, 0, 5);  // quarter of mu~
  REQUIRE(sliced);
  CHECK(sliced->size_in_bytes() < full);
  CHECK(sliced->size_in_bytes() == full / 4);
}

TEST_CASE("dryrun flat result mode_batches tiles the mode extent",
          "[dryrun-result]") {
  auto r = backend_test_regime();
  auto cm = std::make_shared<CostModel const>(r);
  ResultDryRun t{{Index{L"a_3"}, Index{L"i_2"}}, cm};
  Result const& rt = t;
  auto batches = rt.mode_batches(0, 5);  // 20 / 5 = 4 batches
  CHECK(batches.size() == 4);
  CHECK(batches.front().first == 0);
  CHECK(batches.back().second == 20);
}

TEST_CASE("dryrun nested result uses moment-aware inner extent, not extent^k",
          "[dryrun-nested]") {
  auto r = backend_test_regime();
  // Non-trivial second moment: <#PNO^2> != <#PNO>^2 (dispersion inflates it).
  r.csv_pno_moment[1] = 4.0;
  r.csv_pno_moment[2] = 4.0 * 4.0 * 1.5;
  auto cm = std::make_shared<CostModel const>(r);

  Index i1{L"i_1"}, i2{L"i_2"};
  Index a_pno{L"a_1", {i1, i2}};  // proto-indexed (CSV/PNO composite) leg
  container::svector<Index> outer{i1, i2, Index{L"a_3"}};
  container::svector<Index> inner{a_pno};
  ResultDryRunNested c{outer, inner, cm};

  CHECK(c.outer() == outer);
  CHECK(c.inner() == inner);

  container::svector<Index> combined = outer;
  combined.push_back(a_pno);
  CHECK(c.indices() == combined);

  auto const without_composite = cm->memsize(outer);  // just i1*i2*mu~
  auto const with_composite = cm->memsize(combined);  // routes a_pno via k=1
  // The composite contributes the FIRST moment (4.0), not extent(a_pno)^1
  // (a_pno's own "extent" as a bare, non-composite space is never queried
  // here -- only backend_test_regime()'s csv_pno_moment[1] is), so the ratio
  // is exactly the first moment.
  CHECK(with_composite == without_composite * 4);
  Result const& rc = c;
  CHECK(rc.size_in_bytes() == with_composite);
}

TEST_CASE("dryrun make_dryrun_result dispatches flat vs nested by content",
          "[dryrun-nested]") {
  auto r = backend_test_regime();
  auto cm = std::make_shared<CostModel const>(r);
  Index i1{L"i_1"}, i2{L"i_2"};
  Index a_pno{L"a_1", {i1, i2}};

  auto flat = make_dryrun_result({Index{L"i_1"}, Index{L"i_2"}}, cm);
  CHECK(flat->is<ResultDryRun>());

  auto nested = make_dryrun_result({i1, i2, a_pno}, cm);
  CHECK(nested->is<ResultDryRunNested>());
}

TEST_CASE("dryrun leaf yielder builds a sized token from a tensor leaf",
          "[dryrun-leaf]") {
  auto ctx = get_default_context().clone();
  ctx.set_first_dummy_index_ordinal(1000000);
  auto ctx_resetter = set_scoped_default_context(std::move(ctx));

  auto r = backend_test_regime();
  auto cm = std::make_shared<CostModel const>(r);

  auto expr = deserialize<ExprPtr>("g{i_1,i_2;a_4}");
  bool const parsed = static_cast<bool>(expr);
  REQUIRE(parsed);

  SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
  auto node = binarize<EvalExprDryRun>(expr);
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
  REQUIRE(node.leaf());

  DryRunLeafEvaluator yield{cm};
  auto res = yield(node);
  REQUIRE(res);
  bool const is_flat = res->is<ResultDryRun>();
  CHECK(is_flat);

  container::svector<Index> idx{Index{L"i_1"}, Index{L"i_2"}, Index{L"a_4"}};
  CHECK(res->size_in_bytes() == cm->memsize(idx));
}

// ===========================================================================
// Task 6: THE replay harness. Deserializes the real post-transform giant
// term (the FIRST summand of csv_ccsd_doubles_residual_df.txt -- see the
// POST-TRANSFORM VERDICT case above, which already established this is the
// ~13-tensor free-mu~ giant), optimizes it once under the SAME C60-scale
// regime/BatchPolicy the DP verdict used, binarizes with the DP's
// batch-axis annotations, then REPLAYS it through the REAL runtime
// (make_evaluator / evaluate<Trace::On>) against zero-data DryRun tokens --
// witnessing what the runtime actually realizes, not what the DP annotated.
// ===========================================================================

TEST_CASE(
    "dryrun eval backend replays the post-transform giant term through the "
    "real batched runtime",
    "[dryrun-eval]") {
  auto ctx = get_default_context().clone();
  ctx.set_first_dummy_index_ordinal(1000000);
  auto isr = ctx.mutable_index_space_registry();
  REQUIRE(isr != nullptr);
  sequant::mbpt::add_pao_spaces(isr);  // mu~
  sequant::mbpt::add_df_spaces(isr);   // K
  auto ctx_resetter = set_scoped_default_context(std::move(ctx));

  auto const body = slurp(std::string(SEQUANT_UNIT_TESTS_SOURCE_DIR) +
                          "/data/csv_ccsd_doubles_residual_df.txt");
  REQUIRE(!body.empty());
  std::string line = body;
  if (auto nl = line.find('\n'); nl != std::string::npos)
    line = line.substr(0, nl);
  auto expr = deserialize<ExprPtr>(line);
  bool const parsed = static_cast<bool>(expr);
  REQUIRE(parsed);
  bool const is_sum = expr->is<Sum>();
  REQUIRE(is_sum);

  // Focus the replay on the GIANT TERM only -- the batched DP is 20-60s per
  // term (running it on all 55 would take ~20-30 minutes). The [dryrun-df]
  // verdict case above's exhaustive sweep (`for (auto const& s :
  // expr->as<Sum>()) ...`, the SAME deserialize-order iteration as
  // `summands()` below) empirically identified summand index 38 (term 39/55,
  // 1-indexed in that sweep's log) as the giant: the ~13-tensor
  // g.C.g.C.s.C.C.s.C.C.t.t.t chain reporting the 1.2 TB free-mu~
  // intermediate at C60 (mu~=1800, K=4320) scale. A cheap proxy (picking the
  // summand with the most flattened tensor factors) was tried first and
  // picked a DIFFERENT, structurally-similar but much smaller term (14
  // factors, ~0.0005 GB giant) -- factor count alone does not identify the
  // giant, only the DP's actual index-space accounting does. So this uses
  // the exhaustively-verified positional index directly rather than a
  // heuristic that was empirically shown to pick the wrong term.
  auto const& summands = expr->as<Sum>().summands();
  REQUIRE(!summands.empty());
  auto flatten_product = [](ExprPtr const& e) -> ExprPtr {
    if (!e->is<Product>()) return e;
    auto const& p = e->as<Product>();
    return ex<Product>(p.scalar(), p.factors(), Product::Flatten::Yes);
  };
  std::size_t const giant_idx = 38 < summands.size() ? 38 : 0;
  ExprPtr giant = flatten_product(summands[giant_idx]);
  REQUIRE(giant);
  std::size_t const giant_nfactors =
      giant->is<Product>() ? giant->as<Product>().factors().size() : 1;
  std::cerr << "[dryrun-eval] selected giant term (index " << giant_idx
            << "): " << giant_nfactors << " flattened factors (of "
            << summands.size() << " summands)\n";

  // FAITHFUL real C60 config (614336 job log): occ=120, PNO=42, OSV=310,
  // mu~=1800, aux=4320, pao_target_size=256, aux_target_size=72,
  // volatile_weight=20, machine_balance=200, fast_mem_elems=1e6. SAME regime
  // and SAME batchable predicate the [dryrun-df] verdict case uses.
  auto regime = df_regime(/*mu_tilde=*/1800u, /*aux=*/4320u, /*i_occ=*/120u,
                          /*pno=*/42.0, /*osv=*/310.0);
  auto cm = std::make_shared<CostModel const>(regime);

  // ONE BatchPolicy object, reused verbatim for both optimize() and the
  // runtime evaluator factory (make_evaluator) -- the plan's hard
  // constraint, so the DP's and the runtime's notion of "batchable" and
  // "target batch size" cannot drift apart.
  sequant::BatchPolicy policy;
  policy.is_batchable_index = is_df_batchable;
  policy.batch_target_size = [](Index const& ix) -> std::size_t {
    return ix.space().base_key() == L"μ̃" ? std::size_t{256}  // pao_target_size
                                         : std::size_t{72};  // aux_target_size
  };
  policy.is_volatile_leaf = [](Tensor const& t) { return t.label() == L"t"; };
  policy.accumulation_factor = 1.0;
  policy.peak_threshold = 40e9;  // DP-side knob only; the runtime evaluator
                                 // never consults peak_threshold.

  auto axes_map = std::make_shared<std::unordered_map<
      Expr const*, container::vector<container::svector<Index>>>>();
  OptimizeOptions opts;
  opts.objective_function = ObjectiveFunction::DensePeakSizeBatched;
  opts.idx_to_extent = regime.idx_to_extent();
  opts.inner_pow = regime.inner_pow_fn();
  opts.batch_policy = policy;
  opts.volatile_weight = 20.0;
  opts.roofline.machine_balance = 200.0;
  opts.roofline.fast_mem_elems = 1000000.0;
  opts.term_batch_axes = axes_map;

  auto t0 = std::chrono::steady_clock::now();
  auto optimized = optimize(giant, opts);
  auto const opt_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                          std::chrono::steady_clock::now() - t0)
                          .count();
  std::cerr << "[dryrun-eval] optimize(giant) took " << opt_ms << "ms\n";
  bool const have_optimized = static_cast<bool>(optimized);
  REQUIRE(have_optimized);

  auto it = axes_map->find(optimized.get());
  container::vector<container::svector<Index>> node_axes;
  if (it != axes_map->end()) node_axes = it->second;
  REQUIRE(!node_axes.empty());

  BinarizationOptions bopts;
  bopts.node_batch_axes = node_axes;

  SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
  auto node = binarize<EvalExprDryRun>(optimized, {}, bopts);
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_END

  // Locate the giant sub-node (the free-mu~ contraction node whose modeled
  // size dwarfs everything else -- same identification criterion the
  // [dryrun-df] verdict case above used) purely to REPORT what the DP
  // annotated on it before the runtime replay, for side-by-side comparison
  // with what the runtime actually realizes.
  auto memsize_ext = sequant::opt::detail::memsize_counter(
      regime.idx_to_extent(), regime.inner_pow_fn());
  auto has_free_mu_tilde = [](std::vector<Index> const& ixs) {
    for (auto const& ix : ixs)
      if (ix.space().base_key() == L"μ̃") return true;
    return false;
  };
  double giant_nominal_bytes = 0.0;
  std::wstring giant_desc, giant_axes_desc;
  node.visit_internal([&](auto const& n) {
    auto free_ixs = node_free_indices(*n);
    if (!has_free_mu_tilde(free_ixs)) return;
    double const bytes =
        memsize_ext(free_ixs, std::vector<Index>{}, std::vector<Index>{}) * 8.0;
    if (bytes > giant_nominal_bytes) {
      giant_nominal_bytes = bytes;
      giant_desc = describe_indices(free_ixs);
      giant_axes_desc = describe_indices(container::vector<Index>(
          n->batch_axes().begin(), n->batch_axes().end()));
    }
  });
  REQUIRE(giant_nominal_bytes > 0.0);

  auto cache = sequant::cache_manager(std::vector<EvalNodeDryRun>{node});
  cache.set_custom_evaluator(
      sequant::make_evaluator(policy, DryRunLeafEvaluator{cm}));

  // Enable eval tracing (redirected to a private ostringstream, not stdout)
  // so working_set_hwmark() actually accumulates: the engine's per-op hwmark
  // input is gated at RUNTIME on Logger::instance().eval.level > 0
  // (log::printing()), independent of the Trace::On COMPILE-TIME template
  // argument below (which only gates whether the tracing code path exists at
  // all). Restore the previous logger state afterward so this test does not
  // leak global state to others.
  std::ostringstream trace_os;
  auto& logger = Logger::instance();
  auto const prev_level = logger.eval.level;
  auto* const prev_stream = logger.eval.stream;
  logger.eval.level = 2;
  logger.eval.stream = &trace_os;

  std::cerr << "[dryrun-eval] replaying giant term through the batched "
               "runtime evaluator ...\n";
  auto t1 = std::chrono::steady_clock::now();
  ResultPtr result;
  bool threw = false;
  std::string what;
  try {
    result = sequant::evaluate<Trace::On>(node, DryRunLeafEvaluator{cm}, cache);
  } catch (std::exception const& e) {
    threw = true;
    what = e.what();
  }
  auto const eval_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                           std::chrono::steady_clock::now() - t1)
                           .count();

  logger.eval.level = prev_level;
  logger.eval.stream = prev_stream;

  INFO("evaluate() threw: " << what);
  REQUIRE(!threw);
  REQUIRE(result);

  auto const peak = cache.working_set_hwmark();
  auto const root_bytes = result->size_in_bytes();

  // Parse the captured per-op trace for the largest single `result=<N>B`
  // materialized anywhere during the replay -- the most direct witness of
  // whether the giant was EVER realized at (close to) its full un-sliced
  // size, regardless of caching/accounting nuances in working_set_hwmark().
  std::size_t max_single_result_bytes = 0;
  {
    std::string const trace = trace_os.str();
    std::size_t pos = 0;
    while ((pos = trace.find("result=", pos)) != std::string::npos) {
      pos += std::string("result=").size();
      std::size_t end = trace.find('B', pos);
      if (end == std::string::npos) break;
      std::string const num = trace.substr(pos, end - pos);
      if (!num.empty() &&
          num.find_first_not_of("0123456789") == std::string::npos) {
        std::size_t const v = std::stoull(num);
        max_single_result_bytes = std::max(max_single_result_bytes, v);
      }
      pos = end;
    }
  }

  bool const giant_realized_full =
      max_single_result_bytes >= 0.9 * giant_nominal_bytes;

  // Extra corroborating diagnostics: how many ops ran, and how many distinct
  // BatchGroup interceptions fired (a real nested multi-batch replay should
  // show many -- the giant's mu~ (extent 1800, target 100 => ~18 batches) and
  // K (extent 4320, target 100 => ~44 batches) axes, nested, would fire
  // hundreds of small batched ops if genuinely realized).
  std::string const trace = trace_os.str();
  auto count_occurrences = [](std::string const& hay,
                              std::string const& needle) {
    std::size_t n = 0, pos = 0;
    while ((pos = hay.find(needle, pos)) != std::string::npos) {
      ++n;
      pos += needle.size();
    }
    return n;
  };
  std::size_t const n_eval_lines = count_occurrences(trace, "Eval |");
  std::size_t const n_batch_group_begin =
      count_occurrences(trace, "BatchGroup | Begin");
  // Distribution of result= sizes >= 100 MB, to see whether many large
  // (but sub-nominal) intermediates appeared (consistent with a partially-
  // sliced axis) or none did (consistent with the OTHER axis alone already
  // bounding everything well below 100 MB).
  std::size_t n_results_over_100mb = 0;
  {
    std::size_t pos = 0;
    while ((pos = trace.find("result=", pos)) != std::string::npos) {
      pos += std::string("result=").size();
      std::size_t end = trace.find('B', pos);
      if (end == std::string::npos) break;
      std::string const num = trace.substr(pos, end - pos);
      if (!num.empty() &&
          num.find_first_not_of("0123456789") == std::string::npos) {
        if (std::stoull(num) >= 100'000'000ull) ++n_results_over_100mb;
      }
      pos = end;
    }
  }

  std::wcout << L"\n=== [dryrun-eval] GIANT TERM RUNTIME REPLAY ===\n"
             << L"optimize(): " << opt_ms << L"ms, evaluate(): " << eval_ms
             << L"ms\n"
             << L"DP-annotated giant: free={" << giant_desc << L"} nominal="
             << (giant_nominal_bytes / 1e9) << L" GB batch_axes={"
             << giant_axes_desc << L"}\n"
             << L"root result size = " << root_bytes << L" bytes ("
             << (double(root_bytes) / 1e9) << L" GB)\n"
             << L"cache.working_set_hwmark() = " << peak << L" bytes ("
             << (double(peak) / 1e9) << L" GB)\n"
             << L"max single result= observed in trace = "
             << max_single_result_bytes << L" bytes ("
             << (double(max_single_result_bytes) / 1e9) << L" GB)\n"
             << L"trace ops: " << n_eval_lines << L" Eval lines, "
             << n_batch_group_begin << L" BatchGroup interceptions, "
             << n_results_over_100mb << L" results >=100MB\n"
             << L"=> giant realized at (>=90% of) its FULL nominal size during "
                L"replay: "
             << (giant_realized_full
                     ? L"YES (mu~ NOT sliced in practice -- bug "
                       L"reproduced)"
                     : L"NO (sliced down from nominal)")
             << L"\n";

  CHECK(peak > 0);
}

// Task 6 (perf-first validation): optimize the SAME C60 giant term (index 38)
// under BOTH the peak-first (DenseSpaceTimeBatched) and perf-first
// (DenseTimeSpaceBatched) objectives, at the faithful real config, and compare
// the factorization the DP picks. The 4-PAO signature is a contraction node
// carrying >= 4 free mu~ indices (the (mu~ mu~|mu~ mu~) AO integral).
// Peak-first forms it (fully sliceable below the 40 GB threshold, so it
// survives the hard filter despite astronomically higher flops); perf-first is
// flops-primary and must NEVER form it. This is the direct in-harness proof of
// the fix, on ONE term (~seconds), without the [dryrun-df] full-sweep cost.
TEST_CASE("dryrun perf-first vs peak-first factorization of the C60 giant term",
          "[dryrun-perf]") {
  auto ctx = get_default_context().clone();
  ctx.set_first_dummy_index_ordinal(1000000);
  auto isr = ctx.mutable_index_space_registry();
  REQUIRE(isr != nullptr);
  sequant::mbpt::add_pao_spaces(isr);  // mu~
  sequant::mbpt::add_df_spaces(isr);   // K
  auto ctx_resetter = set_scoped_default_context(std::move(ctx));

  auto const body = slurp(std::string(SEQUANT_UNIT_TESTS_SOURCE_DIR) +
                          "/data/csv_ccsd_doubles_residual_df.txt");
  REQUIRE(!body.empty());
  std::string line = body;
  if (auto nl = line.find('\n'); nl != std::string::npos)
    line = line.substr(0, nl);
  auto expr = deserialize<ExprPtr>(line);
  REQUIRE(static_cast<bool>(expr));
  REQUIRE(expr->is<Sum>());
  auto const& summands = expr->as<Sum>().summands();
  REQUIRE(!summands.empty());
  auto flatten_product = [](ExprPtr const& e) -> ExprPtr {
    if (!e->is<Product>()) return e;
    auto const& p = e->as<Product>();
    return ex<Product>(p.scalar(), p.factors(), Product::Flatten::Yes);
  };
  std::size_t const giant_idx = 38 < summands.size() ? 38 : 0;
  ExprPtr giant = flatten_product(summands[giant_idx]);
  REQUIRE(giant);

  // FAITHFUL real C60 config (614336 job log), identical to [dryrun-eval].
  auto regime = df_regime(/*mu_tilde=*/1800u, /*aux=*/4320u, /*i_occ=*/120u,
                          /*pno=*/42.0, /*osv=*/310.0);
  auto memsize = sequant::opt::detail::memsize_counter(regime.idx_to_extent(),
                                                       regime.inner_pow_fn());
  auto ext_of = [](Index const& ix) -> double {
    auto bk = ix.space().base_key();
    if (bk == L"μ̃") return 1800.0;
    if (bk == L"Κ") return 4320.0;
    return 0.0;
  };
  auto tgt_of = [](Index const& ix) -> double {
    return ix.space().base_key() == L"μ̃" ? 256.0 : 72.0;
  };
  auto keyof = [](Index const& ix) { return std::wstring(ix.full_label()); };
  auto cm = std::make_shared<CostModel const>(regime);

  struct Analysis {
    std::size_t max_free_mu = 0;       // >= 4 => 4-PAO AO integral formed
    double largest_realized_gb = 0.0;  // DP-model realized free-mu~ (static)
    std::wstring largest_desc;
    // Runtime-replay metrics (zero-data eval through the REAL eval loop + REAL
    // CacheManager):
    double hwmark_gb =
        0.0;  // cache.working_set_hwmark() = per-op running max
              // of (sum of alive cached entries + in-flight
              // result). "high watermark of cached intermediates"
    double max_single_result_gb = 0.0;  // largest single result= in the trace =
                                        // the largest transient materialized
    std::string largest_result_line;    // full trace line of that result=
                                        // (shows result= and hw= together)
    bool replay_threw = false;
  };

  // Optimize the giant under `obj`, binarize, and walk the tree computing per
  // node (a) its free-mu~ count and (b) its REALIZED free-mu~ size after
  // ancestor slicing (same active-ancestor accounting as the [dryrun-df]
  // verdict case). Then REPLAY the zero-data schedule through the real eval
  // loop + real CacheManager to capture the cached-intermediate high watermark
  // (working_set_hwmark) alongside the largest single transient (result=), so
  // the cache-vs-transient accounting gap is visible per objective.
  auto analyze = [&](ObjectiveFunction obj) -> Analysis {
    sequant::BatchPolicy policy;
    policy.is_batchable_index = is_df_batchable;
    policy.batch_target_size = [](Index const& ix) -> std::size_t {
      return ix.space().base_key() == L"μ̃" ? std::size_t{256} : std::size_t{72};
    };
    policy.is_volatile_leaf = [](Tensor const& t) { return t.label() == L"t"; };
    policy.accumulation_factor = 1.0;
    policy.peak_threshold = 40e9;

    auto axes_map = std::make_shared<std::unordered_map<
        Expr const*, container::vector<container::svector<Index>>>>();
    OptimizeOptions opts;
    opts.objective_function = obj;
    opts.idx_to_extent = regime.idx_to_extent();
    opts.inner_pow = regime.inner_pow_fn();
    opts.batch_policy = policy;
    opts.volatile_weight = 20.0;
    opts.roofline.machine_balance = 200.0;
    opts.roofline.fast_mem_elems = 1000000.0;
    opts.term_batch_axes = axes_map;

    auto t0 = std::chrono::steady_clock::now();
    auto optimized = optimize(giant, opts);
    auto const opt_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                            std::chrono::steady_clock::now() - t0)
                            .count();
    REQUIRE(static_cast<bool>(optimized));
    auto it = axes_map->find(optimized.get());
    container::vector<container::svector<Index>> node_axes;
    if (it != axes_map->end()) node_axes = it->second;
    BinarizationOptions bopts;
    bopts.node_batch_axes = node_axes;
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto node = binarize<EvalExprDryRun>(optimized, {}, bopts);
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END

    // Optional full-tree dump (DRYRUN_PERF_TREE=1): every node's free indices
    // and batch_axes (the indices SLICED, i.e. batched, at that node), in
    // post-order-ish indentation, so the exact schedule is inspectable.
    if (std::getenv("DRYRUN_PERF_TREE")) {
      wchar_t const* on = (obj == ObjectiveFunction::DenseTimeSpaceBatched)
                              ? L"perf-first (DenseTimeSpaceBatched)"
                              : L"peak-first (DenseSpaceTimeBatched)";
      std::wcerr << L"\n[dryrun-perf] TREE for " << on << L":\n";
      std::function<void(std::remove_cvref_t<decltype(node)> const&, int)>
          dump = [&](auto const& n, int depth) {
            std::wstring const pad(2 * depth + 2, L' ');
            auto const free = node_free_indices(*n);
            container::vector<Index> const bax(n->batch_axes().begin(),
                                               n->batch_axes().end());
            std::size_t nmu = 0, nk = 0;
            for (auto const& ix : free) {
              if (ix.space().base_key() == L"μ̃") ++nmu;
              if (ix.space().base_key() == L"Κ") ++nk;
            }
            std::wcerr << pad << (n.leaf() ? L"leaf  " : L"CONTRACT ")
                       << L"free={" << describe_indices(free) << L"} (mu~="
                       << nmu << L" K=" << nk << L")  batch_axes={"
                       << describe_indices(bax) << L"}\n";
            if (!n.leaf()) {
              dump(n.left(), depth + 1);
              dump(n.right(), depth + 1);
            }
          };
      dump(node, 0);
    }

    Analysis a;
    std::function<void(std::remove_cvref_t<decltype(node)> const&,
                       std::map<std::wstring, std::wstring>)>
        walk = [&](auto const& n, std::map<std::wstring, std::wstring> active) {
          auto free_ixs = node_free_indices(*n);
          std::size_t nmu = 0;
          for (auto const& ix : free_ixs)
            if (ix.space().base_key() == L"μ̃") ++nmu;
          if (nmu > a.max_free_mu) a.max_free_mu = nmu;
          if (nmu > 0) {
            double nominal =
                memsize(free_ixs, std::vector<Index>{}, std::vector<Index>{}) *
                8.0;
            double factor = 1.0;
            for (auto const& ix : free_ixs) {
              auto bk = ix.space().base_key();
              if (bk != L"μ̃" && bk != L"Κ") continue;
              if (active.find(keyof(ix)) != active.end())
                factor *= tgt_of(ix) / ext_of(ix);
            }
            double const realized_gb = nominal * factor / 1e9;
            if (realized_gb > a.largest_realized_gb) {
              a.largest_realized_gb = realized_gb;
              a.largest_desc = describe_indices(free_ixs);
            }
          }
          if (!n.leaf()) {
            std::map<std::wstring, std::wstring> child_active = active;
            for (auto const& ax : n->batch_axes())
              child_active[keyof(ax)] = L"y";
            walk(n.left(), child_active);
            walk(n.right(), child_active);
          }
        };
    walk(node, {});

    // ---- runtime replay: predict the cached-intermediate high watermark ----
    // Zero-data eval through the real eval loop + real CacheManager, with the
    // real batched custom evaluator. working_set_hwmark() folds, per op, (sum
    // of alive cached entries + in-flight result) into a running max; the
    // trace's result= field is the largest single materialized transient.
    // (Single-term, default cache_manager: this is the WITHIN-TERM cached
    // working set. The whole-residual / cross-iteration watermark -- one shared
    // cache over the full DAG with the real max_footprint gate -- is the
    // deferred DryRun-as-prefix wiring.)
    auto cache = sequant::cache_manager(std::vector<EvalNodeDryRun>{node});
    cache.set_custom_evaluator(
        sequant::make_evaluator(policy, DryRunLeafEvaluator{cm}));
    std::ostringstream trace_os;
    auto& logger = Logger::instance();
    auto const prev_level = logger.eval.level;
    auto* const prev_stream = logger.eval.stream;
    logger.eval.level = 2;  // gate log::printing() on so hwmark accumulates
    logger.eval.stream = &trace_os;
    try {
      (void)sequant::evaluate<Trace::On>(node, DryRunLeafEvaluator{cm}, cache);
    } catch (std::exception const&) {
      a.replay_threw = true;
    }
    logger.eval.level = prev_level;
    logger.eval.stream = prev_stream;
    a.hwmark_gb = double(cache.working_set_hwmark()) / 1e9;
    // Largest single result= in the trace, and its full line (result= and hw=
    // appear together, so the line exposes whether the giant transient is
    // folded into the cache hwmark).
    {
      std::string const trace = trace_os.str();
      std::size_t pos = 0, best = 0;
      while ((pos = trace.find("result=", pos)) != std::string::npos) {
        std::size_t const np = pos + std::string("result=").size();
        std::size_t const end = trace.find('B', np);
        if (end == std::string::npos) break;
        std::string const num = trace.substr(np, end - np);
        if (!num.empty() &&
            num.find_first_not_of("0123456789") == std::string::npos) {
          std::size_t const v = std::stoull(num);
          if (v > best) {
            best = v;
            std::size_t ls = trace.rfind('\n', pos);
            ls = (ls == std::string::npos) ? 0 : ls + 1;
            std::size_t const le = trace.find('\n', pos);
            a.largest_result_line = trace.substr(
                ls, (le == std::string::npos ? trace.size() : le) - ls);
          }
        }
        pos = end;
      }
      a.max_single_result_gb = double(best) / 1e9;
    }

    wchar_t const* obj_name = (obj == ObjectiveFunction::DenseTimeSpaceBatched)
                                  ? L"perf-first (DenseTimeSpaceBatched)"
                                  : L"peak-first (DenseSpaceTimeBatched)";
    std::wcerr << L"[dryrun-perf] " << obj_name << L": optimize " << opt_ms
               << L"ms  max free-mu~ on a node=" << a.max_free_mu
               << L"  largest realized free-mu~={" << a.largest_desc << L"}="
               << a.largest_realized_gb << L" GB\n               replay: cache "
               << L"working_set_hwmark=" << a.hwmark_gb
               << L" GB  largest single transient(result=)="
               << a.max_single_result_gb << L" GB"
               << (a.replay_threw ? L"  [replay threw]" : L"") << L"\n";
    return a;
  };

  auto peak_first = analyze(ObjectiveFunction::DenseSpaceTimeBatched);
  auto perf_first = analyze(ObjectiveFunction::DenseTimeSpaceBatched);

  auto report = [](wchar_t const* tag, Analysis const& a) {
    std::wcerr << tag << L": 4-PAO node formed = "
               << (a.max_free_mu >= 4 ? L"YES" : L"NO") << L" (max free mu~="
               << a.max_free_mu << L")\n    DP-model largest realized free-mu~="
               << a.largest_realized_gb << L" GB {" << a.largest_desc << L"}\n"
               << L"    replay cache working_set_hwmark=" << a.hwmark_gb
               << L" GB  largest single transient(result=)="
               << a.max_single_result_gb << L" GB\n";
  };
  std::wcerr << L"\n=== [dryrun-perf] VERDICT (C60 giant, index 38) ===\n";
  report(L"peak-first (DenseSpaceTimeBatched)", peak_first);
  report(L"perf-first (DenseTimeSpaceBatched)", perf_first);
  // (b) The cache-vs-transient accounting gap: print the full trace line of the
  // largest single materialized result for each objective. result= (the
  // transient) and hw= (the cache watermark at that op) appear together, so a
  // large result= with a small hw= means the transient is NOT folded into the
  // cache high watermark (it is a streamed/batched result, never a cached
  // co-resident entry) -- i.e. working_set_hwmark tracks CACHED residency, not
  // the transient peak. The two together bound the true peak memory.
  std::wcerr
      << L"\n--- largest transient trace line (result= with its hw=) ---\n";
  {
    std::string const& pl = peak_first.largest_result_line;
    std::string const& tl = perf_first.largest_result_line;
    std::wcerr << L"peak-first: " << std::wstring(pl.begin(), pl.end()) << L"\n"
               << L"perf-first: " << std::wstring(tl.begin(), tl.end())
               << L"\n";
  }

  // INTERPRETATION -- why the outer working_set_hwmark is NOT the whole peak,
  // and why the replay sizes are not directly comparable to the DP-model sizes:
  //
  // (1) Batched transients are invisible to the OUTER cache accessor. The
  //     batched custom evaluator runs each batch against a SEPARATE scratch
  //     cache (detail::make_batched_scratch, eval.hpp:1386), so the per-op hw=
  //     field on the peak-first giant reaches ~38.9 GB (the 34 GB 4-PAO result
  //     + operands) INSIDE that scratch, while the outer
  //     cache.working_set_hwmark() reports only ~0.2 GB (outer-scope ops). So
  //     the outer accessor tracks CACHED (persistent, cross-batch) residency,
  //     not the batched-inner transient peak. Reading the whole peak requires
  //     folding each scratch cache's hwmark into a global accumulator.
  //
  // (2) Twin-PNO composite Results are mis-sized by the runtime backend. The
  //     perf-first largest transient = 358.47 GB = 120^2 * 42^4 * 8 EXACTLY
  //     (occ^2 * PNO^4): the DryRun Result's size_in_bytes() sizes the
  //     twin-PNO result R{a<i,i>,a<i,i>;i,i} as the naive product, NOT the
  //     power-mean moment (occ^2 * PNO^2 ~ 89 GB) the DP cost model uses via
  //     inner_pow. So the perf-first replay hwmark (~358 GB) is a sizing
  //     artifact, ~4x the true ~89 GB; the DP-model realized size (89 GB) is
  //     the trustworthy number. A faithful cached-watermark predictor needs the
  //     runtime Result sizing to be moment-aware too.
  //
  // Bottom line: the moment-aware DP peak model is the reliable predictor
  // today; the runtime replay hwmark is not, until (1) scratch hwmarks
  // propagate and (2) the Result sizing is moment-aware.
  std::wcerr
      << L"\n--- INTERPRETATION ---\n"
      << L"(1) peak-first outer hwmark (" << peak_first.hwmark_gb
      << L" GB) << its batched-inner transient ("
      << peak_first.max_single_result_gb
      << L" GB, hw= inside a make_batched_scratch cache): the OUTER accessor "
         L"misses batched transients.\n"
      << L"(2) perf-first transient " << perf_first.max_single_result_gb
      << L" GB = 120^2*42^4*8 (occ^2*PNO^4) = the twin-PNO size_in_bytes() "
         L"artifact; the moment-aware DP size ("
      << perf_first.largest_realized_gb
      << L" GB) is the real number. Runtime Result sizing is NOT "
         L"moment-aware.\n"
      << L"=> DP-model peak is the reliable predictor; the replay hwmark is "
         L"not "
         L"(yet).\n";

  // Peak-first forms the fully-sliceable 4-PAO AO integral (the C60 pathology).
  CHECK(peak_first.max_free_mu >= 4);
  // Perf-first, being flops-primary, must never form it.
  CHECK(perf_first.max_free_mu < 4);
}

// Whole-residual dry-run trace: optimize+binarize EVERY summand of the
// post-transform PNO-CCSD doubles residual, then replay them all through the
// real eval loop and ONE SHARED CacheManager (so cross-term CSE / persistent
// caching is exercised, as in MPQC's process_for_evaluation), with the per-op
// trace ("Eval | ... result= ... hw= ...") written to a FILE. Hidden ([.])
// because the batched DP runs per summand and the trace can be hundreds of MB;
// select it explicitly:
//   ./unit_tests-sequant "[dryrun-trace]"
// Env knobs (all optional):
//   DRYRUN_OBJ         dense_time_space (default, perf-first) |
//   dense_space_time DRYRUN_TRACE_FILE  output path (default
//   /tmp/dryrun_residual_trace.txt) DRYRUN_MAX_TERMS   cap #summands (default 0
//   = all) -- for a quick smoke run
TEST_CASE("dryrun whole-residual trace (PNO-CCSD doubles)",
          "[.][dryrun-trace]") {
  auto ctx = get_default_context().clone();
  ctx.set_first_dummy_index_ordinal(1000000);
  auto isr = ctx.mutable_index_space_registry();
  REQUIRE(isr != nullptr);
  sequant::mbpt::add_pao_spaces(isr);  // mu~
  sequant::mbpt::add_df_spaces(isr);   // K
  auto ctx_resetter = set_scoped_default_context(std::move(ctx));

  auto const body = slurp(std::string(SEQUANT_UNIT_TESTS_SOURCE_DIR) +
                          "/data/csv_ccsd_doubles_residual_df.txt");
  REQUIRE(!body.empty());
  std::string line = body;
  if (auto nl = line.find('\n'); nl != std::string::npos)
    line = line.substr(0, nl);
  auto expr = deserialize<ExprPtr>(line);
  REQUIRE(static_cast<bool>(expr));
  REQUIRE(expr->is<Sum>());
  auto flatten_product = [](ExprPtr const& e) -> ExprPtr {
    if (!e->is<Product>()) return e;
    auto const& p = e->as<Product>();
    return ex<Product>(p.scalar(), p.factors(), Product::Flatten::Yes);
  };
  std::vector<ExprPtr> terms;
  for (auto const& s : expr->as<Sum>().summands())
    terms.push_back(flatten_product(s));
  REQUIRE(!terms.empty());

  // Objective (perf-first default; dense_space_time / dense_peak_size = peak).
  char const* const obj_env = std::getenv("DRYRUN_OBJ");
  std::string const obj_s = obj_env ? obj_env : "dense_time_space";
  ObjectiveFunction const obj =
      (obj_s == "dense_space_time" || obj_s == "dense_peak_size")
          ? ObjectiveFunction::DenseSpaceTimeBatched
          : ObjectiveFunction::DenseTimeSpaceBatched;
  char const* const max_env = std::getenv("DRYRUN_MAX_TERMS");
  std::size_t const max_terms =
      max_env ? static_cast<std::size_t>(std::atoll(max_env)) : 0u;
  std::size_t const n_terms =
      (max_terms > 0 && max_terms < terms.size()) ? max_terms : terms.size();

  // FAITHFUL real C60 config (614336 job log), identical to [dryrun-eval].
  auto regime = df_regime(1800u, 4320u, 120u, 42.0, 310.0);
  auto cm = std::make_shared<CostModel const>(regime);
  sequant::BatchPolicy policy;
  policy.is_batchable_index = is_df_batchable;
  policy.batch_target_size = [](Index const& ix) -> std::size_t {
    return ix.space().base_key() == L"μ̃" ? std::size_t{256} : std::size_t{72};
  };
  policy.is_volatile_leaf = [](Tensor const& t) { return t.label() == L"t"; };
  policy.accumulation_factor = 1.0;
  policy.peak_threshold = 40e9;

  // Phase 1: optimize + binarize every summand; collect the DryRun trees so one
  // shared cache can count cross-term repeats.
  std::vector<EvalNodeDryRun> nodes;
  nodes.reserve(n_terms);
  for (std::size_t ti = 0; ti < n_terms; ++ti) {
    std::cerr << "[dryrun-trace] optimizing summand " << (ti + 1) << "/"
              << n_terms << " ..." << std::flush;
    auto axes_map = std::make_shared<std::unordered_map<
        Expr const*, container::vector<container::svector<Index>>>>();
    OptimizeOptions opts;
    opts.objective_function = obj;
    opts.idx_to_extent = regime.idx_to_extent();
    opts.inner_pow = regime.inner_pow_fn();
    opts.batch_policy = policy;
    opts.volatile_weight = 20.0;
    opts.roofline.machine_balance = 200.0;
    opts.roofline.fast_mem_elems = 1000000.0;
    opts.term_batch_axes = axes_map;
    auto t0 = std::chrono::steady_clock::now();
    auto optimized = optimize(terms[ti], opts);
    auto const ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                        std::chrono::steady_clock::now() - t0)
                        .count();
    std::cerr << " " << ms << "ms\n";
    if (!static_cast<bool>(optimized)) continue;
    auto it = axes_map->find(optimized.get());
    container::vector<container::svector<Index>> node_axes;
    if (it != axes_map->end()) node_axes = it->second;
    BinarizationOptions bopts;
    bopts.node_batch_axes = node_axes;
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    nodes.push_back(binarize<EvalExprDryRun>(optimized, {}, bopts));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
  }
  REQUIRE(!nodes.empty());

  // Phase 2: one shared cache over ALL trees (cross-term CSE), the real batched
  // evaluator, trace -> file. The cache is NOT reset between summands, so
  // working_set_hwmark() reflects cross-term cached residency (subject to the
  // [dryrun-perf] caveats: batched-inner transients and twin-PNO sizing).
  auto cache = sequant::cache_manager(nodes);
  cache.set_custom_evaluator(
      sequant::make_evaluator(policy, DryRunLeafEvaluator{cm}));

  char const* const file_env = std::getenv("DRYRUN_TRACE_FILE");
  std::string const trace_path =
      file_env ? file_env : "/tmp/dryrun_residual_trace.txt";
  std::ofstream trace_file(trace_path);
  REQUIRE(trace_file.is_open());
  auto& logger = Logger::instance();
  auto const prev_level = logger.eval.level;
  auto* const prev_stream = logger.eval.stream;
  logger.eval.level = 2;  // gate log::printing() on
  logger.eval.stream = &trace_file;

  std::size_t n_ok = 0;
  for (std::size_t i = 0; i < nodes.size(); ++i) {
    try {
      (void)sequant::evaluate<Trace::On>(nodes[i], DryRunLeafEvaluator{cm},
                                         cache);
      ++n_ok;
    } catch (std::exception const& e) {
      trace_file << "  [summand " << i << " threw: " << e.what() << "]\n";
    }
  }

  logger.eval.level = prev_level;
  logger.eval.stream = prev_stream;
  trace_file.flush();
  auto const hwmark_gb = double(cache.working_set_hwmark()) / 1e9;

  std::wcerr << L"\n=== [dryrun-trace] whole-residual trace written ===\n"
             << L"objective    = "
             << (obj == ObjectiveFunction::DenseTimeSpaceBatched
                     ? L"dense_time_space (perf-first)"
                     : L"dense_space_time (peak-first)")
             << L"\nsummands     = " << n_ok << L"/" << nodes.size()
             << L" evaluated\ncache hwmark = " << hwmark_gb
             << L" GB (cross-term cached residency; see [dryrun-perf] caveats)"
             << L"\ntrace file   = "
             << std::wstring(trace_path.begin(), trace_path.end()) << L"\n";
  SUCCEED();
}
