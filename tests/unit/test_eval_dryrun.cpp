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

#include <SeQuant/core/eval/backends/dryrun/size_regime.hpp>

#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/expressions/tensor.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/optimize/optimize.hpp>
#include <SeQuant/core/optimize/options.hpp>
#include <SeQuant/core/optimize/single_term_detail.hpp>
#include <SeQuant/core/runtime.hpp>
#include <SeQuant/core/space.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
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
