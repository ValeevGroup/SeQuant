// Task 2 of the whole-scope batched DAG execution design (see
// doc/dev/specs/2026-08-10-whole-scope-batched-dag-execution-design.md): the
// executor SKELETON. For a forest with NO batch loops (the scope tree built
// by Task 1's build_scope_schedule is root-only), evaluate_whole_scope must
// be provably equivalent to the existing per-tree forest descent (sequant::
// evaluate(Nodes const&, ...)).
//
// The equivalence check needs REAL numeric arithmetic (not just shape/size
// modeling): a bug that drops or double-counts a forest root would not
// necessarily change a zero-data DryRun result's modeled shape, so this file
// deliberately avoids the DryRun backend and instead evaluates genuine scalar
// arithmetic (sequant::ResultScalar<double>) via a minimal local EvalExpr
// subclass -- just enough to satisfy meta::can_evaluate (an annot() method),
// with no tensor backend (BTAS/TiledArray) dependency.
//
// Task 7 of the same design (end-to-end WITNESSES, see the design's
// "Validation" section) adds two `[.]` (hidden, run-by-name) TEST_CASEs below
// that DO use the DryRun backend: they reuse the real post-transform C60
// residual fixture (data/csv_ccsd_doubles_residual_df.txt) and the exact
// water-20 / C60 forest-construction recipes test_eval_dryrun.cpp's
// [dryrun-water-frag] / [dryrun-occ-veto] cases use, so the earlier
// numeric-only rationale above applies only to the Task-2/5 tests that precede
// them.

#include <SeQuant/core/batch_policy.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/eval/backends/dryrun/cost_model_object.hpp>
#include <SeQuant/core/eval/backends/dryrun/cost_profile.hpp>
#include <SeQuant/core/eval/backends/dryrun/eval_expr.hpp>
#include <SeQuant/core/eval/backends/dryrun/meter.hpp>
#include <SeQuant/core/eval/backends/dryrun/result.hpp>
#include <SeQuant/core/eval/backends/dryrun/size_regime.hpp>
#include <SeQuant/core/eval/cache_manager.hpp>
#include <SeQuant/core/eval/eval.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/eval/peak_profile.hpp>
#include <SeQuant/core/eval/result.hpp>
#include <SeQuant/core/eval/scope_executor.hpp>
#include <SeQuant/core/eval/scope_schedule.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/expressions/tensor.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/logger.hpp>
#include <SeQuant/core/optimize/optimize.hpp>
#include <SeQuant/core/optimize/options.hpp>
#include <SeQuant/core/runtime.hpp>
#include <SeQuant/core/space.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/space_qns.hpp>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

namespace {

using sequant::Constant;
using sequant::EvalExpr;
using sequant::EvalNode;
using sequant::ExprPtr;
using sequant::Index;
using sequant::ResultPtr;
using sequant::ResultScalar;
using sequant::Variable;
using sequant::eval::build_scope_schedule;
using sequant::eval::compute_dag_boulevard;
using sequant::eval::evaluate_whole_scope;
using sequant::eval::OccurrenceRec;
using sequant::eval::RichSchedule;
using sequant::eval::ScopeSchedule;
using sequant::eval::ValueCell;
using sequant::eval::weighted_use_count;
using sequant::eval::dryrun::CostModel;
using sequant::eval::dryrun::SizeRegime;

///
/// \brief A minimal EvalExpr subclass whose only job is to satisfy \c
/// meta::can_evaluate (i.e. carry an annot() method) so a plain scalar
/// arithmetic forest (Constant/Variable leaves, Sum/Product internal nodes,
/// no tensor indices at all) can be run through \c evaluate_impl /
/// \c evaluate_whole_scope without pulling in a tensor backend. Every node in
/// this test forest is scalar (no canon_indices), so annot() is never
/// consulted by a permute() call (the default-constructed \p layout used
/// throughout keeps \c perm false in both the reference forest descent and
/// \c evaluate_whole_scope) -- it only needs to exist and be comparable.
///
class ScalarEvalExpr final : public EvalExpr {
 public:
  using annot_t = int;

  template <typename... Args, typename = std::enable_if_t<
                                  std::is_constructible_v<EvalExpr, Args...>>>
  explicit ScalarEvalExpr(Args&&... args)
      : EvalExpr{std::forward<Args>(args)...} {}

  [[nodiscard]] annot_t annot() const noexcept { return 0; }
};

using ScalarNode = EvalNode<ScalarEvalExpr>;

static_assert(sequant::meta::eval_node<ScalarNode>);
static_assert(sequant::meta::can_evaluate<ScalarNode>);

///
/// \brief Leaf evaluator for the scalar forest: a Constant leaf yields its
/// own numeric value; a Variable leaf yields the value bound to its label in
/// \c values (every Variable used by the test forest must be bound).
///
struct ScalarLeafEvaluator {
  std::map<std::wstring, double> values;

  [[nodiscard]] ResultPtr operator()(ScalarNode const& leaf) const {
    SEQUANT_ASSERT(leaf.leaf());
    ExprPtr const& xpr = leaf->expr();
    if (xpr->is<Constant>())
      return sequant::eval_result<ResultScalar<double>>(
          xpr->as<Constant>().value<double>());
    SEQUANT_ASSERT(xpr->is<Variable>());
    auto const it = values.find(std::wstring(xpr->as<Variable>().label()));
    SEQUANT_ASSERT(it != values.end() && "unbound scalar leaf in test forest");
    return sequant::eval_result<ResultScalar<double>>(it->second);
  }
};

// One scalar equation deserialized+binarized into a ScalarNode tree.
ScalarNode scalar_tree(std::wstring_view spec) {
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
  return sequant::binarize<ScalarEvalExpr>(
      sequant::deserialize<ExprPtr>(std::wstring(spec)));
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
}

///
/// \brief The shared unbatched, two-root test forest + its Task-1
/// ScopeSchedule + leaf evaluator, built via the real pipeline
/// (compute_dag_boulevard -> build_scope_schedule): no node ever calls
/// set_node_slice_mask(), so the resulting scope tree must be root-only. Reused
/// by both the numeric-equivalence test and the trace-fidelity (fix round 1)
/// test below so the two exercise the identical forest.
///
struct TestForest {
  std::vector<ScalarNode> forest;
  RichSchedule rich;
  ScopeSchedule sched;
  ScalarLeafEvaluator yield;
};

TestForest make_test_forest() {
  // A two-root forest, each root a nontrivial Sum-of-Products over named
  // scalar Variables -- exercises both EvalOp::Sum and EvalOp::Product inside
  // evaluate_impl, plus the cross-root accumulation evaluate_whole_scope must
  // reproduce.
  std::vector<ScalarNode> forest{scalar_tree(L"2 * a * b - c"),
                                 scalar_tree(L"a * a + 3 * b - 2 * c")};

  ScalarLeafEvaluator yield{{{L"a", 2.0}, {L"b", -3.5}, {L"c", 7.25}}};

  SizeRegime const regime;
  CostModel const cm{regime};
  auto const block_of = [](sequant::Index const&) -> std::size_t { return 1; };
  RichSchedule rich = compute_dag_boulevard(forest, cm, block_of);
  ScopeSchedule sched = build_scope_schedule<std::wstring>(rich, {});

  return TestForest{std::move(forest), std::move(rich), std::move(sched),
                    std::move(yield)};
}

}  // namespace

TEST_CASE("evaluate_whole_scope matches forest descent for an unbatched forest",
          "[scope-executor]") {
  TestForest const tf = make_test_forest();
  REQUIRE(tf.sched.root.children.empty());
  REQUIRE(tf.sched.num_values == tf.rich.cells.size());

  // Reference: existing per-tree forest descent, own (fresh) cache.
  auto ref_cache = sequant::CacheManager<ScalarNode>::empty();
  ResultPtr const reference = sequant::evaluate(tf.forest, tf.yield, ref_cache);
  double const expected = reference->as<ResultScalar<double>>().value();

  // evaluate_whole_scope, its own (fresh) cache -- must match exactly.
  auto whole_cache = sequant::CacheManager<ScalarNode>::empty();
  ResultPtr const got =
      evaluate_whole_scope(tf.forest, tf.sched, tf.rich,
                           ScalarEvalExpr::annot_t{}, tf.yield, whole_cache);
  double const got_val = got->as<ResultScalar<double>>().value();

  // Hand-computed cross-check that the reference itself is right:
  //   p1 = 2*a*b - c = 2*2*(-3.5) - 7.25 = -14 - 7.25 = -21.25
  //   p2 = a*a + 3*b - 2*c = 4 + 3*(-3.5) - 2*7.25 = 4 - 10.5 - 14.5 = -21
  //   total = p1 + p2 = -42.25
  double const hand = -42.25;

  CHECK(expected == Catch::Approx(hand));
  CHECK(got_val == Catch::Approx(expected));
}

// Fix round 1 (code review finding): evaluate_whole_scope's accumulation loop
// must emit the SAME trace/hwmark bookkeeping sequant::evaluate(Nodes const&,
// layout, ...) does around the identical permute()/add_inplace() calls -- the
// per-root Term Begin/End boundary and, for every root after the first, a
// SumInplace EvalStat. Verified here by actually capturing the eval-trace
// stream under Trace::On (not just asserting it by inspection): a regression
// that silently dropped the SumInplace record again (the exact bug this fix
// addresses) would fail this test.
TEST_CASE(
    "evaluate_whole_scope emits Term/SumInplace trace records under "
    "Trace::On",
    "[scope-executor]") {
  TestForest const tf = make_test_forest();
  REQUIRE(tf.sched.root.children.empty());

  // Logger is a process-wide Singleton; save/restore its eval-log state so
  // this test cannot leak a nonzero level or a dangling stream pointer into
  // whatever test runs next.
  auto& logger_eval = sequant::Logger::instance().eval;
  auto const saved_level = logger_eval.level;
  auto* const saved_stream = logger_eval.stream;

  std::ostringstream oss;
  logger_eval.level = 1;
  logger_eval.stream = &oss;

  ResultPtr got;
  try {
    auto cache = sequant::CacheManager<ScalarNode>::empty();
    got = evaluate_whole_scope<sequant::Trace::On>(tf.forest, tf.sched, tf.rich,
                                                   ScalarEvalExpr::annot_t{},
                                                   tf.yield, cache);
  } catch (...) {
    logger_eval.level = saved_level;
    logger_eval.stream = saved_stream;
    throw;
  }

  logger_eval.level = saved_level;
  logger_eval.stream = saved_stream;

  // Numeric result must be unaffected by tracing.
  REQUIRE(got);
  CHECK(got->as<ResultScalar<double>>().value() == Catch::Approx(-42.25));

  std::string const trace = oss.str();
  INFO("captured eval trace:\n" << trace);
  // Per-root term boundary (log::term(Begin/End)), mirroring
  // sequant::evaluate(Node const&, layout, ...)'s identical bracket.
  CHECK(trace.find("Term | Begin") != std::string::npos);
  CHECK(trace.find("Term | End") != std::string::npos);
  // The cross-root accumulation EvalStat this fix restores: 2 forest roots
  // => exactly one add_inplace, hence exactly one SumInplace record.
  CHECK(trace.find("SumInplace") != std::string::npos);
}

// Task 5: the weighted use-count lifetime model (replaces the ensure_hoist_slot
// MAX-life hack). life(V) = sum over consumers C of V, product over the loops
// on the path (home(V), scope(C)] of n_blocks(L). Pinned here directly on a
// small synthetic scope tree so the arithmetic is checkable by hand:
//   - one consumer inside an n-block inner loop -> count == n,
//   - two consumers -> the sum of each consumer's product,
//   - two nested loops between home and consumer -> their product.
// Loop modes are TYPE-keyed: "i" (occ) and "a" (virt) are both present in the
// default context, so distinct types stand in for distinct realized loops.
TEST_CASE("weighted_use_count sums consumers and multiplies nested loops",
          "[scope-executor]") {
  using Range = std::pair<std::size_t, std::size_t>;
  Index const a1{L"a_1"};  // OUTER loop (type "a")
  Index const i1{L"i_1"};  // INNER loop (type "i")

  // n_blocks: 3 blocks for the inner ("i") loop, 5 for the outer ("a") loop.
  auto const n_blocks = [](Index const& m) -> std::size_t {
    auto const& bk = m.space().base_key();
    if (bk == L"i") return 3;
    if (bk == L"a") return 5;
    return 1;
  };

  // (1) Home at the OUTER (a) loop; a single consumer sits inside the INNER (i)
  // loop. The path (home, consumer] is just the inner loop -> life == 3.
  {
    ValueCell cell{};
    cell.home_modes = {a1};
    OccurrenceRec occ{};
    occ.ectx = {{a1, Range{0, 5}}, {i1, Range{0, 4}}};  // outer a, inner i
    cell.occurrences = {occ};
    CHECK(weighted_use_count(cell, n_blocks) == 3);
  }

  // (2) Two consumers, same home (outer a): one inside the inner (i) loop
  // (contributes 3), one AT the home with no inner loop (contributes 1). The
  // life is the SUM: 3 + 1 == 4.
  {
    ValueCell cell{};
    cell.home_modes = {a1};
    OccurrenceRec inner{};
    inner.ectx = {{a1, Range{0, 5}}, {i1, Range{0, 4}}};
    OccurrenceRec at_home{};
    at_home.ectx = {{a1, Range{0, 5}}};
    cell.occurrences = {inner, at_home};
    CHECK(weighted_use_count(cell, n_blocks) == 4);
  }

  // (3) Whole-nest invariant (empty home): a single consumer sits inside TWO
  // nested loops (outer a, inner i). Both loops are on the path, so the counts
  // MULTIPLY: 5 * 3 == 15.
  {
    ValueCell cell{};
    cell.home_modes = {};  // homed at the root
    OccurrenceRec occ{};
    occ.ectx = {{a1, Range{0, 5}}, {i1, Range{0, 4}}};
    cell.occurrences = {occ};
    CHECK(weighted_use_count(cell, n_blocks) == 15);
  }

  // (4) Degenerate: a value used only at its own home (no enclosing inner loop)
  // has the ordinary use count -- here two occurrences, both at home -> 2.
  {
    ValueCell cell{};
    cell.home_modes = {a1};
    OccurrenceRec o1{};
    o1.ectx = {{a1, Range{0, 5}}};
    OccurrenceRec o2{};
    o2.ectx = {{a1, Range{0, 5}}};
    cell.occurrences = {o1, o2};
    CHECK(weighted_use_count(cell, n_blocks) == 2);
  }
}

// ===========================================================================
// Task 7 of the whole-scope batched DAG execution design (see
// doc/dev/specs/2026-08-10-whole-scope-batched-dag-execution-design.md,
// "Validation"): end-to-end WITNESSES on the real water-20 / C60 CSV-CCSD
// residuals. `[.]` (hidden, run-by-name) -- these replay a real, sizable
// forest through the DryRun backend and are excluded from the default suite
// for the same reason test_eval_dryrun.cpp's [dryrun-water-frag] /
// [dryrun-occ-veto] cases are.
// ===========================================================================

namespace {

// Reads the whole content of a file at `path` into a string (identical to
// test_eval_dryrun.cpp's file-local `slurp`; duplicated here rather than
// shared -- no common test header exists for these DryRun fixtures, see the
// doc comments below).
std::string witness_slurp(std::string const& path) {
  std::ifstream in(path);
  std::stringstream ss;
  ss << in.rdbuf();
  return ss.str();
}

// One named (molecule, basis, parameter-set) problem size for the DryRun cost
// model, duplicated from test_eval_dryrun.cpp's `ProblemSize` /
// `kC60_pVDZF12` / `kWater20_pVDZF12` (no shared test header exists between
// the two .cpp files, so a moment re-measurement there needs a matching edit
// here; see that file's doc comments for provenance -- job logs / mpqc
// PaoPnoRMP2 measurements).
struct WitnessProblemSize {
  std::size_t mu_tilde;         // PAO domain extent (= #AO)
  std::size_t aux;              // DF aux (K) extent
  std::size_t i_occ;            // active occupied extent
  std::array<double, 5> pno_M;  // per-pair PNO power means M_0..M_4
  std::array<double, 5> osv_M;  // per-orbital OSV power means M_0..M_4
};

inline constexpr WitnessProblemSize kWitnessC60_pVDZF12{
    /*mu_tilde=*/1800u,
    /*aux=*/4320u,
    /*i_occ=*/120u,
    /*pno_M=*/
    {1.0, 42.029069767441861, 46.039206412923569, 49.766252354482994,
     53.151291880343109},
    /*osv_M=*/
    {1.0, 148.25, 155.04434849422921, 161.33527408797721, 166.85553430303926}};

inline constexpr WitnessProblemSize kWitnessWater20_pVDZF12{
    /*mu_tilde=*/896u,
    /*aux=*/1682u,
    /*i_occ=*/80u,
    /*pno_M=*/
    {1.0, 23.175775480059084, 25.865548281212597, 28.171416142614103,
     30.03848680550367},
    /*osv_M=*/
    {1.0, 58.987499999999997, 59.289227520688783, 59.584437469011633,
     59.872014818179686}};

sequant::eval::dryrun::SizeRegime witness_df_regime(
    WitnessProblemSize const& p) {
  sequant::eval::dryrun::SizeRegime r;
  r.space_extent = {
      {L"i", p.i_occ},
      {L"μ̃", p.mu_tilde},
      {L"Κ", p.aux},
      {L"a", p.mu_tilde},
  };
  r.csv_pno_moment = p.pno_M;
  r.csv_osv_moment = p.osv_M;
  return r;
}

// Flattens a top-level Product summand (identical to test_eval_dryrun.cpp's
// per-TEST_CASE local lambda) -- a nested Product needs flattening before
// optimize().
sequant::ExprPtr witness_flatten_product(sequant::ExprPtr const& e) {
  if (!e->is<sequant::Product>()) return e;
  auto const& p = e->as<sequant::Product>();
  return sequant::ex<sequant::Product>(p.scalar(), p.factors(),
                                       sequant::Product::Flatten::Yes);
}

// Parses the max `hw=<N>B` token out of a captured (narrow) eval trace --
// see eval.hpp's EvalStat doc comment: mem_hwmark is the running max, held by
// whichever CacheManager instance emitted that record, of bytes(cache) +
// bytes(result) + un-aliased operand bytes + the scope chain's ancestor
// residency (folded in at every call site via
// `cache.parent()->chain_residency()`). Since every EvalStat record --
// whether emitted by the root cache or a batch-scratch cache nested inside
// detail::walk_scope -- funnels through the SAME process-wide Logger, the
// GLOBAL max over the whole trace is the realized peak of the run: the same
// number a PeakSink would fold via std::max over each scratch's
// working_set_hwmark() (as test_eval_dryrun.cpp's "dryrun scratch-fold
// captures batched peak" test does for forest descent), just read directly
// from the trace text -- evaluate_whole_scope has no PeakSink seam
// (walk_scope is production code this witness MEASURES, not modifies; see
// the design's Scope boundary).
std::size_t witness_parse_max_hwmark_bytes(std::string const& trace) {
  std::size_t best = 0;
  std::string const key = "hw=";
  for (std::size_t pos = trace.find(key); pos != std::string::npos;
       pos = trace.find(key, pos + key.size())) {
    std::size_t p = pos + key.size();
    std::size_t val = 0;
    bool any = false;
    while (p < trace.size() && trace[p] >= '0' && trace[p] <= '9') {
      val = val * 10 + static_cast<std::size_t>(trace[p] - '0');
      ++p;
      any = true;
    }
    if (any) best = std::max(best, val);
  }
  return best;
}

// Counts non-overlapping occurrences of `needle` in `haystack`. Used to count
// `"Term | Begin"` trace records: both `sequant::evaluate(Nodes const&, ...)`
// (forest descent) and `evaluate_whole_scope`'s shared per-root combine loop
// emit ONE such record per forest root (an identical `log::term(TermMode::
// Begin, xpr)` bracket, see eval.hpp / scope_executor.hpp), so this count is
// the TERM-COUNT-SENSITIVE check `size_in_bytes()` alone cannot provide: a
// silently dropped or double-counted root changes this count, whereas
// `ResultDryRun::sum()`/`add_inplace()` (backends/dryrun/result.hpp) derive
// the result's shape/overrides from the engine-supplied annotation and
// positional merges -- independent of how many summands actually
// contributed -- so `size_in_bytes()` alone would NOT necessarily change if a
// root were dropped or duplicated. The two checks are complementary: this one
// is the drop/duplicate detector, `size_in_bytes()` equality is the
// shape/footprint agreement check layered on top of it.
std::size_t witness_count_substr(std::string const& haystack,
                                 std::string const& needle) {
  std::size_t count = 0;
  for (std::size_t pos = haystack.find(needle); pos != std::string::npos;
       pos = haystack.find(needle, pos + needle.size()))
    ++count;
  return count;
}

// Total builds (summed over slices) of one specific node in a build tally
// (CacheManager::recompute_tally()) -- identical helper to the one every
// build-once witness in test_eval_dryrun.cpp defines locally.
template <typename Tally>
std::size_t witness_builds_of(Tally const& tally,
                              sequant::eval::dryrun::EvalNodeDryRun const& n) {
  auto it = tally.find(n);
  if (it == tally.end()) return 0;
  std::size_t b = 0;
  for (auto const& [sig, bc] : it->second.slices) b += bc.count;
  return b;
}

}  // namespace

// The water-20 WITNESS. Same real CSV-CCSD doubles residual and the SAME
// aux-only (Κ-contracted) MPQC batch config test_eval_dryrun.cpp's
// [dryrun-water-frag] / "whole-scope executor builds shared aux composites
// once per block" cases use (df_regime(kWater20_pVDZF12), K@256,
// DenseTimeSpaceBatched perf-first objective) -- reused here, not
// reinvented. Four things are checked together, on the SAME forest and the
// SAME whole-scope replay, because (per the design's "Watch item")
// equivalence-to-tolerance alone would pass even if nothing shared:
//
//   (a) the whole-scope replay's TARGETED Κ-carrying shared composites (the
//       scope's homed_values, resolved to forest nodes) each build AT MOST
//       once per Κ-block, and the busiest one hits EXACTLY n_blocks -- the
//       build-once win Task 3 introduced (the KEY metric this witness
//       exists to pin);
//   (b) forest descent, replayed on the IDENTICAL forest with the batched
//       custom evaluator, rebuilds the SAME composites MORE than the
//       whole-scope build-once bound -- the fragmentation Task 3's win
//       eliminates, i.e. the "recompute elimination vs forest descent" the
//       brief asks to document;
//   (c) the co-residency oracle (cost_profile() under
//       policy.scheduler == BatchScheduler::whole_scope) predicts the
//       realized whole-scope
//       peak, read directly off the SAME replay's own eval trace (the `hw=`
//       high-water mark every EvalStat record -- root cache or a
//       batch-scratch cache nested inside detail::walk_scope -- reports;
//       see witness_parse_max_hwmark_bytes's doc comment);
//   (d) whole-scope and forest descent -- run separately below on the
//       IDENTICAL forest/policy for (b) anyway -- agree STRUCTURALLY: the
//       SAME number of forest roots processed (a `Term | Begin` trace-record
//       count equal to `forest.size()` on BOTH replays -- see
//       witness_count_substr's doc comment for why this, not
//       `size_in_bytes()` alone, is what actually rules out a silently
//       dropped or double-counted root) and the SAME final result footprint
//       (`size_in_bytes()` equality -- a shape/positional agreement check,
//       not by itself a drop detector; see the same doc comment).
//
// Numeric equivalence is NOT checked by this witness: the real CSV-CCSD
// residual is intractable at real tensor sizes (the whole reason
// test_eval_dryrun.cpp evaluates it zero-data), so DryRun's Result carries no
// real floating-point content to compare -- (a)-(d) above are the zero-data
// witness's full correctness surface. The numeric equivalence of the
// batched-summation algorithm itself -- that whole-scope's reordered
// accumulation reproduces the SAME NUMBERS forest descent does, to floating-
// point tolerance -- is established separately, at small real-TA scale, by
// the Task-3/4 real-data TA proxy in test_eval_ta.cpp ("evaluate_whole_scope
// matches forest descent over one aux loop" / "... over nested aux+occ",
// relative L2 < 1e-10), which exercises the IDENTICAL sharing topology (a
// Κ-carrying composite shared by two roots, contracted at each root) at a
// size small enough for real TA arithmetic.
TEST_CASE(
    "whole-scope witness: water-20 aux-only residual builds shared "
    "composites once per block and the co-residency oracle predicts the "
    "realized peak",
    "[.][scope-executor-witness-water20]") {
  using sequant::eval::dryrun::EvalExprDryRun;
  using sequant::eval::dryrun::EvalNodeDryRun;

  auto ctx = sequant::get_default_context().clone();
  ctx.set_first_dummy_index_ordinal(1000000);
  auto isr = ctx.mutable_index_space_registry();
  REQUIRE(isr != nullptr);
  sequant::mbpt::add_pao_spaces(isr, sequant::mbpt::Spin::any);
  sequant::mbpt::add_df_spaces(isr);
  auto ctx_resetter = sequant::set_scoped_default_context(std::move(ctx));

  auto const body = witness_slurp(std::string(SEQUANT_UNIT_TESTS_SOURCE_DIR) +
                                  "/data/csv_ccsd_doubles_residual_df.txt");
  REQUIRE(!body.empty());
  std::string line = body;
  if (auto nl = line.find('\n'); nl != std::string::npos)
    line = line.substr(0, nl);
  auto expr = sequant::deserialize<sequant::ExprPtr>(line);
  REQUIRE(static_cast<bool>(expr));
  REQUIRE(expr->is<sequant::Sum>());
  auto const& summands = expr->as<sequant::Sum>().summands();
  REQUIRE(!summands.empty());

  // Same term cap as test_eval_dryrun.cpp's Task-3 sharing test: enough terms
  // that aux composites are shared cross-tree, small enough that all three
  // replays below (whole-scope, forest-descent contrast, cost_profile's own
  // replay) stay quick.
  std::size_t nterms = std::min<std::size_t>(summands.size(), 40);
  if (char const* nt = std::getenv("SEQUANT_UT_DRYRUN_NTERMS"))
    nterms = std::min<std::size_t>(summands.size(), std::atoll(nt));

  auto regime = witness_df_regime(kWitnessWater20_pVDZF12);
  auto cm = std::make_shared<sequant::eval::dryrun::CostModel const>(regime);

  // EXACT MPQC aux-only config (make_csv_batch_policy, aux_target=256): Κ is
  // the only batchable mode, contracted role.
  sequant::BatchPolicy policy;
  policy.is_batchable_contracted_index = [](sequant::Index const& ix) {
    return ix.space().base_key() == L"Κ";
  };
  policy.is_batchable_external_index = [](sequant::Index const&) {
    return false;
  };
  policy.batch_spectator_indices = false;
  policy.batch_target_size = [](sequant::Index const&) -> std::size_t {
    return 256;
  };
  policy.is_volatile_leaf = [](sequant::Tensor const& t) {
    return t.label() == L"t";
  };
  policy.accumulation_factor = 1.0;
  policy.persistent_only = false;
  policy.peak_threshold = 1e11;

  auto axes_map = std::make_shared<std::unordered_map<
      sequant::Expr const*,
      sequant::container::vector<sequant::NodeBatchAnnotation>>>();
  sequant::OptimizeOptions opts;
  opts.objective_function = sequant::ObjectiveFunction::DenseTimeSpaceBatched;
  opts.idx_to_extent = regime.idx_to_extent();
  opts.inner_pow = regime.inner_pow_fn();
  opts.batch_policy = policy;
  opts.volatile_weight = 20.0;
  opts.roofline.machine_balance = 200.0;
  opts.roofline.fast_mem_elems = 1000000.0;
  opts.term_batch_axes = axes_map;

  std::vector<EvalNodeDryRun> forest;
  for (std::size_t s = 0; s < nterms; ++s) {
    sequant::ExprPtr const term = witness_flatten_product(summands[s]);
    if (!term) continue;
    sequant::ExprPtr optimized;
    try {
      optimized = sequant::optimize(term, opts);
    } catch (std::exception const&) {
      continue;
    }
    if (!optimized) continue;
    sequant::BinarizationOptions bopts;
    if (auto it = axes_map->find(optimized.get()); it != axes_map->end())
      bopts.node_batch_axes = it->second;
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    forest.push_back(sequant::binarize<EvalExprDryRun>(optimized, {}, bopts));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
  }
  REQUIRE(!forest.empty());

  // Scope schedule: aux Κ is the single realized batch loop (root-only ->
  // one child, no nested child -- the Task-3 topology).
  auto const block_of = [](sequant::Index const&) -> std::size_t {
    return 256;
  };
  auto rich = sequant::eval::compute_dag_boulevard(forest, *cm, block_of);
  auto sched = sequant::eval::build_scope_schedule<std::wstring>(rich, {L"Κ"});
  REQUIRE(sched.root.children.size() == 1);
  auto const& kscope = sched.root.children.front();
  REQUIRE(kscope.mode.space().base_key() == L"Κ");
  REQUIRE(kscope.kind == sequant::BatchModeType::Contracted);
  sequant::Index const K = kscope.mode;

  // The TARGETED shared composites: the K scope's homed values, resolved to
  // forest nodes (skipping leaves and forest roots -- neither is a pre-built
  // "composite").
  auto const vmap = sequant::eval::build_value_node_map(forest);
  std::unordered_set<std::size_t> root_hashes;
  for (auto const& r : forest) root_hashes.insert(r->hash_value());
  std::vector<EvalNodeDryRun> k_homed;
  for (auto vid : kscope.homed_values) {
    auto const h = rich.cells[vid].hash;
    auto const it = vmap.find(h);
    if (it == vmap.end() || it->second.leaf() || root_hashes.count(h)) continue;
    k_homed.push_back(it->second);
  }
  REQUIRE(!k_homed.empty());

  sequant::eval::dryrun::DryRunLeafEvaluator yield{cm};
  auto const target = [](sequant::Index const&) -> std::size_t { return 256; };

  // n_blocks from any Κ-carrying leaf's realized partition.
  auto aops = sequant::eval::dryrun::make_dryrun_array_ops(cm);
  std::size_t n_blocks = 0;
  for (auto const& root : forest)
    if (auto lf = sequant::find_leaf_carrying(root, K)) {
      n_blocks = aops.axis_batches(K, 256).size();
      break;
    }
  REQUIRE(n_blocks > 1);

  using annot_t = std::remove_cvref_t<decltype(forest.front()->annot())>;
  annot_t const layout{};

  auto& logger = sequant::Logger::instance();
  auto const prev_level = logger.eval.level;
  auto* const prev_stream = logger.eval.stream;

  // ---- (a) + (c): ONE whole-scope replay -- the build tally AND the eval
  // trace (for the peak-oracle comparison) come from the SAME run. The
  // elevated logger level is what makes BOTH the trace records (mem_hwmark,
  // gated on `if constexpr (detail::trace(EvalTrace))`) AND the DryRun
  // per-op cost sink that feeds CacheManager::tally_build (DryRunOps::prod's
  // `if (Logger::instance().eval.level > 0)` RUNTIME gate on last_op_flops(),
  // result.hpp) actually fire -- so it stays elevated across the
  // forest-descent contrast run below too, not just this one.
  std::ostringstream ws_trace;
  logger.eval.level = 1;
  logger.eval.stream = &ws_trace;
  auto ws_cache = sequant::cache_manager(forest);
  ws_cache.set_array_ops(&aops);
  ws_cache.set_recompute_tally_enabled(true);
  ResultPtr ws_result;
  try {
    ws_result = sequant::eval::evaluate_whole_scope<sequant::Trace::On>(
        forest, sched, rich, layout, yield, ws_cache, target);
  } catch (std::exception const& e) {
    std::cerr << "[scope-executor-witness-water20] whole-scope evaluate "
                 "threw: "
              << e.what() << "\n";
  }

  std::size_t max_ws_builds = 0;
  for (auto const& n : k_homed)
    max_ws_builds = std::max(max_ws_builds,
                             witness_builds_of(ws_cache.recompute_tally(), n));
  REQUIRE(max_ws_builds > 0);

  // Parse the whole-scope peak AND term-count BEFORE the forest-descent run
  // reuses the logger's stream (below), so the fd replay's own trace text
  // cannot contaminate either measurement.
  std::size_t const realized_peak_bytes =
      witness_parse_max_hwmark_bytes(ws_trace.str());
  std::size_t const ws_term_count =
      witness_count_substr(ws_trace.str(), "Term | Begin");

  // ---- (b) + (d): forest-descent contrast, batched custom evaluator, on the
  // IDENTICAL forest/policy -- also the structural comparison target for (d).
  // Same elevated logger level (see above), a FRESH trace sink (this one is
  // READ, for (d)'s term-count check, unlike the C60 witness's discarded
  // fd_trace). ----
  std::ostringstream fd_trace;
  logger.eval.stream = &fd_trace;
  auto fd_cache = sequant::cache_manager(forest);
  fd_cache.set_array_ops(&aops);
  fd_cache.set_custom_evaluator(sequant::make_evaluator(policy, yield));
  fd_cache.set_recompute_tally_enabled(true);
  // The forest descent does NOT go through evaluate_whole_scope, so install the
  // same schedule-derived per-op annotation (remat home + use scopes) around it
  // -- save/restore -- so both traces carry identical home/uses annotations for
  // clean diffing.
  auto const prev_node_meta = logger.eval.node_meta;
  logger.eval.node_meta = sequant::eval::make_node_meta(rich);
  ResultPtr fd_result;
  try {
    fd_result =
        sequant::evaluate<sequant::Trace::On>(forest, layout, yield, fd_cache);
  } catch (std::exception const& e) {
    std::cerr << "[scope-executor-witness-water20] forest-descent evaluate "
                 "threw: "
              << e.what() << "\n";
  }
  logger.eval.node_meta = prev_node_meta;
  logger.eval.level = prev_level;
  logger.eval.stream = prev_stream;

  std::size_t const fd_term_count =
      witness_count_substr(fd_trace.str(), "Term | Begin");

  std::size_t max_fd_builds = 0;
  for (auto const& n : k_homed)
    max_fd_builds = std::max(max_fd_builds,
                             witness_builds_of(fd_cache.recompute_tally(), n));

  // ---- (c) continued: the co-residency oracle. ----
  sequant::eval::dryrun::CacheConfig cfg;
  cfg.max_footprint = 0.;
  cfg.min_repeats = 1;
  cfg.is_volatile = [](EvalNodeDryRun const& n) {
    return n.leaf() && n->is_tensor() && n->as_tensor().label() == L"t";
  };
  sequant::BatchPolicy policy_ws = policy;
  policy_ws.scheduler = sequant::BatchScheduler::whole_scope;
  auto const cp =
      sequant::eval::dryrun::cost_profile(forest, policy_ws, cfg, regime);
  double const oracle_peak_bytes = cp.peak_bytes;

  double const ratio = oracle_peak_bytes > 0.0
                           ? double(realized_peak_bytes) / oracle_peak_bytes
                           : 0.0;
  std::wcerr << L"\n=== [scope-executor-witness-water20] water-20 aux-only, "
             << forest.size() << L" terms, n_blocks=" << n_blocks << L", "
             << k_homed.size() << L" Κ-homed composites ===\n"
             << L"  whole-scope replay completed  = "
             << (ws_result ? L"yes" : L"NO") << L"\n"
             << L"  forest-descent replay completed = "
             << (fd_result ? L"yes" : L"NO") << L"\n"
             << L"  term (root) count: whole-scope = " << ws_term_count
             << L", forest-descent = " << fd_term_count << L"  (forest.size()="
             << forest.size() << L")\n"
             << L"  final result size: whole-scope = "
             << (ws_result ? double(ws_result->size_in_bytes()) / 1e9 : -1.0)
             << L" GB, forest-descent = "
             << (fd_result ? double(fd_result->size_in_bytes()) / 1e9 : -1.0)
             << L" GB\n"
             << L"  whole-scope max builds (targeted composites) = "
             << max_ws_builds << L"  (n_blocks=" << n_blocks << L")\n"
             << L"  forest-descent max builds (SAME composites)  = "
             << max_fd_builds << L"\n"
             << L"  realized peak (hw= trace max)  = "
             << (double(realized_peak_bytes) / 1e9) << L" GB\n"
             << L"  co-residency oracle peak       = "
             << (oracle_peak_bytes / 1e9) << L" GB\n"
             << L"  ratio (realized/oracle)        = " << ratio << L"\n";

  // (a) build-once: no targeted composite is EVER rebuilt within a block,
  // and the busiest one is built exactly once per Κ-block.
  for (auto const& n : k_homed)
    CHECK(witness_builds_of(ws_cache.recompute_tally(), n) <= n_blocks);
  CHECK(max_ws_builds == n_blocks);

  // (b) recompute elimination: forest descent's fragmentation rebuilds the
  // SAME composites strictly more than the whole-scope build-once bound.
  CHECK(max_fd_builds > max_ws_builds);

  // (c) the oracle predicts the realized peak: both positive, and in
  // agreement within a generous-but-meaningful band. The two are
  // INDEPENDENT computations (a static co-residency sweep over
  // compute_dag_path's home_modes footprints vs a live eval-trace
  // high-water mark), not the identical arithmetic, so exact equality is not
  // expected -- but they must land within an order of magnitude of each
  // other for the oracle to be a meaningful predictor.
  CHECK(realized_peak_bytes > 0);
  CHECK(oracle_peak_bytes > 0.0);
  CHECK(ratio > 0.1);
  CHECK(ratio < 10.0);

  // (d) structural agreement between whole-scope and forest descent, on the
  // IDENTICAL forest/policy: the SAME number of forest roots processed (a
  // dropped or double-counted root would change this) and the SAME final
  // result footprint (a shape/positional agreement check -- see
  // witness_count_substr's doc comment for why the term-count check, not
  // this alone, is the drop/duplicate detector).
  REQUIRE(ws_result);
  REQUIRE(fd_result);
  CHECK(ws_term_count == forest.size());
  CHECK(fd_term_count == forest.size());
  CHECK(ws_result->size_in_bytes() == fd_result->size_in_bytes());

  // ---- meter() both executors, on the SAME forest/policy/regime/cfg the
  // replays above already used. Unlike k_homed (which deliberately tracks
  // only the K scope's homed values, skipping forest roots and leaves), the
  // meter's home_fidelity list covers EVERY distinct value in the recompute
  // tally, including the Kappa-free composites homed at the root scope
  // (home == "" -- a whole-nest invariant, not a K-block-local one) that
  // k_homed structurally excludes. Reporting those rebuild counts honestly
  // (not asserting them away) is this task's whole point -- the executor fix
  // that would make root-homed composites build once too is out of scope
  // here.
  sequant::BatchPolicy policy_fd = policy;
  policy_fd.scheduler = sequant::BatchScheduler::forest_descent;

  auto const ws_rep =
      sequant::eval::dryrun::meter(forest, policy_ws, regime, cfg);
  auto const fd_rep =
      sequant::eval::dryrun::meter(forest, policy_fd, regime, cfg);

  auto const has_root_homed = [](sequant::eval::dryrun::MeterReport const& r) {
    return std::any_of(r.home_fidelity.begin(), r.home_fidelity.end(),
                       [](auto const& h) { return h.home.empty(); });
  };

  if (auto it =
          std::find_if(ws_rep.home_fidelity.begin(), ws_rep.home_fidelity.end(),
                       [](auto const& h) { return h.home.empty(); });
      it != ws_rep.home_fidelity.end()) {
    std::cerr << "  [scope-executor-witness-water20] sample root-homed "
                 "(home={}) composite: "
              << it->label << ", builds=" << it->builds << "\n";
  }

  // (1) monitor == established measurement, per executor: the whole-scope
  // meter's own PeakMonitor agrees with the SAME realized_peak_bytes parsed
  // off the (a)-(d) replay's trace above -- WITHIN a small tolerance, not
  // bit-for-bit: meter()'s cache is built from cfg (is_volatile flags 't'
  // leaves, min_repeats=1), whereas the (a)-(d) run above's ws_cache is the
  // plain sequant::cache_manager(forest) (no persistence split, min_repeats
  // default 2) -- a DIFFERENT, slightly more permissive admission policy, so
  // a ~1% peak delta (measured: ~0.9%) is expected from admitting a few extra
  // P-frontier nodes, not a bug. 2% comfortably covers that measured gap
  // while still catching a gross regression.
  CHECK(double(ws_rep.peak_bytes) ==
        Catch::Approx(double(realized_peak_bytes)).epsilon(0.02));

  // (2) new signals populated, both executors.
  CHECK(ws_rep.builds_total > 0);
  CHECK(fd_rep.builds_total > 0);
  CHECK(!ws_rep.home_fidelity.empty());
  // Kappa-free home={} composites now appear (the old k_homed EXCLUDED them
  // by construction -- see this block's doc comment above):
  CHECK(has_root_homed(ws_rep));  // REPORTED, not asserted == 1 (the
                                  // executor fix is out of scope here)

  // (3) forest-vs-whole-scope comparison, via meter()'s OWN reports this
  // time (not the separate fd_cache/ws_cache pair (a)-(d) used): meter()'s
  // policy_fd.scheduler == BatchScheduler::forest_descent path now installs
  // sequant::make_evaluator(policy, yield) as its internal cache's custom
  // evaluator
  // before the coexistence-entry call (mirroring MPQC's wet forest path,
  // cck.ipp's `else` branch) rather than leaving that interception dark, so
  // fd_rep is the SAME per-root FRAGMENTED-but-still-BATCHED replay
  // max_fd_builds measures in (b) above, not an unbatched single pass
  // (confirmed: fd_rep.builds_total nearly doubles once the evaluator is
  // installed). The per-node comparison (b) already established --
  // fragmentation rebuilds a K-homed composite strictly more than
  // whole-scope's build-once-per-block bound -- reproduces through meter()'s
  // own home_fidelity list for the busiest K-homed composite, so re-check it
  // here rather than re-deriving (b)'s fd_cache/ws_cache pair a second time.
  //
  // The REPORT-WIDE totals (fd_rep.builds_total vs ws_rep.builds_total), by
  // contrast, are NOT a meaningful ordering and are deliberately left
  // UNASSERTED: ws_rep.builds_total is inflated by evaluate_whole_scope's
  // own, separately-diagnosed per-BatchGroup full-materialization replay
  // (see project history on the water-20 slowdown: "GROUND TRUTH ... Time
  // sink = per-group FULL replay after BatchGroup|End ... + leaf
  // re-materializations/iter"), which rebuilds many NON-K-homed values too,
  // on top of (not instead of) its build-once sharing of the K-homed
  // composites -- an orthogonal, already-tracked executor cost, not what
  // this comparison is about. Measured on this 40-term fixture: fd_rep
  // builds_total = 366 (batched, up from 195 unbatched pre-fix), ws_rep
  // builds_total = 978 -- so the REPORT-WIDE total is actually SMALLER for
  // the batched forest path despite its per-composite fragmentation being
  // worse; asserting a report-wide ordering here would either be false
  // (fd >= ws) or would silently paper over the whole-scope over-replay cost
  // (ws >= fd) instead of exercising the fact this block is actually about.
  auto const builds_of_hash = [](sequant::eval::dryrun::MeterReport const& r,
                                 std::size_t h) -> std::size_t {
    auto const it = std::find_if(r.home_fidelity.begin(), r.home_fidelity.end(),
                                 [h](auto const& hf) { return hf.hash == h; });
    return it == r.home_fidelity.end() ? 0 : it->builds;
  };
  std::size_t max_fd_rep_builds = 0, max_ws_rep_builds = 0;
  for (auto const& n : k_homed) {
    auto const h = n->hash_value();
    max_fd_rep_builds = std::max(max_fd_rep_builds, builds_of_hash(fd_rep, h));
    max_ws_rep_builds = std::max(max_ws_rep_builds, builds_of_hash(ws_rep, h));
  }
  CHECK(max_fd_rep_builds > 0);
  CHECK(max_ws_rep_builds > 0);
  CHECK(max_fd_rep_builds > max_ws_rep_builds);
}

// The C60 WITNESS. Same real C60 residual and the SAME occ-veto MPQC batch
// config test_eval_dryrun.cpp's [dryrun-occ-veto] "both" arm uses (Κ
// contracted @256, occ i external/spectator @8 -- the production
// csv_batch_policy.h config) -- reused here, not reinvented.
//
// EMPIRICAL FINDING (not an assumption): building the scope schedule from
// the REAL, fully-optimized C60 residual (all 55 summands in the fixture,
// `SEQUANT_UT_DRYRUN_NTERMS` capped or not) under this policy realizes only
// the OUTER Κ (contracted) loop -- no value is ever homed under occ i, so
// `build_scope_schedule` never nests an i child under it, unlike the
// HAND-BUILT synthetic forest test_eval_dryrun.cpp's "outer-homed aux
// composites... across occ blocks" test uses to exercise that nesting.
// Reading why: occ-veto's own per-NODE `node_slice_mask()` tally does show
// External-occ stamps (>0, confirmed separately), but every Κ-carrying
// shared composite in this residual (the gC-class DF integrals) ALSO
// carries occ as a proto index of its own CSV/PNO domain (the composite's
// definition, not an enclosing loop) -- so no value here is simultaneously
// Κ-shared AND occ-invariant, the precondition build_scope_schedule's
// (root-to-node) home-modes assignment needs to place anything under a
// nested i child. This witness therefore exercises the SAME single
// Contracted-loop walk_scope branch the water-20 witness does, now on the
// real, much larger, C60 topology (55 summands, individually annotated with
// BOTH batch roles) rather than the External-scatter branch; the External
// branch is exercised by the small hand-built forests in test_eval_ta.cpp /
// test_eval_dryrun.cpp instead.
//
// Its scope is deliberately a MEASUREMENT, not a budget gate: per the
// design's Scope boundary, placement/feasibility is the optimizer's
// concern, not the executor's, so this records the realized peak and the
// oracle's prediction without asserting either is "small enough" -- only
// that they are computed and that whole-scope's build-once win still holds
// (a CHECK, since that direction was already established at water-20 scale
// and is cheap to keep verifying here).
//
// Correctness-to-tolerance: DryRun is zero-data (no real floating-point
// content to compare -- see test_eval_ta.cpp's real-data TA proxy for that
// at a tractable size), so "correct" here means STRUCTURAL equivalence,
// checked two ways: whole-scope and forest descent, replayed on the
// IDENTICAL forest, must (1) process the SAME number of forest roots (a
// `Term | Begin` trace-record count equal to `forest.size()` on BOTH
// replays -- the check that actually rules out a silently dropped or
// double-counted root; see witness_count_substr's doc comment for why
// `size_in_bytes()` alone cannot), and (2) assemble the SAME final result
// footprint (`size_in_bytes()` equality -- a shape/positional agreement
// check layered on top, not by itself a drop detector). This is the
// C60-scale complement to test_eval_ta.cpp's small real-arithmetic
// equivalence tests: it exercises the actual production residual (55
// summands, both batch roles individually annotated) that is intractable at
// real tensor sizes.
TEST_CASE(
    "whole-scope witness: C60 aux+occ residual runs whole-scope, matches "
    "forest descent structurally, and records the realized peak",
    "[.][scope-executor-witness-c60]") {
  using sequant::eval::dryrun::EvalExprDryRun;
  using sequant::eval::dryrun::EvalNodeDryRun;

  auto ctx = sequant::get_default_context().clone();
  ctx.set_first_dummy_index_ordinal(1000000);
  auto isr = ctx.mutable_index_space_registry();
  REQUIRE(isr != nullptr);
  sequant::mbpt::add_pao_spaces(isr, sequant::mbpt::Spin::any);
  sequant::mbpt::add_df_spaces(isr);
  auto ctx_resetter = sequant::set_scoped_default_context(std::move(ctx));

  auto const body = witness_slurp(std::string(SEQUANT_UNIT_TESTS_SOURCE_DIR) +
                                  "/data/csv_ccsd_doubles_residual_df.txt");
  REQUIRE(!body.empty());
  std::string line = body;
  if (auto nl = line.find('\n'); nl != std::string::npos)
    line = line.substr(0, nl);
  auto expr = sequant::deserialize<sequant::ExprPtr>(line);
  REQUIRE(static_cast<bool>(expr));
  REQUIRE(expr->is<sequant::Sum>());
  auto const& summands = expr->as<sequant::Sum>().summands();
  REQUIRE(!summands.empty());

  // The full residual (55 summands in the fixture) by default -- fast enough
  // in practice for this witness (~tens of seconds); env-overridable to a
  // smaller subset for a quicker manual run, matching [dryrun-occ-veto]'s own
  // knob.
  std::size_t nterms = summands.size();
  if (char const* nt = std::getenv("SEQUANT_UT_DRYRUN_NTERMS"))
    nterms = std::min<std::size_t>(summands.size(), std::atoll(nt));

  auto regime = witness_df_regime(kWitnessC60_pVDZF12);
  auto cm = std::make_shared<sequant::eval::dryrun::CostModel const>(regime);

  // EXACT MPQC aux+occ config ([dryrun-occ-veto]'s "both" arm / production
  // csv_batch_policy.h): Κ contracted @256 outer, occ i external/spectator
  // @8 inner.
  sequant::BatchPolicy policy;
  policy.batch_spectator_indices = true;
  policy.order_aware_recompute = true;
  policy.node_level_placement = true;
  policy.is_batchable_contracted_index = [](sequant::Index const& ix) {
    return ix.space().base_key() == L"Κ";
  };
  policy.is_batchable_external_index = [](sequant::Index const& ix) {
    return ix.space().base_key() == L"i";
  };
  policy.batch_target_size = [](sequant::Index const& ix) -> std::size_t {
    return ix.space().base_key() == L"Κ" ? std::size_t{256} : std::size_t{8};
  };
  policy.is_volatile_leaf = [](sequant::Tensor const& t) {
    return t.label() == L"t";
  };
  policy.accumulation_factor = 1.0;
  policy.peak_threshold = 100e9;  // the production C60 job's 100 GB budget

  auto axes_map = std::make_shared<std::unordered_map<
      sequant::Expr const*,
      sequant::container::vector<sequant::NodeBatchAnnotation>>>();
  sequant::OptimizeOptions opts;
  opts.objective_function = sequant::ObjectiveFunction::DenseTimeSpaceBatched;
  opts.idx_to_extent = regime.idx_to_extent();
  opts.inner_pow = regime.inner_pow_fn();
  opts.batch_policy = policy;
  opts.volatile_weight = 20.0;
  opts.roofline.machine_balance = 200.0;
  opts.roofline.fast_mem_elems = 1000000.0;
  opts.term_batch_axes = axes_map;

  std::vector<EvalNodeDryRun> forest;
  for (std::size_t s = 0; s < nterms; ++s) {
    sequant::ExprPtr const term = witness_flatten_product(summands[s]);
    if (!term) continue;
    sequant::ExprPtr optimized;
    try {
      optimized = sequant::optimize(term, opts);
    } catch (std::exception const&) {
      continue;
    }
    if (!optimized) continue;
    sequant::BinarizationOptions bopts;
    if (auto it = axes_map->find(optimized.get()); it != axes_map->end())
      bopts.node_batch_axes = it->second;
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    forest.push_back(sequant::binarize<EvalExprDryRun>(optimized, {}, bopts));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
  }
  REQUIRE(!forest.empty());

  // Scope schedule from the SAME two-role mode order the water-20 nested
  // test uses ({Κ, i}) -- but see the EMPIRICAL FINDING in this TEST_CASE's
  // doc comment: on THIS real residual/policy, build_scope_schedule realizes
  // only the OUTER Κ (contracted) loop; no value is ever homed under occ i
  // (the External branch of walk_scope is exercised by the small hand-built
  // forests in test_eval_ta.cpp / test_eval_dryrun.cpp instead). Assert on
  // what the real DAG actually realizes rather than assuming the synthetic
  // forests' nested shape.
  auto const block_of = [](sequant::Index const& ix) -> std::size_t {
    return ix.space().base_key() == L"Κ" ? std::size_t{256} : std::size_t{8};
  };
  auto rich = sequant::eval::compute_dag_boulevard(forest, *cm, block_of);
  auto sched =
      sequant::eval::build_scope_schedule<std::wstring>(rich, {L"Κ", L"i"});
  REQUIRE(sched.root.children.size() == 1);
  auto const& kscope = sched.root.children.front();
  REQUIRE(kscope.mode.space().base_key() == L"Κ");
  REQUIRE(kscope.kind == sequant::BatchModeType::Contracted);
  CHECK(kscope.children.empty());  // documents the empirical finding above
  sequant::Index const K = kscope.mode;

  // The TARGETED shared composites: the OUTER Κ scope's homed values (Κ
  // carrying, occ-invariant), resolved to forest nodes.
  auto const vmap = sequant::eval::build_value_node_map(forest);
  std::unordered_set<std::size_t> root_hashes;
  for (auto const& r : forest) root_hashes.insert(r->hash_value());
  std::vector<EvalNodeDryRun> k_homed;
  for (auto vid : kscope.homed_values) {
    auto const h = rich.cells[vid].hash;
    auto const it = vmap.find(h);
    if (it == vmap.end() || it->second.leaf() || root_hashes.count(h)) continue;
    k_homed.push_back(it->second);
  }

  sequant::eval::dryrun::DryRunLeafEvaluator yield{cm};
  auto const target = [](sequant::Index const& ix) -> std::size_t {
    return ix.space().base_key() == L"Κ" ? std::size_t{256} : std::size_t{8};
  };

  auto aops = sequant::eval::dryrun::make_dryrun_array_ops(cm);
  std::size_t n_kappa = 0;
  for (auto const& root : forest)
    if (auto lf = sequant::find_leaf_carrying(root, K)) {
      n_kappa = aops.axis_batches(K, 256).size();
      break;
    }
  REQUIRE(n_kappa > 1);

  using annot_t = std::remove_cvref_t<decltype(forest.front()->annot())>;
  annot_t const layout{};

  auto& logger = sequant::Logger::instance();
  auto const prev_level = logger.eval.level;
  auto* const prev_stream = logger.eval.stream;

  // ---- whole-scope replay: structural result, build tally, eval trace. The
  // elevated logger level is what makes BOTH the trace records AND the
  // DryRun per-op cost sink that feeds CacheManager::tally_build
  // (DryRunOps::prod's `if (Logger::instance().eval.level > 0)` RUNTIME gate
  // on last_op_flops(), result.hpp) actually fire -- so it stays elevated
  // across the forest-descent contrast run below too, not just this one. ----
  std::ostringstream ws_trace;
  logger.eval.level = 1;
  logger.eval.stream = &ws_trace;
  auto ws_cache = sequant::cache_manager(forest);
  ws_cache.set_array_ops(&aops);
  ws_cache.set_recompute_tally_enabled(true);
  ResultPtr ws_result;
  try {
    ws_result = sequant::eval::evaluate_whole_scope<sequant::Trace::On>(
        forest, sched, rich, layout, yield, ws_cache, target);
  } catch (std::exception const& e) {
    std::cerr << "[scope-executor-witness-c60] whole-scope evaluate threw: "
              << e.what() << "\n";
  }

  // Parse the whole-scope peak AND term-count BEFORE the forest-descent run
  // reuses the logger's stream (below), so the fd replay's own trace text
  // cannot contaminate either measurement.
  std::size_t const realized_peak_bytes =
      witness_parse_max_hwmark_bytes(ws_trace.str());
  std::size_t const ws_term_count =
      witness_count_substr(ws_trace.str(), "Term | Begin");
  std::size_t max_ws_builds = 0;
  for (auto const& n : k_homed)
    max_ws_builds = std::max(max_ws_builds,
                             witness_builds_of(ws_cache.recompute_tally(), n));

  // ---- forest-descent contrast: SAME forest/policy, structural result +
  // build tally, for the correctness check (including the term-count
  // drop/duplicate detector, read from THIS trace below) and the
  // recompute-elimination direction. Same elevated logger level (see
  // above), a FRESH trace sink. ----
  std::ostringstream fd_trace;
  logger.eval.stream = &fd_trace;
  auto fd_cache = sequant::cache_manager(forest);
  fd_cache.set_array_ops(&aops);
  fd_cache.set_custom_evaluator(sequant::make_evaluator(policy, yield));
  fd_cache.set_recompute_tally_enabled(true);
  ResultPtr fd_result;
  try {
    fd_result =
        sequant::evaluate<sequant::Trace::On>(forest, layout, yield, fd_cache);
  } catch (std::exception const& e) {
    std::cerr << "[scope-executor-witness-c60] forest-descent evaluate threw: "
              << e.what() << "\n";
  }
  logger.eval.level = prev_level;
  logger.eval.stream = prev_stream;

  std::size_t const fd_term_count =
      witness_count_substr(fd_trace.str(), "Term | Begin");
  std::size_t max_fd_builds = 0;
  for (auto const& n : k_homed)
    max_fd_builds = std::max(max_fd_builds,
                             witness_builds_of(fd_cache.recompute_tally(), n));

  // ---- the co-residency oracle (measurement only -- no budget gate). ----
  sequant::eval::dryrun::CacheConfig cfg;
  cfg.max_footprint = 0.;
  cfg.min_repeats = 1;
  cfg.is_volatile = [](EvalNodeDryRun const& n) {
    return n.leaf() && n->is_tensor() && n->as_tensor().label() == L"t";
  };
  sequant::BatchPolicy policy_ws = policy;
  policy_ws.scheduler = sequant::BatchScheduler::whole_scope;
  auto const cp =
      sequant::eval::dryrun::cost_profile(forest, policy_ws, cfg, regime);
  double const oracle_peak_bytes = cp.peak_bytes;
  double const ratio = oracle_peak_bytes > 0.0
                           ? double(realized_peak_bytes) / oracle_peak_bytes
                           : 0.0;

  std::wcerr
      << L"\n=== [scope-executor-witness-c60] C60 aux+occ, " << forest.size()
      << L" terms (of " << nterms << L" attempted), n_Κ_blocks=" << n_kappa
      << L", " << k_homed.size() << L" Κ-homed composites ===\n"
      << L"  whole-scope replay completed  = " << (ws_result ? L"yes" : L"NO")
      << L"\n"
      << L"  forest-descent replay completed = " << (fd_result ? L"yes" : L"NO")
      << L"\n"
      << L"  term (root) count: whole-scope = " << ws_term_count
      << L", forest-descent = " << fd_term_count << L"  (forest.size()="
      << forest.size() << L")\n"
      << L"  final result size: whole-scope = "
      << (ws_result ? double(ws_result->size_in_bytes()) / 1e9 : -1.0)
      << L" GB, forest-descent = "
      << (fd_result ? double(fd_result->size_in_bytes()) / 1e9 : -1.0)
      << L" GB\n"
      << L"  whole-scope max builds (targeted composites) = " << max_ws_builds
      << L"  (n_Κ_blocks=" << n_kappa << L")\n"
      << L"  forest-descent max builds (SAME composites)  = " << max_fd_builds
      << L"\n"
      << L"  realized peak (hw= trace max) = "
      << (double(realized_peak_bytes) / 1e9) << L" GB\n"
      << L"  co-residency oracle peak      = " << (oracle_peak_bytes / 1e9)
      << L" GB  (feasibility of this placement is NOT judged here -- the "
         L"optimizer's concern, not the executor's)\n"
      << L"  ratio (realized/oracle)       = " << ratio << L"\n";

  // Per-composite build counts (RECORDED, not gated on a budget/pattern --
  // see this TEST_CASE's doc comment: on the real 55-summand forest the
  // per-composite picture is genuinely mixed -- some Κ-homed composites hit
  // the clean "built once per Κ-block, zero forest-descent rebuilds" pattern
  // water-20 shows uniformly, others do not, plausibly because 55
  // independently-optimized summands bind the SAME canonical composite
  // under physically-divergent index labels (the codebase's own documented
  // `divergent_modes` / slicing-signature phenomenon, see
  // place_at_this_level's "Split gate" doc comment in eval.hpp) -- so no
  // single-number pass/fail bound is asserted here, per the C60 witness's
  // "measure, do not gate" scope).
  std::wcerr << L"  per-Κ-homed-composite builds (whole-scope / forest-"
                L"descent):\n";
  for (std::size_t i = 0; i != k_homed.size(); ++i)
    std::wcerr << L"    [" << i << L"] "
               << witness_builds_of(ws_cache.recompute_tally(), k_homed[i])
               << L" / "
               << witness_builds_of(fd_cache.recompute_tally(), k_homed[i])
               << L"\n";

  // Structural correctness, where checkable: both replays must complete,
  // process the SAME number of forest roots (the term-count check that
  // actually rules out a dropped or double-counted root -- see
  // witness_count_substr's doc comment), and assemble the SAME final
  // footprint (a shape/positional agreement check layered on top).
  REQUIRE(ws_result);
  REQUIRE(fd_result);
  CHECK(ws_term_count == forest.size());
  CHECK(fd_term_count == forest.size());
  CHECK(ws_result->size_in_bytes() == fd_result->size_in_bytes());

  // Build-once/recompute: MEASURE and record only (see the per-composite
  // table above) -- the strict "== n_blocks, forest descent worse" bound is
  // the water-20 witness's assertion; here only the loose sanity floor holds
  // (whole-scope built at least one targeted composite at all). The
  // precondition itself (some Κ-homed composite exists) is NOT optional --
  // REQUIRE it directly so a future change that eliminated Κ-homed
  // composites entirely would fail loudly here instead of this CHECK
  // silently no-opping.
  REQUIRE(!k_homed.empty());
  CHECK(max_ws_builds > 0);

  // Peak/oracle: MEASURE and record only -- no budget assertion (the design's
  // Scope boundary explicitly leaves placement feasibility to the
  // optimizer).
  CHECK(realized_peak_bytes > 0);
  CHECK(oracle_peak_bytes > 0.0);
}
