#include <SeQuant/core/eval/backends/dryrun/cost_model_object.hpp>
#include <SeQuant/core/eval/backends/dryrun/eval_expr.hpp>
#include <SeQuant/core/eval/backends/dryrun/meter.hpp>
#include <SeQuant/core/eval/backends/dryrun/size_regime.hpp>
#include <SeQuant/core/eval/cache_manager.hpp>
#include <SeQuant/core/eval/eval.hpp>
#include <SeQuant/core/eval/peak_monitor.hpp>
#include <SeQuant/core/eval/peak_profile.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/logger.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cstddef>
#include <memory>
#include <sstream>
#include <vector>

TEST_CASE("PeakMonitor tracks hierarchy-wide co-resident high-water",
          "[meter]") {
  sequant::eval::PeakMonitor mon;
  std::vector<std::size_t> peaks;
  mon.on_peak = [&](sequant::eval::PeakEvent const& e) {
    peaks.push_back(e.bytes);
  };
  mon.observe(100, 0xA);
  mon.observe(50, 0xB);   // below hwmark: no advance, no fire
  mon.observe(300, 0xC);  // new peak
  CHECK(mon.hwmark_bytes == 300);
  CHECK(mon.peak.op_hash == 0xC);
  CHECK(peaks == std::vector<std::size_t>{100, 300});
}

TEST_CASE(
    "PeakMonitor wired onto a CacheManager observes the same high-water as "
    "the cache's own working_set_hwmark()",
    "[meter]") {
  using sequant::eval::dryrun::CostModel;
  using sequant::eval::dryrun::DryRunLeafEvaluator;
  using sequant::eval::dryrun::EvalExprDryRun;
  using sequant::eval::dryrun::EvalNodeDryRun;
  using sequant::eval::dryrun::SizeRegime;

  // A tiny, self-consistent regime: two named spaces, no composite (proto-
  // indexed) legs involved, so no CSV/PNO moments are needed.
  SizeRegime regime;
  regime.space_extent = {
      {L"i", 10},
      {L"a", 20},
  };
  auto cm = std::make_shared<CostModel const>(regime);

  // A minimal 3-node forest: two leaves (g, t) contracted into one product,
  // fully contracted (no external indices) so the whole tree is a scalar.
  // Small enough to hand-build directly rather than routing through
  // sequant::optimize (see test-1 brief: a 2-3 node hand-built forest is
  // acceptable for this cache-integration assertion).
  auto expr = sequant::deserialize<sequant::ExprPtr>("g{i_1;a_3} * t{a_3;i_1}");
  REQUIRE(static_cast<bool>(expr));

  SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
  auto node = sequant::binarize<EvalExprDryRun>(expr);
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
  REQUIRE_FALSE(node.leaf());

  auto cache = sequant::cache_manager(std::vector<EvalNodeDryRun>{node});

  sequant::eval::PeakMonitor mon;
  cache.set_peak_monitor(&mon);

  // Redirect the eval trace to a private buffer (rather than stdout) and
  // restore the logger's prior state afterward -- note_working_set()'s
  // per-op hwmark input is only folded when Logger::instance().eval.level >
  // 0 (see cache_manager.hpp / eval.hpp).
  auto& logger = sequant::Logger::instance();
  auto const prev_level = logger.eval.level;
  auto* const prev_stream = logger.eval.stream;
  std::ostringstream trace_os;
  logger.eval.level = 1;
  logger.eval.stream = &trace_os;

  sequant::eval::dryrun::DryRunLeafEvaluator yield{cm};
  sequant::ResultPtr result;
  try {
    result = sequant::evaluate<sequant::Trace::On>(node, yield, cache);
  } catch (...) {
    logger.eval.level = prev_level;
    logger.eval.stream = prev_stream;
    throw;
  }
  logger.eval.level = prev_level;
  logger.eval.stream = prev_stream;

  REQUIRE(result);
  CHECK(mon.hwmark_bytes > 0);
  CHECK(mon.hwmark_bytes == cache.working_set_hwmark());
}

TEST_CASE(
    "assemble_report rolls a metered dry-run replay's per-node build tally "
    "into a MeterReport (peak, persistent/volatile FLOPs+time, build-vs-home)",
    "[meter]") {
  using sequant::eval::compute_dag_boulevard;
  using sequant::eval::PeakMonitor;
  using sequant::eval::dryrun::assemble_report;
  using sequant::eval::dryrun::compute_volatility;
  using sequant::eval::dryrun::CostModel;
  using sequant::eval::dryrun::DryRunLeafEvaluator;
  using sequant::eval::dryrun::EvalExprDryRun;
  using sequant::eval::dryrun::EvalNodeDryRun;
  using sequant::eval::dryrun::SizeRegime;

  SizeRegime regime;
  regime.space_extent = {
      {L"i", 10},
      {L"a", 20},
  };
  auto cm = std::make_shared<CostModel const>(regime);

  // Three-factor product, folded LEFT-TO-RIGHT by binarize (see
  // fold_left_to_node, binary_node.hpp): root == (X * t) with X == (g * h).
  // X is a PERSISTENT internal node (neither g nor h is volatile); the root
  // additionally consumes the volatile leaf t, so the root is volatile --
  // exercising BOTH branches of assemble_report's persistent/volatile split.
  // X keeps i_1 (from g) and i_2 (from h) external after contracting a_3; t
  // contracts i_2 and keeps a_5 external, so the ROOT keeps i_1 and a_5
  // external too (not a bare scalar). That is deliberate: i_1 being a
  // genuine external slot of BOTH X and the root lets the ROOT's own
  // External loop over i_1 (stamped below) give X a NONEMPTY home -- X sits
  // inside that loop but does not own it -- while the root's OWN home stays
  // empty (a node's own realized loop is excluded from its own home; see
  // compute_dag_boulevard's own_modes_union subtraction). That asymmetry is
  // the content check below: a broken cell_by_hash lookup in assemble_report
  // would silently default BOTH to empty, indistinguishable from the root's
  // genuinely-empty case, so only X's nonempty "i" catches a broken lookup.
  auto expr = sequant::deserialize<sequant::ExprPtr>(
      "g{i_1;a_3} * h{a_3;i_2} * t{i_2;a_5}");
  REQUIRE(static_cast<bool>(expr));

  SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
  auto node = sequant::binarize<EvalExprDryRun>(expr);
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
  REQUIRE_FALSE(node.leaf());
  REQUIRE_FALSE(node.left().leaf());  // X == g * h (persistent)
  REQUIRE(node.right().leaf());       // t (volatile)

  sequant::Index const mode{L"i_1"};
  REQUIRE(sequant::index_position(node, mode).has_value());
  REQUIRE(sequant::index_position(node.left(), mode).has_value());
  // Stamp the ROOT's own External loop over i_1. Plain (non-scope)
  // evaluate() ignores batched_here without a custom evaluator (see
  // test_eval_dryrun.cpp's equivalent stamping), so this only feeds
  // compute_dag_boulevard's home/ectx bookkeeping below, not the replay.
  node->set_batched_here({{mode, sequant::BatchModeType::External}});

  std::vector<EvalNodeDryRun> const forest{node};

  auto cache = sequant::cache_manager(forest);
  cache.set_recompute_tally_enabled(true);

  sequant::eval::PeakMonitor mon;
  cache.set_peak_monitor(&mon);

  auto const is_volatile = [](EvalNodeDryRun const& n) {
    return n.leaf() && n->is_tensor() && n->as_tensor().label() == L"t";
  };

  auto const block_of = [](sequant::Index const&) -> std::size_t {
    return 256;
  };
  auto const rich = compute_dag_boulevard(forest, *cm, block_of);

  auto& logger = sequant::Logger::instance();
  auto const prev_level = logger.eval.level;
  auto* const prev_stream = logger.eval.stream;
  std::ostringstream trace_os;
  logger.eval.level = 1;
  logger.eval.stream = &trace_os;

  sequant::eval::dryrun::DryRunLeafEvaluator yield{cm};
  sequant::ResultPtr result;
  try {
    result = sequant::evaluate<sequant::Trace::On>(node, yield, cache);
  } catch (...) {
    logger.eval.level = prev_level;
    logger.eval.stream = prev_stream;
    throw;
  }
  logger.eval.level = prev_level;
  logger.eval.stream = prev_stream;

  REQUIRE(result);

  auto const report = assemble_report(cache, mon, rich, forest, is_volatile,
                                      /*whole_scope=*/false);

  CHECK(report.builds_total > 0);
  CHECK_FALSE(report.home_fidelity.empty());
  CHECK(report.whole_scope == false);

  // EXEC/COST SPLIT: exec is threaded BuildRecord -> tally_build ->
  // assemble_report independently of (and not swapped with) flops -- both
  // the volatile (root) and persistent (X) buckets get a positive exec
  // estimate alongside their flops.
  CHECK(report.flops_volatile > 0.0);
  CHECK(report.cost_volatile > 0.0);
  CHECK(report.flops_persistent > 0.0);
  CHECK(report.cost_persistent > 0.0);

  // HOME-FIDELITY CONTENT: look up X's (persistent, g*h) and the root's
  // entries by hash and check their home/uses STRINGS, not just
  // non-emptiness of the whole list -- see the forest-design comment above
  // for why X's nonempty "i" is the discriminating case.
  auto const find_by_hash = [&](std::size_t h) {
    return std::find_if(report.home_fidelity.begin(),
                        report.home_fidelity.end(),
                        [h](auto const& hf) { return hf.hash == h; });
  };

  auto const x_it = find_by_hash(node.left()->hash_value());
  REQUIRE(x_it != report.home_fidelity.end());
  CHECK(x_it->home == "i");
  CHECK(x_it->uses == "i");

  auto const root_it = find_by_hash(node->hash_value());
  REQUIRE(root_it != report.home_fidelity.end());
  CHECK(root_it->home.empty());
  CHECK(root_it->uses.empty());
}

// Task 3: meter() drives the REAL policy-selected executor (whole-scope or
// forest descent, chosen by BatchPolicy::whole_scope_execution) through the
// sizing backend with its own metered cache, rather than a hand-rolled proxy
// of one. This forest carries a genuine Contracted batch axis (a_3, the
// "aux"-analog): with whole_scope_execution=true, the coexistence entry
// (scope_executor.hpp) rebuilds a scope tree with ONE realized batch loop
// from that stamp and drives the executor's nested batch-scratch caches
// (a Task-1 coverage gap -- no earlier [meter] test exercised a multi-level
// PeakMonitor parent-chain under a real batched walk); with it false, the
// SAME forest runs through today's unbatched per-tree descent. Both modes
// must report a positive peak/build count and the matching `whole_scope`
// flag.
TEST_CASE(
    "meter runs the policy-selected executor (whole-scope and forest) with "
    "the sizing backend and a metered cache",
    "[meter]") {
  using sequant::BatchModeType;
  using sequant::BatchPolicy;
  using sequant::Index;
  using sequant::eval::dryrun::CacheConfig;
  using sequant::eval::dryrun::EvalExprDryRun;
  using sequant::eval::dryrun::EvalNodeDryRun;
  using sequant::eval::dryrun::meter;
  using sequant::eval::dryrun::SizeRegime;

  // Same small, self-consistent regime as the tests above; "a" (extent 20)
  // doubles as the batch-axis space here (target_size 5 => 4 blocks).
  SizeRegime regime;
  regime.space_extent = {
      {L"i", 10},
      {L"a", 20},
  };

  // A single-root, fully-contracted product -- a_3 is a genuine operand
  // index of BOTH factors (contracted away at the root), hand-stamped as the
  // root's one Contracted batch axis (no optimize() involved, matching
  // test_eval_dryrun.cpp's hand-built batch-annotation recipe): the
  // coexistence entry rebuilds its schedule from the forest's OWN
  // batched_here() stamps, not from BatchPolicy predicates, so this alone is
  // enough for it to realize a non-root-only scope tree.
  auto expr =
      sequant::deserialize<sequant::ExprPtr>(L"g{i_1;a_3} * h{a_3;i_1}");
  REQUIRE(static_cast<bool>(expr));

  SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
  auto node = sequant::binarize<EvalExprDryRun>(expr);
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
  REQUIRE_FALSE(node.leaf());

  Index const a3{L"a_3"};
  node->set_batched_here({{a3, BatchModeType::Contracted}});

  std::vector<EvalNodeDryRun> const forest{node};

  CacheConfig const cfg;  // default: no footprint gate, min_repeats=2

  BatchPolicy policy;
  policy.batch_target_size = [](Index const& ix) -> std::size_t {
    return ix.space().base_key() == L"a" ? std::size_t{5} : std::size_t{1};
  };

  // ---- whole-scope: exercises the nested batch-scratch walk. ----
  policy.whole_scope_execution = true;
  auto const ws_report = meter(forest, policy, regime, cfg);
  CHECK(ws_report.peak_bytes > 0.0);
  CHECK(ws_report.builds_total > 0);
  CHECK(ws_report.whole_scope == true);

  // ---- forest descent: the SAME forest/policy, flag off. ----
  policy.whole_scope_execution = false;
  auto const fd_report = meter(forest, policy, regime, cfg);
  CHECK(fd_report.peak_bytes > 0.0);
  CHECK(fd_report.builds_total > 0);
  CHECK(fd_report.whole_scope == false);
}
