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

  // One volatile leaf (t) contracted against one persistent leaf (g); the
  // resulting product is fully contracted to a scalar (no external indices).
  auto expr = sequant::deserialize<sequant::ExprPtr>("g{i_1;a_3} * t{a_3;i_1}");
  REQUIRE(static_cast<bool>(expr));

  SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
  auto node = sequant::binarize<EvalExprDryRun>(expr);
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
  REQUIRE_FALSE(node.leaf());

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
  CHECK(report.flops_volatile > 0.0);
  CHECK_FALSE(report.home_fidelity.empty());
  CHECK(report.whole_scope == false);
}
