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

#include <SeQuant/core/eval/backends/dryrun/cost_model_object.hpp>
#include <SeQuant/core/eval/backends/dryrun/size_regime.hpp>
#include <SeQuant/core/eval/cache_manager.hpp>
#include <SeQuant/core/eval/eval.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/eval/peak_profile.hpp>
#include <SeQuant/core/eval/result.hpp>
#include <SeQuant/core/eval/scope_executor.hpp>
#include <SeQuant/core/eval/scope_schedule.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/logger.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {

using sequant::Constant;
using sequant::EvalExpr;
using sequant::EvalNode;
using sequant::ExprPtr;
using sequant::ResultPtr;
using sequant::ResultScalar;
using sequant::Variable;
using sequant::eval::build_scope_schedule;
using sequant::eval::compute_dag_boulevard;
using sequant::eval::evaluate_whole_scope;
using sequant::eval::RichSchedule;
using sequant::eval::ScopeSchedule;
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
/// set_batched_here(), so the resulting scope tree must be root-only. Reused
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
