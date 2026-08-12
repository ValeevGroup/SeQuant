// SP3 Task 1 of the ordered-scope batched-eval design (the sequel to the
// whole-scope batched DAG execution design, see
// doc/dev/specs/2026-08-10-whole-scope-batched-dag-execution-design.md, and
// ordered_schedule.hpp's own doc comments for the SP2 OrderedSchedule IR):
// the ORDERED executor SKELETON. For a forest with NO batchable index (so
// build_ordered_schedule realizes ordered.root as a flat, topologically
// sorted sequence of BuildSteps -- no nested child ScopeBlock, since no axis
// type is ever realized), evaluate_ordered_schedule must reproduce the SAME
// numeric result the existing per-tree forest descent (sequant::evaluate(
// Nodes const&, ...)) produces.
//
// Mirrors test_scope_executor.cpp's harness exactly: the equivalence check
// needs REAL numeric arithmetic (not just shape/size modeling under the
// zero-data DryRun backend), so this reuses the identical minimal
// ScalarEvalExpr subclass (Constant/Variable leaves, Sum/Product internal
// nodes, no tensor backend) that file introduced -- duplicated here (no
// shared test header exists between the two .cpp files) rather than
// factored out, per that file's own precedent for the DryRun witness
// fixtures.

#include <SeQuant/core/batch_policy.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/eval/backends/dryrun/cost_model_object.hpp>
#include <SeQuant/core/eval/backends/dryrun/size_regime.hpp>
#include <SeQuant/core/eval/cache_manager.hpp>
#include <SeQuant/core/eval/eval.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/eval/legality.hpp>
#include <SeQuant/core/eval/ordered_executor.hpp>
#include <SeQuant/core/eval/ordered_schedule.hpp>
#include <SeQuant/core/eval/peak_profile.hpp>
#include <SeQuant/core/eval/result.hpp>
#include <SeQuant/core/eval/scope_executor.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/expressions/tensor.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/runtime.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <map>
#include <string>
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
using sequant::eval::analyze_legality;
using sequant::eval::build_ordered_schedule;
using sequant::eval::compute_dag_boulevard;
using sequant::eval::evaluate_ordered_schedule;
using sequant::eval::OrderedSchedule;
using sequant::eval::RichSchedule;
using sequant::eval::dryrun::CostModel;
using sequant::eval::dryrun::SizeRegime;

///
/// \brief A minimal EvalExpr subclass whose only job is to satisfy \c
/// meta::can_evaluate (i.e. carry an annot() method) so a plain scalar
/// arithmetic forest (Constant/Variable leaves, Sum/Product internal nodes,
/// no tensor indices at all) can be run through \c evaluate_impl /
/// \c evaluate_ordered_schedule without pulling in a tensor backend.
/// Identical to test_scope_executor.cpp's ScalarEvalExpr.
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
/// Identical to test_scope_executor.cpp's ScalarLeafEvaluator.
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

}  // namespace

TEST_CASE(
    "evaluate_ordered_schedule matches forest descent for an unbatched "
    "forest",
    "[ordered-executor]") {
  // Same two-root forest as test_scope_executor.cpp's TestForest, a
  // nontrivial Sum-of-Products over named scalar Variables sharing a common
  // subexpression pattern (2*a*b appears in both roots' construction path
  // via 'a' and 'b'), exercising both EvalOp::Sum and EvalOp::Product inside
  // evaluate_impl plus the cross-root accumulation the combine step
  // reproduces.
  std::vector<ScalarNode> forest{scalar_tree(L"2 * a * b - c"),
                                 scalar_tree(L"a * a + 3 * b - 2 * c")};

  ScalarLeafEvaluator const yield{{{L"a", 2.0}, {L"b", -3.5}, {L"c", 7.25}}};

  // NO batchable index: build_ordered_schedule must realize a single root
  // block of plain BuildSteps (no child ScopeBlock) -- the shape SP3 Task 1
  // handles. Default-constructed BatchPolicy declines every index in both
  // roles (is_batchable_contracted_index / is_batchable_external_index both
  // default to "false"), so this is the natural no-batching policy, not a
  // hand-suppressed one.
  sequant::BatchPolicy const policy;

  SizeRegime const regime;
  CostModel const cm{regime};
  auto const block_of = [](Index const&) -> std::size_t { return 1; };
  RichSchedule const rich = compute_dag_boulevard(forest, cm, block_of);
  auto const legality = analyze_legality(rich, forest, policy);
  OrderedSchedule const ordered =
      build_ordered_schedule(rich, legality, policy, {});

  // The precondition this test exists to exercise: a flat BuildStep sequence,
  // no nested loop block.
  for (auto const& step : ordered.root.steps)
    REQUIRE(std::holds_alternative<sequant::eval::BuildStep>(step.value));
  REQUIRE(!ordered.root.steps.empty());
  REQUIRE(ordered.num_values == rich.cells.size());

  // Reference: existing per-tree forest descent, own (fresh) cache.
  auto ref_cache = sequant::CacheManager<ScalarNode>::empty();
  ResultPtr const reference = sequant::evaluate(forest, yield, ref_cache);
  double const expected = reference->as<ResultScalar<double>>().value();

  // evaluate_ordered_schedule, its own (fresh) cache -- must match exactly.
  auto ordered_cache = sequant::CacheManager<ScalarNode>::empty();
  std::function<std::size_t(Index const&)> const target =
      [](Index const&) -> std::size_t { return 1; };
  ResultPtr const got = evaluate_ordered_schedule(forest, ordered, rich,
                                                  ScalarEvalExpr::annot_t{},
                                                  yield, ordered_cache, target);
  double const got_val = got->as<ResultScalar<double>>().value();

  // Hand-computed cross-check that the reference itself is right:
  //   p1 = 2*a*b - c = 2*2*(-3.5) - 7.25 = -14 - 7.25 = -21.25
  //   p2 = a*a + 3*b - 2*c = 4 + 3*(-3.5) - 2*7.25 = 4 - 10.5 - 14.5 = -21
  //   total = p1 + p2 = -42.25
  double const hand = -42.25;

  CHECK(expected == Catch::Approx(hand));
  CHECK(got_val == Catch::Approx(expected));
}

// The same equivalence, but reached through the BatchPolicy-gated dispatch
// entry (sequant::evaluate(Nodes const&, BatchPolicy const&, ...),
// scope_executor.hpp) with ordered_schedule_execution set -- the actual
// caller-facing seam SP3 Task 1 wires up, rather than calling
// evaluate_ordered_schedule directly as the test above does.
TEST_CASE(
    "BatchPolicy::ordered_schedule_execution routes through the ordered "
    "executor and matches forest descent",
    "[ordered-executor]") {
  std::vector<ScalarNode> forest{scalar_tree(L"2 * a * b - c"),
                                 scalar_tree(L"a * a + 3 * b - 2 * c")};
  ScalarLeafEvaluator const yield{{{L"a", 2.0}, {L"b", -3.5}, {L"c", 7.25}}};

  auto ref_cache = sequant::CacheManager<ScalarNode>::empty();
  ResultPtr const reference = sequant::evaluate(forest, yield, ref_cache);
  double const expected = reference->as<ResultScalar<double>>().value();

  sequant::BatchPolicy policy;
  policy.ordered_schedule_execution = true;

  auto cache = sequant::CacheManager<ScalarNode>::empty();
  ResultPtr const got = sequant::evaluate(
      forest, policy, ScalarEvalExpr::annot_t{}, yield, cache);
  double const got_val = got->as<ResultScalar<double>>().value();

  CHECK(got_val == Catch::Approx(expected));
  CHECK(got_val == Catch::Approx(-42.25));
}

// Flag-off byte-identical guard: with ordered_schedule_execution left at its
// default (false) and whole_scope_execution also false, the BatchPolicy-
// gated dispatch entry must take the FIRST pre-existing arm (an unconditional
// forward to sequant::evaluate(Nodes const&, layout, leaf_evaluator, cache))
// -- i.e. this task's new branch must not disturb either existing arm when
// its own gating flag is off.
TEST_CASE(
    "ordered_schedule_execution defaults to false and does not disturb the "
    "pre-existing BatchPolicy dispatch arms",
    "[ordered-executor]") {
  std::vector<ScalarNode> forest{scalar_tree(L"2 * a * b - c"),
                                 scalar_tree(L"a * a + 3 * b - 2 * c")};
  ScalarLeafEvaluator const yield{{{L"a", 2.0}, {L"b", -3.5}, {L"c", 7.25}}};

  sequant::BatchPolicy const policy;  // both flags default false
  CHECK_FALSE(policy.ordered_schedule_execution);
  CHECK_FALSE(policy.whole_scope_execution);

  auto cache = sequant::CacheManager<ScalarNode>::empty();
  ResultPtr const got = sequant::evaluate(
      forest, policy, ScalarEvalExpr::annot_t{}, yield, cache);
  CHECK(got->as<ResultScalar<double>>().value() == Catch::Approx(-42.25));
}
