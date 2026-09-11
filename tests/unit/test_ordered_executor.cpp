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
#include <SeQuant/core/eval/backends/dryrun/eval_expr.hpp>
#include <SeQuant/core/eval/backends/dryrun/meter.hpp>
#include <SeQuant/core/eval/backends/dryrun/result.hpp>
#include <SeQuant/core/eval/backends/dryrun/size_regime.hpp>
#include <SeQuant/core/eval/cache_manager.hpp>
#include <SeQuant/core/eval/cell_table_builder.hpp>
#include <SeQuant/core/eval/dag_scope.hpp>
#include <SeQuant/core/eval/eval.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/eval/legality.hpp>
#include <SeQuant/core/eval/node_batch_annotation.hpp>
#include <SeQuant/core/eval/ordered_executor.hpp>
#include <SeQuant/core/eval/ordered_schedule.hpp>
#include <SeQuant/core/eval/peak_monitor.hpp>
#include <SeQuant/core/eval/peak_profile.hpp>
#include <SeQuant/core/eval/result.hpp>
#include <SeQuant/core/eval/scope_executor.hpp>
#include <SeQuant/core/eval/scope_schedule.hpp>
#include <SeQuant/core/eval/slicing_signature.hpp>
#include <SeQuant/core/eval/value_node_map.hpp>
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
#include <array>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
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

// Explicit value cells (SP4 Task 3): evaluate_ordered_schedule now builds and
// statically validates a cell table before evaluating a schedule at all --
// same fixture as the case above, corrupting nothing, so the run must succeed
// and the test-facing diagnostic must report the table's own cell count.
TEST_CASE("ordered executor asserts the cell table before evaluating",
          "[ordered][cell_table]") {
  std::vector<ScalarNode> forest{scalar_tree(L"2 * a * b - c"),
                                 scalar_tree(L"a * a + 3 * b - 2 * c")};

  ScalarLeafEvaluator const yield{{{L"a", 2.0}, {L"b", -3.5}, {L"c", 7.25}}};

  sequant::BatchPolicy const policy;

  SizeRegime const regime;
  CostModel const cm{regime};
  auto const block_of = [](Index const&) -> std::size_t { return 1; };
  RichSchedule const rich = compute_dag_boulevard(forest, cm, block_of);
  auto const legality = analyze_legality(rich, forest, policy);
  OrderedSchedule const ordered =
      build_ordered_schedule(rich, legality, policy, {});

  auto ordered_cache = sequant::CacheManager<ScalarNode>::empty();
  std::function<std::size_t(Index const&)> const target =
      [](Index const&) -> std::size_t { return 1; };
  ResultPtr const got = evaluate_ordered_schedule(forest, ordered, rich,
                                                  ScalarEvalExpr::annot_t{},
                                                  yield, ordered_cache, target);
  REQUIRE(got);
  CHECK(sequant::eval::detail::ordered_last_cell_table_size() > 0);
}

// The same equivalence, but reached through the BatchPolicy-gated dispatch
// entry (sequant::evaluate(Nodes const&, BatchPolicy const&, ...),
// scope_executor.hpp) with policy.scheduler == BatchScheduler::ordered -- the
// actual caller-facing seam SP3 Task 1 wires up, rather than calling
// evaluate_ordered_schedule directly as the test above does.
TEST_CASE(
    "BatchPolicy::scheduler == BatchScheduler::ordered routes through the "
    "ordered executor and matches forest descent",
    "[ordered-executor]") {
  std::vector<ScalarNode> forest{scalar_tree(L"2 * a * b - c"),
                                 scalar_tree(L"a * a + 3 * b - 2 * c")};
  ScalarLeafEvaluator const yield{{{L"a", 2.0}, {L"b", -3.5}, {L"c", 7.25}}};

  auto ref_cache = sequant::CacheManager<ScalarNode>::empty();
  ResultPtr const reference = sequant::evaluate(forest, yield, ref_cache);
  double const expected = reference->as<ResultScalar<double>>().value();

  sequant::BatchPolicy policy;
  policy.scheduler = sequant::BatchScheduler::ordered;

  auto cache = sequant::CacheManager<ScalarNode>::empty();
  ResultPtr const got = sequant::evaluate(
      forest, policy, ScalarEvalExpr::annot_t{}, yield, cache);
  double const got_val = got->as<ResultScalar<double>>().value();

  CHECK(got_val == Catch::Approx(expected));
  CHECK(got_val == Catch::Approx(-42.25));
}

// Flag-off byte-identical guard: with policy.scheduler left at its default
// (BatchScheduler::forest_descent), the BatchPolicy-gated dispatch entry must
// take the FIRST pre-existing arm (an unconditional forward to
// sequant::evaluate(Nodes const&, layout, leaf_evaluator, cache)) -- i.e.
// this task's new branch must not disturb either existing arm when the
// scheduler is not set to ordered or whole_scope.
TEST_CASE(
    "BatchScheduler defaults to forest_descent and does not disturb the "
    "pre-existing BatchPolicy dispatch arms",
    "[ordered-executor]") {
  std::vector<ScalarNode> forest{scalar_tree(L"2 * a * b - c"),
                                 scalar_tree(L"a * a + 3 * b - 2 * c")};
  ScalarLeafEvaluator const yield{{{L"a", 2.0}, {L"b", -3.5}, {L"c", 7.25}}};

  sequant::BatchPolicy const policy;  // scheduler defaults to forest_descent
  CHECK(policy.scheduler == sequant::BatchScheduler::forest_descent);

  auto cache = sequant::CacheManager<ScalarNode>::empty();
  ResultPtr const got = sequant::evaluate(
      forest, policy, ScalarEvalExpr::annot_t{}, yield, cache);
  CHECK(got->as<ResultScalar<double>>().value() == Catch::Approx(-42.25));
}

// ===========================================================================
// SP3 Task 4: the acceptance payoff. On the real water-20 CSV-CCSD doubles
// residual (DF/aux-only Κ batching), all three executors -- forest descent,
// whole-scope, and the ORDERED executor -- must model the SAME final result
// shape (Step 1 / resolution R1: this is a zero-data DryRun forest, so
// "numerical equivalence" here is result-shape/extent equivalence, NOT
// floating-point arithmetic; the real-FP equivalence of the ordered executor
// is already proven at small scale by the [ordered-executor] scalar tests
// above and the [eval][ordered-executor] real-TA tests). AND the Κ-free
// home={} composite (I(i,i;a,a)-shaped: a Reduction that sums Κ away at its
// own node yet is homed at the ROOT scope -- placed by build_ordered_schedule
// as a plain root-level BuildStep, never inside the {Κ} child block) must
// build EXACTLY ONCE under the ordered executor and under forest descent
// (CSE), versus the whole-scope executor's per-block rebuild of the same
// root-homed composite -- the regression the whole ordered-scope effort exists
// to eliminate (Step 2 / R4).
//
// `[.]` hidden (run-by-name): a multi-second DryRun optimize+binarize of ~40
// residual terms, mirroring test_scope_executor.cpp's
// [scope-executor-witness-water20]. The fixture construction below is
// replicated (with an `orderedexec_` prefix) from that witness and from
// test_ordered_schedule.cpp's water-20 fixture rather than shared via a header:
// the sibling test files already duplicate it deliberately under distinct
// prefixes to avoid a CMake UNITY_BUILD anonymous-namespace collision (see
// test_ordered_schedule.cpp's own note at its Task-3 section), and a shared
// header consumed by only this new file would not reduce the existing two-file
// duplication without a risky behavior-changing refactor of both (R5's
// fallback path). The executor-DRIVE logic itself is NOT duplicated -- it lives
// in ordered_executor.hpp / scope_executor.hpp; only the forest/rich/policy
// construction is.
// ===========================================================================

namespace {

std::string orderedexec_witness_slurp(std::string const& path) {
  std::ifstream in(path);
  std::stringstream ss;
  ss << in.rdbuf();
  return ss.str();
}

struct OrderedExecWater20ProblemSize {
  std::size_t mu_tilde;
  std::size_t aux;
  std::size_t i_occ;
  std::array<double, 5> pno_M;
  std::array<double, 5> osv_M;
};

// Same (molecule, basis, parameter-set) size as test_ordered_schedule.cpp's
// kOrderedSchedWater20_pVDZF12 / test_scope_executor.cpp's
// kWitnessWater20_pVDZF12 (job-log / mpqc PaoPnoRMP2 moments) -- duplicated
// per the file-header note above.
inline constexpr OrderedExecWater20ProblemSize kOrderedExecWater20_pVDZF12{
    /*mu_tilde=*/896u,
    /*aux=*/1682u,
    /*i_occ=*/80u,
    /*pno_M=*/
    {1.0, 23.175775480059084, 25.865548281212597, 28.171416142614103,
     30.03848680550367},
    /*osv_M=*/
    {1.0, 58.987499999999997, 59.289227520688783, 59.584437469011633,
     59.872014818179686}};

// C60 pVDZ-F12 problem size (copied from test_eval_dryrun.cpp's kC60_pVDZF12).
inline constexpr OrderedExecWater20ProblemSize kOrderedExecC60_pVDZF12{
    /*mu_tilde=*/1800u,
    /*aux=*/4320u,
    /*i_occ=*/120u,
    /*pno_M=*/
    {1.0, 42.029069767441861, 46.039206412923569, 49.766252354482994,
     53.151291880343109},
    /*osv_M=*/
    {1.0, 148.25, 155.04434849422921, 161.33527408797721, 166.85553430303926}};

sequant::eval::dryrun::SizeRegime orderedexec_witness_df_regime(
    OrderedExecWater20ProblemSize const& p) {
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

sequant::ExprPtr orderedexec_witness_flatten_product(
    sequant::ExprPtr const& e) {
  if (!e->is<sequant::Product>()) return e;
  auto const& p = e->as<sequant::Product>();
  return sequant::ex<sequant::Product>(p.scalar(), p.factors(),
                                       sequant::Product::Flatten::Yes);
}

// Total builds (summed over slices) of one specific node in a build tally
// (CacheManager::recompute_tally()) -- identical helper to the one every
// build-once witness in test_eval_dryrun.cpp / test_scope_executor.cpp defines
// locally.
template <typename Tally>
std::size_t orderedexec_builds_of(
    Tally const& tally, sequant::eval::dryrun::EvalNodeDryRun const& n) {
  auto it = tally.find(n);
  if (it == tally.end()) return 0;
  std::size_t b = 0;
  for (auto const& [sig, bc] : it->second.slices) b += bc.count;
  return b;
}

///
/// \brief A \c CellTableInputs::operands_of callback: the direct operand
/// value ids of a value WITH REPETITION, one per leg of its production tree.
///
/// \details Read off the canonical forest node's two children and resolved
/// back to value ids through their value keys (\c value_key_of). The dependency
/// graph
/// (\c ordered_schedule_dep_graph, and \c OrderedSchedule::operand_vids
/// copied from it) de-duplicates its operand lists, so a value contracted
/// with itself would otherwise contribute ONE read where the runtime performs
/// two home accesses.
/// \note Captures \p rich and \p vmap by reference: both must outlive the
/// returned callable.
///
template <typename NodeT>
std::function<sequant::container::svector<std::size_t>(std::size_t)>
orderedexec_per_leg_operands(
    sequant::eval::RichSchedule const& rich,
    std::unordered_map<std::size_t, NodeT> const& vmap) {
  // Keyed by VALUE id (value_key_of): a child node resolves to the value it
  // is an occurrence of, not to the first value sharing its node hash.
  auto vid_of_key =
      std::make_shared<std::unordered_map<std::size_t, std::size_t>>();
  for (auto const& vc : rich.cells)
    vid_of_key->emplace(sequant::eval::value_key_of(vc), vc.value_id);
  return [&rich, &vmap, vid_of_key](
             std::size_t vid) -> sequant::container::svector<std::size_t> {
    sequant::container::svector<std::size_t> out;
    auto const it = vmap.find(sequant::eval::value_key_of(rich.cells[vid]));
    if (it == vmap.end() || it->second.leaf()) return out;
    auto const add = [&](NodeT const& child) {
      auto const f = vid_of_key->find(sequant::value_key_of(child));
      if (f != vid_of_key->end()) out.push_back(f->second);
    };
    add(it->second.left());
    add(it->second.right());
    return out;
  };
}

std::optional<std::size_t> orderedexec_index_of_build_step(
    sequant::eval::ScopeBlock const& block, std::size_t value_id) {
  for (std::size_t i = 0; i < block.steps.size(); ++i)
    if (auto const* b =
            std::get_if<sequant::eval::BuildStep>(&block.steps[i].value))
      if (b->value_id == value_id) return i;
  return std::nullopt;
}

}  // namespace

TEST_CASE(
    "ordered-executor witness: water-20 aux-only residual -- forest descent, "
    "whole-scope, and the ordered executor model the same result shape, and "
    "the Κ-free home={} composite builds exactly once under the ordered "
    "executor (vs the whole-scope per-block rebuild)",
    "[.][ordered-executor-witness-water20]") {
  using sequant::eval::dryrun::EvalExprDryRun;
  using sequant::eval::dryrun::EvalNodeDryRun;
  using Node = EvalNodeDryRun;

  auto ctx = sequant::get_default_context().clone();
  ctx.set_first_dummy_index_ordinal(1000000);
  auto isr = ctx.mutable_index_space_registry();
  REQUIRE(isr != nullptr);
  sequant::mbpt::add_pao_spaces(isr, sequant::mbpt::Spin::any);
  sequant::mbpt::add_df_spaces(isr);
  auto ctx_resetter = sequant::set_scoped_default_context(std::move(ctx));

  auto const body =
      orderedexec_witness_slurp(std::string(SEQUANT_UNIT_TESTS_SOURCE_DIR) +
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

  std::size_t nterms = std::min<std::size_t>(summands.size(), 40);
  if (char const* nt = std::getenv("SEQUANT_UT_DRYRUN_NTERMS"))
    nterms = std::min<std::size_t>(summands.size(), std::atoll(nt));

  auto regime = orderedexec_witness_df_regime(kOrderedExecWater20_pVDZF12);
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

  std::vector<Node> forest;
  for (std::size_t s = 0; s < nterms; ++s) {
    sequant::ExprPtr const term =
        orderedexec_witness_flatten_product(summands[s]);
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

  // ONE forest + ONE rich shared by all three pipelines (R2).
  auto const block_of = [](sequant::Index const&) -> std::size_t {
    return 256;
  };
  auto const rich = sequant::eval::compute_dag_boulevard(forest, *cm, block_of);
  REQUIRE(!rich.cells.empty());

  // ---- Identify the Κ-free home={} composite (I(i,i;a,a)-shaped), R3:
  // the composite that CONSUMES a Κ-reduction result (a non-leaf that
  // contracts Κ at its own node but does not carry Κ free) but does not
  // itself carry or contract Κ -- placed by build_ordered_schedule as a plain
  // root-level BuildStep. Same identification as test_ordered_schedule.cpp's
  // water-20 acceptance test (Target 1 -> its structural parent).
  auto const is_K = [](sequant::Index const& ix) {
    return ix.space().base_key() == L"Κ";
  };
  auto const carries_type =
      [&](sequant::container::svector<sequant::Index> const& v,
          auto const& pred) { return std::any_of(v.begin(), v.end(), pred); };
  auto const vmap = sequant::eval::build_value_node_map(forest);

  std::optional<std::size_t> mu_mu_hash;
  {
    auto const is_mu_mu_pair = [](auto const& carried) {
      return carried.size() == 2 &&
             std::all_of(carried.begin(), carried.end(),
                         [](sequant::Index const& ix) {
                           return ix.space().base_key() == L"μ̃";
                         });
    };
    for (bool const require_mu_mu : {true, false}) {
      if (mu_mu_hash) break;
      for (auto const& vc : rich.cells) {
        auto const it = vmap.find(sequant::eval::value_key_of(vc));
        if (it == vmap.end() || it->second.leaf()) continue;
        if (carries_type(vc.carried, is_K)) continue;
        auto const contracted = sequant::contracted_indices(it->second);
        if (std::find_if(contracted.begin(), contracted.end(), is_K) ==
            contracted.end())
          continue;
        if (require_mu_mu && !is_mu_mu_pair(vc.carried)) continue;
        mu_mu_hash = vc.hash;
        break;
      }
    }
  }
  REQUIRE(mu_mu_hash.has_value());

  std::optional<Node> parent;  // the Κ-free composite consuming the reduction
  {
    std::function<void(Node const&)> find_parent = [&](Node const& n) {
      if (parent || n.leaf()) return;
      if (n.left()->hash_value() == *mu_mu_hash ||
          n.right()->hash_value() == *mu_mu_hash) {
        parent = n;
        return;
      }
      find_parent(n.left());
      find_parent(n.right());
    };
    for (auto const& tree : forest) {
      find_parent(tree);
      if (parent) break;
    }
  }
  REQUIRE(parent.has_value());
  Node const I_node = *parent;
  auto const parent_hash = I_node->hash_value();
  auto const parent_cell =
      std::find_if(rich.cells.begin(), rich.cells.end(),
                   [&](auto const& vc) { return vc.hash == parent_hash; });
  REQUIRE(parent_cell != rich.cells.end());
  std::size_t const parent_value_id = parent_cell->value_id;

  // The composite must be Κ-free at its own node (no Κ carried, no Κ
  // contracted) -- i.e. genuinely home={} (root-homed), not a Κ-loop-local.
  CHECK_FALSE(carries_type(parent_cell->carried, is_K));
  {
    auto const pc = sequant::contracted_indices(I_node);
    CHECK(std::find_if(pc.begin(), pc.end(), is_K) == pc.end());
  }

  using annot_t = std::remove_cvref_t<decltype(forest.front()->annot())>;
  annot_t const layout{};
  sequant::eval::dryrun::DryRunLeafEvaluator const yield{cm};
  std::function<std::size_t(sequant::Index const&)> const target =
      [](sequant::Index const&) -> std::size_t { return 256; };

  auto& logger = sequant::Logger::instance();
  auto const prev_level = logger.eval.level;
  auto* const prev_stream = logger.eval.stream;
  logger.eval.level = 1;  // arms tally_build (DryRunOps::prod's runtime gate)

  auto aops = sequant::eval::dryrun::make_dryrun_array_ops(cm);

  // ---- (1) ordered executor. ----
  auto const legality = sequant::eval::analyze_legality(rich, forest, policy);
  REQUIRE(legality.cells.size() == rich.cells.size());
  auto const ordered =
      sequant::eval::build_ordered_schedule(rich, legality, policy, {L"Κ"});
  REQUIRE(sequant::eval::well_formed(ordered));
  // R3: confirm the composite is a plain ROOT-level BuildStep (home={}), never
  // inside the {Κ} child block.
  REQUIRE(orderedexec_index_of_build_step(ordered.root, parent_value_id)
              .has_value());

  std::ostringstream ord_trace;
  logger.eval.stream = &ord_trace;
  auto ordered_cache = sequant::cache_manager(forest);
  ordered_cache.set_array_ops(&aops);
  ordered_cache.set_recompute_tally_enabled(true);
  // R1/R2: install the hierarchy-wide PeakMonitor on the ROOT cache only; it
  // propagates to every per-batch scratch via peak_monitor()'s parent_
  // fallthrough (the scratch caches set_parent to this root), so note_working_
  // set() calls anywhere in the ordered walk fold into ONE high-water mark. No
  // executor wiring -- the install lives entirely here.
  sequant::eval::PeakMonitor ord_mon;
  ordered_cache.set_peak_monitor(&ord_mon);
  // Task 4: the NODE-level lift of policy.is_volatile_leaf (mirrors
  // make_evaluator's is_volatile_node lift, eval.hpp): a leaf tensor labeled
  // "t" is volatile, every internal node non-volatile. Threaded into the
  // ordered executor so it classifies each root-homed composite
  // volatile-vs-persistent and eagerly releases the volatile ones -- the
  // reclaim this witness pins.
  std::function<bool(Node const&)> const is_volatile_node =
      [p = policy.is_volatile_leaf](Node const& n) -> bool {
    if (!n.leaf() || !n->is_tensor()) return false;
    return p && p(n->as_tensor());
  };
  ResultPtr ord_result;
  try {
    ord_result = sequant::eval::evaluate_ordered_schedule<sequant::Trace::On>(
        forest, ordered, rich, layout, yield, ordered_cache, target, {},
        is_volatile_node);
  } catch (std::exception const& e) {
    std::cerr << "[ordered-executor-witness-water20] ordered evaluate threw: "
              << e.what() << "\n";
  }

  // ---- (2) whole-scope executor, SAME forest/rich. ----
  std::ostringstream ws_trace;
  logger.eval.stream = &ws_trace;
  auto const sched =
      sequant::eval::build_scope_schedule<std::wstring>(rich, {L"Κ"});
  auto ws_cache = sequant::cache_manager(forest);
  ws_cache.set_array_ops(&aops);
  ws_cache.set_recompute_tally_enabled(true);
  // R2: same measurement on the whole-scope executor's own root cache, so the
  // ordered-vs-whole-scope peak comparison below is apples-to-apples (both
  // realized peaks, both via a root PeakMonitor, both under Trace::On).
  sequant::eval::PeakMonitor ws_mon;
  ws_cache.set_peak_monitor(&ws_mon);
  ResultPtr ws_result;
  try {
    ws_result = sequant::eval::evaluate_whole_scope<sequant::Trace::On>(
        forest, sched, rich, layout, yield, ws_cache, target);
  } catch (std::exception const& e) {
    std::cerr << "[ordered-executor-witness-water20] whole-scope evaluate "
                 "threw: "
              << e.what() << "\n";
  }

  // ---- (3) forest descent (reference), SAME forest, plain CSE cache. ----
  std::ostringstream fd_trace;
  logger.eval.stream = &fd_trace;
  auto fd_cache = sequant::cache_manager(forest);
  fd_cache.set_recompute_tally_enabled(true);
  // Reference peak too: forest descent also builds-once (CSE) and keeps shared
  // composites resident, so its realized peak is the natural build-once
  // baseline the ordered executor's peak should track (both keep the Kappa-free
  // composite resident once) -- reported alongside for context.
  sequant::eval::PeakMonitor fd_mon;
  fd_cache.set_peak_monitor(&fd_mon);
  ResultPtr fd_result;
  try {
    fd_result =
        sequant::evaluate<sequant::Trace::On>(forest, layout, yield, fd_cache);
  } catch (std::exception const& e) {
    std::cerr << "[ordered-executor-witness-water20] forest-descent evaluate "
                 "threw: "
              << e.what() << "\n";
  }

  logger.eval.level = prev_level;
  logger.eval.stream = prev_stream;

  std::size_t const ord_builds =
      orderedexec_builds_of(ordered_cache.recompute_tally(), I_node);
  std::size_t const ws_builds =
      orderedexec_builds_of(ws_cache.recompute_tally(), I_node);
  std::size_t const fd_builds =
      orderedexec_builds_of(fd_cache.recompute_tally(), I_node);

  std::wcerr << L"\n=== [ordered-executor-witness-water20] water-20 aux-only, "
             << forest.size() << L" terms ===\n"
             << L"  ordered replay completed     = "
             << (ord_result ? L"yes" : L"NO") << L"\n"
             << L"  whole-scope replay completed = "
             << (ws_result ? L"yes" : L"NO") << L"\n"
             << L"  forest-descent replay done   = "
             << (fd_result ? L"yes" : L"NO") << L"\n"
             << L"  final result size (bytes): ordered = "
             << (ord_result ? double(ord_result->size_in_bytes()) : -1.0)
             << L", whole-scope = "
             << (ws_result ? double(ws_result->size_in_bytes()) : -1.0)
             << L", forest-descent = "
             << (fd_result ? double(fd_result->size_in_bytes()) : -1.0) << L"\n"
             << L"  Κ-free home={} composite builds: ordered = " << ord_builds
             << L", whole-scope = " << ws_builds << L", forest-descent = "
             << fd_builds << L"\n";

  // Systematic contrast over EVERY root-homed (root-level BuildStep), Κ-free
  // composite -- not just the single I_node above: the ordered executor homes
  // each at the root cache and builds it exactly once, while the whole-scope
  // executor rebuilds the busiest one more than once (its per-block
  // regression).
  std::size_t worst_ord = 0, worst_ws = 0;
  for (auto const& vc : rich.cells) {
    auto const vit = vmap.find(sequant::eval::value_key_of(vc));
    if (vit == vmap.end() || vit->second.leaf()) continue;
    if (carries_type(vc.carried, is_K)) continue;
    auto const pc = sequant::contracted_indices(vit->second);
    if (std::find_if(pc.begin(), pc.end(), is_K) != pc.end()) continue;
    if (!orderedexec_index_of_build_step(ordered.root, vc.value_id).has_value())
      continue;
    worst_ord = std::max(
        worst_ord,
        orderedexec_builds_of(ordered_cache.recompute_tally(), vit->second));
    worst_ws = std::max(worst_ws, orderedexec_builds_of(
                                      ws_cache.recompute_tally(), vit->second));
  }
  std::wcerr << L"  worst build count over ALL root-homed Κ-free composites: "
             << L"ordered = " << worst_ord << L", whole-scope = " << worst_ws
             << L"\n";

  // ---- Step 1 (R1): all three model the SAME final result shape. ----
  REQUIRE(ord_result);
  REQUIRE(ws_result);
  REQUIRE(fd_result);
  CHECK(ord_result->size_in_bytes() == fd_result->size_in_bytes());
  CHECK(ws_result->size_in_bytes() == fd_result->size_in_bytes());

  // ---- Step 2 (R4): the build-count acceptance. The ordered executor
  // builds the Κ-free root-homed composite EXACTLY ONCE (home store/lookup
  // de-aliases it to the root scope), and forest descent builds it once too
  // (CSE) -- the sanity anchor. The whole-scope executor rebuilds the SAME
  // root-homed composite MORE than once (its per-block regression), which the
  // ordered executor is designed to eliminate. If ord_builds != 1, that is a
  // REAL ordered_executor.hpp bug to fix, NOT an assertion to weaken.
  CHECK(fd_builds == 1);   // reference: CSE build-once
  CHECK(ord_builds == 1);  // THE acceptance: ordered build-once
  CHECK(ws_builds > 1);    // contrast: whole-scope per-block rebuild

  // The build-once property holds for EVERY root-homed Κ-free composite (the
  // systematic form of the acceptance), and the whole-scope per-block rebuild
  // contrast holds for the busiest of them too.
  CHECK(worst_ord == 1);
  CHECK(worst_ws > 1);

  // ---- Step 1 (R1/R2): realized-peak monitoring. All three peaks are measured
  // the SAME way (a PeakMonitor on each executor's root cache, every run under
  // Trace::On, all on THIS shared forest), emitted for context.
  //
  // The peak assertion is vs the BUILD-ONCE baseline (forest descent, pure
  // CSE), NOT vs whole-scope. Whole-scope's per-block REBUILD is a recompute /
  // peak-MINIMIZATION strategy: it drops each rebuild transiently, so fewer
  // values are co-resident and its peak is a LOWER bound (the recompute floor),
  // NOT an upper bound a build-once executor should be expected to meet. Both
  // build-once executors (forest descent and ordered) instead keep shared
  // composites resident across the whole walk -- the space/time tradeoff. The
  // right invariant is that ordered's build-once homing does not raise peak
  // over the build-once REFERENCE (forest descent), which it satisfies (and
  // here strictly improves on). Per-composite home-vs-rebuild peak control
  // (choosing which composites to home vs recompute to cap peak) is SP4's
  // cost-driven job; SP3 executes whatever build_ordered_schedule produces.
  std::wcerr << L"  realized peak (bytes): ordered = " << ord_mon.hwmark_bytes
             << L", whole-scope = " << ws_mon.hwmark_bytes
             << L", forest-descent = " << fd_mon.hwmark_bytes << L"\n";
  INFO("realized peak (bytes): ordered = "
       << ord_mon.hwmark_bytes << ", whole-scope = " << ws_mon.hwmark_bytes
       << ", forest-descent = " << fd_mon.hwmark_bytes);
  CAPTURE(ord_mon.hwmark_bytes, ws_mon.hwmark_bytes, fd_mon.hwmark_bytes);
  // R1: if the ordered high-water is 0, the ordered build path never called
  // note_working_set() where the reference does -- a real executor gap, not a
  // test artifact. Under Trace::On every evaluate_impl / combine step notes its
  // working set, so a nonzero mark confirms the monitor observes the ordered
  // walk.
  CHECK(ord_mon.hwmark_bytes > 0);
  CHECK(ws_mon.hwmark_bytes > 0);
  CHECK(fd_mon.hwmark_bytes > 0);
  // R2 (ruling: assert vs the build-once baseline): the ordered executor's
  // build-once homing must not raise peak over forest descent's build-once CSE.
  CHECK(ord_mon.hwmark_bytes <= fd_mon.hwmark_bytes);

  // ---- Task 4 (THE PAYOFF): eager release of volatile homed values.
  //
  // With is_volatile_node threaded in, every root-homed composite whose subtree
  // carries a volatile ("t"-labeled) leaf is homed NON-persistent -- released
  // at its genuine last cross-block use instead of pinned resident for the
  // whole walk. That reclaim drops the realized ordered peak strictly below the
  // pinned baseline (988'732'399'293 == the peak WITHOUT eager release, the
  // number this same witness reports when is_volatile is NOT threaded), while
  // build-once (ord_builds == 1, worst_ord == 1, asserted above) and numerical/
  // shape equivalence (ord vs forest-descent, ord vs whole-scope, above) hold:
  // a value is released only after its true last use, never rebuilt.
  //
  // Attribute the reclaim: split the root-homed composites into the volatile
  // floor (released) and the persistent floor (held) by the SAME subtree_any
  // classification the executor applies at its homing sites, and sum each
  // group's modeled home footprint (cell_footprint, the peak-composition
  // metric this file already uses elsewhere).
  auto const foot = [&](sequant::eval::ValueCell const& vc) -> std::size_t {
    return sequant::eval::detail::cell_footprint(vc.carried, vc.home_modes, *cm,
                                                 block_of);
  };
  std::size_t vol_bytes = 0, persist_bytes = 0;
  std::size_t vol_count = 0, persist_count = 0;
  for (auto const& vc : rich.cells) {
    auto const vit = vmap.find(sequant::eval::value_key_of(vc));
    if (vit == vmap.end() || vit->second.leaf()) continue;
    // only the composites the ROOT walk actually homes (root-level BuildSteps)
    if (!orderedexec_index_of_build_step(ordered.root, vc.value_id).has_value())
      continue;
    if (sequant::subtree_any(vit->second, is_volatile_node)) {
      ++vol_count;
      vol_bytes += foot(vc);
    } else {
      ++persist_count;
      persist_bytes += foot(vc);
    }
  }
  std::wcerr << L"  root-homed composites: volatile = " << vol_count << L" ("
             << vol_bytes << L" B), persistent = " << persist_count << L" ("
             << persist_bytes << L" B)\n";
  INFO("root-homed composites: volatile = "
       << vol_count << " (" << vol_bytes << " B); persistent = "
       << persist_count << " (" << persist_bytes << " B)");
  CAPTURE(vol_count, vol_bytes, persist_count, persist_bytes);
  // At least one composite must actually be classified volatile, else the
  // reclaim below is vacuous and the seed is not taking effect.
  CHECK(vol_count > 0);
  // THE acceptance: eager release drops the realized peak strictly below the
  // pinned no-eager-release baseline. If this FAILS at the baseline value, the
  // seed is not wired; if peak drops but ord_builds/worst_ord break, a home
  // value was released before its last use (a count bug to fix, NOT an
  // assertion to weaken).
  CHECK(ord_mon.hwmark_bytes < 988732399293ull);

  // ---- Step 1 (R3): test-side run-completeness cross-check, complementing
  // the SEQUANT_ASSERT inside evaluate_ordered_schedule (a no-op in this
  // IGNORE-assert release build; ACTIVE in cmake-build-debug, where the
  // executor-side `built` ledger -- set at every real production site, product
  // or sum -- is the AUTHORITATIVE completeness gate; it is proven to hold on
  // this witness in debug). This proxy cross-checks completeness from OUTSIDE
  // the executor against the forest-descent (pure-CSE) reference tally: the
  // ordered executor must build, at least once, every PRODUCT the reference
  // run builds at least once. Comparing to the reference tally (rather than
  // asserting >=1 outright) is the correct calibration -- a scheduled product
  // that CSE folds to a home/cache alias is built ONCE under some other cell's
  // identity and served as a dedup hit thereafter, so recompute_tally() shows
  // 0 for THAT vid under BOTH runs; that is not a gap. A vid the reference
  // builds but ordered does not (ord 0, fd >=1) WOULD be a real skip. (Only
  // products are compared: tally_build records only the contraction branch of
  // evaluate_impl -- a Sum-typed node is summed, never tallied.)
  {
    sequant::container::vector<std::size_t> prod_ids;
    sequant::eval::detail::collect_production_ids(ordered.root, prod_ids);
    std::size_t checked = 0, ord_missing_vs_fd = 0, dedup_alias = 0;
    for (std::size_t const vid : prod_ids) {
      REQUIRE(vid < ordered.num_values);
      auto const it = vmap.find(sequant::eval::value_key_of(rich.cells[vid]));
      REQUIRE(it != vmap.end());
      if (it->second.leaf()) continue;  // leaves are handed back, never built
      if (!it->second->is_product()) continue;  // only products are tallied
      auto const ord_b =
          orderedexec_builds_of(ordered_cache.recompute_tally(), it->second);
      auto const fd_b =
          orderedexec_builds_of(fd_cache.recompute_tally(), it->second);
      if (fd_b >= 1) {
        ++checked;
        if (ord_b < 1) ++ord_missing_vs_fd;
      } else if (ord_b < 1) {
        ++dedup_alias;  // built by neither as a distinct product (CSE alias)
      }
    }
    std::wcerr << L"  completeness cross-check: " << checked
               << L" reference-built products, ordered-missing = "
               << ord_missing_vs_fd << L", CSE-alias (exempt) = " << dedup_alias
               << L"\n";
    CHECK(checked > 0);
    CHECK(ord_missing_vs_fd == 0);
  }

  // ---- R4: confirm SP2's build-time non-innermost forced-split tripwire
  // stands (ordered_schedule.hpp:839-860). That SEQUANT_ASSERT in
  // build_ordered_schedule prevents a non-innermost forced-split
  // OrderedSchedule from EVER being constructed, so evaluate_ordered_schedule
  // (which takes no `legality` and carries no runtime split marker) need not --
  // and cannot -- re-derive it; the executor's own loud guards (exhaustive
  // Step-variant, OutputKind, and BatchModeType dispatch, all
  // SEQUANT_ASSERT(false, ...) in ordered_executor.hpp) refuse any OTHER
  // unsupported construct. The schedule built here is well_formed (asserted
  // above) and executed to completion.
  CHECK(sequant::eval::well_formed(ordered));
}

// ===========================================================================
// AUX+OCC dry-run WALK reproducer (KNOWN-FAILING, opt-in [.] tag): the MPQC w20
// CSV-CCk residual with occ batched as an EXTERNAL index (plus Κ aux
// contracted, spectator batching, and node-level placement -- the exact
// make_csv_batch_policy(occ_target>0) config) builds a well_formed schedule,
// but WALKING it with the zero-data DryRun backend trips the SAME failure class
// the wet w20 run aborts on. The dry-run walk exercises the identical
// run_ordered_contracted_block / evaluate_impl home-read + scatter logic as the
// wet run -- only the TA tile math is stubbed -- so it reproduces the
// schedule/slice defect LOCALLY in seconds (no cluster round-trip). This is the
// vehicle for the frame-correct-slicing / multi-level-escape design work; see
// doc/dev/specs/2026-08-31-occ-use-induced-slicing-and-escape-chain-design.md.
//
// It FAILS today by design (REQUIRE_NOTHROW reproduces the open bug), so it is
// tagged [.] (hidden -- run by name, e.g. `unit_tests-sequant
// [w20-auxocc-walk]`) to keep it out of the default suite until the design fix
// lands, at which point the [.] is dropped and REQUIRE_NOTHROW becomes the
// acceptance. The distinct sliced_modes-duplication defect this reproducer
// first surfaced IS already fixed (lifetime_mask.hpp acc-dedup, w8-lossless);
// what remains here is the use-induced slicing of whole-produced shared
// operands and the multi-level escape chain -- both left to the design pass.
// ===========================================================================
namespace {
// SEQUANT_UT_SCHED_TREE: the block tree of an ordered schedule (depth, axis,
// kind, step count, escaped outputs), independent of whether the cell table
// later validates.
void ut_dump_sched_tree(sequant::eval::OrderedSchedule const& sched) {
  std::function<void(sequant::eval::ScopeBlock const&, int)> dump =
      [&](sequant::eval::ScopeBlock const& b, int depth) {
        std::size_t builds = 0, children = 0;
        for (auto const& st : b.steps) {
          if (std::get_if<sequant::eval::BuildStep>(&st.value)) ++builds;
          if (std::get_if<sequant::eval::ScopeBlock>(&st.value)) ++children;
        }
        std::wcerr << L"[sched-tree] " << std::wstring(2 * depth, L' ')
                   << L"depth=" << depth << L" axis="
                   << (b.axis ? b.axis.full_label() : L"<root>") << L" kind="
                   << (b.kind == sequant::BatchModeType::External
                           ? L"external"
                           : L"contracted")
                   << L" builds=" << builds << L" children=" << children
                   << L" outs=";
        for (auto const& [ovid, okind] : b.outputs)
          std::wcerr << ovid << L":"
                     << (okind == sequant::eval::OutputKind::AccumulateSum
                             ? L"sum"
                         : okind == sequant::eval::OutputKind::AccumulateScatter
                             ? L"scatter"
                             : L"other")
                     << L" ";
        std::wcerr << L"\n";
        for (auto const& st : b.steps)
          if (auto const* child =
                  std::get_if<sequant::eval::ScopeBlock>(&st.value))
            dump(*child, depth + 1);
      };
  dump(sched.root, 0);
}

// SEQUANT_UT_ROLE_DIAG=<vid>[,<vid>...]: per-axis legality roles and the loop
// identities (carried positions' loop_slot, reduced_slot, modes contracted in
// batches, loops opened), the operands and the consumer of every occurrence
// of those values.
void ut_dump_role_diag(sequant::eval::RichSchedule const& rich,
                       sequant::eval::LegalitySchedule const& legality,
                       char const* rd) {
  std::set<std::size_t> want;
  std::istringstream toks{rd};
  for (std::string tok; std::getline(toks, tok, ',');)
    if (!tok.empty()) want.insert(std::stoul(tok));
  auto role_name = [](sequant::eval::LoopRole r) -> wchar_t const* {
    switch (r) {
      case sequant::eval::LoopRole::LoopLocal:
        return L"LoopLocal";
      case sequant::eval::LoopRole::Reduction:
        return L"Reduction";
      case sequant::eval::LoopRole::LoopCarried:
        return L"LoopCarried";
      case sequant::eval::LoopRole::LoopInvariant:
        return L"LoopInvariant";
    }
    return L"?";
  };
  for (std::size_t v : want) {
    if (v >= rich.cells.size()) continue;
    auto const& cell = rich.cells[v];
    std::wcerr << L"[role-diag] vid=" << v << L" h=" << (cell.hash % 100000)
               << L" key=" << (sequant::eval::value_key_of(cell) % 100000)
               << L" same_hash_vids=";
    for (auto const& o : rich.cells)
      if (o.hash == cell.hash && o.value_id != v)
        std::wcerr << o.value_id << L",";
    std::wcerr << L" per_axis={";
    for (auto const& lc : legality.cells)
      if (lc.hash == sequant::eval::value_key_of(cell))
        for (auto const& ax : lc.per_axis)
          std::wcerr << ax.axis.full_label() << L":" << role_name(ax.role)
                     << L" ";
    std::wcerr << L"}\n";
    for (auto const& occ : cell.occurrences) {
      std::wcerr << L"[role-diag]   occ carried={";
      for (std::size_t p = 0; p < occ.carried.size(); ++p)
        std::wcerr << occ.carried[p].full_label() << L"#"
                   << (p < occ.loop_slot.size() ? occ.loop_slot[p] : -9)
                   << L" ";
      std::wcerr << L"} reduced_slot={";
      for (auto const& [ix, sl] : occ.reduced_slot)
        std::wcerr << ix.full_label() << L"#" << sl << L" ";
      std::wcerr << L"} contracted_batched={";
      for (auto const& ix : occ.contracted_batched)
        std::wcerr << ix.full_label() << L" ";
      std::wcerr << L"} opens={";
      for (auto const& [ix, k] : occ.opens)
        std::wcerr << ix.full_label()
                   << (k == sequant::BatchModeType::External ? L":E" : L":C")
                   << L" ";
      std::wcerr << L"} ectx={";
      for (auto const& [ix, rng] : occ.ectx)
        std::wcerr << ix.full_label() << L" ";
      std::wcerr << L"} home={";
      for (auto const& ix : occ.home) std::wcerr << ix.full_label() << L" ";
      std::wcerr << L"} point=" << occ.point << L" consumer_vid=";
      {
        bool found = false;
        for (std::size_t u = 0; u < rich.cells.size() && !found; ++u)
          for (auto const& co : rich.cells[u].occurrences)
            if (co.point == occ.consumer_point && co.point != occ.point) {
              std::wcerr << u << L"(h=" << (rich.cells[u].hash % 100000)
                         << L")";
              found = true;
              break;
            }
        if (!found) std::wcerr << L"<root>";
      }
      std::wcerr << L"\n";
      // Its operands: every occurrence (of any value) whose structural
      // consumer is this occurrence.
      for (std::size_t u = 0; u < rich.cells.size(); ++u)
        for (auto const& ch : rich.cells[u].occurrences) {
          if (ch.consumer_point != occ.point || ch.point == occ.point) continue;
          std::wcerr << L"[role-diag]     operand vid=" << u << L" h="
                     << (rich.cells[u].hash % 100000)
                     << (rich.cells[u].is_leaf ? L" leaf" : L"")
                     << L" carried={";
          for (std::size_t p = 0; p < ch.carried.size(); ++p)
            std::wcerr << ch.carried[p].full_label() << L"#"
                       << (p < ch.loop_slot.size() ? ch.loop_slot[p] : -9)
                       << L" ";
          std::wcerr << L"} home={";
          for (auto const& ix : ch.home) std::wcerr << ix.full_label() << L" ";
          std::wcerr << L"} ectx={";
          for (auto const& [ix, rng] : ch.ectx)
            std::wcerr << ix.full_label() << L" ";
          std::wcerr << L"}\n";
        }
    }
  }
}
}  // namespace

TEST_CASE(
    "ordered executor: water-20 aux+occ residual dry-run walk completes "
    "without "
    "a vanished home value",
    "[ordered][w20-auxocc-walk]") {
  using sequant::eval::dryrun::EvalExprDryRun;
  using sequant::eval::dryrun::EvalNodeDryRun;
  using Node = EvalNodeDryRun;

  auto ctx = sequant::get_default_context().clone();
  ctx.set_first_dummy_index_ordinal(1000000);
  auto isr = ctx.mutable_index_space_registry();
  REQUIRE(isr != nullptr);
  sequant::mbpt::add_pao_spaces(isr, sequant::mbpt::Spin::any);
  sequant::mbpt::add_df_spaces(isr);
  auto ctx_resetter = sequant::set_scoped_default_context(std::move(ctx));

  auto const body =
      orderedexec_witness_slurp(std::string(SEQUANT_UNIT_TESTS_SOURCE_DIR) +
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

  // FULL residual (all summands) to match the MPQC w20 run; overridable for
  // bisecting which term first breaks the walk.
  std::size_t nterms = summands.size();
  if (char const* nt = std::getenv("SEQUANT_UT_DRYRUN_NTERMS"))
    nterms = std::min<std::size_t>(summands.size(), std::atoll(nt));

  auto regime = orderedexec_witness_df_regime(kOrderedExecWater20_pVDZF12);
  auto cm = std::make_shared<sequant::eval::dryrun::CostModel const>(regime);

  // AUX+OCC: Κ batchable-contracted (aux), occ batchable-EXTERNAL; spectator
  // batching + node-level placement ON (make_csv_batch_policy, occ_target>0).
  sequant::BatchPolicy policy;
  policy.is_batchable_contracted_index = [](sequant::Index const& ix) {
    return ix.space().base_key() == L"Κ";
  };
  policy.is_batchable_external_index = [](sequant::Index const& ix) {
    auto const reg = sequant::get_default_context().index_space_registry();
    return reg && ix.space() && reg->is_pure_occupied(ix.space());
  };
  policy.batch_spectator_indices = true;
  policy.node_level_placement = true;
  policy.batch_target_size = [](sequant::Index const& ix) -> std::size_t {
    return ix.space().base_key() == L"Κ" ? 256 : 16;
  };
  policy.is_volatile_leaf = [](sequant::Tensor const& t) {
    return t.label() == L"t";
  };
  policy.accumulation_factor = 1.0;
  policy.persistent_only = false;
  policy.peak_threshold = 1e11;
  // Mirror a specific MPQC input (z820 w20 csv-cck-diag.json: batch
  // peak_threshold 25e9, optimize objective dense_time_space) via env.
  if (char const* pt = std::getenv("SEQUANT_UT_PEAK_THRESHOLD"))
    policy.peak_threshold = std::atof(pt);

  auto axes_map = std::make_shared<std::unordered_map<
      sequant::Expr const*,
      sequant::container::vector<sequant::NodeBatchAnnotation>>>();
  sequant::OptimizeOptions opts;
  opts.objective_function = sequant::ObjectiveFunction::DenseTimeSpaceBatched;
  if (char const* ob = std::getenv("SEQUANT_UT_OBJECTIVE")) {
    std::string const o{ob};
    if (o == "dense_time_space")
      opts.objective_function = sequant::ObjectiveFunction::DenseTimeSpace;
    else if (o == "dense_time_space_batched")
      opts.objective_function =
          sequant::ObjectiveFunction::DenseTimeSpaceBatched;
    else if (o == "dense_space_time")
      opts.objective_function = sequant::ObjectiveFunction::DenseSpaceTime;
    else if (o == "dense_space_time_batched")
      opts.objective_function =
          sequant::ObjectiveFunction::DenseSpaceTimeBatched;
  }
  opts.idx_to_extent = regime.idx_to_extent();
  opts.inner_pow = regime.inner_pow_fn();
  opts.batch_policy = policy;
  opts.volatile_weight = 20.0;
  opts.roofline.machine_balance = 200.0;
  opts.roofline.fast_mem_elems = 1000000.0;
  opts.term_batch_axes = axes_map;

  std::vector<Node> forest;
  for (std::size_t s = 0; s < nterms; ++s) {
    sequant::ExprPtr const term =
        orderedexec_witness_flatten_product(summands[s]);
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

  auto const block_of = [](sequant::Index const& ix) -> std::size_t {
    return ix.space().base_key() == L"Κ" ? 256 : 16;
  };
  auto const rich = sequant::eval::compute_dag_boulevard(forest, *cm, block_of);
  REQUIRE(!rich.cells.empty());

  auto const legality = sequant::eval::analyze_legality(rich, forest, policy);
  REQUIRE(legality.cells.size() == rich.cells.size());

  // MPQC passes an EMPTY mode_order (build_ordered_schedule derives the forced
  // split axes from the legality) -- match that.
  auto const ordered = sequant::eval::build_ordered_schedule(
      rich, legality, policy, std::initializer_list<std::wstring>{});
  if (std::getenv("SEQUANT_UT_SCHED_TREE")) ut_dump_sched_tree(ordered);
  if (char const* rd = std::getenv("SEQUANT_UT_ROLE_DIAG"))
    ut_dump_role_diag(rich, legality, rd);
  REQUIRE(sequant::eval::well_formed(ordered));

  if (std::getenv("SEQUANT_UT_SCHED_TREE")) {
    auto const vmap_dump = sequant::eval::build_value_node_map(forest);
    auto carried_str = [&](std::size_t vid) {
      std::string s;
      if (vid < rich.cells.size())
        for (auto const& x : rich.cells[vid].carried)
          s += sequant::toUtf8(x.full_label()) + " ";
      return s;
    };
    auto sliced_str = [&](std::size_t vid) {
      std::string s;
      if (vid < rich.cells.size()) {
        auto const it =
            vmap_dump.find(sequant::eval::value_key_of(rich.cells[vid]));
        if (it != vmap_dump.end())
          for (auto const& x : it->second->sliced_modes())
            s += sequant::toUtf8(x.full_label()) + " ";
      }
      return s;
    };
    auto home_str = [&](std::size_t vid) {
      std::string s;
      if (vid < rich.cells.size())
        for (auto const& x : rich.cells[vid].home_modes)
          s += sequant::toUtf8(x.full_label()) + " ";
      return s;
    };
    std::function<void(sequant::eval::ScopeBlock const&, int)> dump =
        [&](sequant::eval::ScopeBlock const& b, int d) {
          std::string ind(2 * d, ' ');
          if (d > 0)
            std::cerr << ind
                      << "BLOCK axis=" << sequant::toUtf8(b.axis.full_label())
                      << " kind="
                      << (b.kind == sequant::BatchModeType::External ? "EXT"
                                                                     : "CON")
                      << " depth=" << b.level.depth
                      << " slot=" << b.level.loop_slot
                      << " lat=" << b.latitude_ordinal << "\n";
          for (auto const& s : b.steps) {
            if (auto const* bs =
                    std::get_if<sequant::eval::BuildStep>(&s.value))
              std::cerr << ind << "  build vid=" << bs->value_id
                        << " h=" << (rich.cells[bs->value_id].hash % 100000)
                        << " carried=[" << carried_str(bs->value_id)
                        << "] sliced=[" << sliced_str(bs->value_id) << "]\n";
            else
              dump(std::get<sequant::eval::ScopeBlock>(s.value), d + 1);
          }
          for (auto const& [ovid, ok] : b.outputs)
            std::cerr << ind << "  OUT vid=" << ovid
                      << " h=" << (rich.cells[ovid].hash % 100000) << " kind="
                      << (ok == sequant::eval::OutputKind::AccumulateScatter
                              ? "SCATTER"
                          : ok == sequant::eval::OutputKind::AccumulateSum
                              ? "SUM"
                              : "?")
                      << " carried=[" << carried_str(ovid) << "] sliced=["
                      << sliced_str(ovid) << "] home=[" << home_str(ovid)
                      << "]\n";
        };
    std::cerr << "=== ORDERED SCHEDULE TREE ===\n";
    dump(ordered.root, 0);
    std::cerr << "=== END TREE ===\n";

    // Children dump: SEQUANT_UT_CHILDREN=<hash%100000> prints the node's two
    // operands (hash, carried, leaf?) straight from the forest.
    if (char const* ch = std::getenv("SEQUANT_UT_CHILDREN")) {
      auto const want = std::strtoul(ch, nullptr, 10);
      for (auto const& [h, nd] : vmap_dump) {
        if ((h % 100000u) != want || nd.leaf()) continue;
        auto pr = [&](char const* tag, auto const& c) {
          std::cerr << "[children] " << tag
                    << " h=" << (c->hash_value() % 100000u)
                    << " leaf=" << (int)c.leaf() << " carried=[";
          for (auto const& x : c->canon_indices())
            std::cerr << sequant::toUtf8(x.full_label()) << " ";
          std::cerr << "] sliced=[";
          for (auto const& x : c->sliced_modes())
            std::cerr << sequant::toUtf8(x.full_label()) << " ";
          std::cerr << "]\n";
        };
        std::cerr << "[children] parent h=" << (h % 100000u) << "\n";
        pr("left ", nd.left());
        pr("right", nd.right());
      }
    }

    // Residency-consistency check (the sanity invariant): a value homed
    // inside a loop over physical mode m, that CARRIES m, must be SLICED on m.
    // A build is homed inside its enclosing block chain (incl. its own block);
    // an escape output is homed ONE LEVEL OUT (its own block excluded). A
    // carried-but-not-sliced mode means the value is placed inside a loop its
    // own residency says it is invariant to -- the vanish/scatter bug.
    std::vector<std::string> viol;
    auto label_set = [&](std::size_t vid, bool sliced) {
      std::vector<std::string> v;
      if (sliced) {
        auto const it =
            vmap_dump.find(sequant::eval::value_key_of(rich.cells[vid]));
        if (it != vmap_dump.end())
          for (auto const& x : it->second->sliced_modes())
            v.push_back(sequant::toUtf8(x.full_label()));
      } else if (vid < rich.cells.size()) {
        for (auto const& x : rich.cells[vid].carried)
          v.push_back(sequant::toUtf8(x.full_label()));
      }
      return v;
    };
    auto check_value = [&](std::size_t vid,
                           std::vector<std::string> const& encl,
                           const char* what) {
      auto const carried = label_set(vid, /*sliced=*/false);
      auto const sliced = label_set(vid, /*sliced=*/true);
      for (auto const& m : encl) {
        bool const carries =
            std::find(carried.begin(), carried.end(), m) != carried.end();
        bool const is_sliced =
            std::find(sliced.begin(), sliced.end(), m) != sliced.end();
        if (carries && !is_sliced)
          viol.push_back(std::string(what) + " vid=" + std::to_string(vid) +
                         " h=" + std::to_string(rich.cells[vid].hash % 100000) +
                         " carries " + m +
                         " (enclosing loop) but is NOT sliced on it");
      }
    };
    std::function<void(sequant::eval::ScopeBlock const&,
                       std::vector<std::string>, int)>
        chk = [&](sequant::eval::ScopeBlock const& b,
                  std::vector<std::string> encl, int d) {
          // outputs escape one level out -> checked against encl (b excluded)
          for (auto const& [ovid, ok] : b.outputs)
            check_value(ovid, encl, "OUT");
          // builds & nested blocks live inside b -> b's own axis encloses them
          // (the root scope, d==0, opens no loop)
          if (d > 0) encl.push_back(sequant::toUtf8(b.axis.full_label()));
          for (auto const& s : b.steps) {
            if (auto const* bs =
                    std::get_if<sequant::eval::BuildStep>(&s.value))
              check_value(bs->value_id, encl, "build");
            else
              chk(std::get<sequant::eval::ScopeBlock>(s.value), encl, d + 1);
          }
        };
    chk(ordered.root, {}, 0);
    std::cerr << "=== RESIDENCY-CONSISTENCY: " << viol.size()
              << " violation(s) ===\n";
    for (auto const& v : viol) std::cerr << "  " << v << "\n";
    std::cerr << "=== END RESIDENCY-CONSISTENCY ===\n";

    // Occurrence trace: the meet stamps sliced=i_1 on a canonical node iff
    // EVERY forest occurrence has i_1 in its node_modes = (i_1 opened at/above
    // it) AND (it carries i_1). Record, per Κ-block member, each occurrence's
    // (i1_open_above, carries_i1) so a cross-occurrence intersection that kills
    // i_1 (some occurrence lacks the open) is visible directly.
    {
      std::unordered_map<std::size_t, std::string> tag;
      for (auto const& [vid, t] :
           std::initializer_list<std::pair<int, const char*>>{{43, "43"},
                                                              {44, "44"},
                                                              {80, "80"},
                                                              {81, "81"},
                                                              {106, "106"},
                                                              {107, "107"},
                                                              {45, "45=64060"},
                                                              {82, "82"},
                                                              {108, "108"},
                                                              {192, "192"}})
        tag[rich.cells[vid].hash] = t;
      // per target hash: list of "(open,carry)" occurrence strings
      std::unordered_map<std::size_t, std::vector<std::string>> occ;
      using NodeT = std::remove_cvref_t<decltype(forest.front())>;
      std::function<void(NodeT const&, bool)> ow = [&](NodeT const& n,
                                                       bool i1_open_above) {
        if (n.leaf()) return;
        bool opens_i1 = false;
        for (auto const& [ix, k] : n->batch_loops_opened_here())
          if (sequant::toUtf8(ix.full_label()) == "i_1") opens_i1 = true;
        bool const i1_here = i1_open_above || opens_i1;
        if (auto it = tag.find(n->hash_value()); it != tag.end()) {
          bool carries = false;
          for (auto const& x : n->canon_indices())
            if (sequant::toUtf8(x.full_label()) == "i_1") carries = true;
          occ[n->hash_value()].push_back(
              std::string("(open_above=") + (i1_open_above ? "Y" : "n") +
              " opens_here=" + (opens_i1 ? "Y" : "n") +
              " carries=" + (carries ? "Y" : "n") + " => node_modes has i_1: " +
              ((i1_here && carries) ? "YES" : "no") + ")");
        }
        ow(n.left(), i1_here);
        ow(n.right(), i1_here);
      };
      for (auto const& t : forest) ow(t, false);
      std::cerr << "=== OCCURRENCE TRACE (i_1) ===\n";
      for (auto const& [vid, t] :
           std::initializer_list<std::pair<int, const char*>>{{43, "43"},
                                                              {44, "44"},
                                                              {80, "80"},
                                                              {81, "81"},
                                                              {106, "106"},
                                                              {107, "107"},
                                                              {45, "45=64060"},
                                                              {82, "82"},
                                                              {108, "108"},
                                                              {192, "192"}}) {
        auto const& list = occ[rich.cells[vid].hash];
        std::cerr << "  " << t << " (h=" << (rich.cells[vid].hash % 100000)
                  << "): " << list.size() << " occurrence(s)\n";
        for (auto const& s : list) std::cerr << "      " << s << "\n";
      }
      std::cerr << "=== END OCCURRENCE TRACE ===\n";
    }

    // per_axis role + home_floor for the violators: confirm residency home.
    {
      auto role_name = [](sequant::eval::LoopRole r) -> const char* {
        switch (r) {
          case sequant::eval::LoopRole::LoopLocal:
            return "Local";
          case sequant::eval::LoopRole::Reduction:
            return "Reduction";
          case sequant::eval::LoopRole::LoopCarried:
            return "Carried";
          case sequant::eval::LoopRole::LoopInvariant:
            return "Invariant";
        }
        return "?";
      };
      std::unordered_map<std::size_t, std::size_t> cell_of;
      for (std::size_t i = 0; i < legality.cells.size(); ++i)
        cell_of[legality.cells[i].hash] = i;
      std::cerr << "=== PER_AXIS / HOME_FLOOR (violators + peers) ===\n";
      for (int vid : {43, 44, 80, 81, 106, 107, 45, 108, 82, 192}) {
        auto it = cell_of.find(rich.cells[vid].hash);
        if (it == cell_of.end()) continue;
        auto const& cl = legality.cells[it->second];
        std::cerr << "  vid=" << vid << " h=" << (rich.cells[vid].hash % 100000)
                  << " per_axis={";
        for (auto const& ac : cl.per_axis)
          std::cerr << sequant::toUtf8(ac.axis.full_label()) << ":"
                    << role_name(ac.role) << " ";
        std::cerr << "} home_floor={";
        for (auto const& x : cl.home_floor)
          std::cerr << sequant::toUtf8(x.full_label()) << " ";
        std::cerr << "}\n";
      }
      std::cerr << "=== END PER_AXIS / HOME_FLOOR ===\n";
    }
  }

  using annot_t = std::remove_cvref_t<decltype(forest.front()->annot())>;
  annot_t const layout{};
  sequant::eval::dryrun::DryRunLeafEvaluator const yield{cm};
  std::function<std::size_t(sequant::Index const&)> const target =
      [](sequant::Index const& ix) -> std::size_t {
    return ix.space().base_key() == L"Κ" ? 256 : 16;
  };
  std::function<bool(Node const&)> const is_volatile_node =
      [p = policy.is_volatile_leaf](Node const& n) -> bool {
    if (!n.leaf() || !n->is_tensor()) return false;
    return p && p(n->as_tensor());
  };

  auto& logger = sequant::Logger::instance();
  auto const prev_level = logger.eval.level;
  logger.eval.level = 0;
  auto aops = sequant::eval::dryrun::make_dryrun_array_ops(cm);
  auto ordered_cache = sequant::cache_manager(forest);
  ordered_cache.set_array_ops(&aops);

  // THE reproducer: the dry-run walk of the aux+occ schedule must run to
  // completion. It exercises the full ordered-executor runtime on the real w20
  // schedule AS-IS (no schedule rewrite): consumer-aware residency homing,
  // non-decrementing reuse, per-occurrence use-induced slicing with explicit
  // invariant facts, per-level scatter axis selection, and inner-scatter
  // destination sizing. Any regression in those surfaces here as a throw.
  // STRICT dry-run test-drive (2026-09-02): (1) the dry-run backend now
  // carries slice LOBOUNDS and throws on any shared-label range mismatch in
  // prod/sum/add_inplace/scatter (the analogue of TA's index-map and
  // is_range_set_congruent asserts, which Release elides); (2) cache-fill-once
  // is a hard error (a value cell stored twice without a reset = a consumer
  // missed the resident cell and rebuilt it, e.g. a frame-sensitive key);
  // (3) strict read-from-home is on in the ordered executor. Together these
  // make the walk catch, in seconds, the classes of defect that took
  // RelWithDebInfo wet runs to find on w20.
  {
    auto const sma =
        sequant::eval::compute_sliced_mode_assignment(ordered, rich);
    auto const vmap = sequant::eval::build_value_node_map(forest);
    sequant::eval::CellTableInputs in;
    in.ordered = &ordered;
    in.rich = &rich;
    in.sliced = &sma;
    in.sliced_modes_of = [&](std::size_t vid) {
      auto const it = vmap.find(sequant::eval::value_key_of(rich.cells[vid]));
      REQUIRE(it != vmap.end());
      return sequant::eval::detail::home_modes_in_cell_frame(rich, vid,
                                                             it->second);
    };
    in.volatile_of = [&](std::size_t vid) {
      auto const it = vmap.find(sequant::eval::value_key_of(rich.cells[vid]));
      return it != vmap.end() &&
             sequant::subtree_any(it->second, is_volatile_node);
    };
    // SP4 Task 4 fix1 item 4: the static gate must validate the same lives
    // the executor enforces at runtime -- real per-loop-instance batch
    // counts, not the constant-1 stub (which under-counts a value read from
    // outside n*m nested real batches, exactly the gap that made
    // CellRegistry::read throw "read past its life" before ordered_executor
    // wired the real function into run_ordered_schedule_pre_results).
    in.n_batches_of = sequant::eval::detail::ordered_n_batches_by_loop(
        ordered, target, &aops);
    in.operands_of = orderedexec_per_leg_operands(rich, vmap);
    auto const table = sequant::eval::build_cell_table(in);
    auto const violations = sequant::eval::validate_cell_table(
        table, ordered.root, in.n_batches_of);
    // The report goes to stderr as well as to Catch2. MEASURED (Catch2 v3.9.1,
    // this build, console/compact reporters and --verbosity high): an
    // UNSCOPED_INFO emitted immediately before a failing REQUIRE in this
    // fixture is NOT printed with the failure, so the Catch2 messages alone
    // would leave the gate failing without naming a single cell. Both loops
    // are silent when the table is clean.
    for (auto const& v : violations) {
      UNSCOPED_INFO("[" << v.rule << "] " << v.what);
      std::cerr << "[" << v.rule << "] " << v.what << "\n";
    }
    // A sliced mode that matched no enclosing loop instance is recorded whole
    // -- a form the schedule may not produce; that list must be empty, exactly
    // like a violation.
    for (auto const& [cid, pos] : table.unresolved) {
      UNSCOPED_INFO("[unresolved] cell#" << cid << " position " << pos
                                         << " (value "
                                         << table.cells[cid].value_id << ")");
      std::cerr << "[unresolved] cell#" << cid << " position " << pos
                << " (value " << table.cells[cid].value_id << ")\n";
    }
    // SEQUANT_UT_READS_OF=<vid>: every table read whose consumer is a cell
    // of that value (source cell, source value, its scope depth and slices).
    if (char const* rv = std::getenv("SEQUANT_UT_READS_OF")) {
      std::size_t const want = std::strtoul(rv, nullptr, 10);
      for (auto const& r : table.reads) {
        if (table.cells[r.consumer].value_id != want) continue;
        auto const& sc = table.cells[r.source];
        std::cerr << "[reads-of] consumer cell#" << r.consumer << " (value "
                  << want << ") reads cell#" << r.source << " (value "
                  << sc.value_id << ", scope depth " << sc.scope.path.size()
                  << ", sliced=" << sc.sliced.size()
                  << ") slices=" << r.slice.size() << "\n";
      }
    }
    // The static gate, unconditional: EVERY configuration this fixture is run
    // under (the default one, or one mirrored from an input through the
    // environment overrides above) must derive a valid table before the walk
    // executes it.
    REQUIRE(violations.empty());
    REQUIRE(table.unresolved.empty());
  }

  setenv("SEQUANT_UT_STRICT_FILL_ONCE", "1", 1);
  REQUIRE_NOTHROW(sequant::eval::evaluate_ordered_schedule<sequant::Trace::Off>(
      forest, ordered, rich, layout, yield, ordered_cache, target, {},
      is_volatile_node));

  logger.eval.level = prev_level;
}

// ===========================================================================
// Explicit value cells (SP4 Task 3): derive the cell table's productions
// (Build/Assemble/Leaf) from an ordered schedule already built above by the
// [w20-auxocc-walk] fixture. This checks structure only (one Build per
// BuildStep, one Assemble per block output entry, at least one Leaf, every
// Assemble strictly enclosing its source, every Build cell's sliced entries
// naming an enclosing loop instance of its own scope) -- reads and lives are
// Task 4.
// ===========================================================================
TEST_CASE("cell table: cells derived from the w20 default schedule",
          "[cell_table][ordered]") {
  using sequant::eval::dryrun::EvalExprDryRun;
  using sequant::eval::dryrun::EvalNodeDryRun;
  using Node = EvalNodeDryRun;
  // Same construction as the [w20-auxocc-walk] case up to the schedule.
  auto ctx = sequant::get_default_context().clone();
  ctx.set_first_dummy_index_ordinal(1000000);
  auto isr = ctx.mutable_index_space_registry();
  REQUIRE(isr != nullptr);
  sequant::mbpt::add_pao_spaces(isr, sequant::mbpt::Spin::any);
  sequant::mbpt::add_df_spaces(isr);
  auto ctx_resetter = sequant::set_scoped_default_context(std::move(ctx));

  auto const body =
      orderedexec_witness_slurp(std::string(SEQUANT_UNIT_TESTS_SOURCE_DIR) +
                                "/data/csv_ccsd_doubles_residual_df.txt");
  REQUIRE(!body.empty());
  std::string line = body;
  if (auto nl = line.find('\n'); nl != std::string::npos)
    line = line.substr(0, nl);
  auto expr = sequant::deserialize<sequant::ExprPtr>(line);
  REQUIRE(expr->is<sequant::Sum>());
  auto const& summands = expr->as<sequant::Sum>().summands();
  auto regime = orderedexec_witness_df_regime(kOrderedExecWater20_pVDZF12);
  auto cm = std::make_shared<sequant::eval::dryrun::CostModel const>(regime);
  sequant::BatchPolicy policy;
  policy.is_batchable_contracted_index = [](sequant::Index const& ix) {
    return ix.space().base_key() == L"Κ";
  };
  policy.is_batchable_external_index = [](sequant::Index const& ix) {
    auto const reg = sequant::get_default_context().index_space_registry();
    return reg && ix.space() && reg->is_pure_occupied(ix.space());
  };
  policy.batch_spectator_indices = true;
  policy.node_level_placement = true;
  policy.batch_target_size = [](sequant::Index const& ix) -> std::size_t {
    return ix.space().base_key() == L"Κ" ? 256 : 16;
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
  std::vector<Node> forest;
  for (auto const& s : summands) {
    sequant::ExprPtr const term = orderedexec_witness_flatten_product(s);
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
  auto const block_of = [](sequant::Index const& ix) -> std::size_t {
    return ix.space().base_key() == L"Κ" ? 256 : 16;
  };
  auto const rich = sequant::eval::compute_dag_boulevard(forest, *cm, block_of);
  auto const legality = sequant::eval::analyze_legality(rich, forest, policy);
  auto const ordered = sequant::eval::build_ordered_schedule(
      rich, legality, policy, std::initializer_list<std::wstring>{});
  REQUIRE(sequant::eval::well_formed(ordered));
  auto const sma = sequant::eval::compute_sliced_mode_assignment(ordered, rich);
  auto const vmap = sequant::eval::build_value_node_map(forest);
  // SP4 Task 4 fix1 item 4: real per-loop batch counts for n_batches_of (see
  // the walk-gate case's own note); aops is otherwise unused in this
  // static-only fixture.
  auto aops = sequant::eval::dryrun::make_dryrun_array_ops(cm);

  sequant::eval::CellTableInputs in;
  in.ordered = &ordered;
  in.rich = &rich;
  in.sliced = &sma;
  in.sliced_modes_of = [&](std::size_t vid) {
    auto const it = vmap.find(sequant::eval::value_key_of(rich.cells[vid]));
    REQUIRE(it != vmap.end());
    return sequant::eval::detail::home_modes_in_cell_frame(rich, vid,
                                                           it->second);
  };
  in.volatile_of = [&](std::size_t vid) {
    auto const it = vmap.find(sequant::eval::value_key_of(rich.cells[vid]));
    return it != vmap.end() &&
           sequant::subtree_any(it->second, [&](auto const& n) {
             return n.leaf() && n->is_tensor() &&
                    n->as_tensor().label() == L"t";
           });
  };
  in.n_batches_of = sequant::eval::detail::ordered_n_batches_by_loop(
      ordered, policy.batch_target_size, &aops);
  in.operands_of = orderedexec_per_leg_operands(rich, vmap);
  auto const table = sequant::eval::build_cell_table(in);

  // one Build cell per BuildStep, plus one implicit Build cell per output
  // entry (b, ovid) whose value has no BuildStep of its own in b and no
  // direct child block of b that already lists ovid among its own outputs
  // (a nested escape, whose own Build/Assemble already covers it); one
  // Assemble per output entry; Leaf cells.
  std::size_t builds = 0, outputs = 0, implicit = 0;
  std::size_t scatter_outputs = 0, blocks_with_sum_output = 0;
  std::function<void(sequant::eval::ScopeBlock const&)> count =
      [&](sequant::eval::ScopeBlock const& b) {
        std::unordered_set<std::size_t> own_builds, child_outputs;
        for (auto const& st : b.steps) {
          if (auto const* bs =
                  std::get_if<sequant::eval::BuildStep>(&st.value)) {
            ++builds;
            own_builds.insert(bs->value_id);
          } else {
            auto const& child = std::get<sequant::eval::ScopeBlock>(st.value);
            for (auto const& [cvid, ckind] : child.outputs)
              child_outputs.insert(cvid);
            count(child);
          }
        }
        outputs += b.outputs.size();
        bool has_sum_output = false;
        for (auto const& [ovid, okind] : b.outputs) {
          if (!own_builds.count(ovid) && !child_outputs.count(ovid)) ++implicit;
          if (okind == sequant::eval::OutputKind::AccumulateScatter)
            ++scatter_outputs;
          if (okind == sequant::eval::OutputKind::AccumulateSum)
            has_sum_output = true;
        }
        if (has_sum_output) ++blocks_with_sum_output;
      };
  count(ordered.root);
  std::size_t nb = 0, na = 0, nl = 0;
  for (auto const& c : table.cells) {
    if (c.production.kind == sequant::eval::ProductionKind::Build) ++nb;
    if (c.production.kind == sequant::eval::ProductionKind::Assemble) ++na;
    if (c.production.kind == sequant::eval::ProductionKind::Leaf) ++nl;
  }
  CHECK(nb == builds + implicit);
  CHECK(na == outputs);
  CHECK(nl > 0);
  // non-vacuous: the walk actually exercises sliced Build cells, Scatter
  // assembly with a real scatter map (when the tree has any Scatter output
  // at all), and partial-sum Build cells (at least one per block with a Sum
  // output).
  std::size_t n_sliced_builds = 0, n_scatter_with_map = 0,
              n_partial_over_builds = 0;
  for (auto const& c : table.cells) {
    if (c.production.kind == sequant::eval::ProductionKind::Build) {
      if (!c.sliced.empty()) ++n_sliced_builds;
      if (!c.partial_over.empty()) ++n_partial_over_builds;
    }
    if (c.production.kind == sequant::eval::ProductionKind::Assemble &&
        c.production.assemble == sequant::eval::AssembleKind::Scatter &&
        !c.production.scatter_map.empty())
      ++n_scatter_with_map;
  }
  CHECK(n_sliced_builds > 0);
  if (scatter_outputs > 0) CHECK(n_scatter_with_map > 0);
  CHECK(n_partial_over_builds >= blocks_with_sum_output);
  CHECK(n_partial_over_builds > 0);
  // every Assemble encloses its source strictly, assembles a form of its
  // OWN value, and (for a Sum) the source records the closing instance as
  // partial-summed over; every Build cell's sliced entries name enclosing
  // instances of its own scope.
  for (auto const& c : table.cells) {
    if (c.production.kind == sequant::eval::ProductionKind::Assemble) {
      auto const& s = table.cells[c.production.source];
      CHECK(c.scope.encloses(s.scope));
      // STRICT enclosure, not "exactly one level deeper": an escape chain may
      // SKIP a level the value is invariant to, so an Assemble's source can
      // sit more than one level inside it.
      CHECK(s.scope.path.size() > c.scope.path.size());
      CHECK(s.value_id == c.value_id);
      if (c.production.assemble == sequant::eval::AssembleKind::Sum) {
        auto const& closing = s.scope.path.back().first;
        bool found = false;
        for (auto const& k : s.partial_over)
          if (k.depth == closing.depth && k.loop_slot == closing.loop_slot)
            found = true;
        CHECK(found);
      }
    }
    if (c.production.kind == sequant::eval::ProductionKind::Build)
      for (auto const& [pos, k] : c.sliced) {
        bool enclosing = false;
        for (auto const& [pk, lat] : c.scope.path)
          if (pk.depth == k.depth && pk.loop_slot == k.loop_slot)
            enclosing = true;
        CHECK(enclosing);
      }
  }
  // every sliced mode resolved to an enclosing loop instance
  for (auto const& [cid, pos] : table.unresolved)
    UNSCOPED_INFO("[unresolved] cell#" << cid << " position " << pos
                                       << " (value "
                                       << table.cells[cid].value_id << ")");
  CHECK(table.unresolved.empty());

  // reads: exactly one per LEG of a consumer Build cell's production tree --
  // counted from the SAME per-leg source the builder consumed, so a consumer
  // whose two legs read one value contributes two. The de-duplicated
  // dependency graph is the lower bound: one read per distinct (consumer,
  // operand) pair.
  std::size_t edges = 0, distinct_pairs = 0, double_leg_consumers = 0;
  for (auto const& c : table.cells)
    if (c.production.kind == sequant::eval::ProductionKind::Build) {
      auto const ops = in.operands_of(c.value_id);
      edges += ops.size();
      std::unordered_set<std::size_t> distinct(ops.begin(), ops.end());
      distinct_pairs += distinct.size();
      if (distinct.size() < ops.size()) ++double_leg_consumers;
    }
  CHECK(table.reads.size() == edges);
  CHECK(table.reads.size() >= distinct_pairs);
  UNSCOPED_INFO("reads " << table.reads.size()
                         << ", distinct (consumer, "
                            "operand) pairs "
                         << distinct_pairs
                         << ", consumers reading one value "
                            "on both legs "
                         << double_leg_consumers);
  CHECK(double_leg_consumers > 0);
  for (auto const& r : table.reads) {
    REQUIRE(r.source < table.cells.size());
    CHECK(table.cells[r.source].value_id == r.operand_value_id);
    // a declared slice never binds an instance the source is already sliced by
    for (auto const& [p, k] : r.slice)
      for (auto const& [sp, sk] : table.cells[r.source].sliced)
        CHECK_FALSE(
            (sp == p && sk.depth == k.depth && sk.loop_slot == k.loop_slot));
  }
  // residency flags: both kinds occur, and each means what it says.
  //
  // Persistence is the FRONTIER of the invariant region (TableCell::
  // persistent, detail::apply_persistence_frontier): the CANDIDATES are the
  // cells carrying no volatile leaf and bound to no enclosing loop, and a
  // candidate stays persistent only if some CONSUMER of it is volatile, or it
  // has no consumer at all (a forest root). Consumers are counted by SOURCE
  // CELL -- the reads sourcing the cell, plus the Assembles sourcing it --
  // exactly the edge set the runtime's cache-halt closure walks.
  std::size_t const n_cells = table.cells.size();
  std::vector<char> has_consumer(n_cells, 0), has_volatile_consumer(n_cells, 0);
  {
    auto const note = [&](std::size_t source, std::size_t consumer) {
      has_consumer[source] = 1;
      if (in.volatile_of(table.cells[consumer].value_id))
        has_volatile_consumer[source] = 1;
    };
    for (auto const& r : table.reads) note(r.source, r.consumer);
    for (std::size_t c = 0; c < n_cells; ++c)
      if (table.cells[c].production.kind ==
          sequant::eval::ProductionKind::Assemble)
        note(table.cells[c].production.source, c);
  }
  std::size_t n_persistent_candidates = 0;
  for (auto const& c : table.cells)
    if (!in.volatile_of(c.value_id) &&
        sequant::eval::detail::bound_instances(c).empty())
      ++n_persistent_candidates;

  std::size_t n_persistent = 0, n_produce_if_absent = 0;
  for (std::size_t cid = 0; cid < n_cells; ++cid) {
    auto const& c = table.cells[cid];
    if (c.persistent) {
      ++n_persistent;
      // cross-evaluation invariance: bound to no loop instance, no volatile
      // leaf underneath
      CHECK(sequant::eval::detail::bound_instances(c).empty());
      CHECK_FALSE(in.volatile_of(c.value_id));
      // ... and a reader on a LATER evaluation: a volatile consumer, or no
      // consumer at all. Without this, a whole invariant sub-DAG whose every
      // consumer is itself skipped on the warm evaluation would be held in
      // the persistent value store for the life of the cache handle.
      CHECK((has_volatile_consumer[cid] || !has_consumer[cid]));
    }
    if (c.produce_if_absent) {
      ++n_produce_if_absent;
      // first-visit production only makes sense inside a loop
      CHECK_FALSE(c.scope.path.empty());
    }
    // a cell bound to the innermost loop of its own scope is rebuilt every
    // batch of that loop, never reused across visits
    if (c.production.kind == sequant::eval::ProductionKind::Build &&
        !c.scope.path.empty()) {
      bool bound_innermost = false;
      for (auto const& k : sequant::eval::detail::bound_instances(c))
        if (k.depth == c.scope.path.back().first.depth &&
            k.loop_slot == c.scope.path.back().first.loop_slot)
          bound_innermost = true;
      if (bound_innermost) CHECK_FALSE(c.produce_if_absent);
    }
  }
  CHECK(n_persistent > 0);
  CHECK(n_produce_if_absent > 0);
  // The frontier is a strict subset of the candidates on this fixture: the
  // invariant region really does run deeper than one level here, so the
  // demotion has something to do (and the numbers are reported so a change in
  // either is visible).
  WARN("persistent cells " << n_persistent << " of " << n_persistent_candidates
                           << " candidates (no volatile leaf, bound to no "
                              "enclosing loop), out of "
                           << n_cells << " cells");
  CHECK(n_persistent <= n_persistent_candidates);
  CHECK(n_persistent < n_persistent_candidates);

  // Validate the derived table once; reuse the same violations both for the
  // life-rule count and the final full-table check.
  auto const violations =
      sequant::eval::validate_cell_table(table, ordered.root, in.n_batches_of);
  std::size_t life_violations = 0;
  for (auto const& v : violations)
    if (v.rule == "life") ++life_violations;
  CHECK(life_violations == 0);
  for (auto const& v : violations)
    UNSCOPED_INFO("[" << v.rule << "] " << v.what);
  REQUIRE(violations.empty());
}

// ===========================================================================
// The SAME derivation on the configuration mirrored from a real input, set
// DIRECTLY here (no environment variable, so this runs in every suite): the
// only difference from the default configuration above is the finite peak
// budget, which is what makes the optimizer produce the two schedule shapes
// this stage had to learn to describe --
//   (1) an escape chain that SKIPS a level the value is invariant to (the
//       Assemble's source sits more than one level deeper), and
//   (2) a member MATERIALIZED across the forced loop split: built in its home
//       block for its in-nest readers AND escaped out of it for the other
//       pass, so the value has a Build cell and an Assemble cell at a
//       strictly shallower scope.
// The objective is the same enumerator the [w20-auxocc-walk] fixture's
// SEQUANT_UT_OBJECTIVE=dense_time_space_batched override selects
// (DenseTimeSpaceBatched -- also its default). The executor is NOT run here;
// the [w20-auxocc-walk] case does that.
// ===========================================================================
TEST_CASE("cell table: the input-mirrored configuration derives a valid table",
          "[cell_table][ordered]") {
  using sequant::eval::dryrun::EvalExprDryRun;
  using sequant::eval::dryrun::EvalNodeDryRun;
  using Node = EvalNodeDryRun;
  auto ctx = sequant::get_default_context().clone();
  ctx.set_first_dummy_index_ordinal(1000000);
  auto isr = ctx.mutable_index_space_registry();
  REQUIRE(isr != nullptr);
  sequant::mbpt::add_pao_spaces(isr, sequant::mbpt::Spin::any);
  sequant::mbpt::add_df_spaces(isr);
  auto ctx_resetter = sequant::set_scoped_default_context(std::move(ctx));

  auto const body =
      orderedexec_witness_slurp(std::string(SEQUANT_UNIT_TESTS_SOURCE_DIR) +
                                "/data/csv_ccsd_doubles_residual_df.txt");
  REQUIRE(!body.empty());
  std::string line = body;
  if (auto nl = line.find('\n'); nl != std::string::npos)
    line = line.substr(0, nl);
  auto expr = sequant::deserialize<sequant::ExprPtr>(line);
  REQUIRE(expr->is<sequant::Sum>());
  auto const& summands = expr->as<sequant::Sum>().summands();
  auto regime = orderedexec_witness_df_regime(kOrderedExecWater20_pVDZF12);
  auto cm = std::make_shared<sequant::eval::dryrun::CostModel const>(regime);
  sequant::BatchPolicy policy;
  policy.is_batchable_contracted_index = [](sequant::Index const& ix) {
    return ix.space().base_key() == L"Κ";
  };
  policy.is_batchable_external_index = [](sequant::Index const& ix) {
    auto const reg = sequant::get_default_context().index_space_registry();
    return reg && ix.space() && reg->is_pure_occupied(ix.space());
  };
  policy.batch_spectator_indices = true;
  policy.node_level_placement = true;
  policy.batch_target_size = [](sequant::Index const& ix) -> std::size_t {
    return ix.space().base_key() == L"Κ" ? 256 : 16;
  };
  policy.is_volatile_leaf = [](sequant::Tensor const& t) {
    return t.label() == L"t";
  };
  policy.accumulation_factor = 1.0;
  policy.persistent_only = false;
  // THE mirrored setting (the walk fixture reads it from
  // SEQUANT_UT_PEAK_THRESHOLD; here it is set directly).
  policy.peak_threshold = 25e9;
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
  std::vector<Node> forest;
  for (auto const& s : summands) {
    sequant::ExprPtr const term = orderedexec_witness_flatten_product(s);
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
  auto const block_of = [](sequant::Index const& ix) -> std::size_t {
    return ix.space().base_key() == L"Κ" ? 256 : 16;
  };
  auto const rich = sequant::eval::compute_dag_boulevard(forest, *cm, block_of);
  auto const legality = sequant::eval::analyze_legality(rich, forest, policy);
  auto const ordered = sequant::eval::build_ordered_schedule(
      rich, legality, policy, std::initializer_list<std::wstring>{});
  REQUIRE(sequant::eval::well_formed(ordered));
  auto const sma = sequant::eval::compute_sliced_mode_assignment(ordered, rich);
  auto const vmap = sequant::eval::build_value_node_map(forest);
  // SP4 Task 4 fix1 item 4: real per-loop batch counts for n_batches_of (see
  // the walk-gate case's own note); aops is otherwise unused in this
  // static-only fixture.
  auto aops = sequant::eval::dryrun::make_dryrun_array_ops(cm);

  sequant::eval::CellTableInputs in;
  in.ordered = &ordered;
  in.rich = &rich;
  in.sliced = &sma;
  in.sliced_modes_of = [&](std::size_t vid) {
    auto const it = vmap.find(sequant::eval::value_key_of(rich.cells[vid]));
    REQUIRE(it != vmap.end());
    return sequant::eval::detail::home_modes_in_cell_frame(rich, vid,
                                                           it->second);
  };
  in.volatile_of = [&](std::size_t vid) {
    auto const it = vmap.find(sequant::eval::value_key_of(rich.cells[vid]));
    return it != vmap.end() &&
           sequant::subtree_any(it->second, [&](auto const& n) {
             return n.leaf() && n->is_tensor() &&
                    n->as_tensor().label() == L"t";
           });
  };
  in.n_batches_of = sequant::eval::detail::ordered_n_batches_by_loop(
      ordered, policy.batch_target_size, &aops);
  in.operands_of = orderedexec_per_leg_operands(rich, vmap);
  auto const table = sequant::eval::build_cell_table(in);

  auto const violations2 =
      sequant::eval::validate_cell_table(table, ordered.root, in.n_batches_of);
  for (auto const& v : violations2)
    UNSCOPED_INFO("[" << v.rule << "] " << v.what);
  for (auto const& [cid, pos] : table.unresolved)
    UNSCOPED_INFO("[unresolved] cell#" << cid << " position " << pos
                                       << " (value "
                                       << table.cells[cid].value_id << ")");
  REQUIRE(violations2.empty());
  REQUIRE(table.unresolved.empty());

  // (1) at least one escape chain SKIPS a level: the Assemble's source sits
  // two or more levels deeper than the Assemble itself.
  std::size_t level_skipping_chains = 0;
  for (auto const& c : table.cells)
    if (c.production.kind == sequant::eval::ProductionKind::Assemble &&
        table.cells[c.production.source].scope.path.size() >=
            c.scope.path.size() + 2)
      ++level_skipping_chains;
  CHECK(level_skipping_chains > 0);

  // (2) at least one value is MATERIALIZED across the forced split. In the
  // TABLE that shows up as a Build cell whose value also has an Assemble cell
  // at a strictly shallower scope -- necessary but not sufficient, since the
  // implicit per-batch Build the builder synthesizes for an ordinary escape
  // matches it too.
  std::size_t build_with_shallower_assemble = 0;
  for (auto const& b : table.cells) {
    if (b.production.kind != sequant::eval::ProductionKind::Build) continue;
    for (auto const& a : table.cells)
      if (a.production.kind == sequant::eval::ProductionKind::Assemble &&
          a.value_id == b.value_id &&
          a.scope.path.size() < b.scope.path.size()) {
        ++build_with_shallower_assemble;
        break;
      }
  }
  CHECK(build_with_shallower_assemble > 0);

  // The unambiguous witness is in the SCHEDULE: one block that BOTH holds a
  // BuildStep for a value AND lists that same value among its outputs. Only
  // the mixed-pass materialization emits that shape.
  std::size_t built_and_escaped_here = 0;
  std::function<void(sequant::eval::ScopeBlock const&)> scan =
      [&](sequant::eval::ScopeBlock const& b) {
        std::unordered_set<std::size_t> own;
        for (auto const& st : b.steps) {
          if (auto const* bs = std::get_if<sequant::eval::BuildStep>(&st.value))
            own.insert(bs->value_id);
          else
            scan(std::get<sequant::eval::ScopeBlock>(st.value));
        }
        for (auto const& [ovid, kind] : b.outputs) {
          (void)kind;
          if (own.count(ovid)) ++built_and_escaped_here;
        }
      };
  scan(ordered.root);
  // Under value identity (explicit-cells design section 11) the water-20
  // default schedule no longer materializes a member across a split: the
  // later-pass reader of a nest-homed value was a merged-frame occurrence,
  // now its own value in its own nest. The mixed-pass shape stays pinned by
  // the [per-nest-split] fixtures; here it is reported, not required.
  WARN("built_and_escaped_here = " << built_and_escaped_here);
}

// ===========================================================================
// Diagnostic-only dump (Task 1 of the explicit-cells stage-2 plan):
// characterize, on the mirrored w20 configuration, the values that the
// stage-1 table builder and schedule builder mis-handled -- where each is
// built or escaped, and who consumes it and where. Those gaps are FIXED: the
// mirrored configuration now derives a table that validates clean (see the
// unconditional gate block in the [w20-auxocc-walk] case, and the
// "cell table: the input-mirrored configuration derives a valid table" case
// which asserts it without any environment variable). This case is kept as
// the shape characterization behind the fix. Hidden ("[.]") so it never runs
// in the default or CI suites; run explicitly with the mirrored environment
// (SEQUANT_UT_PEAK_THRESHOLD=25e9,
// SEQUANT_UT_OBJECTIVE=dense_time_space_batched) to see the report on
// stderr. The watched value ids are only meaningful under that mirrored
// configuration.
// ===========================================================================
TEST_CASE(
    "cell table: dump the mixed-pass members of the mirrored w20 schedule",
    "[cell_table][ordered][.][mixed-pass-dump]") {
  using sequant::eval::dryrun::EvalExprDryRun;
  using sequant::eval::dryrun::EvalNodeDryRun;
  using Node = EvalNodeDryRun;

  auto ctx = sequant::get_default_context().clone();
  ctx.set_first_dummy_index_ordinal(1000000);
  auto isr = ctx.mutable_index_space_registry();
  REQUIRE(isr != nullptr);
  sequant::mbpt::add_pao_spaces(isr, sequant::mbpt::Spin::any);
  sequant::mbpt::add_df_spaces(isr);
  auto ctx_resetter = sequant::set_scoped_default_context(std::move(ctx));

  auto const body =
      orderedexec_witness_slurp(std::string(SEQUANT_UNIT_TESTS_SOURCE_DIR) +
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

  // FULL residual (all summands) to match the MPQC w20 run; overridable for
  // bisecting which term first breaks the walk.
  std::size_t nterms = summands.size();
  if (char const* nt = std::getenv("SEQUANT_UT_DRYRUN_NTERMS"))
    nterms = std::min<std::size_t>(summands.size(), std::atoll(nt));

  auto regime = orderedexec_witness_df_regime(kOrderedExecWater20_pVDZF12);
  auto cm = std::make_shared<sequant::eval::dryrun::CostModel const>(regime);

  // AUX+OCC: Κ batchable-contracted (aux), occ batchable-EXTERNAL; spectator
  // batching + node-level placement ON (make_csv_batch_policy, occ_target>0).
  sequant::BatchPolicy policy;
  policy.is_batchable_contracted_index = [](sequant::Index const& ix) {
    return ix.space().base_key() == L"Κ";
  };
  policy.is_batchable_external_index = [](sequant::Index const& ix) {
    auto const reg = sequant::get_default_context().index_space_registry();
    return reg && ix.space() && reg->is_pure_occupied(ix.space());
  };
  policy.batch_spectator_indices = true;
  policy.node_level_placement = true;
  policy.batch_target_size = [](sequant::Index const& ix) -> std::size_t {
    return ix.space().base_key() == L"Κ" ? 256 : 16;
  };
  policy.is_volatile_leaf = [](sequant::Tensor const& t) {
    return t.label() == L"t";
  };
  policy.accumulation_factor = 1.0;
  policy.persistent_only = false;
  policy.peak_threshold = 1e11;
  // Mirror a specific MPQC input (z820 w20 csv-cck-diag.json: batch
  // peak_threshold 25e9, optimize objective dense_time_space) via env.
  if (char const* pt = std::getenv("SEQUANT_UT_PEAK_THRESHOLD"))
    policy.peak_threshold = std::atof(pt);

  auto axes_map = std::make_shared<std::unordered_map<
      sequant::Expr const*,
      sequant::container::vector<sequant::NodeBatchAnnotation>>>();
  sequant::OptimizeOptions opts;
  opts.objective_function = sequant::ObjectiveFunction::DenseTimeSpaceBatched;
  if (char const* ob = std::getenv("SEQUANT_UT_OBJECTIVE")) {
    std::string const o{ob};
    if (o == "dense_time_space")
      opts.objective_function = sequant::ObjectiveFunction::DenseTimeSpace;
    else if (o == "dense_time_space_batched")
      opts.objective_function =
          sequant::ObjectiveFunction::DenseTimeSpaceBatched;
    else if (o == "dense_space_time")
      opts.objective_function = sequant::ObjectiveFunction::DenseSpaceTime;
    else if (o == "dense_space_time_batched")
      opts.objective_function =
          sequant::ObjectiveFunction::DenseSpaceTimeBatched;
  }
  opts.idx_to_extent = regime.idx_to_extent();
  opts.inner_pow = regime.inner_pow_fn();
  opts.batch_policy = policy;
  opts.volatile_weight = 20.0;
  opts.roofline.machine_balance = 200.0;
  opts.roofline.fast_mem_elems = 1000000.0;
  opts.term_batch_axes = axes_map;

  std::vector<Node> forest;
  for (std::size_t s = 0; s < nterms; ++s) {
    sequant::ExprPtr const term =
        orderedexec_witness_flatten_product(summands[s]);
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

  auto const block_of = [](sequant::Index const& ix) -> std::size_t {
    return ix.space().base_key() == L"Κ" ? 256 : 16;
  };
  auto const rich = sequant::eval::compute_dag_boulevard(forest, *cm, block_of);
  REQUIRE(!rich.cells.empty());

  auto const legality = sequant::eval::analyze_legality(rich, forest, policy);
  REQUIRE(legality.cells.size() == rich.cells.size());

  // MPQC passes an EMPTY mode_order (build_ordered_schedule derives the forced
  // split axes from the legality) -- match that.
  auto const ordered = sequant::eval::build_ordered_schedule(
      rich, legality, policy, std::initializer_list<std::wstring>{});
  REQUIRE(sequant::eval::well_formed(ordered));

  auto const g = sequant::eval::detail::ordered_schedule_dep_graph(rich);
  std::set<std::size_t> const watch = {45, 48, 197, 198};
  // (1) where each watched value is built or escaped
  std::function<void(sequant::eval::ScopeBlock const&, std::string)> walk =
      [&](sequant::eval::ScopeBlock const& b, std::string path) {
        path += "[d" + std::to_string(b.level.depth) + " s" +
                std::to_string(b.level.loop_slot) + " lat" +
                std::to_string(b.latitude_ordinal) + "]";
        for (auto const& st : b.steps) {
          if (auto const* bs =
                  std::get_if<sequant::eval::BuildStep>(&st.value)) {
            if (watch.count(bs->value_id))
              std::cerr << "BUILD v" << bs->value_id
                        << " hash=" << rich.cells[bs->value_id].hash % 100000u
                        << " at " << path << "\n";
          } else {
            walk(std::get<sequant::eval::ScopeBlock>(st.value), path);
          }
        }
        for (auto const& [ovid, kind] : b.outputs)
          if (watch.count(ovid))
            std::cerr << "OUTPUT v" << ovid
                      << " kind=" << static_cast<int>(kind) << " of " << path
                      << "\n";
      };
  walk(ordered.root, "");
  // (2) consumers of each watched value and where they are built
  std::unordered_map<std::size_t, std::string> built_at;
  std::function<void(sequant::eval::ScopeBlock const&, std::string)> index =
      [&](sequant::eval::ScopeBlock const& b, std::string path) {
        path += "[d" + std::to_string(b.level.depth) + " s" +
                std::to_string(b.level.loop_slot) + " lat" +
                std::to_string(b.latitude_ordinal) + "]";
        for (auto const& st : b.steps) {
          if (auto const* bs = std::get_if<sequant::eval::BuildStep>(&st.value))
            built_at[bs->value_id] = path;
          else
            index(std::get<sequant::eval::ScopeBlock>(st.value), path);
        }
        for (auto const& [ovid, kind] : b.outputs)
          built_at.try_emplace(ovid, path + "(output)");
      };
  index(ordered.root, "");
  for (std::size_t v : watch) {
    std::cerr << "v" << v << " hash=" << rich.cells[v].hash % 100000u
              << " leaf=" << rich.cells[v].is_leaf << " consumers:";
    if (auto it = g.consumers_of.find(v); it != g.consumers_of.end())
      for (std::size_t c : it->second)
        std::cerr << " v" << c << "@" << built_at[c];
    std::cerr << "\n  operands:";
    if (auto it = g.depends_on.find(v); it != g.depends_on.end())
      for (std::size_t o : it->second)
        std::cerr << " v" << o << "@" << built_at[o];
    std::cerr << "\n";
  }
  SUCCEED("characterization printed to stderr");
}

// ===========================================================================
// Cache-halt across CC iterations: a Κ-free PERSISTENT composite (I(i,i;a,a),
// e.g. the 4-PNO-2-occ integral) built by contracting Κ between Κ-carrying
// prerequisites must be built ONCE (iteration 1) and reused thereafter, AND
// its Κ-carrying prerequisites (loop-local Transients of the {Κ} batch block
// that forms it) must NOT be re-formed on later iterations -- their only
// consumer, the composite, is already resident. This mirrors forest descent's
// "descent halts at a cache hit". The ordered executor's top-level needed_build
// gate already halts at resident nodes (excluding the prerequisites), but it
// was not threaded into run_ordered_contracted_block, so the batch block re-ran
// its Transients every iteration. Reuses the SAME water-20 fixture as the
// witness above.
// ===========================================================================
TEST_CASE(
    "ordered executor: cache-halt gate does not re-form a resident persistent "
    "composite, and skips its dead batch-block prerequisites when present",
    "[ordered][cache-halt]") {
  using sequant::eval::dryrun::EvalExprDryRun;
  using sequant::eval::dryrun::EvalNodeDryRun;
  using Node = EvalNodeDryRun;

  auto ctx = sequant::get_default_context().clone();
  ctx.set_first_dummy_index_ordinal(1000000);
  auto isr = ctx.mutable_index_space_registry();
  REQUIRE(isr != nullptr);
  sequant::mbpt::add_pao_spaces(isr, sequant::mbpt::Spin::any);
  sequant::mbpt::add_df_spaces(isr);
  auto ctx_resetter = sequant::set_scoped_default_context(std::move(ctx));

  auto const body =
      orderedexec_witness_slurp(std::string(SEQUANT_UNIT_TESTS_SOURCE_DIR) +
                                "/data/csv_ccsd_doubles_residual_df.txt");
  REQUIRE(!body.empty());
  std::string firstline = body;
  if (auto nl = firstline.find('\n'); nl != std::string::npos)
    firstline = firstline.substr(0, nl);
  auto expr = sequant::deserialize<sequant::ExprPtr>(firstline);
  REQUIRE(static_cast<bool>(expr));
  REQUIRE(expr->is<sequant::Sum>());
  auto const& summands = expr->as<sequant::Sum>().summands();
  REQUIRE(!summands.empty());
  std::size_t nterms = std::min<std::size_t>(summands.size(), 40);
  if (char const* nt = std::getenv("SEQUANT_UT_DRYRUN_NTERMS"))
    nterms = std::min<std::size_t>(summands.size(), std::atoll(nt));

  auto regime = orderedexec_witness_df_regime(kOrderedExecWater20_pVDZF12);
  auto cm = std::make_shared<sequant::eval::dryrun::CostModel const>(regime);

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

  std::vector<Node> forest;
  for (std::size_t s = 0; s < nterms; ++s) {
    sequant::ExprPtr const term =
        orderedexec_witness_flatten_product(summands[s]);
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

  auto const block_of = [](sequant::Index const&) -> std::size_t {
    return 256;
  };
  auto const rich = sequant::eval::compute_dag_boulevard(forest, *cm, block_of);
  REQUIRE(!rich.cells.empty());
  auto const vmap = sequant::eval::build_value_node_map(forest);

  using annot_t = std::remove_cvref_t<decltype(forest.front()->annot())>;
  annot_t const layout{};
  sequant::eval::dryrun::DryRunLeafEvaluator const yield{cm};
  std::function<std::size_t(sequant::Index const&)> const target =
      [](sequant::Index const&) -> std::size_t { return 256; };

  auto& logger = sequant::Logger::instance();
  auto const prev_level = logger.eval.level;
  auto* const prev_stream = logger.eval.stream;
  logger.eval.level = 1;  // arms tally_build
  std::ostringstream sink;
  logger.eval.stream = &sink;

  auto const legality = sequant::eval::analyze_legality(rich, forest, policy);
  auto const ordered =
      sequant::eval::build_ordered_schedule(rich, legality, policy, {L"Κ"});
  REQUIRE(sequant::eval::well_formed(ordered));

  std::function<bool(Node const&)> const is_volatile_node =
      [p = policy.is_volatile_leaf](Node const& n) -> bool {
    if (!n.leaf() || !n->is_tensor()) return false;
    return p && p(n->as_tensor());
  };

  // The volatility-aware overload stamps the persistence classification the
  // executor reads when homing (MPQC's build_cache_manager path).
  auto cache = sequant::cache_manager(forest, is_volatile_node, 2);
  auto aops = sequant::eval::dryrun::make_dryrun_array_ops(cm);
  cache.set_array_ops(&aops);
  cache.set_recompute_tally_enabled(true);

  auto run_iter = [&]() {
    try {
      (void)sequant::eval::evaluate_ordered_schedule<sequant::Trace::On>(
          forest, ordered, rich, layout, yield, cache, target, {},
          is_volatile_node);
    } catch (std::exception const& e) {
      FAIL("ordered evaluate threw: " << e.what());
    }
  };

  // Iteration 1 builds everything and homes persistent (Κ-invariant)
  // composites.
  run_iter();

  // Stage 3: persistence is the TABLE's own classification (no volatile leaf,
  // bound to no loop instance) and lives in the cache handle's cross-call
  // PersistentValueStore, which the cell registry publishes to on every
  // production of such a cell and seeds itself from at the next call's entry.
  // So the composites that must not be re-formed are exactly the ones the
  // store now holds -- collect them and their iteration-1 build counts.
  auto const& store = cache.persistent_values();
  REQUIRE(store.size() > 0);  // iteration 1 actually filled the store
  std::vector<std::pair<Node, std::size_t>> persistent_b1;
  for (auto const& vc : rich.cells) {
    auto const it = vmap.find(sequant::eval::value_key_of(vc));
    if (it == vmap.end() || it->second.leaf()) continue;
    if (!store.holds(vc.hash)) continue;
    persistent_b1.emplace_back(
        it->second, orderedexec_builds_of(cache.recompute_tally(), it->second));
  }
  REQUIRE(
      !persistent_b1.empty());  // fixture actually has persistent composites
  {
    // One specific known-persistent value, held by its own canonical hash:
    // the composite is genuinely batch-invariant (carries no volatile leaf),
    // which is what makes it eligible to survive between evaluations at all.
    Node const& known = persistent_b1.front().first;
    CHECK(store.holds(known->hash_value()));
    CHECK_FALSE(sequant::subtree_any(known, is_volatile_node));
  }

  cache.reset();
  // reset() leaves the persistent store alone by design -- that is what makes
  // the next evaluation's seeding (and so the cache-halt below) possible.
  CHECK(store.size() > 0);

  // Replicate the executor's cache-halt skip set (ordered_cache_halt_skip:
  // a cell is skipped when it is persistent and already held, or when every
  // consumer of its value is skipped) projected onto the value DAG, which is
  // what this fixture can observe: descend from the roots and halt at any
  // value the store holds -- a value the descent never reaches is one whose
  // every consumer chain was cut by a held persistent, i.e. skipped. A
  // NON-LEAF BuildStep inside a {Κ} ScopeBlock whose node is NOT in this set
  // is a DEAD prerequisite -- it must not be re-formed in iteration 2.
  sequant::container::set<std::size_t> needed;
  {
    sequant::container::svector<Node> stack;
    sequant::container::set<std::size_t> visited;
    for (auto&& n : forest) stack.push_back(n);
    while (!stack.empty()) {
      Node const n = stack.back();
      stack.pop_back();
      if (n.leaf()) continue;
      if (!visited.insert(n->hash_value()).second) continue;
      if (store.holds(n->hash_value())) continue;  // held: read, do not descend
      needed.insert(n->hash_value());
      stack.push_back(n.left());
      stack.push_back(n.right());
    }
  }
  std::optional<Node> dead_transient;  // a dead {Κ}-block Transient, if any
  std::size_t n_Kblocks = 0;
  {
    std::function<void(sequant::eval::ScopeBlock const&, bool)> scan =
        [&](sequant::eval::ScopeBlock const& b, bool inK) {
          bool const here = inK || b.axis.space().base_key() == L"Κ";
          if (b.axis.space().base_key() == L"Κ") ++n_Kblocks;
          if (here)
            for (auto const& s : b.steps)
              if (auto const* bs =
                      std::get_if<sequant::eval::BuildStep>(&s.value)) {
                auto const it2 = vmap.find(
                    sequant::eval::value_key_of(rich.cells[bs->value_id]));
                if (it2 != vmap.end() && !it2->second.leaf() &&
                    !needed.count(rich.cells[bs->value_id].hash) &&
                    !dead_transient)
                  dead_transient = it2->second;
              }
          for (auto const& s : b.steps)
            if (auto const* ch =
                    std::get_if<sequant::eval::ScopeBlock>(&s.value))
              scan(*ch, here);
        };
    scan(ordered.root, false);
  }
  std::size_t const b1_dead =
      dead_transient
          ? orderedexec_builds_of(cache.recompute_tally(), *dead_transient)
          : 0;

  // Iteration 2: resident persistents are reused; dead prerequisites skipped.
  run_iter();

  logger.eval.level = prev_level;
  logger.eval.stream = prev_stream;

  // SAFETY invariant (guards the fix): no persistent composite is re-formed in
  // iteration 2 -- its build count is frozen at the iteration-1 value.
  for (auto const& [node, b1] : persistent_b1) {
    std::size_t const b2 = orderedexec_builds_of(cache.recompute_tally(), node);
    CHECK(b2 == b1);
  }

  // Skipped consumers spend the reads they will not perform (CellRegistry::
  // forgo, called at every skip site -- see ordered_forgo_reads): so when the
  // walk is over, the ONLY things the registry is still holding are the cells
  // that are held ON PURPOSE -- the persistent ones (they survive into the
  // next evaluation) and the forest roots (just handed to the caller as
  // pre_results, read by nobody in the table). Anything else would be a
  // non-persistent intermediate whose declared life never ran out: it would
  // pin its memory for the rest of the evaluation and keep looking shared to
  // every later reader, which is what disables in-place accumulation.
  {
    auto const res = sequant::eval::detail::ordered_last_registry_residency();
    INFO("registry residency after iteration 2: live="
         << res.live << " persistent=" << res.persistent
         << " roots=" << res.roots);
    CHECK(res.persistent > 0);  // the fixture really does have persistents
    CHECK(res.live == res.persistent + res.roots);
  }

  WARN("Kblocks=" << n_Kblocks
                  << " dead_transient_found=" << dead_transient.has_value());
  if (dead_transient) {
    // IMPROVEMENT invariant: a batch-block prerequisite that feeds only a
    // now-resident persistent composite is NOT re-formed in iteration 2.
    std::size_t const b2_dead =
        orderedexec_builds_of(cache.recompute_tally(), *dead_transient);
    CHECK(b2_dead == b1_dead);
  }
}

// ===========================================================================
// [.][dryrun-2iter-report] (hidden report, not a strict gate): a 2-ITERATION
// dry-run cost report -- iter 1 (cold cache) vs iter 2 (warm cache: persistent
// composites resident + the needed_build cache-halt active) -- of builds /
// FLOPs / peak, for BOTH forest descent and the ordered/DAG executor, at w20
// residual scale (ALL terms by default). Locally predicts the cache-halt fix's
// steady-state benefit without an Owl run.
// ===========================================================================
TEST_CASE(
    "dry-run 2-iteration report: forest vs ordered, cold vs warm cache "
    "(builds/FLOPs/peak) at w20 residual scale",
    "[.][dryrun-2iter-report]") {
  using sequant::eval::dryrun::EvalExprDryRun;
  using sequant::eval::dryrun::EvalNodeDryRun;
  using Node = EvalNodeDryRun;

  auto ctx = sequant::get_default_context().clone();
  ctx.set_first_dummy_index_ordinal(1000000);
  auto isr = ctx.mutable_index_space_registry();
  REQUIRE(isr != nullptr);
  sequant::mbpt::add_pao_spaces(isr, sequant::mbpt::Spin::any);
  sequant::mbpt::add_df_spaces(isr);
  auto ctx_resetter = sequant::set_scoped_default_context(std::move(ctx));

  auto const body =
      orderedexec_witness_slurp(std::string(SEQUANT_UNIT_TESTS_SOURCE_DIR) +
                                "/data/csv_ccsd_doubles_residual_df.txt");
  REQUIRE(!body.empty());
  std::string firstline = body;
  if (auto nl = firstline.find('\n'); nl != std::string::npos)
    firstline = firstline.substr(0, nl);
  auto expr = sequant::deserialize<sequant::ExprPtr>(firstline);
  REQUIRE(static_cast<bool>(expr));
  REQUIRE(expr->is<sequant::Sum>());
  auto const& summands = expr->as<sequant::Sum>().summands();
  REQUIRE(!summands.empty());
  // DEFAULT: ALL terms (no 40 cap) so the batched persistent-composite
  // structure has the best chance to appear.
  std::size_t nterms = summands.size();
  if (char const* nt = std::getenv("SEQUANT_UT_DRYRUN_NTERMS"))
    nterms = std::min<std::size_t>(summands.size(), std::atoll(nt));

  // SYSTEM selection (SEQUANT_UT_DRYRUN_SYSTEM = water20 [default] | c60).
  std::string system = "water20";
  if (char const* s = std::getenv("SEQUANT_UT_DRYRUN_SYSTEM")) system = s;
  auto const regime = orderedexec_witness_df_regime(
      system == "c60" ? kOrderedExecC60_pVDZF12 : kOrderedExecWater20_pVDZF12);
  // The dry-run cost model carries the SAME roofline parameters the optimizer
  // ran with (opts.roofline above), so the sink's `exec` is the realized
  // roofline cost of the schedule (max(flops, beta*Q) per executed op at its
  // sliced extents), comparable with the DP's chosen_flops objective.
  auto cm = std::make_shared<sequant::eval::dryrun::CostModel const>(
      regime, sequant::RooflineParams{.machine_balance = 200.0,
                                      .fast_mem_elems = 1000000.0});

  // BATCH mode (SEQUANT_UT_DRYRUN_BATCH = none | aux [default] | aux_occ |
  // pao | aux_pao | aux_pao_occ). "pao" adds the μ̃ (PAO) contracted axis --
  // the lever that lets the factorizer contract K in the compact PAO basis and
  // slice μ̃, avoiding the aux-free 4-occ/2-PNO integrals entirely.
  std::string batch = "aux";
  if (char const* b = std::getenv("SEQUANT_UT_DRYRUN_BATCH")) batch = b;
  bool const batch_aux = (batch == "aux" || batch == "aux_occ" ||
                          batch == "aux_pao" || batch == "aux_pao_occ");
  bool const batch_pao =
      (batch == "pao" || batch == "aux_pao" || batch == "aux_pao_occ");
  bool const batch_occ = (batch == "aux_occ" || batch == "aux_pao_occ");
  // Block targets (SEQUANT_UT_DRYRUN_AUX_TS / OCC_TS / PAO_TS override the
  // defaults, e.g. to mirror a production input's batch:*_target_size).
  auto env_size = [](char const* k, std::size_t dflt) -> std::size_t {
    if (char const* v = std::getenv(k))
      return static_cast<std::size_t>(std::atoll(v));
    return dflt;
  };
  std::size_t const kAuxBlock = env_size("SEQUANT_UT_DRYRUN_AUX_TS", 256),
                    kOccBlock = env_size("SEQUANT_UT_DRYRUN_OCC_TS", 8),
                    kPaoBlock = env_size("SEQUANT_UT_DRYRUN_PAO_TS", 256);

  sequant::BatchPolicy policy;
  // SEQUANT_UT_DRYRUN_OCC_CONTRACTED=1: the occupied space is ALSO a batchable
  // CONTRACTED index (MPQC wires occupied batching for external indices only),
  // so a consumer that contracts an occupied pair can loop over it in batches
  // instead of reading the producer's whole assembled value.
  bool const occ_contracted =
      batch_occ && std::getenv("SEQUANT_UT_DRYRUN_OCC_CONTRACTED") != nullptr;
  policy.is_batchable_contracted_index =
      [batch_aux, batch_pao, occ_contracted](sequant::Index const& ix) {
        return (batch_aux && ix.space().base_key() == L"Κ") ||
               (batch_pao && ix.space().base_key() == L"μ̃") ||
               (occ_contracted && ix.space().base_key() == L"i");
      };
  policy.is_batchable_external_index = [batch_occ](sequant::Index const& ix) {
    return batch_occ && ix.space().base_key() == L"i";
  };
  policy.batch_spectator_indices = batch_occ;
  policy.node_level_placement = batch_occ;  // occ external placement needs it
  policy.batch_target_size =
      [batch_aux, batch_pao, batch_occ, kAuxBlock, kOccBlock,
       kPaoBlock](sequant::Index const& ix) -> std::size_t {
    if (batch_aux && ix.space().base_key() == L"Κ") return kAuxBlock;
    if (batch_pao && ix.space().base_key() == L"μ̃") return kPaoBlock;
    if (batch_occ && ix.space().base_key() == L"i") return kOccBlock;
    return 1;
  };
  policy.is_volatile_leaf = [](sequant::Tensor const& t) {
    return t.label() == L"t";
  };
  policy.accumulation_factor = 1.0;
  policy.persistent_only = false;
  // Peak budget (SEQUANT_UT_DRYRUN_PEAK_THR_GB overrides; default 100 GB).
  policy.peak_threshold =
      (std::getenv("SEQUANT_UT_DRYRUN_PEAK_THR_GB")
           ? std::atof(std::getenv("SEQUANT_UT_DRYRUN_PEAK_THR_GB"))
           : 100.0) *
      1e9;

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

  std::vector<Node> forest;
  for (std::size_t s = 0; s < nterms; ++s) {
    sequant::ExprPtr const term =
        orderedexec_witness_flatten_product(summands[s]);
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

  auto const block_of = [](sequant::Index const&) -> std::size_t {
    return 256;
  };
  auto const rich = sequant::eval::compute_dag_boulevard(forest, *cm, block_of);
  REQUIRE(!rich.cells.empty());
  auto const vmap = sequant::eval::build_value_node_map(forest);

  using annot_t = std::remove_cvref_t<decltype(forest.front()->annot())>;
  annot_t const layout{};
  sequant::eval::dryrun::DryRunLeafEvaluator const yield{cm};
  // The executor slices each block's axis at THIS target; it must be the
  // policy's own per-axis target (the one the optimizer planned with) -- a
  // flat 256 left the 80-orbital occupied axis as ONE batch, so the occupied
  // loops were degenerate: loop structure (and lost persistence) without any
  // slicing.
  std::function<std::size_t(sequant::Index const&)> const target =
      policy.batch_target_size;

  auto& logger = sequant::Logger::instance();
  auto const prev_level = logger.eval.level;
  auto* const prev_stream = logger.eval.stream;
  std::ostringstream sink_os;
  logger.eval.level = 2;  // arms tally_build AND working_set_hwmark/PeakMonitor
  logger.eval.stream = &sink_os;

  auto const legality = sequant::eval::analyze_legality(rich, forest, policy);
  // Seed batch axes: {} none, {Κ} aux, {Κ,i} aux_occ. The ordered-schedule
  // builder may throw on some batch shapes (e.g. the occ+aux nested-batch-group
  // path has a known topo-sort assert); tolerate it so the FOREST rows and the
  // note still print.
  std::optional<sequant::eval::OrderedSchedule> ordered_opt;
  try {
    if (batch_aux && !batch_pao && !batch_occ)
      // Pure aux: Κ is the sole (hence innermost) axis.
      ordered_opt =
          sequant::eval::build_ordered_schedule(rich, legality, policy, {L"Κ"});
    else
      // Empty mode_order matches production (cck.ipp); base-key sort nests the
      // realized axes by base_key. Contracted axes (Κ, μ̃) nest cleanly; an
      // external occ (i) nests OUTERMOST and forces a non-innermost split that
      // is still unimplemented -- so any *_occ build asserts, caught below.
      ordered_opt =
          sequant::eval::build_ordered_schedule(rich, legality, policy, {});
  } catch (std::exception const& e) {
    WARN("build_ordered_schedule failed for BATCH=" << batch << ": "
                                                    << e.what());
  }
  bool const ordered_ok = ordered_opt.has_value();
  if (ordered_ok) REQUIRE(sequant::eval::well_formed(*ordered_opt));
  if (ordered_ok && std::getenv("SEQUANT_UT_SCHED_TREE"))
    ut_dump_sched_tree(*ordered_opt);
  if (char const* rd = std::getenv("SEQUANT_UT_ROLE_DIAG"))
    ut_dump_role_diag(rich, legality, rd);

  std::function<bool(Node const&)> const is_volatile_node =
      [p = policy.is_volatile_leaf](Node const& n) -> bool {
    if (!n.leaf() || !n->is_tensor()) return false;
    return p && p(n->as_tensor());
  };

  // # of {Κ} ScopeBlocks the ordered schedule realized.
  std::size_t n_Kblocks = 0;
  {
    std::function<void(sequant::eval::ScopeBlock const&)> cnt =
        [&](sequant::eval::ScopeBlock const& b) {
          if (b.axis.space().base_key() == L"Κ") ++n_Kblocks;
          for (auto const& s : b.steps)
            if (auto const* ch =
                    std::get_if<sequant::eval::ScopeBlock>(&s.value))
              cnt(*ch);
        };
    if (ordered_ok) cnt(ordered_opt->root);
  }

  auto const total_builds = [](auto const& tally) -> std::size_t {
    std::size_t b = 0;
    for (auto const& [n, t] : tally)
      for (auto const& [sig, bc] : t.slices) b += bc.count;
    return b;
  };

  // ONE CostSink on the shared model; per-iteration FLOPs = its delta.
  sequant::eval::dryrun::CostSink sink;
  cm->set_cost_sink(&sink);

  struct Row {
    std::size_t builds = 0;
    double flops = 0;
    double exec = 0;        // realized roofline cost (sink.exec)
    std::size_t n_ops = 0;  // executed product ops (per-batch, sink.n_ops)
    std::size_t peak = 0;
  };
  Row f1, f2, o1, o2;

  // Optional peak-liveset dump (SEQUANT_UT_PEAK_LIVESET) of the ordered COLD
  // run: what set of co-resident values realizes the peak, labeled by node-kind
  // and carried-index signature, flagging whether each carries Κ. Answers
  // whether the aux-only peak is dominated by aux-FREE (unsliceable) values.
  bool const dump_liveset = std::getenv("SEQUANT_UT_PEAK_LIVESET") != nullptr;
  // hash -> rich.cells index, for labeling a live entry by its carried
  // signature.
  std::unordered_map<std::size_t, std::size_t> hash_to_cell;
  for (std::size_t i = 0; i < rich.cells.size(); ++i)
    hash_to_cell.emplace(rich.cells[i].hash, i);
  auto const space_sig =
      [](sequant::container::svector<sequant::Index> const& v) -> std::wstring {
    std::wstring s;
    for (std::size_t i = 0; i < v.size(); ++i) {
      if (i) s += L",";
      s += std::wstring(v[i].space().base_key());
    }
    return s;
  };
  auto const node_kind = [&](std::size_t hash) -> std::wstring {
    auto it = vmap.find(hash);
    if (it == vmap.end()) return L"<?>";
    if (it->second.leaf()) return L"leaf";
    return L"I";
  };
  // SEQUANT_UT_VALUE_LABELS: name every value of the ordered schedule
  // (value id, truncated hash, kind{carried-space signature}), so the
  // executor's [READ]/[BLOCK] diagnostics (which speak value ids) can be
  // read by label.
  if (std::getenv("SEQUANT_UT_VALUE_LABELS")) {
    std::wcerr << L"--- [value-labels] vid  h  label\n";
    for (std::size_t v = 0; v < rich.cells.size(); ++v)
      std::wcerr << L"  vid=" << v << L" h=" << (rich.cells[v].hash % 100000)
                 << L"  " << node_kind(rich.cells[v].hash) << L"{"
                 << space_sig(rich.cells[v].carried) << L"}\n";
  }

  // ---------- ORDERED / DAG executor, 2 iterations ----------
  if (ordered_ok) {
    auto cache = sequant::cache_manager(forest, is_volatile_node, 2);
    auto aops = sequant::eval::dryrun::make_dryrun_array_ops(cm);
    cache.set_array_ops(&aops);
    cache.set_recompute_tally_enabled(true);
    sequant::eval::PeakMonitor mon;
    std::vector<sequant::eval::PeakLiveEntry> peak_live;
    std::size_t peak_total = 0;
    mon.on_peak_liveset =
        [&](std::size_t total,
            std::vector<sequant::eval::PeakLiveEntry> const& live) {
          if (total >= peak_total) {  // keep the largest-total co-resident set
            peak_total = total;
            peak_live = live;
          }
        };
    cache.set_peak_monitor(&mon);
    auto run = [&](Row& r) {
      double const f0 = sink.flops.load();
      double const e0 = sink.exec.load();
      std::size_t const n0 = sink.n_ops.load();
      std::size_t const b0 = total_builds(cache.recompute_tally());
      mon.hwmark_bytes = 0;
      try {
        (void)sequant::eval::evaluate_ordered_schedule<sequant::Trace::On>(
            forest, *ordered_opt, rich, layout, yield, cache, target, {},
            is_volatile_node);
      } catch (std::exception const& e) {
        WARN("ordered evaluate threw: " << e.what());
      }
      r.flops = sink.flops.load() - f0;
      r.exec = sink.exec.load() - e0;
      r.n_ops = sink.n_ops.load() - n0;
      r.builds = total_builds(cache.recompute_tally()) - b0;
      r.peak = mon.hwmark_bytes;
    };
    run(o1);  // COLD

    if (dump_liveset) {
      // Label each co-resident entry; flag carriesΚ via its carried signature.
      struct LE {
        std::size_t bytes;
        bool carriesK;
        std::wstring label;
        std::size_t hash;
      };
      std::vector<LE> rows;
      std::size_t all_sum = 0, auxfree_sum = 0;
      for (auto const& e : peak_live) {
        std::wstring sig, kind = node_kind(e.hash);
        bool carriesK = false;
        if (auto hc = hash_to_cell.find(e.hash); hc != hash_to_cell.end()) {
          auto const& carried = rich.cells[hc->second].carried;
          sig = space_sig(carried);
          carriesK = sig.find(L"Κ") != std::wstring::npos;
        }
        all_sum += e.bytes;
        if (!carriesK) auxfree_sum += e.bytes;
        rows.push_back({e.bytes, carriesK, kind + L"{" + sig + L"}", e.hash});
      }
      std::sort(rows.begin(), rows.end(),
                [](LE const& a, LE const& b) { return a.bytes > b.bytes; });
      auto const GB = [](std::size_t b) { return double(b) / 1e9; };
      std::wcerr << L"\n--- [peak-liveset] ordered COLD peak co-resident set, "
                 << L"entries > 0.1 GB (peak_total=" << GB(peak_total)
                 << L" GB) ---\n";
      for (auto const& r : rows)
        if (GB(r.bytes) > 0.1) {
          std::wcerr << L"  " << GB(r.bytes) << L" GB  carriesΚ="
                     << (r.carriesK ? L"yes" : L"no ") << L"  h="
                     << (r.hash % 100000) << L"  " << r.label << L"  vids=";
          for (auto const& vc : rich.cells)
            if (vc.hash == r.hash) std::wcerr << vc.value_id << L",";
          std::wcerr << L"\n";
        }
      std::wcerr << L"  TOTAL co-resident = " << GB(all_sum)
                 << L" GB;  aux-FREE (no-Κ, aux-batching-immune) floor = "
                 << GB(auxfree_sum) << L" GB\n";
    }

    cache.reset();
    run(o2);  // WARM
    // SEQUANT_UT_TOP_BUILDS=<N>: the N values built most often over the two
    // iterations (the recompute tally), with their labels and occurrence
    // counts -- where a fused chain's interposed loops rebuild a value.
    if (char const* tb = std::getenv("SEQUANT_UT_TOP_BUILDS")) {
      std::size_t const n_top = std::strtoul(tb, nullptr, 10);
      std::vector<std::tuple<std::size_t, std::size_t>> rows;  // builds, vid
      for (auto const& vc : rich.cells) {
        if (vc.is_leaf) continue;
        auto const it = vmap.find(sequant::eval::value_key_of(vc));
        if (it == vmap.end()) continue;
        rows.emplace_back(
            orderedexec_builds_of(cache.recompute_tally(), it->second),
            vc.value_id);
      }
      std::sort(rows.begin(), rows.end(), [](auto const& a, auto const& b) {
        return std::get<0>(a) > std::get<0>(b);
      });
      std::wcerr << L"--- [top-builds] builds vid h label occurrences\n";
      for (std::size_t k = 0; k < rows.size() && k < n_top; ++k) {
        auto const [b, v] = rows[k];
        std::wcerr << L"  " << b << L"  vid=" << v << L" h="
                   << (rich.cells[v].hash % 100000) << L"  "
                   << node_kind(rich.cells[v].hash) << L"{"
                   << space_sig(rich.cells[v].carried) << L"}  occ="
                   << rich.cells[v].occurrences.size() << L"\n";
      }
    }
  }

  // Warm-iter needed_build skip count: replicate the executor's gate (BFS from
  // volatile roots halting at cache-alive) on a fresh warm cache, then count
  // the non-leaf BuildSteps inside {Κ} blocks that the gate would skip.
  std::size_t warm_skipped = 0;
  if (ordered_ok) {
    auto cache = sequant::cache_manager(forest, is_volatile_node, 2);
    auto aops = sequant::eval::dryrun::make_dryrun_array_ops(cm);
    cache.set_array_ops(&aops);
    // prime persistents by one cold run, then reset (persistents stay resident)
    logger.eval.stream = &sink_os;
    try {
      (void)sequant::eval::evaluate_ordered_schedule<sequant::Trace::On>(
          forest, *ordered_opt, rich, layout, yield, cache, target, {},
          is_volatile_node);
    } catch (std::exception const&) {
    }
    cache.reset();
    sequant::container::set<std::size_t> needed;
    sequant::container::svector<Node> stack;
    sequant::container::set<std::size_t> visited;
    for (auto&& n : forest) stack.push_back(n);
    while (!stack.empty()) {
      Node const n = stack.back();
      stack.pop_back();
      if (n.leaf()) continue;
      if (!visited.insert(n->hash_value()).second) continue;
      if (cache.alive(n)) continue;
      needed.insert(n->hash_value());
      stack.push_back(n.left());
      stack.push_back(n.right());
    }
    std::function<void(sequant::eval::ScopeBlock const&, bool)> sc =
        [&](sequant::eval::ScopeBlock const& b, bool inK) {
          bool const here = inK || b.axis.space().base_key() == L"Κ";
          if (here)
            for (auto const& s : b.steps)
              if (auto const* bs =
                      std::get_if<sequant::eval::BuildStep>(&s.value)) {
                auto const it = vmap.find(
                    sequant::eval::value_key_of(rich.cells[bs->value_id]));
                if (it != vmap.end() && !it->second.leaf() &&
                    !needed.count(rich.cells[bs->value_id].hash))
                  ++warm_skipped;
              }
          for (auto const& s : b.steps)
            if (auto const* ch =
                    std::get_if<sequant::eval::ScopeBlock>(&s.value))
              sc(*ch, here);
        };
    sc(ordered_opt->root, false);
  }

  // ---------- FOREST descent, 2 iterations ----------
  {
    auto cache = sequant::cache_manager(forest, is_volatile_node, 2);
    auto aops = sequant::eval::dryrun::make_dryrun_array_ops(cm);
    cache.set_array_ops(&aops);
    cache.set_recompute_tally_enabled(true);
    sequant::eval::PeakMonitor mon;
    cache.set_peak_monitor(&mon);
    auto run = [&](Row& r) {
      double const f0 = sink.flops.load();
      double const e0 = sink.exec.load();
      std::size_t const n0 = sink.n_ops.load();
      std::size_t const b0 = total_builds(cache.recompute_tally());
      mon.hwmark_bytes = 0;
      std::atomic<double> peak{0.0};
      for (auto const& root : forest) {
        cache.set_custom_evaluator(sequant::make_evaluator(
            policy, yield, sequant::make_no_scope_guard{}, &peak));
        try {
          (void)sequant::evaluate<sequant::Trace::On>(root, yield, cache);
        } catch (std::exception const&) {
        }
      }
      r.flops = sink.flops.load() - f0;
      r.exec = sink.exec.load() - e0;
      r.n_ops = sink.n_ops.load() - n0;
      r.builds = total_builds(cache.recompute_tally()) - b0;
      r.peak = std::max<std::size_t>(
          mon.hwmark_bytes,
          std::max<std::size_t>(cache.working_set_hwmark(),
                                static_cast<std::size_t>(peak.load())));
      // DIAG: which of the three peak sources dominates
      std::wcerr << L"  [peak-components] monitor_hwmark=" << mon.hwmark_bytes
                 << L" cache_working_set_hwmark=" << cache.working_set_hwmark()
                 << L" sink_peak=" << static_cast<std::size_t>(peak.load())
                 << L"\n";
    };
    run(f1);
    cache.reset();
    run(f2);
    // SEQUANT_UT_TOP_BUILDS=<N>, forest side: the same tally over the
    // forest-descent evaluation, for the per-value forest-vs-ordered
    // comparison.
    if (char const* tb = std::getenv("SEQUANT_UT_TOP_BUILDS")) {
      std::size_t const n_top = std::strtoul(tb, nullptr, 10);
      std::vector<std::tuple<std::size_t, std::size_t>> rows;
      for (auto const& vc : rich.cells) {
        if (vc.is_leaf) continue;
        auto const it = vmap.find(sequant::eval::value_key_of(vc));
        if (it == vmap.end()) continue;
        rows.emplace_back(
            orderedexec_builds_of(cache.recompute_tally(), it->second),
            vc.value_id);
      }
      std::sort(rows.begin(), rows.end(), [](auto const& a, auto const& b) {
        return std::get<0>(a) > std::get<0>(b);
      });
      std::wcerr << L"--- [top-builds-forest] builds vid h label occurrences\n";
      for (std::size_t k = 0; k < rows.size() && k < n_top; ++k) {
        auto const [b, v] = rows[k];
        std::wcerr << L"  " << b << L"  vid=" << v << L" h="
                   << (rich.cells[v].hash % 100000) << L"  "
                   << node_kind(rich.cells[v].hash) << L"{"
                   << space_sig(rich.cells[v].carried) << L"}  occ="
                   << rich.cells[v].occurrences.size() << L"\n";
      }
    }
  }

  cm->set_cost_sink(nullptr);
  logger.eval.level = prev_level;
  logger.eval.stream = prev_stream;

  auto pr = [](wchar_t const* tag, Row const& r) {
    std::wcerr << L"  " << tag << L"  builds=" << r.builds << L"  ops="
               << r.n_ops << L"  exec=" << r.exec << L"  FLOPs="
               << std::scientific << r.flops << L"  peak_bytes=" << r.peak
               << L"\n";
  };
  std::wcerr << L"\n=== [dryrun-2iter-report] SYSTEM="
             << std::wstring(system.begin(), system.end()) << L" BATCH="
             << std::wstring(batch.begin(), batch.end()) << L", "
             << forest.size() << L" terms, n_Kblocks=" << n_Kblocks
             << L", warm-iter needed_build-skipped {Κ}-block BuildSteps="
             << warm_skipped << L", peak_threshold=" << std::scientific
             << policy.peak_threshold << L" B, blocks aux/occ/pao=" << kAuxBlock
             << L"/" << kOccBlock << L"/" << kPaoBlock
             << (occ_contracted ? L", occ also CONTRACTED-batchable" : L"")
             << L" ===\n";
  pr(L"forest  iter1 (cold)", f1);
  pr(L"forest  iter2 (warm)", f2);
  if (ordered_ok) {
    pr(L"ordered iter1 (cold)", o1);
    pr(L"ordered iter2 (warm)", o2);
  } else {
    std::wcerr << L"  ordered: SCHEDULE BUILD FAILED for this batch mode "
                  L"(pre-existing build_ordered_schedule assert) -- forest "
                  L"rows only\n";
  }

  // Light sanity: warm iter never builds MORE than the cold iter.
  CHECK(f2.builds <= f1.builds);
  if (ordered_ok) CHECK(o2.builds <= o1.builds);
}

// ===========================================================================
// b3 (2026-08-12 eager-home-release plan, Task b): dry == wet schedule
// EQUIVALENCE for the ordered executor. b1 fixed a real dry-run/wet-run
// fidelity bug in meter.hpp: the install `if (!policy.whole_scope_execution)
// cache.set_custom_evaluator(...)` used to fire for BatchScheduler::ordered
// too (ordered is neither whole_scope nor, under the OLD two-bool encoding,
// distinguishable from forest_descent by that single negated check), routing
// the dry ordered replay's root-level BuildSteps through the FOREST custom
// evaluator instead of evaluate_ordered_schedule's own run_ordered_
// contracted_block -- a real WET/DRY divergence, since MPQC's wet ordered
// path (cck.ipp's `cache.set_whole_scope_driver`) installs NO custom
// evaluator at all. This test proves the fix: meter()'s ordered replay
// (routed through the SAME sequant::evaluate(Nodes const&, BatchPolicy
// const&, ...) coexistence entry MPQC's wet path drives) realizes the exact
// SAME peak as a directly-invoked evaluate_ordered_schedule call wired with
// its own PeakMonitor -- the "wet-style" invocation the witness TEST_CASE
// above uses. Reuses the SAME water-20 fixture (forest/policy/rich
// construction) as that witness -- SP3 Task 4's expensive optimize+binarize
// scaffold is not rebuilt from scratch elsewhere in this file, so it is
// duplicated here as-is, matching the file's own stated precedent for these
// fixtures (see the file-header note above the witness TEST_CASE).
// ===========================================================================
TEST_CASE(
    "b3: meter()'s ordered replay peak equals a direct evaluate_ordered_"
    "schedule (PeakMonitor) replay on the water-20 aux-only residual forest",
    "[.][ordered-executor][meter]") {
  using sequant::eval::dryrun::CacheConfig;
  using sequant::eval::dryrun::EvalExprDryRun;
  using sequant::eval::dryrun::EvalNodeDryRun;
  using sequant::eval::dryrun::meter;
  using Node = EvalNodeDryRun;

  auto ctx = sequant::get_default_context().clone();
  ctx.set_first_dummy_index_ordinal(1000000);
  auto isr = ctx.mutable_index_space_registry();
  REQUIRE(isr != nullptr);
  sequant::mbpt::add_pao_spaces(isr, sequant::mbpt::Spin::any);
  sequant::mbpt::add_df_spaces(isr);
  auto ctx_resetter = sequant::set_scoped_default_context(std::move(ctx));

  auto const body =
      orderedexec_witness_slurp(std::string(SEQUANT_UNIT_TESTS_SOURCE_DIR) +
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

  std::size_t nterms = std::min<std::size_t>(summands.size(), 40);
  if (char const* nt = std::getenv("SEQUANT_UT_DRYRUN_NTERMS"))
    nterms = std::min<std::size_t>(summands.size(), std::atoll(nt));

  auto regime = orderedexec_witness_df_regime(kOrderedExecWater20_pVDZF12);
  auto cm = std::make_shared<sequant::eval::dryrun::CostModel const>(regime);

  // EXACT MPQC aux-only config (make_csv_batch_policy, aux_target=256), plus
  // the b0 enum-based scheduler selection: BatchScheduler::ordered.
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
  policy.scheduler = sequant::BatchScheduler::ordered;

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

  std::vector<Node> forest;
  for (std::size_t s = 0; s < nterms; ++s) {
    sequant::ExprPtr const term =
        orderedexec_witness_flatten_product(summands[s]);
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

  // ONE forest + ONE rich shared by both replays (mirrors the witness's own
  // R2 invariant -- meter() below builds its OWN internal rich from the SAME
  // forest/block_of, so this local `rich` is only used to build the direct
  // ordered schedule, not shared code with meter()).
  auto const block_of = [](sequant::Index const&) -> std::size_t {
    return 256;
  };
  auto const rich = sequant::eval::compute_dag_boulevard(forest, *cm, block_of);
  REQUIRE(!rich.cells.empty());

  using annot_t = std::remove_cvref_t<decltype(forest.front()->annot())>;
  annot_t const layout{};
  sequant::eval::dryrun::DryRunLeafEvaluator const yield{cm};
  std::function<std::size_t(sequant::Index const&)> const target =
      [](sequant::Index const&) -> std::size_t { return 256; };

  std::function<bool(Node const&)> const is_volatile_node =
      [p = policy.is_volatile_leaf](Node const& n) -> bool {
    if (!n.leaf() || !n->is_tensor()) return false;
    return p && p(n->as_tensor());
  };

  auto& logger = sequant::Logger::instance();
  auto const prev_level = logger.eval.level;
  auto* const prev_stream = logger.eval.stream;
  logger.eval.level = 1;  // arms tally_build (DryRunOps::prod's runtime gate)

  // ---- (A) direct evaluate_ordered_schedule replay, "wet-style": a fresh
  // cache + PeakMonitor, no custom evaluator installed (MPQC's wet ordered
  // path -- cck.ipp's set_whole_scope_driver -- never installs one either).
  auto const legality = sequant::eval::analyze_legality(rich, forest, policy);
  REQUIRE(legality.cells.size() == rich.cells.size());
  auto const ordered =
      sequant::eval::build_ordered_schedule(rich, legality, policy, {L"Κ"});
  REQUIRE(sequant::eval::well_formed(ordered));

  std::ostringstream direct_trace;
  logger.eval.stream = &direct_trace;
  auto direct_cache = sequant::cache_manager(forest);
  auto aops = sequant::eval::dryrun::make_dryrun_array_ops(cm);
  direct_cache.set_array_ops(&aops);
  direct_cache.set_recompute_tally_enabled(true);
  sequant::eval::PeakMonitor direct_mon;
  direct_cache.set_peak_monitor(&direct_mon);
  ResultPtr direct_result;
  try {
    direct_result =
        sequant::eval::evaluate_ordered_schedule<sequant::Trace::On>(
            forest, ordered, rich, layout, yield, direct_cache, target, {},
            is_volatile_node);
  } catch (std::exception const& e) {
    std::cerr << "[b3-ordered-dry-wet-equivalence] direct evaluate threw: "
              << e.what() << "\n";
  }
  REQUIRE(direct_result);

  // ---- (B) meter()'s own ordered replay, on the SAME forest/policy/regime.
  // meter() builds its own internal cache (persistence-aware, is_volatile
  // matching policy.is_volatile_leaf via cfg.is_volatile below) and drives
  // the SAME coexistence entry (sequant::evaluate(Nodes const&, BatchPolicy
  // const&, ...), scope_executor.hpp) that MPQC's wet whole_scope/ordered
  // path (cck.ipp's set_whole_scope_driver) drives -- the point of this test.
  std::ostringstream meter_trace;
  logger.eval.stream = &meter_trace;
  CacheConfig cfg;
  cfg.is_volatile = is_volatile_node;
  cfg.min_repeats = 2;
  cfg.max_footprint = 0.;  // no footprint gate, matching the direct replay
  auto const meter_report = meter(forest, policy, regime, cfg);

  logger.eval.level = prev_level;
  logger.eval.stream = prev_stream;

  std::wcerr << L"\n=== [b3-ordered-dry-wet-equivalence] water-20 aux-only, "
             << forest.size() << L" terms ===\n"
             << L"  direct (wet-style) evaluate_ordered_schedule peak = "
             << direct_mon.hwmark_bytes << L" B ("
             << (double(direct_mon.hwmark_bytes) / 1e9) << L" GB)\n"
             << L"  meter() ordered replay peak                      = "
             << meter_report.peak_bytes << L" B ("
             << (meter_report.peak_bytes / 1e9) << L" GB)\n";
  INFO("direct (wet-style) peak = " << direct_mon.hwmark_bytes
                                    << " B; meter() peak = "
                                    << meter_report.peak_bytes << " B");
  CAPTURE(direct_mon.hwmark_bytes, meter_report.peak_bytes);

  CHECK(direct_mon.hwmark_bytes > 0);
  CHECK(meter_report.peak_bytes > 0.0);
  CHECK(meter_report.scheduler == sequant::BatchScheduler::ordered);

  // THE payoff: dry (meter) == wet-style (direct) peak, exactly. If these
  // disagree, that is a REAL remaining dry/wet fidelity gap -- do not loosen
  // this to a tolerance to force a pass; a mismatch means b1's fix is
  // incomplete or another divergence exists.
  CHECK(meter_report.peak_bytes ==
        Catch::Approx(double(direct_mon.hwmark_bytes)).margin(1.0));
}

// ===========================================================================
// ANALYSIS PROBE (uncommitted diagnostic): decompose the water-20 ORDERED
// executor dense-model peak into tier A (root-homed Κ-free composites) vs
// tier B (Κ-carrying aux-loop working set), and snapshot the co-resident set
// at the instant of realized high-water. Reuses the same fixture as the
// witness above. Run by name: [.][w20-peak-composition].
// ===========================================================================
TEST_CASE("w20 peak composition: tier-A/tier-B decomposition at realized peak",
          "[.][w20-peak-composition]") {
  using sequant::eval::dryrun::EvalExprDryRun;
  using sequant::eval::dryrun::EvalNodeDryRun;
  using Node = EvalNodeDryRun;

  auto ctx = sequant::get_default_context().clone();
  ctx.set_first_dummy_index_ordinal(1000000);
  auto isr = ctx.mutable_index_space_registry();
  REQUIRE(isr != nullptr);
  sequant::mbpt::add_pao_spaces(isr, sequant::mbpt::Spin::any);
  sequant::mbpt::add_df_spaces(isr);
  auto ctx_resetter = sequant::set_scoped_default_context(std::move(ctx));

  auto const body =
      orderedexec_witness_slurp(std::string(SEQUANT_UNIT_TESTS_SOURCE_DIR) +
                                "/data/csv_ccsd_doubles_residual_df.txt");
  REQUIRE(!body.empty());
  std::string line = body;
  if (auto nl = line.find('\n'); nl != std::string::npos)
    line = line.substr(0, nl);
  auto expr = sequant::deserialize<sequant::ExprPtr>(line);
  REQUIRE(expr->is<sequant::Sum>());
  auto const& summands = expr->as<sequant::Sum>().summands();

  std::size_t nterms = std::min<std::size_t>(summands.size(), 40);
  if (char const* nt = std::getenv("SEQUANT_UT_DRYRUN_NTERMS"))
    nterms = std::min<std::size_t>(summands.size(), std::atoll(nt));

  auto regime = orderedexec_witness_df_regime(kOrderedExecWater20_pVDZF12);
  auto cm = std::make_shared<sequant::eval::dryrun::CostModel const>(regime);

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

  std::vector<Node> forest;
  for (std::size_t s = 0; s < nterms; ++s) {
    sequant::ExprPtr const term =
        orderedexec_witness_flatten_product(summands[s]);
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

  auto const block_of = [](sequant::Index const&) -> std::size_t {
    return 256;
  };
  auto const rich = sequant::eval::compute_dag_boulevard(forest, *cm, block_of);
  REQUIRE(!rich.cells.empty());

  auto const legality = sequant::eval::analyze_legality(rich, forest, policy);
  auto const ordered =
      sequant::eval::build_ordered_schedule(rich, legality, policy, {L"Κ"});
  REQUIRE(sequant::eval::well_formed(ordered));

  auto const vmap = sequant::eval::build_value_node_map(forest);

  // hash -> ValueCell index
  std::unordered_map<std::size_t, std::size_t> hash_to_cell;
  for (auto const& vc : rich.cells) hash_to_cell.emplace(vc.hash, vc.value_id);
  // value_id -> cell index (identity here, but explicit)
  auto const cell_of_vid =
      [&](std::size_t vid) -> sequant::eval::ValueCell const& {
    return rich.cells[vid];
  };

  auto const is_K = [](sequant::Index const& ix) {
    return ix.space().base_key() == L"Κ";
  };
  auto const space_sig =
      [](sequant::container::svector<sequant::Index> const& v) -> std::wstring {
    std::wstring s;
    for (std::size_t i = 0; i < v.size(); ++i) {
      if (i) s += L",";
      s += std::wstring(v[i].space().base_key());
    }
    return s;
  };
  auto const foot = [&](sequant::eval::ValueCell const& vc) -> std::size_t {
    return sequant::eval::detail::cell_footprint(vc.carried, vc.home_modes, *cm,
                                                 block_of);
  };
  auto const node_kind = [&](std::size_t hash) -> std::wstring {
    auto it = vmap.find(hash);
    if (it == vmap.end()) return L"<?>";
    if (it->second.leaf()) return L"leaf:" + it->second->to_latex();
    return L"I";
  };

  // ---------------- (1) TIER A: root-level BuildSteps ----------------
  std::vector<std::size_t> tierA_vids;
  for (auto const& step : ordered.root.steps)
    if (auto const* b = std::get_if<sequant::eval::BuildStep>(&step.value))
      tierA_vids.push_back(b->value_id);

  struct Row {
    std::size_t vid, bytes;
    std::wstring label;
  };
  std::vector<Row> tierA;
  std::size_t tierA_sum = 0;
  for (auto vid : tierA_vids) {
    auto const& vc = cell_of_vid(vid);
    std::size_t const b = foot(vc);
    tierA_sum += b;
    tierA.push_back(
        {vid, b, node_kind(vc.hash) + L"{" + space_sig(vc.carried) + L"}"});
  }
  std::sort(tierA.begin(), tierA.end(),
            [](Row const& a, Row const& b) { return a.bytes > b.bytes; });

  // ---------------- (2) TIER B: aux ScopeBlock nested steps ----------------
  std::vector<Row> tierB;
  std::size_t tierB_sum = 0;
  std::optional<Row> a1_row;
  std::function<void(sequant::eval::ScopeBlock const&)> walk_block =
      [&](sequant::eval::ScopeBlock const& blk) {
        auto record = [&](std::size_t vid, bool is_output,
                          sequant::eval::OutputKind kind) {
          auto const& vc = cell_of_vid(vid);
          std::size_t const b = foot(vc);
          tierB_sum += b;
          std::wstring tag =
              is_output ? (kind == sequant::eval::OutputKind::AccumulateSum
                               ? L" [out:Sum]"
                               : L" [out:Scatter]")
                        : L" [build]";
          Row r{vid, b,
                node_kind(vc.hash) + L"{" + space_sig(vc.carried) + L"}" + tag};
          tierB.push_back(r);
          // a1{i1,i2;mutilde1;K1}: the K1-carrying PPL prerequisite. Match any
          // Κ-carrying LoopLocal [build] cell (NOT an escape output -- those
          // are Κ-reductions that carry Κ only BELOW), pick the LARGEST
          // per-block.
          bool const carries_K =
              std::any_of(vc.carried.begin(), vc.carried.end(), is_K);
          if (!is_output && carries_K && (!a1_row || b > a1_row->bytes))
            a1_row = r;
        };
        for (auto const& step : blk.steps)
          if (auto const* b =
                  std::get_if<sequant::eval::BuildStep>(&step.value))
            record(b->value_id, false,
                   sequant::eval::OutputKind::AccumulateSum);
        for (auto const& [vid, kind] : blk.outputs) record(vid, true, kind);
        for (auto const& step : blk.steps)
          if (auto const* child =
                  std::get_if<sequant::eval::ScopeBlock>(&step.value))
            walk_block(*child);
      };
  // enumerate every root-level ScopeBlock (the aux Κ loop(s))
  for (auto const& step : ordered.root.steps)
    if (auto const* child = std::get_if<sequant::eval::ScopeBlock>(&step.value))
      walk_block(*child);
  std::sort(tierB.begin(), tierB.end(),
            [](Row const& a, Row const& b) { return a.bytes > b.bytes; });

  // ---------------- run the ordered executor with peak monitor ----------
  using annot_t = std::remove_cvref_t<decltype(forest.front()->annot())>;
  annot_t const layout{};
  sequant::eval::dryrun::DryRunLeafEvaluator const yield{cm};
  std::function<std::size_t(sequant::Index const&)> const target =
      [](sequant::Index const&) -> std::size_t { return 256; };

  auto& logger = sequant::Logger::instance();
  auto const prev_level = logger.eval.level;
  auto* const prev_stream = logger.eval.stream;
  std::ostringstream sink;
  logger.eval.level = 1;
  logger.eval.stream = &sink;

  auto ordered_cache = sequant::cache_manager(forest);
  auto aops = sequant::eval::dryrun::make_dryrun_array_ops(cm);
  ordered_cache.set_array_ops(&aops);
  sequant::eval::PeakMonitor mon;
  std::size_t peak_total_snapshot = 0;
  std::size_t peak_clock = 0;
  std::vector<sequant::eval::PeakLiveEntry> peak_liveset;
  mon.on_peak_liveset =
      [&](std::size_t total,
          std::vector<sequant::eval::PeakLiveEntry> const& live) {
        peak_total_snapshot = total;  // last (== max) advance wins
        peak_liveset = live;
        // DIAGNOSTIC: timestamp THIS high-water in the access-clock timeline.
        // on_peak_liveset fires (in note_working_set) BEFORE observe advances
        // the mark, and only on a genuine new global high-water, so the FINAL
        // firing (max peak) captures the clock value at the peak instant. The
        // peak op's own operands were read just before this note_working_set
        // call, so their stamps are <= peak_clock.
        peak_clock = sequant::eval::AccessClock::now();
      };
  ordered_cache.set_peak_monitor(&mon);
  // DIAGNOSTIC: zero the access clock + last-read record for THIS measured run.
  sequant::eval::AccessClock::reset();
  ResultPtr ord_result;
  try {
    ord_result = sequant::eval::evaluate_ordered_schedule<sequant::Trace::On>(
        forest, ordered, rich, layout, yield, ordered_cache, target);
  } catch (std::exception const& e) {
    std::cerr << "[w20-peak-composition] ordered evaluate threw: " << e.what()
              << "\n";
  }
  logger.eval.level = prev_level;
  logger.eval.stream = prev_stream;

  // DIAGNOSTIC: hashes still ALIVE on the root cache AFTER the run == pinned
  // home entries (tier-A composites + block escape outputs the executor homes
  // at root). A genuinely use-counted (LoopLocal transient) value was
  // released at its last use and is NOT alive here. This
  // is the true discriminator for the tier-B sanity check: "dead-but-retained"
  // is EXPECTED for a pinned entry, but a NON-pinned (use-counted) entry that
  // looks dead-but-retained at peak signals a missed read path.
  std::unordered_set<std::size_t> pinned_now;
  ordered_cache.for_each_key([&](Node const& k) {
    if (ordered_cache.alive(k)) pinned_now.insert(k->hash_value());
  });

  auto const GB = [](std::size_t b) { return double(b) / 1e9; };
  auto const TB = [](std::size_t b) { return double(b) / 1e12; };

  std::wcerr << L"\n================= [w20-peak-composition] " << forest.size()
             << L" terms =================\n";
  std::wcerr << L"realized ordered peak (hwmark) = " << mon.hwmark_bytes
             << L" bytes (" << TB(mon.hwmark_bytes) << L" TB)\n";

  // ---- (1) tier A ----
  std::wcerr << L"\n--- TIER A: root-homed Κ-free composites (" << tierA.size()
             << L" BuildSteps), SUM = " << tierA_sum << L" bytes ("
             << TB(tierA_sum) << L" TB) ---\n";
  for (std::size_t i = 0; i < tierA.size() && i < 20; ++i)
    std::wcerr << L"  vid=" << tierA[i].vid << L"  " << GB(tierA[i].bytes)
               << L" GB  " << tierA[i].label << L"\n";

  // ---- (2) tier B ----
  std::wcerr << L"\n--- TIER B: Κ-carrying aux-loop working set ("
             << tierB.size() << L" steps), per-Κ-block SUM = " << tierB_sum
             << L" bytes (" << GB(tierB_sum) << L" GB) ---\n";
  for (std::size_t i = 0; i < tierB.size() && i < 25; ++i)
    std::wcerr << L"  vid=" << tierB[i].vid << L"  " << GB(tierB[i].bytes)
               << L" GB  " << tierB[i].label << L"\n";
  if (a1_row)
    std::wcerr << L"\n  a1{i1,i2;mutilde1;K1} LOCATED in aux block (tier B): "
               << L"vid=" << a1_row->vid << L"  per-block " << GB(a1_row->bytes)
               << L" GB  " << a1_row->label << L"\n";
  else
    std::wcerr << L"\n  a1{i1,i2;mutilde1;K1}: NOT FOUND among tier-B cells\n";

  // ---- (3) peak op ----
  std::wcerr << L"\n--- (3) PEAK OP: op_hash=" << mon.peak.op_hash
             << L"  bytes=" << mon.peak.bytes << L" (" << TB(mon.peak.bytes)
             << L" TB) ---\n";
  {
    auto it = vmap.find(mon.peak.op_hash);
    std::wstring where = L"unknown";
    bool in_tierA = std::any_of(tierA.begin(), tierA.end(), [&](Row const& r) {
      return cell_of_vid(r.vid).hash == mon.peak.op_hash;
    });
    bool in_tierB = std::any_of(tierB.begin(), tierB.end(), [&](Row const& r) {
      return cell_of_vid(r.vid).hash == mon.peak.op_hash;
    });
    if (in_tierA)
      where = L"TIER A (root-homed)";
    else if (in_tierB)
      where = L"TIER B (aux loop)";
    if (it != vmap.end()) {
      auto const hc = hash_to_cell.find(mon.peak.op_hash);
      std::wstring lab = hc != hash_to_cell.end()
                             ? node_kind(mon.peak.op_hash) + L"{" +
                                   space_sig(rich.cells[hc->second].carried) +
                                   L"}"
                             : it->second->to_latex();
      std::wcerr << L"  peak op node = " << lab << L"   [" << where << L"]\n";
    } else {
      std::wcerr << L"  peak op_hash not in vmap (transient/leaf) [" << where
                 << L"]\n";
    }
  }

  // ---- (4) co-resident set at peak ----
  std::wstring const kA = L"A", kB = L"B", kLeaf = L"leaf/input", kOther = L"?";
  std::unordered_map<std::size_t, std::wstring> tierA_hashes, tierB_hashes;
  for (auto const& r : tierA) tierA_hashes[cell_of_vid(r.vid).hash] = r.label;
  for (auto const& r : tierB) tierB_hashes[cell_of_vid(r.vid).hash] = r.label;

  std::sort(peak_liveset.begin(), peak_liveset.end(),
            [](auto const& a, auto const& b) { return a.bytes > b.bytes; });
  std::size_t liveA = 0, liveB = 0, liveLeaf = 0, liveOther = 0, liveTotal = 0;
  std::wcerr << L"\n--- (4) CO-RESIDENT SET AT PEAK (total snapshot = "
             << peak_total_snapshot << L" bytes, " << TB(peak_total_snapshot)
             << L" TB; " << peak_liveset.size() << L" alive entries) ---\n";
  for (auto const& e : peak_liveset) {
    liveTotal += e.bytes;
    std::wstring tier, lab;
    if (auto it = tierA_hashes.find(e.hash); it != tierA_hashes.end()) {
      tier = kA;
      lab = it->second;
      liveA += e.bytes;
    } else if (auto it2 = tierB_hashes.find(e.hash);
               it2 != tierB_hashes.end()) {
      tier = kB;
      lab = it2->second;
      liveB += e.bytes;
    } else {
      auto vit = vmap.find(e.hash);
      if (vit != vmap.end() && vit->second.leaf()) {
        tier = kLeaf;
        lab = L"leaf:" + vit->second->to_latex();
        liveLeaf += e.bytes;
      } else {
        tier = kOther;
        auto hc = hash_to_cell.find(e.hash);
        lab = hc != hash_to_cell.end()
                  ? node_kind(e.hash) + L"{" +
                        space_sig(rich.cells[hc->second].carried) + L"}"
                  : L"<not-in-rich>";
        liveOther += e.bytes;
      }
    }
    if (e.bytes > 1e8)  // only print entries > 0.1 GB to keep it readable
      std::wcerr << L"  [" << tier << L"] " << GB(e.bytes) << L" GB  " << lab
                 << L"\n";
  }
  std::wcerr << L"\n  co-resident totals: tierA = " << GB(liveA)
             << L" GB, tierB = " << GB(liveB) << L" GB, leaf/input = "
             << GB(liveLeaf) << L" GB, other = " << GB(liveOther)
             << L" GB ; SUM = " << GB(liveTotal) << L" GB\n";

  // ---- (5) hypothesis: ordered - whole_scope ~= tierA_sum? ----
  double const whole_scope_TB = 0.36;
  double const ordered_TB = TB(mon.hwmark_bytes);
  std::wcerr << L"\n--- (5) HYPOTHESIS: ordered - whole_scope (=" << ordered_TB
             << L" - 0.36 = " << (ordered_TB - whole_scope_TB)
             << L" TB) ~= tierA home-floor sum (" << TB(tierA_sum)
             << L" TB)? ---\n";

  // ==================================================================
  // (6) EAGER-RELEASE RECLAIM: of the tier-A home floor alive at the
  //     realized peak, how much is DEAD-BUT-RETAINED -- its genuine last
  //     cache READ (AccessClock last_access) already happened BEFORE the
  //     peak instant (peak_clock) yet it stays pinned resident-until-reset.
  //     That sum is exactly the memory an eager release-at-last-use would
  //     reclaim on THIS peak.
  // ==================================================================
  auto const& lam = sequant::eval::AccessClock::last_access_map();
  std::size_t const total_ticks = sequant::eval::AccessClock::now();
  auto const last_access_of = [&](std::size_t h) -> std::size_t {
    auto it = lam.find(h);
    return it == lam.end() ? 0 : it->second;  // 0 == never cache-read
  };

  // forest-root hashes: their result buffer is ALSO pinned by the executor's
  // value_results vector until the final combine (after peak), so releasing
  // only the cache entry would not actually free them here -- reported as a
  // caveat, not subtracted from the primary metric (which is defined purely on
  // last_access < peak).
  std::unordered_set<std::size_t> root_hashes;
  for (auto&& n : forest) root_hashes.insert(n->hash_value());

  // which hashes are alive at the realized peak (the snapshot set).
  std::unordered_set<std::size_t> alive_at_peak;
  for (auto const& e : peak_liveset) alive_at_peak.insert(e.hash);

  // peak-liveset bytes per hash (the modeled size AT peak).
  std::unordered_map<std::size_t, std::size_t> peak_bytes_of;
  for (auto const& e : peak_liveset) peak_bytes_of[e.hash] = e.bytes;

  struct DR {
    std::size_t vid, bytes, last;
    bool is_root;
    std::wstring label;
  };
  std::vector<DR> dead_rows, live_rows;
  std::size_t reclaimable = 0, reclaimable_nonroot = 0;
  std::size_t neverread_bytes = 0;
  for (auto const& r : tierA) {
    std::size_t const h = cell_of_vid(r.vid).hash;
    if (!alive_at_peak.count(h)) continue;  // only floor ALIVE at peak
    std::size_t const la = last_access_of(h);
    std::size_t const bytes =
        peak_bytes_of.count(h) ? peak_bytes_of[h] : r.bytes;
    bool const is_root = root_hashes.count(h) != 0;
    DR row{r.vid, bytes, la, is_root, r.label};
    if (la < peak_clock) {  // dead-but-retained
      reclaimable += bytes;
      if (!is_root) reclaimable_nonroot += bytes;
      if (la == 0) neverread_bytes += bytes;
      dead_rows.push_back(row);
    } else {
      live_rows.push_back(row);
    }
  }
  std::sort(dead_rows.begin(), dead_rows.end(),
            [](DR const& a, DR const& b) { return a.bytes > b.bytes; });

  std::size_t const floor = tierA_sum;
  std::size_t const peak = mon.hwmark_bytes;
  std::wcerr << L"\n============ (6) EAGER-RELEASE RECLAIM ============\n";
  std::wcerr << L"  peak_clock = " << peak_clock << L" of " << total_ticks
             << L" total reads  ("
             << (total_ticks ? 100.0 * peak_clock / total_ticks : 0.0)
             << L"% through the walk => peak is "
             << (peak_clock * 2 < total_ticks ? L"EARLY" : L"LATE") << L")\n";
  std::wcerr << L"  peak op last cache-read clock = "
             << last_access_of(mon.peak.op_hash) << L" (op_hash="
             << mon.peak.op_hash << L")\n";
  std::wcerr << L"  tier-A home floor              = " << floor << L" B ("
             << TB(floor) << L" TB)\n";
  std::wcerr << L"  realized peak                  = " << peak << L" B ("
             << TB(peak) << L" TB)\n";
  std::wcerr << L"  DEAD-BUT-RETAINED (reclaimable)= " << reclaimable << L" B ("
             << TB(reclaimable) << L" TB)\n";
  std::wcerr << L"     = " << (floor ? 100.0 * reclaimable / floor : 0.0)
             << L"% of the 0.712 TB floor, "
             << (peak ? 100.0 * reclaimable / peak : 0.0)
             << L"% of the 0.99 TB peak\n";
  std::wcerr << L"     of which NEVER cache-read (last_access==0) = "
             << neverread_bytes << L" B (" << TB(neverread_bytes) << L" TB)\n";
  std::wcerr << L"     of which value_results-pinned forest ROOTS = "
             << (reclaimable - reclaimable_nonroot) << L" B; NON-root "
             << L"(truly free to drop here) = " << reclaimable_nonroot
             << L" B (" << TB(reclaimable_nonroot) << L" TB)\n";
  std::wcerr << L"  projected post-eager-release peak = "
             << (peak - reclaimable) << L" B (" << TB(peak - reclaimable)
             << L" TB)\n";
  std::wcerr << L"  (non-root-only projection         = "
             << (peak - reclaimable_nonroot) << L" B ("
             << TB(peak - reclaimable_nonroot) << L" TB))\n";

  std::wcerr
      << L"\n  -- top DEAD-BUT-RETAINED tier-A composites (by size) --\n";
  for (std::size_t i = 0; i < dead_rows.size() && i < 20; ++i)
    std::wcerr << L"    vid=" << dead_rows[i].vid << L"  "
               << GB(dead_rows[i].bytes) << L" GB  last_access="
               << dead_rows[i].last << L" < peak_clock=" << peak_clock
               << (dead_rows[i].is_root ? L"  [ROOT]" : L"") << L"  "
               << dead_rows[i].label << L"\n";

  // The specific big composites the analysis flags.
  std::wcerr << L"\n  -- fate of the BIG composites --\n";
  for (std::size_t vid :
       {std::size_t(88), std::size_t(163), std::size_t(115), std::size_t(87)}) {
    if (vid >= rich.cells.size()) continue;
    std::size_t const h = cell_of_vid(vid).hash;
    std::size_t const la = last_access_of(h);
    bool const alive = alive_at_peak.count(h) != 0;
    std::size_t const bytes =
        peak_bytes_of.count(h) ? peak_bytes_of[h] : foot(cell_of_vid(vid));
    std::wstring const verdict =
        !alive ? L"NOT alive at peak"
               : (la < peak_clock ? L"RECLAIMABLE (dead-but-retained)"
                                  : L"PENDING at peak (read >= peak_clock)");
    std::wcerr << L"    vid=" << vid << L"  " << GB(bytes)
               << L" GB  last_access=" << la << L" vs peak_clock=" << peak_clock
               << L"  " << node_kind(h) << L"{"
               << space_sig(cell_of_vid(vid).carried) << L"}  => " << verdict
               << L"\n";
  }

  // (7) SANITY CHECK (corrected classification): the probe's "tier-B" set
  //     (values enumerated inside aux ScopeBlocks) MIXES two kinds:
  //       (a) genuine LoopLocal transients ([build]) -- use-counted, released
  //           at last use, NOT alive post-run; and
  //       (b) block ESCAPE outputs ([out:Sum]/[out:Scatter]) -- which the
  //           executor HOMES at the root cache (pinned) on block close, so
  //           they behave exactly like tier-A.
  //     Only (a) is use-counted, so only (a) must have last_access >=
  //     peak_clock when alive at peak. Discriminate with pinned_now (alive on
  //     root post- run): a pinned entry that looks dead-but-retained is
  //     EXPECTED (it is just more pinned floor); a NON-pinned entry that looks
  //     dead is the real missed-read-path signal.
  std::size_t tierB_alive = 0, tierB_pending = 0;
  std::size_t tierB_pinned_dead = 0, tierB_pinned_dead_bytes = 0;
  std::size_t tierB_transient_dead = 0, tierB_transient_dead_bytes = 0;
  for (auto const& r : tierB) {
    std::size_t const h = cell_of_vid(r.vid).hash;
    if (!alive_at_peak.count(h)) continue;
    ++tierB_alive;
    std::size_t const la = last_access_of(h);
    std::size_t const bytes =
        peak_bytes_of.count(h) ? peak_bytes_of[h] : r.bytes;
    if (la >= peak_clock) {
      ++tierB_pending;
    } else if (pinned_now.count(h)) {
      ++tierB_pinned_dead;  // homed escape output: dead-but-retained EXPECTED
      tierB_pinned_dead_bytes += bytes;
    } else {
      ++tierB_transient_dead;  // genuine transient looking dead: the RED FLAG
      tierB_transient_dead_bytes += bytes;
    }
  }
  std::wcerr << L"\n============ (7) SANITY CHECK (tier-B, corrected) ======\n";
  std::wcerr << L"  tier-B alive at peak = " << tierB_alive << L"\n";
  std::wcerr << L"    pending (last>=peak, use-counted still-needed) = "
             << tierB_pending << L"\n";
  std::wcerr << L"    PINNED homed-output dead-but-retained (EXPECTED) = "
             << tierB_pinned_dead << L" (" << TB(tierB_pinned_dead_bytes)
             << L" TB)\n";
  std::wcerr << L"    NON-pinned transient looking dead (RED FLAG if >0) = "
             << tierB_transient_dead << L" (" << TB(tierB_transient_dead_bytes)
             << L" TB)\n";
  std::wcerr << L"  => "
             << (tierB_transient_dead == 0
                     ? L"PASS: every use-counted transient alive at peak is "
                       L"still-pending; instrumentation catches all reads"
                     : L"SUSPECT: a use-counted transient looks dead -- missed "
                       L"read path")
             << L"\n";

  // (8) The FULL pinned peak: tier-A floor + homed escape outputs are BOTH
  //     pinned-until-reset and BOTH mostly dead-but-retained at the peak.
  //     Report the combined reclaimable so the caller sees the whole pinned
  //     overhang.
  std::size_t const pinned_dead_total = reclaimable + tierB_pinned_dead_bytes;
  std::wcerr << L"\n============ (8) TOTAL PINNED DEAD-BUT-RETAINED ========\n";
  std::wcerr << L"  tier-A floor reclaimable          = " << reclaimable
             << L" B (" << TB(reclaimable) << L" TB)\n";
  std::wcerr << L"  + homed escape-output reclaimable = "
             << tierB_pinned_dead_bytes << L" B ("
             << TB(tierB_pinned_dead_bytes) << L" TB)\n";
  std::wcerr << L"  = TOTAL pinned dead-but-retained  = " << pinned_dead_total
             << L" B (" << TB(pinned_dead_total) << L" TB), "
             << (peak ? 100.0 * pinned_dead_total / peak : 0.0)
             << L"% of the realized peak\n";
  std::wcerr << L"  residual live at peak instant     = "
             << (peak - pinned_dead_total) << L" B ("
             << TB(peak - pinned_dead_total) << L" TB)\n";
  std::wcerr << L"  NOTE: releasing at last-use MOVES the global peak; this "
                L"residual is a LOWER BOUND on the post-release peak, and the "
                L"reclaim an UPPER BOUND on the saving at THIS instant.\n";

  std::wcerr << L"====================================================\n";
  CHECK(mon.hwmark_bytes > 0);
}

// ===========================================================================
// A value MATERIALIZED across a forced loop split is built in its home block
// (its in-nest readers take the per-batch cell) and escapes from that same
// block outward (the other pass takes the assembled form). well_formed admits
// exactly that: every block listing the value in `outputs` either holds its
// BuildStep or is an ancestor of the block that does. Any other combination
// of build and escape sites is still duplicate production.
// ===========================================================================
TEST_CASE("well_formed accepts a value built and escaped in its home block",
          "[ordered][escape-chain]") {
  using sequant::eval::BuildStep;
  using sequant::eval::OrderedSchedule;
  using sequant::eval::OutputKind;
  using sequant::eval::ScopeBlock;
  using sequant::eval::Step;
  using sequant::eval::well_formed;

  // outer loop { inner loop { build 1; escape 1 } escape 1 }
  auto const make = [](bool build_in_inner) {
    ScopeBlock inner;
    inner.axis = sequant::Index{L"i_2"};
    inner.latitude_ordinal = 0;
    if (build_in_inner) inner.steps.push_back(Step{BuildStep{1}});
    inner.outputs.push_back({1, OutputKind::AccumulateScatter});

    ScopeBlock sibling;
    sibling.axis = sequant::Index{L"i_2"};
    sibling.latitude_ordinal = 1;  // distinct ordinal: sibling, not a chain
    if (!build_in_inner) sibling.steps.push_back(Step{BuildStep{1}});

    ScopeBlock outer;
    outer.axis = sequant::Index{L"i_1"};
    outer.latitude_ordinal = 0;
    outer.steps.push_back(Step{std::move(inner)});
    outer.steps.push_back(Step{std::move(sibling)});
    outer.outputs.push_back({1, OutputKind::AccumulateScatter});

    OrderedSchedule sched;
    sched.root.steps.push_back(Step{std::move(outer)});
    sched.num_values = 2;
    return sched;
  };

  // built where it escapes; the outer escape is an ancestor of that block
  CHECK(well_formed(make(/*build_in_inner=*/true)));
  // the same escape chain, but the BuildStep moved to a SIBLING block: the
  // inner escape is neither the home block nor an ancestor of it
  CHECK_FALSE(well_formed(make(/*build_in_inner=*/false)));
}

// A plain (non-template) `requires(...) { c.f(); }` at namespace scope has no
// substitution to be SFINAE-safe over, so an absent member is a hard compile
// error there, not a `false` result -- exactly the failure this test wants to
// witness (a name lookup failure), so a template wrapper is needed to make
// the check itself well-formed.
namespace {
template <typename Cache>
concept has_loop_colored_slice_seam =
    requires(Cache& c) { c.loop_colored_slice_seam(); };
}  // namespace

TEST_CASE("ordered executor has no slice seam", "[ordered][cell_table]") {
  static_assert(!has_loop_colored_slice_seam<sequant::CacheManager<ScalarNode>>,
                "the loop-colored slice seam must be gone from the cache");
  SUCCEED();
}

// Task 5 (Stage 3): the value-keyed cache machinery (home slots, coloring,
// the release bridge) that only the ordered path used is gone now that the
// executor runs entirely on cells. Same SFINAE-friendly-lookup rationale as
// the slice-seam concept above: a plain out-of-line `requires` expression
// would be a hard compile error on a genuinely absent member, not a `false`
// result, so each check is wrapped in its own concept.
namespace {
template <typename Cache>
concept has_ensure_home_slot_1 =
    requires(Cache& c, typename Cache::cache_key_type const& k) {
      c.ensure_home_slot(k);
    };
template <typename Cache>
concept has_ensure_home_slot_3 =
    requires(Cache& c, typename Cache::cache_key_type const& k) {
      c.ensure_home_slot(k, std::size_t{1}, false);
    };
template <typename Cache>
concept has_release_at = requires(
    Cache& c, typename Cache::cache_key_type const& k) { c.release_at(k); };
template <typename Cache>
concept has_recolor = requires(
    Cache& c, typename Cache::cache_key_type const& k) { c.recolor(k); };
template <typename Schedule>
concept has_home_mode_depth = requires(Schedule& s) { s.home_mode_depth; };
// Stage 4 (explicit-cells): the residency/role inference machinery that no
// production path consults any more (verified zero-callers). Same
// SFINAE-friendly-lookup rationale as above.
template <typename Cache>
concept has_peek_at = requires(
    Cache& c, typename Cache::cache_key_type const& k) { c.peek_at(k); };
template <typename Cache>
concept has_entry_is_persistent =
    requires(Cache const& c, typename Cache::cache_key_type const& k) {
      c.entry_is_persistent(k);
    };
template <typename Cache>
concept has_dump_entries_for_hash =
    requires(Cache const& c) { c.dump_entries_for_hash(std::size_t{0}); };
template <typename Cache>
concept has_stored_this_eval =
    requires(Cache const& c, typename Cache::cache_key_type const& k) {
      c.stored_this_eval(k);
    };
template <typename M2L>
concept has_mode_of =
    requires(M2L const& m, sequant::DagScopeLevel const& l) { m.mode_of(l); };
// Task 2 (explicit cells, cache key): the value-id slice coloring is gone --
// CachedValue holds only the node.
template <typename CV>
concept has_coloring_member = requires(CV const& v) { v.coloring; };
}  // namespace

TEST_CASE(
    "cache manager and ordered schedule have no legacy value-keyed cache "
    "machinery",
    "[ordered][cell_table]") {
  static_assert(!has_ensure_home_slot_1<sequant::CacheManager<ScalarNode>>,
                "ensure_home_slot(key) must be gone from the cache");
  static_assert(
      !has_ensure_home_slot_3<sequant::CacheManager<ScalarNode>>,
      "ensure_home_slot(key, use_count, persistent) must be gone from the "
      "cache");
  static_assert(!has_release_at<sequant::CacheManager<ScalarNode>>,
                "release_at must be gone from the cache");
  static_assert(!has_recolor<sequant::CacheManager<ScalarNode>>,
                "recolor must be gone from the cache");
  static_assert(!has_home_mode_depth<OrderedSchedule>,
                "home_mode_depth must be gone from the ordered schedule");
  static_assert(!has_peek_at<sequant::CacheManager<ScalarNode>>,
                "peek_at must be gone from the cache");
  static_assert(!has_entry_is_persistent<sequant::CacheManager<ScalarNode>>,
                "entry_is_persistent must be gone from the cache");
  static_assert(!has_dump_entries_for_hash<sequant::CacheManager<ScalarNode>>,
                "dump_entries_for_hash must be gone from the cache");
  static_assert(!has_stored_this_eval<sequant::CacheManager<ScalarNode>>,
                "stored_this_eval(key) must be gone from the cache");
  static_assert(!has_mode_of<sequant::ModeToLevel>,
                "ModeToLevel::mode_of must be gone");
  static_assert(
      !has_coloring_member<sequant::eval::CachedValue<ScalarNode>>,
      "CachedValue must hold only the node -- the value-id slice coloring "
      "must be gone");
  SUCCEED();
}

// ===========================================================================
// Explicit value cells, final round (C1): a LEAF's FIRST touch must be served
// through its Read, sliced. The read resolver defers a leaf whose cell has no
// result yet (nothing has recorded it), the leaf evaluator runs, and the leaf
// is recorded -- but the value handed to the consumer used to be the WHOLE
// leaf: the only slicing left on that path is `slice_to_use`, whose ordered
// arm is a no-op (it fires on an `exact_axis`, which no table-driven fetch
// sets). So the first consumer of a leaf inside a batch loop got a whole
// operand where the schedule declares a batch slice, while every LATER
// consumer of the same leaf got the declared slice -- a silent
// whole-against-sliced pairing.
//
// Driven at the read path directly (a hand-built table + registry + resolver
// on a cache carrying one batch-loop context) rather than through a whole
// derived schedule: the defect is entirely in which value the leaf branch
// finalizes, and this pins it without depending on a schedule shape that
// happens to touch some leaf first inside a loop.
// ===========================================================================
TEST_CASE(
    "ordered executor: a leaf's first touch inside a batch loop is served "
    "through its Read, with the declared slice",
    "[ordered][cell_registry]") {
  using sequant::eval::dryrun::EvalExprDryRun;
  using sequant::eval::dryrun::EvalNodeDryRun;

  auto ctx0 = sequant::get_default_context().clone();
  ctx0.set_first_dummy_index_ordinal(1000000);
  auto ctx_resetter = sequant::set_scoped_default_context(std::move(ctx0));

  sequant::eval::dryrun::SizeRegime regime;
  regime.space_extent = {{L"i", 8}, {L"a", 4}};
  auto const cm =
      std::make_shared<sequant::eval::dryrun::CostModel const>(regime);

  // C = X(i_1) * P(a_1): X is the leaf whose Read declares a slice on the
  // enclosing loop; P is deliberately NOT a value of the table (a transient,
  // which the resolver reports as such and the caller evaluates in place).
  auto const mk = [](std::wstring const& label, sequant::Index const& ix) {
    return sequant::ex<sequant::Tensor>(
        label, sequant::bra(sequant::container::svector<sequant::Index>{ix}),
        sequant::ket{}, sequant::Symmetry::Nonsymm,
        sequant::BraKetSymmetry::Symm, sequant::ColumnSymmetry::Nonsymm);
  };
  auto const prod = sequant::ex<sequant::Product>(
      1,
      sequant::ExprPtrList{mk(L"X", sequant::Index{L"i_1"}),
                           mk(L"P", sequant::Index{L"a_1"})},
      sequant::Product::Flatten::No);
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
  auto const node = sequant::binarize<EvalExprDryRun>(prod);
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
  REQUIRE_FALSE(node.leaf());
  // Whichever leg carries i_1 is the table's leaf value; the other is the
  // transient.
  auto const carries_i = [](EvalNodeDryRun const& n) {
    for (auto const& ix : n->canon_indices())
      if (ix.label() == L"i_1") return true;
    return false;
  };
  REQUIRE(carries_i(node.left()) != carries_i(node.right()));
  std::size_t const x_hash = carries_i(node.left())
                                 ? node.left()->hash_value()
                                 : node.right()->hash_value();

  // The table: cell 0 = the Leaf (value 0); cell 1 = the consumer's Build
  // (value 1) at the loop's scope; one Read of value 0 declaring the slice of
  // position 0 by loop instance (1,0).
  sequant::eval::CellTable table;
  {
    sequant::eval::TableCell leaf;
    leaf.value_id = 0;
    leaf.production.kind = sequant::eval::ProductionKind::Leaf;
    leaf.life = 1;
    table.cells.push_back(leaf);
    sequant::eval::TableCell build;
    build.value_id = 1;
    build.production.kind = sequant::eval::ProductionKind::Build;
    build.scope.path = {{sequant::eval::LoopKey{1, 0}, 0}};
    build.sliced = {{0, sequant::eval::LoopKey{1, 0}}};
    table.cells.push_back(build);
    table.reads.push_back(
        sequant::eval::Read{1, 0, 0, {{0, sequant::eval::LoopKey{1, 0}}}, {}});
  }
  sequant::eval::CellRegistry registry(table);
  sequant::eval::CellReadResolver resolver(
      registry, [x_hash](std::size_t h) -> std::optional<std::size_t> {
        if (h == x_hash) return std::size_t{0};
        return std::nullopt;  // the partner leg is a transient
      });
  resolver.begin_consumer(1);

  auto cache = sequant::CacheManager<EvalNodeDryRun>::empty();
  auto aops = sequant::eval::dryrun::make_dryrun_array_ops(cm);
  cache.set_array_ops(&aops);
  sequant::eval::BatchContext bctx;
  bctx.push_back({sequant::Index{L"i_1"},
                  sequant::eval::DagScopeLevel{1, L"i", 0, 0, 0},
                  {2, 4},
                  std::nullopt});
  cache.set_batch_context(bctx);
  cache.set_cell_read_resolver(&resolver);

  sequant::eval::dryrun::DryRunLeafEvaluator const yield{cm};
  sequant::ResultPtr const got =
      sequant::evaluate_impl<sequant::Trace::Off>(node, yield, cache);
  REQUIRE(got);

  // The leaf reached the contraction SLICED to [2,4): the product's own i_1
  // mode carries the batch extent and the absolute lower bound. Served whole
  // (the pre-fix behavior) the result has no override on that position at
  // all, and lobound 0.
  auto const idx = sequant::eval::dryrun::detail::indices_of(*got);
  auto const ov = sequant::eval::dryrun::detail::overrides_of(*got);
  auto const lob = sequant::eval::dryrun::detail::lobounds_of(*got);
  std::optional<std::size_t> i_pos;
  for (std::size_t p = 0; p < idx.size(); ++p)
    if (idx[p].label() == L"i_1") i_pos = p;
  REQUIRE(i_pos.has_value());
  REQUIRE(ov.count(*i_pos) == 1);
  CHECK(ov.at(*i_pos) == 2);
  REQUIRE(lob.count(*i_pos) == 1);
  CHECK(lob.at(*i_pos) == 2);

  // The Read was consumed exactly once by the two fetches the leg makes (the
  // deferring first touch left it unconsumed; the second served it), so the
  // consumer has no read of that value left.
  CHECK_THROWS(resolver.fetch(x_hash, bctx));
}
// ===========================================================================
// Explicit value cells, final round (I2): the input-mirrored configuration is
// EXECUTED, not merely validated. The static gate above ("cell table: the
// input-mirrored configuration derives a valid table") proves the table is
// well formed under that configuration; only a run proves the executor can
// walk it -- the level-skipping escape chain and the materialized member the
// mirrored schedule contains are runtime shapes the default configuration
// never produces. Same construction as the static case (peak threshold 25e9,
// the batched time-then-space objective, set directly rather than through the
// environment), then the strict dry-run walk the default fixture uses:
// range/lobound checks in the backend plus the cache-fill-once tripwire.
// The [w20-auxocc-walk] fixture keeps its environment-driven mirroring
// unchanged; this case makes the mirrored RUN part of the default suite.
// ===========================================================================
TEST_CASE(
    "ordered executor: the input-mirrored configuration RUNS the strict "
    "dry-run walk to completion",
    "[ordered][cell_table]") {
  using sequant::eval::dryrun::EvalExprDryRun;
  using sequant::eval::dryrun::EvalNodeDryRun;
  using Node = EvalNodeDryRun;
  auto ctx = sequant::get_default_context().clone();
  ctx.set_first_dummy_index_ordinal(1000000);
  auto isr = ctx.mutable_index_space_registry();
  REQUIRE(isr != nullptr);
  sequant::mbpt::add_pao_spaces(isr, sequant::mbpt::Spin::any);
  sequant::mbpt::add_df_spaces(isr);
  auto ctx_resetter = sequant::set_scoped_default_context(std::move(ctx));

  auto const body =
      orderedexec_witness_slurp(std::string(SEQUANT_UNIT_TESTS_SOURCE_DIR) +
                                "/data/csv_ccsd_doubles_residual_df.txt");
  REQUIRE(!body.empty());
  std::string line = body;
  if (auto nl = line.find('\n'); nl != std::string::npos)
    line = line.substr(0, nl);
  auto expr = sequant::deserialize<sequant::ExprPtr>(line);
  REQUIRE(expr->is<sequant::Sum>());
  auto const& summands = expr->as<sequant::Sum>().summands();
  auto regime = orderedexec_witness_df_regime(kOrderedExecWater20_pVDZF12);
  auto cm = std::make_shared<sequant::eval::dryrun::CostModel const>(regime);
  sequant::BatchPolicy policy;
  policy.is_batchable_contracted_index = [](sequant::Index const& ix) {
    return ix.space().base_key() == L"Κ";
  };
  policy.is_batchable_external_index = [](sequant::Index const& ix) {
    auto const reg = sequant::get_default_context().index_space_registry();
    return reg && ix.space() && reg->is_pure_occupied(ix.space());
  };
  policy.batch_spectator_indices = true;
  policy.node_level_placement = true;
  policy.batch_target_size = [](sequant::Index const& ix) -> std::size_t {
    return ix.space().base_key() == L"Κ" ? 256 : 16;
  };
  policy.is_volatile_leaf = [](sequant::Tensor const& t) {
    return t.label() == L"t";
  };
  policy.accumulation_factor = 1.0;
  policy.persistent_only = false;
  policy.peak_threshold = 25e9;  // THE mirrored setting
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
  std::vector<Node> forest;
  for (auto const& s : summands) {
    sequant::ExprPtr const term = orderedexec_witness_flatten_product(s);
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
  auto const block_of = [](sequant::Index const& ix) -> std::size_t {
    return ix.space().base_key() == L"Κ" ? 256 : 16;
  };
  auto const rich = sequant::eval::compute_dag_boulevard(forest, *cm, block_of);
  auto const legality = sequant::eval::analyze_legality(rich, forest, policy);
  auto const ordered = sequant::eval::build_ordered_schedule(
      rich, legality, policy, std::initializer_list<std::wstring>{});
  REQUIRE(sequant::eval::well_formed(ordered));

  using annot_t = std::remove_cvref_t<decltype(forest.front()->annot())>;
  annot_t const layout{};
  sequant::eval::dryrun::DryRunLeafEvaluator const yield{cm};
  std::function<std::size_t(sequant::Index const&)> const target =
      [](sequant::Index const& ix) -> std::size_t {
    return ix.space().base_key() == L"Κ" ? 256 : 16;
  };
  std::function<bool(Node const&)> const is_volatile_node =
      [p = policy.is_volatile_leaf](Node const& n) -> bool {
    if (!n.leaf() || !n->is_tensor()) return false;
    return p && p(n->as_tensor());
  };

  auto& logger = sequant::Logger::instance();
  auto const prev_level = logger.eval.level;
  logger.eval.level = 0;
  auto aops = sequant::eval::dryrun::make_dryrun_array_ops(cm);
  auto ordered_cache = sequant::cache_manager(forest);
  ordered_cache.set_array_ops(&aops);

  // The strict tripwires the default walk runs under, set and RESTORED (this
  // case does not own the process environment).
  char const* const prev_strict = std::getenv("SEQUANT_UT_STRICT_FILL_ONCE");
  std::string const prev_strict_val = prev_strict ? prev_strict : "";
  setenv("SEQUANT_UT_STRICT_FILL_ONCE", "1", 1);
  REQUIRE_NOTHROW(sequant::eval::evaluate_ordered_schedule<sequant::Trace::Off>(
      forest, ordered, rich, layout, yield, ordered_cache, target, {},
      is_volatile_node));
  if (prev_strict)
    setenv("SEQUANT_UT_STRICT_FILL_ONCE", prev_strict_val.c_str(), 1);
  else
    unsetenv("SEQUANT_UT_STRICT_FILL_ONCE");

  logger.eval.level = prev_level;
}
// ===========================================================================
// An Assemble step's destination is sized from ITS OWN CELL'S FORM (the
// table), never inferred from where the escaped axis happens to sit on the
// node: the descriptor is the value's own canonical index list, narrowed by
// exactly the (position, loop instance) pairs the Assemble cell declares in
// its `sliced` -- which, at the root scope, is none at all, so the
// destination is the value's FULL extent with lobound 0 on every mode. The
// deleted alternative walked the enclosing batch context and narrowed by
// whatever positions matched a loop's space and fusion slot, which is a
// different (inferred) answer whenever the two disagree.
//
// The fixture is the aux+occ water-20 residual: aux (Κ) is batchable-
// contracted and occ is batchable-external, so the schedule nests loops of
// both kinds and a value reduced on an inner loop is carried (scattered) on
// an outer one -- the two-level shape this rule is about.
// ===========================================================================
TEST_CASE(
    "ordered executor: an Assemble step's scatter destination is sized from "
    "its own cell form",
    "[ordered][assemble-dest]") {
  using sequant::eval::dryrun::EvalExprDryRun;
  using sequant::eval::dryrun::EvalNodeDryRun;
  using Node = EvalNodeDryRun;

  auto ctx = sequant::get_default_context().clone();
  ctx.set_first_dummy_index_ordinal(1000000);
  auto isr = ctx.mutable_index_space_registry();
  REQUIRE(isr != nullptr);
  sequant::mbpt::add_pao_spaces(isr, sequant::mbpt::Spin::any);
  sequant::mbpt::add_df_spaces(isr);
  auto ctx_resetter = sequant::set_scoped_default_context(std::move(ctx));

  auto const body =
      orderedexec_witness_slurp(std::string(SEQUANT_UNIT_TESTS_SOURCE_DIR) +
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

  auto regime = orderedexec_witness_df_regime(kOrderedExecWater20_pVDZF12);
  auto cm = std::make_shared<sequant::eval::dryrun::CostModel const>(regime);

  sequant::BatchPolicy policy;
  policy.is_batchable_contracted_index = [](sequant::Index const& ix) {
    return ix.space().base_key() == L"Κ";
  };
  policy.is_batchable_external_index = [](sequant::Index const& ix) {
    auto const reg = sequant::get_default_context().index_space_registry();
    return reg && ix.space() && reg->is_pure_occupied(ix.space());
  };
  policy.batch_spectator_indices = true;
  policy.node_level_placement = true;
  policy.batch_target_size = [](sequant::Index const& ix) -> std::size_t {
    return ix.space().base_key() == L"Κ" ? 256 : 16;
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

  std::vector<Node> forest;
  for (auto const& s : summands) {
    sequant::ExprPtr const term = orderedexec_witness_flatten_product(s);
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

  auto const block_of = [](sequant::Index const& ix) -> std::size_t {
    return ix.space().base_key() == L"Κ" ? 256 : 16;
  };
  auto const rich = sequant::eval::compute_dag_boulevard(forest, *cm, block_of);
  auto const legality = sequant::eval::analyze_legality(rich, forest, policy);
  auto const ordered = sequant::eval::build_ordered_schedule(
      rich, legality, policy, std::initializer_list<std::wstring>{});
  REQUIRE(sequant::eval::well_formed(ordered));

  using annot_t = std::remove_cvref_t<decltype(forest.front()->annot())>;
  annot_t const layout{};
  sequant::eval::dryrun::DryRunLeafEvaluator const yield{cm};
  std::function<std::size_t(sequant::Index const&)> const target =
      [](sequant::Index const& ix) -> std::size_t {
    return ix.space().base_key() == L"Κ" ? 256 : 16;
  };
  std::function<bool(Node const&)> const is_volatile_node =
      [p = policy.is_volatile_leaf](Node const& n) -> bool {
    if (!n.leaf() || !n->is_tensor()) return false;
    return p && p(n->as_tensor());
  };

  auto& logger = sequant::Logger::instance();
  auto const prev_level = logger.eval.level;
  logger.eval.level = 0;
  auto aops = sequant::eval::dryrun::make_dryrun_array_ops(cm);

  // ---- The table the executor will run on, derived exactly as it derives
  // it, so the structural claims below are about the very same cells.
  auto const vmap = sequant::eval::build_value_node_map(forest);
  auto const sma = sequant::eval::compute_sliced_mode_assignment(ordered, rich);
  sequant::eval::CellTableInputs in;
  in.ordered = &ordered;
  in.rich = &rich;
  in.sliced = &sma;
  in.sliced_modes_of = [&](std::size_t vid) {
    auto const it = vmap.find(sequant::eval::value_key_of(rich.cells[vid]));
    REQUIRE(it != vmap.end());
    return sequant::eval::detail::home_modes_in_cell_frame(rich, vid,
                                                           it->second);
  };
  in.volatile_of = [&](std::size_t vid) {
    auto const it = vmap.find(sequant::eval::value_key_of(rich.cells[vid]));
    return it != vmap.end() &&
           sequant::subtree_any(it->second, is_volatile_node);
  };
  in.n_batches_of =
      sequant::eval::detail::ordered_n_batches_by_loop(ordered, target, &aops);
  in.operands_of = orderedexec_per_leg_operands(rich, vmap);
  auto const table = sequant::eval::build_cell_table(in);
  REQUIRE(
      sequant::eval::validate_cell_table(table, ordered.root, in.n_batches_of)
          .empty());

  // (1) The two-level shape is present: a Scatter Assemble at the ROOT scope
  // whose per-batch source is itself a partial (reduced on an inner loop, or
  // assembled one level in), and at least one Assemble living INSIDE a loop
  // (whose destination therefore IS narrowed, by its own declared `sliced`).
  std::size_t n_root_scatter = 0, n_nested_assemble = 0, n_two_level = 0;
  for (sequant::eval::CellId c = 0; c < table.cells.size(); ++c) {
    auto const& a = table.cells[c];
    if (a.production.kind != sequant::eval::ProductionKind::Assemble) continue;
    if (!a.scope.path.empty()) ++n_nested_assemble;
    if (a.production.assemble != sequant::eval::AssembleKind::Scatter) continue;
    if (!a.scope.path.empty()) continue;
    ++n_root_scatter;
    auto const& src = table.cells[a.production.source];
    if (!src.partial_over.empty() ||
        src.production.kind == sequant::eval::ProductionKind::Assemble)
      ++n_two_level;
  }
  REQUIRE(n_root_scatter > 0);
  REQUIRE(n_nested_assemble > 0);
  REQUIRE(n_two_level > 0);

  // (1b) The NARROWED case is present too, and the narrowing is a real one:
  // a Scatter Assemble whose own `sliced` is non-empty is a destination that
  // lives INSIDE an enclosing loop, so it covers only that loop's current
  // batch of the sliced position. Assert the executor's own two steps --
  // make_zeros on the value's full index list, then slice_mode(pos, range) --
  // land exactly on that batch: extent = the batch's width, lobound = its
  // start. (Reproduced here on the very descriptor and range the executor
  // uses; the destination object itself is internal to the run.)
  {
    std::size_t n_narrowed = 0;
    for (sequant::eval::CellId c = 0; c < table.cells.size(); ++c) {
      auto const& a = table.cells[c];
      if (a.production.kind != sequant::eval::ProductionKind::Assemble ||
          a.production.assemble != sequant::eval::AssembleKind::Scatter ||
          a.sliced.empty())
        continue;
      auto const nit =
          vmap.find(sequant::eval::value_key_of(rich.cells[a.value_id]));
      if (nit == vmap.end()) continue;
      auto const [pos, key] = a.sliced.front();
      // The enclosing loop instance the position is narrowed by, and its own
      // batch partition (the same one the executor iterates).
      std::optional<sequant::Index> axis;
      std::function<void(sequant::eval::ScopeBlock const&)> find =
          [&](sequant::eval::ScopeBlock const& b) {
            if (!axis && sequant::eval::detail::same_key(b.level.key(), key))
              axis = b.axis;
            for (auto const& st : b.steps)
              if (auto const* ch =
                      std::get_if<sequant::eval::ScopeBlock>(&st.value))
                find(*ch);
          };
      find(ordered.root);
      REQUIRE(axis.has_value());  // the table names a loop the schedule has
      auto const batches = aops.axis_batches(*axis, target(*axis));
      REQUIRE(batches.size() > 1);  // a genuine narrowing, not the whole axis
      auto const [lo, hi] = batches.front();
      auto const full = aops.make_zeros(nit->second->canon_indices());
      auto const narrowed = full->slice_mode(pos, lo, hi);
      auto const ov = sequant::eval::dryrun::detail::overrides_of(*narrowed);
      auto const lb = sequant::eval::dryrun::detail::lobounds_of(*narrowed);
      REQUIRE(ov.count(pos) == 1);
      CHECK(ov.at(pos) == hi - lo);
      REQUIRE(lb.count(pos) == 1);
      CHECK(lb.at(pos) == lo);
      ++n_narrowed;
    }
    INFO("Scatter Assembles with a declared narrowing: " << n_narrowed);
    REQUIRE(n_narrowed > 0);
  }

  // (2) The sizing rule's own precondition, stated on the table: every
  // position an Assemble cell declares sliced names a loop instance that is
  // OPEN at that cell's scope (so the executor can bind it to a batch range),
  // and a root-scope Assemble declares none -- its destination is the value's
  // full extent.
  for (sequant::eval::CellId c = 0; c < table.cells.size(); ++c) {
    auto const& a = table.cells[c];
    if (a.production.kind != sequant::eval::ProductionKind::Assemble) continue;
    for (auto const& [pos, key] : a.sliced) {
      bool on_path = false;
      for (auto const& [pk, lat] : a.scope.path)
        if (sequant::eval::detail::same_key(pk, key)) on_path = true;
      INFO("cell#" << c << " (value " << a.value_id << ") slices position "
                   << pos << " on a loop instance not open at its scope");
      CHECK(on_path);
    }
    if (a.scope.path.empty()) {
      INFO("root-scope Assemble cell#" << c << " (value " << a.value_id
                                       << ") declares a narrowed form");
      CHECK(a.sliced.empty());
    }
  }

  // (3) Every destination the run allocates is created at the value's OWN
  // full index list, with no extent override and no lobound -- the narrowing,
  // where the table declares one, is applied on top of that afterwards and is
  // never baked into what make_zeros is asked for.
  // Snapshotted AT CREATION: the destination object itself is then narrowed
  // (where the table declares a narrowing) and scattered into, and
  // write_into_slice records that coverage ON it, so its state after the run
  // says nothing about how it was SIZED.
  struct MadeZeros {
    sequant::container::vector<sequant::Index> descriptor;
    sequant::eval::dryrun::ExtentOverrides overrides, lobounds;
  };
  std::vector<MadeZeros> zeros;
  auto const make_zeros_inner = aops.make_zeros;
  aops.make_zeros =
      [&](sequant::container::vector<sequant::Index> const& d) -> ResultPtr {
    auto r = make_zeros_inner(d);
    zeros.push_back(MadeZeros{d,
                              sequant::eval::dryrun::detail::overrides_of(*r),
                              sequant::eval::dryrun::detail::lobounds_of(*r)});
    return r;
  };

  auto ordered_cache = sequant::cache_manager(forest);
  ordered_cache.set_array_ops(&aops);
  char const* const prev_strict = std::getenv("SEQUANT_UT_STRICT_FILL_ONCE");
  std::string const prev_strict_val = prev_strict ? prev_strict : "";
  setenv("SEQUANT_UT_STRICT_FILL_ONCE", "1", 1);
  REQUIRE_NOTHROW(sequant::eval::evaluate_ordered_schedule<sequant::Trace::Off>(
      forest, ordered, rich, layout, yield, ordered_cache, target, {},
      is_volatile_node));
  if (prev_strict)
    setenv("SEQUANT_UT_STRICT_FILL_ONCE", prev_strict_val.c_str(), 1);
  else
    unsetenv("SEQUANT_UT_STRICT_FILL_ONCE");
  logger.eval.level = prev_level;

  REQUIRE(!zeros.empty());  // the run really did assemble by scattering
  for (auto const& z : zeros) {
    // The descriptor is some escaping value's own canonical index list.
    bool matches_a_value = false;
    for (auto const& a : table.cells) {
      if (a.production.kind != sequant::eval::ProductionKind::Assemble ||
          a.production.assemble != sequant::eval::AssembleKind::Scatter)
        continue;
      auto const it =
          vmap.find(sequant::eval::value_key_of(rich.cells[a.value_id]));
      if (it == vmap.end()) continue;
      auto const& ci = it->second->canon_indices();
      if (sequant::container::vector<sequant::Index>(ci.begin(), ci.end()) ==
          z.descriptor)
        matches_a_value = true;
    }
    CHECK(matches_a_value);
    // Full extent, lobound 0: the destination is created whole, on the
    // value's own index list -- no extent override and no lobound anywhere.
    CHECK(z.overrides.empty());
    CHECK(z.lobounds.empty());
  }
}

// ===========================================================================
// A loop-invariant escape is not re-formed on later batches.
//
// The case: a block whose escape is loop-invariant to the loop its Assemble
// cell sits in -- the table marks that cell `produce_if_absent`, so it is
// assembled on the first visit and REUSED afterwards; the steps that feed it
// are then dead on every later visit. What this pins is the fill-once
// property this reuse depends on, not the whole-block skip machinery that
// (on THIS fixture) happens to implement it: with strictness enabled on the
// cache handle (`set_strict_fill_once`, the per-manager knob -- NOT the
// `SEQUANT_UT_STRICT_FILL_ONCE` environment variable, which is read only
// once, on the first `eval::strict_fill_once()` call anywhere in the
// process, so a later setenv in this translation unit has no effect on it),
// `CellRegistry::set` throws if a non-persistent cell is filled a second
// time while its prior fill is still live (see cell_registry.hpp's
// `CellRegistryHooks::strict_fill_once`, which
// `run_ordered_schedule_pre_results` now reads off the cache handle), so
// REQUIRE_NOTHROW over the whole evaluation is a direct proof that the
// invariant escape's cell is filled EXACTLY ONCE across the loop's batches, not
// once per batch. (An earlier version of this case instead checked
// `ordered_last_block_skips() > 0` -- the vestigial signal of the OLD builder's
// separate consumer-pass block, which the per-nest forced-split design no
// longer emits; the loop-invariant escape this case is actually about is
// rebuilt on every batch in both old and new code and needs the fill-once
// property above to be pinned at all.)
// ===========================================================================
TEST_CASE(
    "ordered executor: a loop-invariant escape is not re-formed on later "
    "batches",
    "[ordered][block-skip]") {
  using sequant::eval::dryrun::EvalExprDryRun;
  using sequant::eval::dryrun::EvalNodeDryRun;
  using Node = EvalNodeDryRun;

  auto ctx = sequant::get_default_context().clone();
  ctx.set_first_dummy_index_ordinal(1000000);
  auto isr = ctx.mutable_index_space_registry();
  REQUIRE(isr != nullptr);
  sequant::mbpt::add_pao_spaces(isr, sequant::mbpt::Spin::any);
  sequant::mbpt::add_df_spaces(isr);
  auto ctx_resetter = sequant::set_scoped_default_context(std::move(ctx));

  auto const body =
      orderedexec_witness_slurp(std::string(SEQUANT_UNIT_TESTS_SOURCE_DIR) +
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

  auto regime = orderedexec_witness_df_regime(kOrderedExecWater20_pVDZF12);
  auto cm = std::make_shared<sequant::eval::dryrun::CostModel const>(regime);

  // aux+occ: Κ batchable-contracted, occ batchable-external -- nested loops,
  // which is what makes a block-inside-a-block (and so a loop-invariant
  // escape) possible at all.
  sequant::BatchPolicy policy;
  policy.is_batchable_contracted_index = [](sequant::Index const& ix) {
    return ix.space().base_key() == L"Κ";
  };
  policy.is_batchable_external_index = [](sequant::Index const& ix) {
    auto const reg = sequant::get_default_context().index_space_registry();
    return reg && ix.space() && reg->is_pure_occupied(ix.space());
  };
  policy.batch_spectator_indices = true;
  policy.node_level_placement = true;
  policy.batch_target_size = [](sequant::Index const& ix) -> std::size_t {
    return ix.space().base_key() == L"Κ" ? 256 : 16;
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

  std::vector<Node> forest;
  for (auto const& s : summands) {
    sequant::ExprPtr const term = orderedexec_witness_flatten_product(s);
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

  auto const block_of = [](sequant::Index const& ix) -> std::size_t {
    return ix.space().base_key() == L"Κ" ? 256 : 16;
  };
  auto const rich = sequant::eval::compute_dag_boulevard(forest, *cm, block_of);
  auto const legality = sequant::eval::analyze_legality(rich, forest, policy);
  auto const ordered = sequant::eval::build_ordered_schedule(
      rich, legality, policy, std::initializer_list<std::wstring>{});
  REQUIRE(sequant::eval::well_formed(ordered));

  using annot_t = std::remove_cvref_t<decltype(forest.front()->annot())>;
  annot_t const layout{};
  sequant::eval::dryrun::DryRunLeafEvaluator const yield{cm};
  std::function<std::size_t(sequant::Index const&)> const target =
      [](sequant::Index const& ix) -> std::size_t {
    return ix.space().base_key() == L"Κ" ? 256 : 16;
  };
  std::function<bool(Node const&)> const is_volatile_node =
      [p = policy.is_volatile_leaf](Node const& n) -> bool {
    if (!n.leaf() || !n->is_tensor()) return false;
    return p && p(n->as_tensor());
  };

  auto const vmap = sequant::eval::build_value_node_map(forest);
  auto aops = sequant::eval::dryrun::make_dryrun_array_ops(cm);

  // Precondition: the table really does have a loop-invariant escape -- an
  // Assemble cell INSIDE a loop (non-root scope) that the builder marked
  // `produce_if_absent`, i.e. one whose value does not vary with the loop its
  // own scope sits in.
  {
    auto const sma =
        sequant::eval::compute_sliced_mode_assignment(ordered, rich);
    sequant::eval::CellTableInputs in;
    in.ordered = &ordered;
    in.rich = &rich;
    in.sliced = &sma;
    in.sliced_modes_of = [&](std::size_t vid) {
      auto const it = vmap.find(sequant::eval::value_key_of(rich.cells[vid]));
      REQUIRE(it != vmap.end());
      return sequant::eval::detail::home_modes_in_cell_frame(rich, vid,
                                                             it->second);
    };
    in.volatile_of = [&](std::size_t vid) {
      auto const it = vmap.find(sequant::eval::value_key_of(rich.cells[vid]));
      return it != vmap.end() &&
             sequant::subtree_any(it->second, is_volatile_node);
    };
    in.n_batches_of = sequant::eval::detail::ordered_n_batches_by_loop(
        ordered, target, &aops);
    in.operands_of = orderedexec_per_leg_operands(rich, vmap);
    auto const table = sequant::eval::build_cell_table(in);
    std::size_t n_invariant_escapes = 0;
    for (auto const& c : table.cells)
      if (c.production.kind == sequant::eval::ProductionKind::Assemble &&
          !c.scope.path.empty() && c.produce_if_absent)
        ++n_invariant_escapes;
    INFO("loop-invariant escapes in the table: " << n_invariant_escapes);
    REQUIRE(n_invariant_escapes > 0);
  }

  auto ordered_cache = sequant::cache_manager(forest);
  ordered_cache.set_array_ops(&aops);
  // Deterministic per-instance override (cache_manager.hpp's
  // set_strict_fill_once), not the environment: SEQUANT_UT_STRICT_FILL_ONCE
  // is read once, via a function-local static latched on the first
  // eval::strict_fill_once() call anywhere in the process (already long
  // past by this point in the file), so a setenv here would have no effect
  // on it. This override is robust to that latch and to the assert
  // configuration alike.
  ordered_cache.set_strict_fill_once(true);
  // THE pinned property: the whole evaluation runs to completion under
  // strict fill-once, so the invariant escape's cell (and every other cell)
  // is filled EXACTLY ONCE while live -- a re-formed escape on a later batch
  // would throw here. (The recompute tally -- `set_recompute_tally_enabled`
  // -- was tried as a more targeted signal and dropped: MEASURED on this
  // fixture, the block that would be skipped whole holds no value whose
  // ONLY production site is inside it -- every one of them is also built by
  // a step of another block or another pass -- so no per-value build count
  // in the tally moves when the whole-block skip fires or not; it cannot
  // express this case's claim on this fixture, and REQUIRE_NOTHROW is the
  // pinned property instead.)
  REQUIRE_NOTHROW(sequant::eval::evaluate_ordered_schedule<sequant::Trace::Off>(
      forest, ordered, rich, layout, yield, ordered_cache, target, {},
      is_volatile_node));
}

// ===========================================================================
// A `produce_if_absent` cell that is BOUND to an enclosing loop instance is
// re-produced on every batch of that loop.
//
// Such a cell is invariant to the loop it is homed in (that is what the flag
// says) but NOT to an outer one, so the outer loop's per-batch reset
// (CellRegistry::clear_bound_to) empties it, and the next visit has to build
// it again. The hazard is the per-visit skip set: it is computed once at a
// block's entry from the cells the registry then holds and is consulted
// across every batch of that block and inside every nested block, so a cell
// seeded there and cleared afterwards would be skipped while unproduced --
// its consumers would then read a value that is not there. Only cells that no
// clear can reach (bound to no loop instance at all) are seeded.
//
// The INPUT-MIRRORED configuration is the one whose schedule has such a cell
// (MEASURED: the default configuration's produce_if_absent cells are all
// unbound), so this case sets that configuration directly.
// ===========================================================================
TEST_CASE(
    "ordered executor: a produce-if-absent cell bound to an enclosing loop is "
    "re-produced on that loop's batches",
    "[ordered][pia-rebind]") {
  using sequant::eval::dryrun::EvalExprDryRun;
  using sequant::eval::dryrun::EvalNodeDryRun;
  using Node = EvalNodeDryRun;

  auto ctx = sequant::get_default_context().clone();
  ctx.set_first_dummy_index_ordinal(1000000);
  auto isr = ctx.mutable_index_space_registry();
  REQUIRE(isr != nullptr);
  sequant::mbpt::add_pao_spaces(isr, sequant::mbpt::Spin::any);
  sequant::mbpt::add_df_spaces(isr);
  auto ctx_resetter = sequant::set_scoped_default_context(std::move(ctx));

  auto const body =
      orderedexec_witness_slurp(std::string(SEQUANT_UNIT_TESTS_SOURCE_DIR) +
                                "/data/csv_ccsd_doubles_residual_df.txt");
  REQUIRE(!body.empty());
  std::string line = body;
  if (auto nl = line.find('\n'); nl != std::string::npos)
    line = line.substr(0, nl);
  auto expr = sequant::deserialize<sequant::ExprPtr>(line);
  REQUIRE(expr->is<sequant::Sum>());
  auto const& summands = expr->as<sequant::Sum>().summands();
  auto regime = orderedexec_witness_df_regime(kOrderedExecWater20_pVDZF12);
  auto cm = std::make_shared<sequant::eval::dryrun::CostModel const>(regime);
  sequant::BatchPolicy policy;
  policy.is_batchable_contracted_index = [](sequant::Index const& ix) {
    return ix.space().base_key() == L"Κ";
  };
  policy.is_batchable_external_index = [](sequant::Index const& ix) {
    auto const reg = sequant::get_default_context().index_space_registry();
    return reg && ix.space() && reg->is_pure_occupied(ix.space());
  };
  policy.batch_spectator_indices = true;
  policy.node_level_placement = true;
  policy.batch_target_size = [](sequant::Index const& ix) -> std::size_t {
    return ix.space().base_key() == L"Κ" ? 256 : 16;
  };
  policy.is_volatile_leaf = [](sequant::Tensor const& t) {
    return t.label() == L"t";
  };
  policy.accumulation_factor = 1.0;
  policy.persistent_only = false;
  policy.peak_threshold = 25e9;  // THE mirrored setting
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
  std::vector<Node> forest;
  for (auto const& s : summands) {
    sequant::ExprPtr const term = orderedexec_witness_flatten_product(s);
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

  auto const block_of = [](sequant::Index const& ix) -> std::size_t {
    return ix.space().base_key() == L"Κ" ? 256 : 16;
  };
  auto const rich = sequant::eval::compute_dag_boulevard(forest, *cm, block_of);
  auto const legality = sequant::eval::analyze_legality(rich, forest, policy);
  auto const ordered = sequant::eval::build_ordered_schedule(
      rich, legality, policy, std::initializer_list<std::wstring>{});
  REQUIRE(sequant::eval::well_formed(ordered));

  using annot_t = std::remove_cvref_t<decltype(forest.front()->annot())>;
  annot_t const layout{};
  sequant::eval::dryrun::DryRunLeafEvaluator const yield{cm};
  std::function<std::size_t(sequant::Index const&)> const target =
      [](sequant::Index const& ix) -> std::size_t {
    return ix.space().base_key() == L"Κ" ? 256 : 16;
  };
  std::function<bool(Node const&)> const is_volatile_node =
      [p = policy.is_volatile_leaf](Node const& n) -> bool {
    if (!n.leaf() || !n->is_tensor()) return false;
    return p && p(n->as_tensor());
  };
  auto const vmap = sequant::eval::build_value_node_map(forest);
  auto aops = sequant::eval::dryrun::make_dryrun_array_ops(cm);

  // The cell this case is about, from the schedule's own derived table.
  auto const sma = sequant::eval::compute_sliced_mode_assignment(ordered, rich);
  sequant::eval::CellTableInputs in;
  in.ordered = &ordered;
  in.rich = &rich;
  in.sliced = &sma;
  in.sliced_modes_of = [&](std::size_t vid) {
    auto const it = vmap.find(sequant::eval::value_key_of(rich.cells[vid]));
    REQUIRE(it != vmap.end());
    return sequant::eval::detail::home_modes_in_cell_frame(rich, vid,
                                                           it->second);
  };
  in.volatile_of = [&](std::size_t vid) {
    auto const it = vmap.find(sequant::eval::value_key_of(rich.cells[vid]));
    return it != vmap.end() &&
           sequant::subtree_any(it->second, is_volatile_node);
  };
  in.n_batches_of =
      sequant::eval::detail::ordered_n_batches_by_loop(ordered, target, &aops);
  in.operands_of = orderedexec_per_leg_operands(rich, vmap);
  auto const table = sequant::eval::build_cell_table(in);

  std::optional<sequant::eval::CellId> bound_pia;
  for (sequant::eval::CellId c = 0; c < table.cells.size() && !bound_pia; ++c) {
    auto const& tc = table.cells[c];
    if (!tc.produce_if_absent) continue;
    if (sequant::eval::detail::bound_instances(tc).empty()) continue;
    bound_pia = c;
  }
  REQUIRE(bound_pia.has_value());  // this configuration really has one
  auto const& pia = table.cells[*bound_pia];
  auto const bound_key = sequant::eval::detail::bound_instances(pia).front();
  // The loop it is bound to is a real, MULTI-batch loop of this schedule --
  // otherwise "once per batch of the enclosing loop" and "once" coincide and
  // the assertion below would be vacuous.
  std::size_t const n_bound_batches = in.n_batches_of(bound_key);
  INFO("produce_if_absent cell#"
       << *bound_pia << " (value " << pia.value_id << ", scope depth "
       << pia.scope.path.size() << ") is bound to a loop of " << n_bound_batches
       << " batches");
  REQUIRE(n_bound_batches > 1);

  auto ordered_cache = sequant::cache_manager(forest);
  ordered_cache.set_array_ops(&aops);
  ordered_cache.set_recompute_tally_enabled(true);
  auto& logger = sequant::Logger::instance();
  auto const prev_level = logger.eval.level;
  logger.eval.level = 1;  // arms the build tally
  std::ostringstream sink;
  auto* const prev_stream = logger.eval.stream;
  logger.eval.stream = &sink;
  char const* const prev_strict = std::getenv("SEQUANT_UT_STRICT_FILL_ONCE");
  std::string const prev_strict_val = prev_strict ? prev_strict : "";
  setenv("SEQUANT_UT_STRICT_FILL_ONCE", "1", 1);
  REQUIRE_NOTHROW(sequant::eval::evaluate_ordered_schedule<sequant::Trace::On>(
      forest, ordered, rich, layout, yield, ordered_cache, target, {},
      is_volatile_node));
  if (prev_strict)
    setenv("SEQUANT_UT_STRICT_FILL_ONCE", prev_strict_val.c_str(), 1);
  else
    unsetenv("SEQUANT_UT_STRICT_FILL_ONCE");
  logger.eval.level = prev_level;
  logger.eval.stream = prev_stream;

  // The cell's own value is produced once per batch of the loop it is bound
  // to -- NOT once. (Its production is the per-batch form the Assemble folds,
  // which the executor evaluates through the same evaluate_impl the build
  // tally counts.)
  auto const nit =
      vmap.find(sequant::eval::value_key_of(rich.cells[pia.value_id]));
  REQUIRE(nit != vmap.end());
  std::size_t const builds =
      orderedexec_builds_of(ordered_cache.recompute_tally(), nit->second);
  INFO("value " << pia.value_id << " built " << builds << " times; the loop it "
                << "is bound to has " << n_bound_batches << " batches");
  CHECK(builds >= n_bound_batches);
}

// The rule that closes the hazard the case above characterizes, on a
// hand-built table so the two halves of it are visible side by side: a
// `produce_if_absent` cell BOUND to a loop instance loses its value when that
// loop advances, so it must never be marked skipped for a whole visit (the
// mark outlives the clear); an UNBOUND one is reachable by no clear and so may
// be. Before the fix the seeding rule was "produce_if_absent and currently
// held", which admits the bound cell and elides its re-production.
TEST_CASE(
    "ordered executor: only an unbound produce-if-absent cell may be seeded "
    "into a per-visit skip set",
    "[ordered][pia-rebind]") {
  using sequant::eval::CellRegistry;
  using sequant::eval::CellTable;
  using sequant::eval::LoopKey;
  using sequant::eval::ProductionKind;
  using sequant::eval::TableCell;

  LoopKey const outer{1, 0}, inner{2, 0};
  CellTable t;
  // cell 0: homed inside the INNER loop, sliced by the OUTER one -- invariant
  // to its own loop (so produce_if_absent) but not to the outer one.
  TableCell bound;
  bound.value_id = 0;
  bound.production.kind = ProductionKind::Build;
  bound.scope.path = {{outer, 0}, {inner, 0}};
  bound.sliced = {{0, outer}};
  bound.produce_if_absent = true;
  bound.life = 1;
  t.cells.push_back(bound);
  // cell 1: same home, bound to nothing at all.
  TableCell unbound;
  unbound.value_id = 1;
  unbound.production.kind = ProductionKind::Build;
  unbound.scope.path = {{outer, 0}, {inner, 0}};
  unbound.produce_if_absent = true;
  unbound.life = 1;
  t.cells.push_back(unbound);

  CellRegistry reg(t);
  sequant::eval::dryrun::SizeRegime regime;
  regime.space_extent = {{L"i", 8}};
  auto cm = std::make_shared<sequant::eval::dryrun::CostModel const>(regime);
  sequant::ResultPtr const r =
      std::make_shared<sequant::eval::dryrun::ResultDryRun>(
          sequant::container::svector<sequant::Index>{sequant::Index{L"i_1"}},
          cm);
  reg.set(0, r);
  reg.set(1, r);

  // The inner loop's own batch boundary keeps both (that is what
  // produce_if_absent means); the OUTER loop's boundary empties the bound one.
  reg.clear_bound_to(inner);
  CHECK(reg.peek(0) == r);
  CHECK(reg.peek(1) == r);
  reg.clear_bound_to(outer);
  CHECK_FALSE(reg.peek(0));  // cleared: it must be produced again
  CHECK(reg.peek(1) == r);   // no clear can reach it

  // ... which is exactly what the seeding rule has to encode, since a
  // per-visit mark is decided before those clears and consulted after them.
  CHECK_FALSE(sequant::eval::detail::ordered_visit_skip_seedable(t.cells[0]));
  CHECK(sequant::eval::detail::ordered_visit_skip_seedable(t.cells[1]));
}
