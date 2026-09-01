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
        auto const it = vmap.find(vc.hash);
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

  // Task 2 (eager-home-release plan): detail::ordered_n_blocks reports > 1
  // blocks for the batched Κ (aux) axis at target 256 (matches this test's
  // own batch_target_size policy above), and exactly 1 for an unbatched
  // index TYPE (never realized as a ScopeBlock axis anywhere in `ordered`).
  {
    std::optional<sequant::Index> kappa_index;
    std::function<void(sequant::eval::ScopeBlock const&)> find_kappa_axis =
        [&](sequant::eval::ScopeBlock const& b) {
          if (kappa_index) return;
          if (b.axis.space().base_key() == L"Κ") {
            kappa_index = b.axis;
            return;
          }
          for (auto const& s : b.steps)
            if (auto const* child =
                    std::get_if<sequant::eval::ScopeBlock>(&s.value))
              find_kappa_axis(*child);
        };
    for (auto const& s : ordered.root.steps)
      if (auto const* child = std::get_if<sequant::eval::ScopeBlock>(&s.value))
        find_kappa_axis(*child);
    REQUIRE(kappa_index.has_value());

    std::optional<sequant::Index> unbatched_index;
    for (auto const& vc : rich.cells) {
      for (sequant::Index const& ix : vc.carried) {
        if (ix.space().base_key() != L"Κ") {
          unbatched_index = ix;
          break;
        }
      }
      if (unbatched_index) break;
    }
    REQUIRE(unbatched_index.has_value());

    auto const nb = sequant::eval::detail::ordered_n_blocks(
        ordered, rich, vmap, yield, target, /*aops=*/&aops);
    CHECK(nb(*kappa_index) > 1);
    CHECK(nb(*unbatched_index) == 1);
  }

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

  // ---- Task 4 exactness gate: predicted (ordered_home_reads) == measured.
  // A second ordered run with a NON-EVICTING home life lets us read the ACTUAL
  // number of times each homed value's home entry is accessed (max_life - life)
  // and compare it to the static home_reads prediction. If they differ, the
  // static model is WRONG (an undercount rebuilds; an overcount over-retains)
  // -- this is the ground truth the payoff must satisfy.
  {
    auto meas_cache = sequant::cache_manager(forest);
    meas_cache.set_array_ops(&aops);
    std::function<std::size_t(std::size_t)> const huge =
        [](std::size_t) -> std::size_t { return 1000000000ull; };
    try {
      (void)sequant::eval::evaluate_ordered_schedule<sequant::Trace::On>(
          forest, ordered, rich, layout, yield, meas_cache, target, {},
          is_volatile_node, huge);
    } catch (std::exception const&) {
    }
    auto const n_blocks_m = sequant::eval::detail::ordered_n_blocks<Node>(
        ordered, rich, vmap, yield, target, /*aops=*/&aops);
    auto const predicted = sequant::eval::detail::ordered_home_reads<Node>(
        ordered, rich, vmap, n_blocks_m);
    std::size_t n_checked = 0, n_mismatch = 0;
    for (auto const& vc : rich.cells) {
      auto const vit = vmap.find(vc.hash);
      if (vit == vmap.end() || vit->second.leaf()) continue;
      if (!sequant::subtree_any(vit->second, is_volatile_node)) continue;
      // Only ROOT-homed composites live in the root meas_cache and so are
      // measurable post-run; a block-homed value's entry is in a transient
      // scratch (gone by now), so measuring it against the root cache is a
      // false 0. Restrict the exactness gate to root-level BuildSteps.
      if (!orderedexec_index_of_build_step(ordered.root, vc.value_id)
               .has_value())
        continue;
      int const ml = meas_cache.max_life(vit->second);
      int const lf = meas_cache.life(vit->second);
      if (ml <= 0) continue;  // never homed non-persistently (persistent path)
      std::size_t const measured = std::size_t(ml - lf);
      // Root-level BuildStep home scope is ROOT (empty key); per-cell
      // ordered_home_reads takes the home-scope key.
      std::size_t const pred = predicted(vc.value_id, {});
      ++n_checked;
      if (pred != measured) {
        ++n_mismatch;
        if (std::getenv("SEQUANT_UT_T4_DIAG"))
          std::wcerr << L"    [HRDIAG] vid=" << vc.value_id << L" predicted="
                     << pred << L" measured=" << measured << L"\n";
      }
    }
    std::wcerr << L"  home_reads exactness: checked " << n_checked
               << L" volatile homed composites, " << n_mismatch
               << L" mismatches\n";
    INFO("home_reads predicted vs measured: " << n_mismatch << "/" << n_checked
                                              << " mismatch");
    CHECK(n_mismatch == 0);
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
    auto const vit = vmap.find(vc.hash);
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
    auto const vit = vmap.find(vc.hash);
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
      auto const it = vmap.find(rich.cells[vid].hash);
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
TEST_CASE(
    "ordered executor: water-20 aux+occ residual dry-run walk completes "
    "without "
    "a vanished home value",
    "[.][ordered][w20-auxocc-walk]") {
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
        auto const it = vmap_dump.find(rich.cells[vid].hash);
        if (it != vmap_dump.end())
          for (auto const& x : it->second->sliced_modes())
            s += sequant::toUtf8(x.full_label()) + " ";
      }
      return s;
    };
    std::function<void(sequant::eval::ScopeBlock const&, int)> dump =
        [&](sequant::eval::ScopeBlock const& b, int d) {
          std::string ind(2 * d, ' ');
          if (d > 0)
            std::cerr << ind << "BLOCK axis="
                      << sequant::toUtf8(b.axis.space().base_key()) << " kind="
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
                        << " carried=[" << carried_str(bs->value_id) << "]\n";
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
                      << sliced_str(ovid) << "]\n";
        };
    std::cerr << "=== ORDERED SCHEDULE TREE ===\n";
    dump(ordered.root, 0);
    std::cerr << "=== END TREE ===\n";
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
  // completion. It throws today ("read-from-home value vanished ... 98403") --
  // the over-homing defect; fixing it makes this pass.
  REQUIRE_NOTHROW(sequant::eval::evaluate_ordered_schedule<sequant::Trace::Off>(
      forest, ordered, rich, layout, yield, ordered_cache, target, {},
      is_volatile_node));

  logger.eval.level = prev_level;
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

  // Collect the composites the executor homed PERSISTENT and their iter-1 build
  // counts. The cache-halt gate must not cause any of them to be rebuilt in
  // iteration 2 (they survive reset() and are read, not re-formed).
  std::vector<std::pair<Node, std::size_t>> persistent_b1;
  for (auto const& vc : rich.cells) {
    auto const it = vmap.find(vc.hash);
    if (it == vmap.end() || it->second.leaf()) continue;
    if (cache.entry_is_persistent(it->second))
      persistent_b1.emplace_back(
          it->second,
          orderedexec_builds_of(cache.recompute_tally(), it->second));
  }
  REQUIRE(
      !persistent_b1.empty());  // fixture actually has persistent composites

  cache.reset();

  // Replicate the executor's needed_build gate: BFS from the volatile roots,
  // halting at cache-alive (resident persistent) nodes. A NON-LEAF BuildStep
  // inside a {Κ} ScopeBlock whose node is NOT in this set is a DEAD
  // prerequisite
  // -- the cache-halt fix must not re-form it in iteration 2.
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
      if (cache.alive(n)) continue;  // resident: read, do not descend
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
                auto const it2 = vmap.find(rich.cells[bs->value_id].hash);
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
  auto cm = std::make_shared<sequant::eval::dryrun::CostModel const>(regime);

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
  constexpr std::size_t kAuxBlock = 256, kOccBlock = 8, kPaoBlock = 256;

  sequant::BatchPolicy policy;
  policy.is_batchable_contracted_index = [batch_aux,
                                          batch_pao](sequant::Index const& ix) {
    return (batch_aux && ix.space().base_key() == L"Κ") ||
           (batch_pao && ix.space().base_key() == L"μ̃");
  };
  policy.is_batchable_external_index = [batch_occ](sequant::Index const& ix) {
    return batch_occ && ix.space().base_key() == L"i";
  };
  policy.batch_spectator_indices = batch_occ;
  policy.node_level_placement = batch_occ;  // occ external placement needs it
  policy.batch_target_size = [batch_aux, batch_pao, batch_occ](
                                 sequant::Index const& ix) -> std::size_t {
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
        rows.push_back({e.bytes, carriesK, kind + L"{" + sig + L"}"});
      }
      std::sort(rows.begin(), rows.end(),
                [](LE const& a, LE const& b) { return a.bytes > b.bytes; });
      auto const GB = [](std::size_t b) { return double(b) / 1e9; };
      std::wcerr << L"\n--- [peak-liveset] ordered COLD peak co-resident set, "
                 << L"entries > 0.1 GB (peak_total=" << GB(peak_total)
                 << L" GB) ---\n";
      for (auto const& r : rows)
        if (GB(r.bytes) > 0.1)
          std::wcerr << L"  " << GB(r.bytes) << L" GB  carriesΚ="
                     << (r.carriesK ? L"yes" : L"no ") << L"  " << r.label
                     << L"\n";
      std::wcerr << L"  TOTAL co-resident = " << GB(all_sum)
                 << L" GB;  aux-FREE (no-Κ, aux-batching-immune) floor = "
                 << GB(auxfree_sum) << L" GB\n";
    }

    cache.reset();
    run(o2);  // WARM
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
                auto const it = vmap.find(rich.cells[bs->value_id].hash);
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
      r.builds = total_builds(cache.recompute_tally()) - b0;
      r.peak = std::max<std::size_t>(
          mon.hwmark_bytes,
          std::max<std::size_t>(cache.working_set_hwmark(),
                                static_cast<std::size_t>(peak.load())));
    };
    run(f1);
    cache.reset();
    run(f2);
  }

  cm->set_cost_sink(nullptr);
  logger.eval.level = prev_level;
  logger.eval.stream = prev_stream;

  auto pr = [](wchar_t const* tag, Row const& r) {
    std::wcerr << L"  " << tag << L"  builds=" << r.builds << L"  FLOPs="
               << std::scientific << r.flops << L"  peak_bytes=" << r.peak
               << L"\n";
  };
  std::wcerr << L"\n=== [dryrun-2iter-report] SYSTEM="
             << std::wstring(system.begin(), system.end()) << L" BATCH="
             << std::wstring(batch.begin(), batch.end()) << L", "
             << forest.size() << L" terms, n_Kblocks=" << n_Kblocks
             << L", warm-iter needed_build-skipped {Κ}-block BuildSteps="
             << warm_skipped << L" ===\n";
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
  // at root via ensure_home_slot). A genuinely use-counted (LoopLocal
  // transient) value was released at its last use and is NOT alive here. This
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
  //           executor HOMES at the root cache (ensure_home_slot, pinned) on
  //           block close, so they behave exactly like tier-A.
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

// Task 1 of the eager-home-release plan (the sequel to this file's own
// ordered-executor design): the bounded/persistent CacheManager::
// ensure_home_slot(key, use_count, persistent) overload, which later tasks
// use to replace the unconditional SIZE_MAX pin the zero-arg
// ensure_home_slot(key) installs (see its call site's doc comment above,
// "at root via ensure_home_slot"). Proves both lifetimes directly against
// entry::access()/store()/reset(), independent of the ordered executor
// itself: a volatile (non-persistent) slot is released at its genuine
// use_count-th access, and a persistent slot is never drained and survives
// reset().
TEST_CASE(
    "ensure_home_slot bounds life for volatile and persists for invariant",
    "[ordered-executor]") {
  // A realistic ScalarNode cache key: the root of a scalar-forest tree (see
  // scalar_tree above), a non-leaf node like the ones the ordered executor
  // actually homes.
  ScalarNode const n = scalar_tree(L"2 * a * b - c");
  REQUIRE_FALSE(n.leaf());

  // Volatile: bounded life == use_count, non-persistent, released at its
  // genuine last use. CacheManager::store() itself performs an implicit
  // access() (see its doc comment: "Implictly accesses the stored data,
  // hence, decays the lifetime"), so the store below IS use 1 of 2; one more
  // explicit access_at() is use 2 of 2 and drains the entry.
  auto cache = sequant::CacheManager<ScalarNode>::empty();
  cache.ensure_home_slot(n, /*use_count=*/2, /*persistent=*/false);
  (void)cache.store(
      n, sequant::eval_result<ResultScalar<double>>(1.0));  // use 1 of 2
  CHECK(cache.access_at(n).ptr);        // use 2 of 2 -> drains
  CHECK_FALSE(cache.access_at(n).ptr);  // released after last use

  // Persistent: never drained regardless of use_count, and survives reset().
  auto cache2 = sequant::CacheManager<ScalarNode>::empty();
  cache2.ensure_home_slot(n, /*use_count=*/1, /*persistent=*/true);
  (void)cache2.store(n, sequant::eval_result<ResultScalar<double>>(2.0));
  CHECK(cache2.access_at(n).ptr);
  CHECK(
      cache2.access_at(n).ptr);  // still alive despite use_count=1 (persistent)
  cache2.reset();
  CHECK(cache2.access_at(n).ptr);  // survives reset()

  // Upgrade-existing-entry branch (the `!inserted` arm): a SECOND
  // ensure_home_slot call on the SAME key in the SAME cache must re-arm the
  // already-present entry via set_life()/make_persistent(), not leave its
  // original (possibly already-exhausted) count in place. Task 4 relies on
  // this branch to re-home a node the base cache_manager may have already
  // registered, so it must be exercised, not merely inspected.

  // (a) insert bounded life=1, then UPGRADE to bounded life=3 on the same
  // key/cache. If the upgrade were a no-op (stale life=1 left in place), the
  // store() below -- which itself performs an implicit access() -- would
  // already be the 1-of-1 last use and drain the entry on the spot, so the
  // very next access_at() would immediately return nullptr. Instead it must
  // survive two more explicit reads and drain only at the upgraded 3rd.
  auto cache3 = sequant::CacheManager<ScalarNode>::empty();
  cache3.ensure_home_slot(n, /*use_count=*/1, /*persistent=*/false);  // insert
  cache3.ensure_home_slot(n, /*use_count=*/3,
                          /*persistent=*/false);  // upgrade: set_life(3)
  (void)cache3.store(
      n, sequant::eval_result<ResultScalar<double>>(3.0));  // use 1 of 3
  CHECK(cache3.access_at(n).ptr);        // use 2 of 3 (would be dead if stale)
  CHECK(cache3.access_at(n).ptr);        // use 3 of 3 -> drains
  CHECK_FALSE(cache3.access_at(n).ptr);  // released at the UPGRADED bound

  // (b) insert bounded life=1, then UPGRADE to persistent on the same
  // key/cache. If the upgrade were a no-op (stale non-persistent life=1),
  // store()'s implicit access() would again already be the last use and
  // drain the entry, so the very next access_at() would return nullptr.
  // Instead it must survive repeated reads past the original count=1 AND a
  // reset(), proving make_persistent() latched on the pre-existing entry.
  auto cache4 = sequant::CacheManager<ScalarNode>::empty();
  cache4.ensure_home_slot(n, /*use_count=*/1, /*persistent=*/false);  // insert
  cache4.ensure_home_slot(n, /*use_count=*/1,
                          /*persistent=*/true);  // upgrade: make_persistent()
  (void)cache4.store(n, sequant::eval_result<ResultScalar<double>>(4.0));
  CHECK(cache4.access_at(n).ptr);
  CHECK(cache4.access_at(n).ptr);  // past the original count=1 (would be dead
                                   // if stale)
  cache4.reset();
  CHECK(cache4.access_at(n).ptr);  // survives reset() too
}

// ===========================================================================
// PILLAR 2 (sliced-value canonical-layout / loop-coloring design, task 7
// part 2): the SYMMETRIC-shared-value end-to-end proof of the water-8 fix.
//
// ONE homed-full intermediate S{;i_1,i_2} is fetched by TWO member roots under
// the SAME external-occ loop, each binding the loop to a DIFFERENT free occ
// mode (R1 keeps i_1 / contracts i_2; R2 keeps i_2 / contracts i_1). The
// consumer-blind seam returned S's FIRST stamped mode to BOTH consumers (the
// bug: "pos0 here, pos0 there"); the consumer-aware slice_to_use binds each
// consumer to its OWN mode. A slice_observer witnesses both fetches and proves
// the two consumers slice DIFFERENT physical modes of the shared value.
// ===========================================================================
namespace {

sequant::eval::dryrun::EvalNodeDryRun symslice_leaf(std::string_view tensor) {
  auto expr = sequant::deserialize<sequant::ExprPtr>(std::string(tensor));
  REQUIRE(static_cast<bool>(expr));
  return sequant::eval::dryrun::EvalNodeDryRun{
      sequant::eval::dryrun::EvalExprDryRun{expr->as<sequant::Tensor>()}};
}

sequant::eval::dryrun::EvalNodeDryRun symslice_inode(
    std::string_view result, sequant::eval::dryrun::EvalNodeDryRun l,
    sequant::eval::dryrun::EvalNodeDryRun r) {
  auto expr = sequant::deserialize<sequant::ExprPtr>(std::string(result));
  REQUIRE(static_cast<bool>(expr));
  sequant::eval::dryrun::EvalExprDryRun data{expr->as<sequant::Tensor>()};
  sequant::EvalOpSetter{}.set(data, sequant::EvalOp::Product);
  return sequant::eval::dryrun::EvalNodeDryRun{std::move(data), std::move(l),
                                               std::move(r)};
}

// Stamp every occ ("i") index carried by n as an EXTERNAL batched axis.
void symslice_stamp_ext_occ(sequant::eval::dryrun::EvalNodeDryRun& n) {
  using sequant::BatchModeType;
  sequant::container::svector<std::pair<Index, BatchModeType>> stamps;
  for (auto const& ix : n->canon_indices())
    if (ix.space().base_key() == L"i")
      stamps.push_back({ix, BatchModeType::External});
  if (!stamps.empty()) n->set_node_slice_mask(stamps);
}

}  // namespace

// TEMPORARILY DISABLED (2026-08-25 loop-open-vs-sliced-mask plan, Task 4): this
// case builds the seam manually and asserts the transitional position seam
// (mode_of by consumer hash); Task 4 replaces the seam with the per-occurrence
// positional resolution and rewrites this case against the final semantics.
// Guarded out until then so the suite links.
#if 0
TEST_CASE(
    "consumer-aware slice_to_use: two member roots slice one homed symmetric "
    "shared value on DIFFERENT free modes under the same occ loop (w8 fix)",
    "[ordered-executor][sliced-layout]") {
  using sequant::eval::dryrun::EvalExprDryRun;
  using sequant::eval::dryrun::EvalNodeDryRun;
  using Node = EvalNodeDryRun;

  auto ctx = sequant::get_default_context().clone();
  auto isr = ctx.mutable_index_space_registry();
  REQUIRE(isr != nullptr);
  ctx.set(sequant::AssertStrictBraKetSymmetry::No);
  auto ctx_resetter = sequant::set_scoped_default_context(std::move(ctx));

  // S{;i_1,i_2} = X{a_1;i_1} * Y{a_1;i_2} (contract virtual a_1) -- an occ-occ
  // intermediate, homed full (no node_slice_mask stamp). Two roots share it.
  auto X = symslice_leaf("X{a_1;i_1}");
  auto Y = symslice_leaf("Y{a_1;i_2}");
  auto S = symslice_inode("S{;i_1,i_2}", X, Y);
  auto T = symslice_leaf("T{;i_2}");
  auto U = symslice_leaf("U{;i_1}");
  auto R1 = symslice_inode("R1{;i_1}", S, T);  // contract i_2, external i_1
  auto R2 = symslice_inode("R2{;i_2}", S, U);  // contract i_1, external i_2
  symslice_stamp_ext_occ(R1);
  symslice_stamp_ext_occ(R2);

  std::vector<Node> forest{R1, R2};
  sequant::stamp_lifetime_masks(forest);

  sequant::BatchPolicy policy;
  policy.is_batchable_external_index = [](Index const& ix) {
    return ix.space().base_key() == L"i";
  };
  policy.is_batchable_contracted_index = [](Index const&) { return false; };
  policy.batch_spectator_indices = true;
  policy.node_level_placement = true;

  sequant::eval::dryrun::SizeRegime regime;
  regime.space_extent = {{L"i", 8u}, {L"a", 4u}};
  auto cm = std::make_shared<sequant::eval::dryrun::CostModel const>(regime);

  auto const block_of = [](Index const&) -> std::size_t { return 4; };
  auto rich = compute_dag_boulevard(forest, *cm, block_of);
  REQUIRE(!rich.cells.empty());

  auto const legality = analyze_legality(rich, forest, policy);
  auto const ordered = build_ordered_schedule(rich, legality, policy, {L"i"});
  REQUIRE(sequant::eval::well_formed(ordered));

  auto const assignment =
      sequant::eval::compute_sliced_mode_assignment(ordered, rich);

  // Locate S's cell + hash; both its occ modes are stamped under the ONE occ
  // loop -- the symmetric ambiguity the fix resolves.
  auto const s_hash = S->hash_value();
  auto const s_it =
      std::find_if(rich.cells.begin(), rich.cells.end(),
                   [&](auto const& vc) { return vc.hash == s_hash; });
  REQUIRE(s_it != rich.cells.end());
  REQUIRE(s_it->occurrences.size() == 2);  // shared by R1 and R2
  REQUIRE(s_it->carried.size() == 2);      // carries both i_1 and i_2
  REQUIRE(assignment.levels.size() == 1);  // exactly ONE realized occ loop
  auto const by = assignment.by_value.find(s_it->value_id);
  REQUIRE(by != assignment.by_value.end());
  // Both of S's occ modes are sliced by the SAME loop 0 -- the consumer-blind
  // "one value, one loop, two modes" hole this fix closes.
  REQUIRE(by->second.size() == 2);
  CHECK(by->second[0].second == by->second[1].second);  // same LoopId
  CHECK(by->second[0].first != by->second[1].first);    // different modes
  // The per-occurrence facts carry the consumer attribution.
  REQUIRE(assignment.occ_facts.size() == 2);

  auto const r1h = R1->hash_value();
  auto const r2h = R2->hash_value();

  // ---- The seam (built exactly as the executor builds it): consumer-blind
  // vs consumer-aware resolution. ----
  sequant::LoopColoredSliceSeam seam;
  seam.levels = assignment.levels;
  for (auto const& c : rich.cells) {
    auto const it = assignment.by_value.find(c.value_id);
    if (it != assignment.by_value.end())
      seam.by_hash.emplace(c.hash, it->second);
  }
  for (auto const& [vid, mode, loop, cvid] : assignment.occ_facts)
    seam.by_hash_consumer[rich.cells[vid].hash].push_back(
        std::make_tuple(mode, loop, rich.cells[cvid].hash));

  std::size_t const occ_loop = 0;
  // Consumer-BLIND (2-arg) returns the SAME mode regardless of consumer -- the
  // bug: it would serve one occurrence its sibling's slice.
  auto const blind = seam.mode_of(s_hash, occ_loop);
  REQUIRE(blind.has_value());
  // Consumer-AWARE: each use-site gets its OWN mode, and they DIFFER.
  auto const for_r1 = seam.mode_of(s_hash, occ_loop, r1h);
  auto const for_r2 = seam.mode_of(s_hash, occ_loop, r2h);
  REQUIRE(for_r1.has_value());
  REQUIRE(for_r2.has_value());
  CHECK(*for_r1 != *for_r2);  // THE FIX: the two consumers bind different modes
  // The blind resolution collapses onto exactly one of the two -- proving the
  // consumer-aware step is load-bearing (without it, one consumer mis-slices).
  bool const blind_is_r1 = (blind == for_r1);
  bool const blind_is_r2 = (blind == for_r2);
  CHECK((blind_is_r1 || blind_is_r2));
  CHECK((blind_is_r1 != blind_is_r2));

  // ---- END-TO-END: evaluate through the ordered executor's shared core (the
  // full batched schedule walk -- homing + per-read slicing) with a slice
  // observer, and confirm each consumer's FETCH of S slices its own mode. ----
  sequant::eval::dryrun::DryRunLeafEvaluator const yield{cm};
  std::function<std::size_t(Index const&)> const target =
      [](Index const&) -> std::size_t { return 4; };

  auto& logger = sequant::Logger::instance();
  auto const prev_level = logger.eval.level;
  logger.eval.level = 1;
  auto aops = sequant::eval::dryrun::make_dryrun_array_ops(cm);
  auto cache = sequant::cache_manager(forest);
  cache.set_array_ops(&aops);

  // (node_hash, consumer) -> set of physical slice positions witnessed.
  std::map<std::pair<std::size_t, std::size_t>, std::set<std::size_t>> seen;
  cache.set_slice_observer([&](std::size_t node_hash,
                               std::optional<std::size_t> consumer,
                               std::size_t pmode, int) {
    if (consumer) seen[{node_hash, *consumer}].insert(pmode);
  });

  auto const pre_results =
      sequant::eval::detail::run_ordered_schedule_pre_results<
          sequant::Trace::On>(forest, ordered, rich, yield, cache, target);
  logger.eval.level = prev_level;
  REQUIRE(pre_results.size() == 2);

  // Each consumer sliced S at exactly ONE physical mode, and the two consumers
  // sliced DIFFERENT modes (the runtime witness of the w8 fix). The consumer-
  // blind executor would have sliced the same physical mode for both.
  auto const& s_by_r1 = seen[{s_hash, r1h}];
  auto const& s_by_r2 = seen[{s_hash, r2h}];
  REQUIRE(s_by_r1.size() == 1);
  REQUIRE(s_by_r2.size() == 1);
  CHECK(*s_by_r1.begin() != *s_by_r2.begin());
  // And each matches its seam-resolved mode's physical slot on S.
  CHECK(*s_by_r1.begin() == sequant::index_position(S, *for_r1).value());
  CHECK(*s_by_r2.begin() == sequant::index_position(S, *for_r2).value());
}
#endif  // loop-open-vs-sliced-mask Task 4 rewrite pending
