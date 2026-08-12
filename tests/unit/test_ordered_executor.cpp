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
#include <SeQuant/core/eval/backends/dryrun/result.hpp>
#include <SeQuant/core/eval/backends/dryrun/size_regime.hpp>
#include <SeQuant/core/eval/cache_manager.hpp>
#include <SeQuant/core/eval/eval.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/eval/legality.hpp>
#include <SeQuant/core/eval/node_batch_annotation.hpp>
#include <SeQuant/core/eval/ordered_executor.hpp>
#include <SeQuant/core/eval/ordered_schedule.hpp>
#include <SeQuant/core/eval/peak_profile.hpp>
#include <SeQuant/core/eval/result.hpp>
#include <SeQuant/core/eval/scope_executor.hpp>
#include <SeQuant/core/eval/scope_schedule.hpp>
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
  ordered_cache.set_recompute_tally_enabled(true);
  ResultPtr ord_result;
  try {
    ord_result = sequant::eval::evaluate_ordered_schedule<sequant::Trace::On>(
        forest, ordered, rich, layout, yield, ordered_cache, target);
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
  ws_cache.set_recompute_tally_enabled(true);
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
}
