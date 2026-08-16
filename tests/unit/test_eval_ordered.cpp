// Task 4 of the multi-root single-DAG eval plan (see
// doc/dev/specs/2026-08-15-multiroot-single-dag-eval-design.md in the
// consuming MPQC repo): SeQuant-side MULTI-ROOT ordered evaluate -- one
// schedule built over SEVERAL INDEPENDENT root trees, returning one result
// PER ROOT (a map, no cross-root summation), so a subexpression shared
// across roots is built exactly once.
//
// Mirrors test_ordered_executor.cpp's minimal harness exactly: the
// equivalence/build-count checks need REAL numeric arithmetic (not just
// shape/size modeling under the zero-data DryRun backend), so this reuses
// the identical minimal ScalarEvalExpr subclass (Constant/Variable leaves,
// Sum/Product internal nodes, no tensor backend) that file introduced --
// duplicated here (no shared test header exists between the two .cpp
// files), per that file's own precedent for its DryRun witness fixtures.

#include <SeQuant/core/batch_policy.hpp>
#include <SeQuant/core/container.hpp>
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
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/expressions/tensor.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cstdlib>
#include <functional>
#include <map>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <variant>

namespace {

using sequant::Constant;
using sequant::EvalExpr;
using sequant::EvalNode;
using sequant::ExprPtr;
using sequant::ResultPtr;
using sequant::ResultScalar;
using sequant::Variable;
using sequant::container::svector;
using sequant::eval::analyze_legality;
using sequant::eval::build_ordered_schedule;
using sequant::eval::BuildStep;
using sequant::eval::compute_dag_boulevard;
using sequant::eval::evaluate_ordered_multiroot;
using sequant::eval::evaluate_ordered_schedule;
using sequant::eval::OrderedSchedule;
using sequant::eval::RichSchedule;
using sequant::eval::ScopeBlock;
using sequant::eval::dryrun::CostModel;
using sequant::eval::dryrun::SizeRegime;

///
/// \brief A minimal EvalExpr subclass whose only job is to satisfy \c
/// meta::can_evaluate (i.e. carry an annot() method) so a plain scalar
/// arithmetic forest (Constant/Variable leaves, Sum/Product internal nodes,
/// no tensor indices at all) can be run through \c evaluate_impl / \c
/// evaluate_ordered_schedule / \c evaluate_ordered_multiroot without pulling
/// in a tensor backend. Identical to test_ordered_executor.cpp's (and
/// test_scope_executor.cpp's) ScalarEvalExpr.
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
/// Identical to test_ordered_executor.cpp's ScalarLeafEvaluator.
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

// Number of ValueCell's in `rich` whose hash equals `hash` -- the CSE
// acceptance: a value shared structurally across the roots that produced
// `rich` (via compute_dag_boulevard) must fold to exactly ONE cell, no
// matter how many roots reference it. Mirrors the counting
// evaluate_ordered_schedule's own SEQUANT_UT_CELL_CHECK diagnostic does
// internally (ordered_executor.hpp), specialized to one target hash.
std::size_t cells_with_hash(RichSchedule const& rich, std::size_t hash) {
  std::size_t n = 0;
  for (auto const& c : rich.cells)
    if (c.hash == hash) ++n;
  return n;
}

// Number of BuildStep/output production sites in `ordered`'s WHOLE
// ScopeBlock tree whose value_id resolves (via rich.cells) to `hash` -- the
// schedule-structural "built once" proof: well_formed()'s static single-
// producer invariant (ordered_schedule.hpp) plus evaluate_ordered_schedule's
// own run-completeness assert (SEQUANT_ASSERT(built[vid])) together
// guarantee this count is exactly the number of times the value is actually
// realized at runtime -- neither more (duplicate production is a static
// well-formedness violation) nor fewer (an unrealized promise aborts the
// run). Walks the same shape (recursing into every nested ScopeBlock)
// ordered_executor.hpp's SEQUANT_UT_CELL_CHECK diagnostic already does.
std::size_t production_sites_with_hash(OrderedSchedule const& ordered,
                                       RichSchedule const& rich,
                                       std::size_t hash) {
  std::unordered_map<std::size_t, std::size_t> vid_to_hash;
  vid_to_hash.reserve(rich.cells.size());
  for (auto const& c : rich.cells) vid_to_hash[c.value_id] = c.hash;

  std::size_t n = 0;
  auto const at_hash = [&](std::size_t vid) {
    auto const it = vid_to_hash.find(vid);
    return it != vid_to_hash.end() && it->second == hash;
  };
  auto walk = [&](auto&& self, ScopeBlock const& b) -> void {
    for (auto const& s : b.steps) {
      if (auto const* bs = std::get_if<BuildStep>(&s.value)) {
        if (at_hash(bs->value_id)) ++n;
      } else if (auto const* ch = std::get_if<ScopeBlock>(&s.value)) {
        self(self, *ch);
      }
    }
    for (auto const& [vid, kind] : b.outputs) {
      (void)kind;
      if (at_hash(vid)) ++n;
    }
  };
  walk(walk, ordered.root);
  return n;
}

}  // namespace

TEST_CASE(
    "evaluate_ordered_multiroot builds a subexpression shared across two "
    "independent roots exactly once, and matches per-root forest descent",
    "[eval][ordered]") {
  // Two INDEPENDENT roots (not summands of one equation) that structurally
  // share the SAME subexpression A = a * b -- root1 = c + (a*b), root2 =
  // d + (a*b). Deserialize+binarize each root SEPARATELY (as any two
  // genuinely independent equations would be), so any folding to one A is
  // entirely compute_dag_boulevard's structural-hash CSE, not object
  // identity.
  //
  // A is deliberately placed as the RIGHT summand in both roots (never the
  // left): every binary Sum node binarize() produces marks itself
  // accumulate_in_place (EvalExpr::accumulate_in_place, eval_expr.cpp's
  // binarize(Sum const&, ...)), which unconditionally mutates its LEFT
  // operand's own Result object in place (result = move(f.left); result->
  // add_inplace(*f.right)) rather than allocating a fresh one. That is safe
  // when the left operand is used exactly once (a per-root leaf like c/d
  // here) but would corrupt a value with a SECOND independent consumer if it
  // were ever the left operand -- exactly A's situation once it is CSE'd
  // across two roots. Putting A on the right sidesteps this pre-existing,
  // orthogonal hazard (guarded upstream by a chain_holds() assert that is
  // compiled out under this Release build's SEQUANT_ASSERT_BEHAVIOR_=IGNORE)
  // without weakening what THIS test verifies: A is still the identical
  // shared subexpression, built once and read (never mutated) by both roots.
  ScalarNode const root1 = scalar_tree(L"c + (a * b)");
  ScalarNode const root2 = scalar_tree(L"d + (a * b)");
  // A probed independently (never inserted into either root), purely to
  // learn its canonical hash for the cell/production-site counts below.
  ScalarNode const a_probe = scalar_tree(L"a * b");
  std::size_t const a_hash = a_probe->hash_value();

  ScalarLeafEvaluator const yield{
      {{L"a", 2.0}, {L"b", -3.5}, {L"c", 7.25}, {L"d", -1.5}}};

  // NO batchable index -- same no-batching BatchPolicy as
  // test_ordered_executor.cpp's unbatched-forest test.
  sequant::BatchPolicy const policy;
  SizeRegime const regime;
  CostModel const cm{regime};
  auto const block_of = [](sequant::Index const&) -> std::size_t { return 1; };
  std::function<std::size_t(sequant::Index const&)> const target =
      [](sequant::Index const&) -> std::size_t { return 1; };

  svector<ScalarNode> const roots{root1, root2};
  RichSchedule const rich = compute_dag_boulevard(roots, cm, block_of);
  auto const legality = analyze_legality(rich, roots, policy);
  OrderedSchedule const ordered =
      build_ordered_schedule(rich, legality, policy, {});

  // ---- Structural CSE proof: exactly one cell, one production site. ----
  CHECK(cells_with_hash(rich, a_hash) == 1);
  CHECK(production_sites_with_hash(ordered, rich, a_hash) == 1);

  // ---- Dynamic build-count proof: a custom_evaluator interception (see
  // CacheManager::custom_evaluator_type) fires only on a genuine build
  // attempt (AFTER the Enter-stage cache probe misses), so counting its
  // calls at A's hash is the actual runtime "how many times was A built"
  // tally -- the counting-evaluator mechanism test_ordered_executor.cpp's
  // own build-tally helpers are the dryrun-backend analog of. ----
  std::size_t a_build_count = 0;
  auto cache = sequant::CacheManager<ScalarNode>::empty();
  cache.set_custom_evaluator(
      [&](ScalarNode const& node,
          sequant::CacheManager<ScalarNode>&) -> ResultPtr {
        if (!node.leaf() && node->hash_value() == a_hash) ++a_build_count;
        return nullptr;  // decline: let the standard scheme build it
      });

  svector<ResultPtr> const res = evaluate_ordered_multiroot(
      roots, ordered, rich, ScalarEvalExpr::annot_t{}, yield, cache, target);

  REQUIRE(res.size() == 2);
  CHECK(a_build_count == 1);

  // ---- Per-root numeric correctness: each result matches evaluating that
  // root ALONE (its own fresh cache, plain per-tree forest descent). ----
  auto ref_cache1 = sequant::CacheManager<ScalarNode>::empty();
  ResultPtr const ref1 = sequant::evaluate(root1, yield, ref_cache1);
  auto ref_cache2 = sequant::CacheManager<ScalarNode>::empty();
  ResultPtr const ref2 = sequant::evaluate(root2, yield, ref_cache2);

  double const got1 = res[0]->as<ResultScalar<double>>().value();
  double const got2 = res[1]->as<ResultScalar<double>>().value();

  // Hand-computed cross-check that the references themselves are right:
  //   A = a*b = 2*(-3.5) = -7
  //   root1 = c + A = 7.25 + (-7) = 0.25
  //   root2 = d + A = -1.5 + (-7) = -8.5
  CHECK(ref1->as<ResultScalar<double>>().value() == Catch::Approx(0.25));
  CHECK(ref2->as<ResultScalar<double>>().value() == Catch::Approx(-8.5));

  CHECK(got1 == Catch::Approx(ref1->as<ResultScalar<double>>().value()));
  CHECK(got2 == Catch::Approx(ref2->as<ResultScalar<double>>().value()));

  // NOT summed: res[0] + res[1] would equal -14, the (wrong, for this API)
  // whole-scope-style forest-descent sum. Assert the map contract instead.
  CHECK(got1 != Catch::Approx(got1 + got2));
}

TEST_CASE(
    "evaluate_ordered_multiroot on a single root is a regression anchor: "
    "identical to the existing single-forest ordered path for one root",
    "[eval][ordered]") {
  ScalarNode const root = scalar_tree(L"2 * a * b - c");
  ScalarLeafEvaluator const yield{{{L"a", 2.0}, {L"b", -3.5}, {L"c", 7.25}}};

  sequant::BatchPolicy const policy;
  SizeRegime const regime;
  CostModel const cm{regime};
  auto const block_of = [](sequant::Index const&) -> std::size_t { return 1; };
  std::function<std::size_t(sequant::Index const&)> const target =
      [](sequant::Index const&) -> std::size_t { return 1; };

  svector<ScalarNode> const roots{root};
  RichSchedule const rich = compute_dag_boulevard(roots, cm, block_of);
  auto const legality = analyze_legality(rich, roots, policy);
  OrderedSchedule const ordered =
      build_ordered_schedule(rich, legality, policy, {});

  // The existing single-forest ordered path.
  auto ordered_cache = sequant::CacheManager<ScalarNode>::empty();
  ResultPtr const single =
      evaluate_ordered_schedule(roots, ordered, rich, ScalarEvalExpr::annot_t{},
                                yield, ordered_cache, target);

  // The new multi-root entry, called with exactly one root.
  auto multiroot_cache = sequant::CacheManager<ScalarNode>::empty();
  svector<ResultPtr> const multi = evaluate_ordered_multiroot(
      roots, ordered, rich, ScalarEvalExpr::annot_t{}, yield, multiroot_cache,
      target);

  REQUIRE(multi.size() == 1);
  double const single_val = single->as<ResultScalar<double>>().value();
  double const multi_val = multi[0]->as<ResultScalar<double>>().value();
  CHECK(multi_val == Catch::Approx(single_val));
  CHECK(single_val == Catch::Approx(-21.25));
}

TEST_CASE(
    "sequant::evaluate_multiroot dispatches through cache.multiroot_driver() "
    "and throws when unset (no per-root fallback)",
    "[eval][ordered]") {
  // Two INDEPENDENT roots that structurally SHARE the subexpression A = a * b.
  //
  // root1 = (a*b + c) * (a*b): A appears BOTH as the LEFT operand of the
  // inner Sum (a*b + c) AND as the right factor of the outer Product. Every
  // binary Sum binarize() marks accumulate_in_place (EvalExpr::binarize(Sum)),
  // and the in-place branch (eval.hpp) would, unless it detects the hazard at
  // runtime, MOVE-then-add_inplace its left operand's own Result buffer.
  // Because A is CSE'd (used three times across the two roots) it is homed
  // RESIDENT by the ordered executor (ordered_executor.hpp's ensure_home_slot),
  // so that one buffer is read by every consumer. Accumulating (a*b + c) in
  // place therefore MUTATES the shared A buffer, and the SAME root's outer
  // Product then reads the corrupted A -- a deterministic, order-independent
  // corruption (the inner Sum is evaluated before the outer Product's own
  // read of A, within one post-order descent). root2 = (a*b) * d likewise
  // reads the shared A. The pre-Release guard against this was a
  // SEQUANT_ASSERT(!chain_holds...) that compiles out under this build
  // (SEQUANT_ASSERT_ENABLED undefined), so the safety must be a RUNTIME
  // condition: eval.hpp gates in-place on !cache.chain_holds_shared(f.left),
  // falling back to the allocating sum() when the accumulator is a shared
  // (multi-consumer, resident) buffer -- see CacheManager::chain_holds_shared.
  //
  // NB the two roots are chosen so their top nodes do NOT collide under the
  // CSE hash (a Product and a differently-shaped Product), so A is genuinely a
  // shared child of two distinct root cells rather than the two roots folding
  // into one.
  //
  // Hand-computed (A = a*b = 2*(-3.5) = -7):
  //   root1 = (A + c) * A = (-7 + 7.25) * (-7) = 0.25 * (-7) = -1.75
  //   root2 = A * d       = (-7) * 1.5         = -10.5
  // With the in-place hazard UNGUARDED, (a*b + c) mutates A to 0.25, so root1
  // would wrongly yield 0.25 * 0.25 = 0.0625 (and root2 0.25 * 1.5 = 0.375).
  ScalarNode const root1 = scalar_tree(L"((a * b) + c) * (a * b)");
  ScalarNode const root2 = scalar_tree(L"(a * b) * d");
  ScalarLeafEvaluator const yield{
      {{L"a", 2.0}, {L"b", -3.5}, {L"c", 7.25}, {L"d", 1.5}}};

  sequant::BatchPolicy const policy;
  SizeRegime const regime;
  CostModel const cm{regime};
  auto const block_of = [](sequant::Index const&) -> std::size_t { return 1; };
  std::function<std::size_t(sequant::Index const&)> const target =
      [](sequant::Index const&) -> std::size_t { return 1; };

  svector<ScalarNode> const roots{root1, root2};
  RichSchedule const rich = compute_dag_boulevard(roots, cm, block_of);
  auto const legality = analyze_legality(rich, roots, policy);
  OrderedSchedule const ordered =
      build_ordered_schedule(rich, legality, policy, {});

  SECTION("no driver installed: throws, no per-root fallback") {
    auto cache = sequant::CacheManager<ScalarNode>::empty();
    CHECK_THROWS_AS(
        sequant::evaluate_multiroot(roots, std::string{}, yield, cache),
        std::logic_error);
  }

  SECTION("driver installed: routes to evaluate_ordered_multiroot") {
    auto cache = sequant::CacheManager<ScalarNode>::empty();
    cache.set_multiroot_driver(
        [&](svector<ScalarNode> const& rs, std::string const& lay,
            sequant::CacheManager<ScalarNode>& c) -> svector<ResultPtr> {
          return evaluate_ordered_multiroot(rs, ordered, rich, lay, yield, c,
                                            target);
        });

    svector<ResultPtr> const res =
        sequant::evaluate_multiroot(roots, std::string{}, yield, cache);
    REQUIRE(res.size() == 2);
    // Each root reads the shared, homed A = a*b; the in-place guard keeps
    // (a*b + c)'s accumulation from mutating that shared buffer. Unguarded,
    // root1 would be 0.0625 (and root2 0.375); see the header comment.
    CHECK(res[0]->as<ResultScalar<double>>().value() == Catch::Approx(-1.75));
    CHECK(res[1]->as<ResultScalar<double>>().value() == Catch::Approx(-10.5));

    // Cross-check the references (independent per-tree forest descent, each
    // with its own fresh cache -- no cross-root homing, so no in-place
    // hazard) agree with the multi-root results.
    auto rc1 = sequant::CacheManager<ScalarNode>::empty();
    auto rc2 = sequant::CacheManager<ScalarNode>::empty();
    CHECK(sequant::evaluate(root1, yield, rc1)
              ->as<ResultScalar<double>>()
              .value() == Catch::Approx(-1.75));
    CHECK(sequant::evaluate(root2, yield, rc2)
              ->as<ResultScalar<double>>()
              .value() == Catch::Approx(-10.5));
  }
}
