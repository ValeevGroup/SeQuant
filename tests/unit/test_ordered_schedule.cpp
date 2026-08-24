// Task 2 of the ordered-scope batched-eval design (SP2): pins the
// OrderedSchedule IR (SeQuant/core/eval/ordered_schedule.hpp) -- an ORDERED
// tree of loop blocks and build steps -- plus its well_formed structural
// sanity check. No sequencer/executor here.
//
// Task 3 (below, "[ordered-schedule]" water-20 acceptance test) pins
// build_ordered_schedule -- the deterministic sequencer that lowers SP1's
// LegalitySchedule + the RichSchedule into an OrderedSchedule for the
// NON-SPLIT case.

#include <SeQuant/core/batch_policy.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/eval/backends/dryrun/cost_model_object.hpp>
#include <SeQuant/core/eval/backends/dryrun/eval_expr.hpp>
#include <SeQuant/core/eval/backends/dryrun/size_regime.hpp>
#include <SeQuant/core/eval/eval.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/eval/legality.hpp>
#include <SeQuant/core/eval/lifetime_mask.hpp>
#include <SeQuant/core/eval/ordered_schedule.hpp>
#include <SeQuant/core/eval/peak_profile.hpp>
#include <SeQuant/core/eval/scope_executor.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/optimize/optimize.hpp>
#include <SeQuant/core/optimize/options.hpp>
#include <SeQuant/core/runtime.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/space_qns.hpp>  // mbpt::Spin

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <numeric>
#include <optional>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <variant>
#include <vector>

using sequant::Index;
using sequant::eval::BuildStep;
using sequant::eval::OrderedSchedule;
using sequant::eval::OutputKind;
using sequant::eval::ScopeBlock;
using sequant::eval::Step;
using sequant::eval::well_formed;

TEST_CASE(
    "well_formed accepts a root block with one build step and one child "
    "loop block carrying an AccumulateSum output",
    "[ordered-schedule]") {
  Index const i1{L"i_1"};

  // value_id 1 is the per-iteration transient built inside the loop; value_id
  // 2 is the DISTINCT accumulator that the block's own AccumulateSum output
  // entry produces on close (its single production site IS that output
  // entry, not a BuildStep) -- giving both the same id would double-produce
  // value_id 1 (once as a BuildStep, again as an output), which the
  // single-producer invariant below correctly rejects.
  ScopeBlock child;
  child.axis = i1;
  child.ordinal = 0;
  child.steps.push_back(Step{BuildStep{1}});
  child.outputs.push_back({2, OutputKind::AccumulateSum});

  OrderedSchedule sched;
  sched.root.steps.push_back(Step{BuildStep{0}});
  sched.root.steps.push_back(Step{std::move(child)});
  sched.num_values = 3;

  CHECK(well_formed(sched));
}

// SP2 non-innermost forced split: fork_subchain partitions an already-built
// inner sub-chain into a producer-side and a consumer-side copy by an
// in_consumer(value_id) predicate. A BuildStep goes wholly to one side; a
// nested loop block is forked (duplicated across both sides when its steps
// straddle the partition), with its steps and escape outputs partitioned and a
// side kept only when it has surviving steps.
TEST_CASE("fork_subchain forks an inner sub-chain by consumer-pass membership",
          "[ordered-schedule][fork]") {
  Index const ax{L"i_1"};

  // Inner sub-chain: BuildStep{0}, a nested Κ-loop whose two builds straddle
  // the split (1 producer, 2 consumer) and whose two AccumulateSum outputs
  // likewise straddle it (3 producer, 4 consumer), then BuildStep{5}.
  ScopeBlock inner;
  inner.axis = ax;
  inner.ordinal = 0;
  inner.kind = sequant::BatchModeType::Contracted;
  inner.steps.push_back(Step{BuildStep{1}});
  inner.steps.push_back(Step{BuildStep{2}});
  inner.outputs.push_back({3, OutputKind::AccumulateSum});
  inner.outputs.push_back({4, OutputKind::AccumulateSum});

  sequant::container::vector<Step> steps;
  steps.push_back(Step{BuildStep{0}});
  steps.push_back(Step{std::move(inner)});
  steps.push_back(Step{BuildStep{5}});

  // consumer side = {2, 4, 5}; producer side = {0, 1, 3}.
  std::function<bool(std::size_t)> const in_consumer = [](std::size_t v) {
    return v == 2 || v == 4 || v == 5;
  };

  auto const forked = sequant::eval::detail::fork_subchain(steps, in_consumer);

  auto const build_id = [](Step const& s) {
    return std::get<BuildStep>(s.value).value_id;
  };

  // Producer: BuildStep{0}, then a Κ-loop copy holding only BuildStep{1} and
  // its AccumulateSum output {3}.
  REQUIRE(forked.producer.size() == 2);
  CHECK(build_id(forked.producer[0]) == 0);
  auto const& p_loop = std::get<ScopeBlock>(forked.producer[1].value);
  CHECK(p_loop.axis == ax);
  CHECK(p_loop.ordinal == 0);
  REQUIRE(p_loop.steps.size() == 1);
  CHECK(build_id(p_loop.steps[0]) == 1);
  REQUIRE(p_loop.outputs.size() == 1);
  CHECK(p_loop.outputs[0].first == 3);

  // Consumer: a Κ-loop copy holding only BuildStep{2} and its output {4}, then
  // BuildStep{5}. Order among the surviving steps is preserved.
  REQUIRE(forked.consumer.size() == 2);
  auto const& c_loop = std::get<ScopeBlock>(forked.consumer[0].value);
  CHECK(c_loop.axis == ax);
  REQUIRE(c_loop.steps.size() == 1);
  CHECK(build_id(c_loop.steps[0]) == 2);
  REQUIRE(c_loop.outputs.size() == 1);
  CHECK(c_loop.outputs[0].first == 4);
  CHECK(build_id(forked.consumer[1]) == 5);
}

// A nested loop that lands entirely on one side is copied whole to that side
// and NOT emitted (empty) on the other.
TEST_CASE("fork_subchain drops the empty side of a one-sided nested loop",
          "[ordered-schedule][fork]") {
  Index const ax{L"i_1"};

  ScopeBlock inner;
  inner.axis = ax;
  inner.steps.push_back(Step{BuildStep{1}});
  inner.steps.push_back(Step{BuildStep{2}});
  inner.outputs.push_back({3, OutputKind::AccumulateSum});

  sequant::container::vector<Step> steps;
  steps.push_back(Step{std::move(inner)});

  // Whole nested loop is producer-side.
  std::function<bool(std::size_t)> const in_consumer = [](std::size_t) {
    return false;
  };

  auto const forked = sequant::eval::detail::fork_subchain(steps, in_consumer);

  REQUIRE(forked.producer.size() == 1);
  auto const& p_loop = std::get<ScopeBlock>(forked.producer[0].value);
  CHECK(p_loop.steps.size() == 2);
  CHECK(p_loop.outputs.size() == 1);
  CHECK(forked.consumer.empty());  // empty side dropped, no stranded output
}

// SP2 multi-level escape chain: a value that reduces an inner axis AND is
// carried on an outer one escapes at BOTH -- AccumulateSum at the inner block,
// AccumulateScatter at the outer block. well_formed accepts the same value_id
// escaping at two blocks WHEN they nest (inner is a descendant of outer).
TEST_CASE("well_formed accepts a nested multi-level escape chain",
          "[ordered-schedule][escape-chain]") {
  Index const outer{L"i_1"};
  Index const inner{L"i_2"};

  // value_id 2 escapes: AccumulateSum on the inner block's close (partial ->
  // accumulator), then AccumulateScatter on the outer block's close (-> full).
  ScopeBlock inner_block;
  inner_block.axis = inner;
  inner_block.ordinal = 0;
  inner_block.steps.push_back(Step{BuildStep{1}});  // per-iteration partial
  inner_block.outputs.push_back({2, OutputKind::AccumulateSum});

  ScopeBlock outer_block;
  outer_block.axis = outer;
  outer_block.ordinal = 0;
  outer_block.steps.push_back(Step{std::move(inner_block)});
  outer_block.outputs.push_back({2, OutputKind::AccumulateScatter});

  OrderedSchedule sched;
  sched.root.steps.push_back(Step{std::move(outer_block)});
  sched.num_values = 3;

  CHECK(well_formed(sched));
}

// The same value escaping at two UNRELATED (sibling) blocks is NOT a chain --
// it is duplicate production, rejected.
TEST_CASE("well_formed rejects the same escape in two sibling blocks",
          "[ordered-schedule][escape-chain]") {
  Index const ax{L"i_2"};

  ScopeBlock a;
  a.axis = ax;
  a.ordinal = 0;
  a.steps.push_back(Step{BuildStep{1}});
  a.outputs.push_back({2, OutputKind::AccumulateSum});

  ScopeBlock b;
  b.axis = ax;
  b.ordinal =
      1;  // distinct ordinal: passes the same-axis-sibling ordinal check
  b.steps.push_back(Step{BuildStep{3}});
  b.outputs.push_back({2, OutputKind::AccumulateScatter});  // same value_id 2

  OrderedSchedule sched;
  sched.root.steps.push_back(Step{std::move(a)});
  sched.root.steps.push_back(Step{std::move(b)});
  sched.num_values = 4;

  CHECK_FALSE(well_formed(sched));
}

TEST_CASE("well_formed rejects an out-of-range BuildStep::value_id",
          "[ordered-schedule]") {
  OrderedSchedule sched;
  sched.root.steps.push_back(Step{BuildStep{5}});
  sched.num_values = 1;  // value_id 5 is out of range

  CHECK_FALSE(well_formed(sched));
}

TEST_CASE(
    "well_formed rejects duplicate ordinals among same-axis sibling blocks",
    "[ordered-schedule]") {
  Index const i1{L"i_1"};
  Index const i2{L"i_2"};  // same TYPE ("i") as i1, different physical label

  ScopeBlock child_a;
  child_a.axis = i1;
  child_a.ordinal = 0;
  child_a.steps.push_back(Step{BuildStep{0}});

  ScopeBlock child_b;
  child_b.axis = i2;
  child_b.ordinal = 0;  // duplicate ordinal at the same axis TYPE as child_a
  child_b.steps.push_back(Step{BuildStep{1}});

  OrderedSchedule sched;
  sched.root.steps.push_back(Step{std::move(child_a)});
  sched.root.steps.push_back(Step{std::move(child_b)});
  sched.num_values = 2;

  CHECK_FALSE(well_formed(sched));
}

TEST_CASE("well_formed rejects an out-of-range output value_id",
          "[ordered-schedule]") {
  ScopeBlock child;
  child.axis = Index{L"i_1"};
  child.outputs.push_back({7, OutputKind::AccumulateScatter});

  OrderedSchedule sched;
  sched.root.steps.push_back(Step{std::move(child)});
  sched.num_values = 1;  // output value_id 7 is out of range

  CHECK_FALSE(well_formed(sched));
}

TEST_CASE(
    "well_formed rejects a value_id produced twice: once as a root "
    "BuildStep and again as a child block's AccumulateSum output",
    "[ordered-schedule]") {
  // Single-producer (SSA-like) invariant, checked WHOLE-SCHEDULE (not just
  // within one block): value_id 0 is built directly at the root AND is also
  // claimed as the accumulated output of an unrelated child loop -- two
  // production sites for the same value_id, which is never legal regardless
  // of how far apart in the tree they sit.
  ScopeBlock child;
  child.axis = Index{L"i_1"};
  child.outputs.push_back({0, OutputKind::AccumulateSum});

  OrderedSchedule sched;
  sched.root.steps.push_back(Step{BuildStep{0}});
  sched.root.steps.push_back(Step{std::move(child)});
  sched.num_values = 1;

  CHECK_FALSE(well_formed(sched));
}

// ===========================================================================
// Task 3: build_ordered_schedule, validated on the real water-20 CSV-CCSD
// doubles residual (DF/aux-only batching) -- the exact fixture test_legality.
// cpp's "classify_axis / analyze_legality: four-way axis classification on
// the water-20 aux-only residual" test already exercises (same recipe,
// duplicated under an `orderedsched_` prefix per that file's own convention:
// no shared test header exists for these DryRun fixtures, and same-named
// anonymous-namespace helpers would collide under CMake UNITY_BUILD grouping).
// ===========================================================================

namespace {

std::string orderedsched_witness_slurp(std::string const& path) {
  std::ifstream in(path);
  std::stringstream ss;
  ss << in.rdbuf();
  return ss.str();
}

struct OrderedSchedWater20ProblemSize {
  std::size_t mu_tilde;
  std::size_t aux;
  std::size_t i_occ;
  std::array<double, 5> pno_M;
  std::array<double, 5> osv_M;
};

inline constexpr OrderedSchedWater20ProblemSize kOrderedSchedWater20_pVDZF12{
    /*mu_tilde=*/896u,
    /*aux=*/1682u,
    /*i_occ=*/80u,
    /*pno_M=*/
    {1.0, 23.175775480059084, 25.865548281212597, 28.171416142614103,
     30.03848680550367},
    /*osv_M=*/
    {1.0, 58.987499999999997, 59.289227520688783, 59.584437469011633,
     59.872014818179686}};

sequant::eval::dryrun::SizeRegime orderedsched_witness_df_regime(
    OrderedSchedWater20ProblemSize const& p) {
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

sequant::ExprPtr orderedsched_witness_flatten_product(
    sequant::ExprPtr const& e) {
  if (!e->is<sequant::Product>()) return e;
  auto const& p = e->as<sequant::Product>();
  return sequant::ex<sequant::Product>(p.scalar(), p.factors(),
                                       sequant::Product::Flatten::Yes);
}

// Recursively count how many BuildStep/child-ScopeBlock steps sit inside
// `block`'s OWN steps list before the first step whose ScopeBlock axis TYPE
// is `axis_key` -- used to confirm relative ORDER (not just presence) among
// a block's steps.
std::optional<std::size_t> orderedsched_index_of_child_block(
    ScopeBlock const& block, std::wstring const& axis_key) {
  for (std::size_t i = 0; i < block.steps.size(); ++i) {
    if (auto const* child = std::get_if<ScopeBlock>(&block.steps[i].value))
      if (child->axis.space().base_key() == axis_key) return i;
  }
  return std::nullopt;
}

std::optional<std::size_t> orderedsched_index_of_build_step(
    ScopeBlock const& block, std::size_t value_id) {
  for (std::size_t i = 0; i < block.steps.size(); ++i) {
    if (auto const* b = std::get_if<BuildStep>(&block.steps[i].value))
      if (b->value_id == value_id) return i;
  }
  return std::nullopt;
}

}  // namespace

TEST_CASE(
    "build_ordered_schedule: water-20 aux-only residual places the "
    "Κ-contraction result as an AccumulateSum output of the {Κ} block, "
    "ordered before the root-level composite that reads it",
    "[ordered-schedule]") {
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
      orderedsched_witness_slurp(std::string(SEQUANT_UNIT_TESTS_SOURCE_DIR) +
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

  auto regime = orderedsched_witness_df_regime(kOrderedSchedWater20_pVDZF12);
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
        orderedsched_witness_flatten_product(summands[s]);
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
  auto rich = sequant::eval::compute_dag_boulevard(forest, *cm, block_of);
  REQUIRE(!rich.cells.empty());

  auto const legality = sequant::eval::analyze_legality(rich, forest, policy);
  REQUIRE(legality.cells.size() == rich.cells.size());

  auto const sched =
      sequant::eval::build_ordered_schedule(rich, legality, policy, {L"Κ"});
  REQUIRE(well_formed(sched));
  CHECK(sched.num_values == rich.cells.size());

  // The root's steps must contain exactly one {Κ} child ScopeBlock.
  auto const k_child_idx = orderedsched_index_of_child_block(sched.root, L"Κ");
  REQUIRE(k_child_idx.has_value());
  ScopeBlock const& k_block =
      std::get<ScopeBlock>(sched.root.steps[*k_child_idx].value);
  CHECK(k_block.axis.space().base_key() == L"Κ");

  // The Κ-contraction RESULT (a non-leaf value that does not itself carry Κ
  // but reduces it at its own node -- classify_axis Reduction, same target
  // identification as test_legality.cpp's water-20 test) must be an
  // AccumulateSum output of the {Κ} block, with NO BuildStep anywhere in the
  // whole schedule (single-producer: its only production site is that
  // output entry).
  auto const is_K = [](sequant::Index const& ix) {
    return ix.space().base_key() == L"Κ";
  };
  auto const carries_type =
      [&](sequant::container::svector<sequant::Index> const& v,
          auto const& pred) { return std::any_of(v.begin(), v.end(), pred); };
  auto const vmap = sequant::eval::build_value_node_map(forest);

  std::optional<std::size_t> mu_mu_value_id;
  std::optional<std::size_t> parent_value_id;  // I(i,i;a,a)-shaped consumer
  {
    std::optional<std::size_t> mu_mu_hash;
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
        auto const k_it =
            std::find_if(contracted.begin(), contracted.end(), is_K);
        if (k_it == contracted.end()) continue;
        if (require_mu_mu && !is_mu_mu_pair(vc.carried)) continue;
        mu_mu_hash = vc.hash;
        mu_mu_value_id = vc.value_id;
        break;
      }
    }
    REQUIRE(mu_mu_hash.has_value());

    // Its structural parent (the Κ-free composite consuming it) --
    // identical search to test_legality.cpp's Target 2.
    std::optional<Node> parent;
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
    REQUIRE(parent.has_value());
    auto const parent_hash = (*parent)->hash_value();
    auto const cell_it =
        std::find_if(rich.cells.begin(), rich.cells.end(),
                     [&](auto const& vc) { return vc.hash == parent_hash; });
    REQUIRE(cell_it != rich.cells.end());
    parent_value_id = cell_it->value_id;
  }
  REQUIRE(mu_mu_value_id.has_value());
  REQUIRE(parent_value_id.has_value());

  // AccumulateSum output of the {Κ} block, not a BuildStep anywhere.
  auto const out_it =
      std::find_if(k_block.outputs.begin(), k_block.outputs.end(),
                   [&](auto const& p) { return p.first == *mu_mu_value_id; });
  REQUIRE(out_it != k_block.outputs.end());
  CHECK(out_it->second == OutputKind::AccumulateSum);
  CHECK_FALSE(
      orderedsched_index_of_build_step(k_block, *mu_mu_value_id).has_value());
  CHECK_FALSE(orderedsched_index_of_build_step(sched.root, *mu_mu_value_id)
                  .has_value());

  // I(i,i;a,a) (the parent) is a root-level BuildStep, ordered AFTER the
  // {Κ} child block in the root's steps.
  auto const parent_idx =
      orderedsched_index_of_build_step(sched.root, *parent_value_id);
  REQUIRE(parent_idx.has_value());
  CHECK(*parent_idx > *k_child_idx);

  // At least one Κ-carrying LoopLocal intermediate (test_legality.cpp's
  // Target 3) is Transient: a BuildStep INSIDE the {Κ} block, with no
  // outputs entry anywhere (not an accumulate output).
  bool found_k_local_transient = false;
  for (auto const& cl : legality.cells) {
    bool const k_local = std::any_of(
        cl.per_axis.begin(), cl.per_axis.end(), [&](auto const& ac) {
          return is_K(ac.axis) && ac.role == sequant::eval::LoopRole::LoopLocal;
        });
    if (!k_local) continue;
    auto const vid_it =
        std::find_if(rich.cells.begin(), rich.cells.end(),
                     [&](auto const& vc) { return vc.hash == cl.hash; });
    REQUIRE(vid_it != rich.cells.end());
    std::size_t const vid = vid_it->value_id;
    if (orderedsched_index_of_build_step(k_block, vid).has_value()) {
      found_k_local_transient = true;
      CHECK_FALSE(std::any_of(k_block.outputs.begin(), k_block.outputs.end(),
                              [&](auto const& p) { return p.first == vid; }));
    }
  }
  REQUIRE(found_k_local_transient);
}

// ===========================================================================
// Fix round 1 (design review): a single scalar sort key can place a child
// block BEFORE every value that reads its output (the direction water-20
// above already exercises), but it has NO corresponding guarantee that the
// block sorts AFTER every value ITS OWN content reads as an input -- water-
// 20's {Κ} block happens to be leaf-only, so that direction was never
// stressed. This is a small HAND-BUILT fixture (mirrors test_scope_schedule.
// cpp's own scope_eval_tensor/scope_leaf/scope_inode helpers, prefixed
// distinctly to avoid a UNITY_BUILD anonymous-namespace collision) exercising
// BOTH directions in one forest:
//   R{a_1;a_2}          -- Κ-INDEPENDENT leaf: carries/contracts no Κ at its
//                          own node, so it is ROOT-homed even though it sits
//                          structurally under a realized Κ loop.
//   G{Κ_1;a_1}, H{Κ_1;a_3} -- Κ-carrying leaves, lockstep with the enclosing
//                          Κ loop -> LoopLocal (BuildStep's inside {Κ}).
//   V{Κ_1;a_2} = G * R  -- Κ-local (LoopLocal): DIRECTLY READS R, the
//                          root-homed leaf above -- the "block follows its
//                          own input" direction.
//   W{a_2;a_3} = V * H, with Κ realized AT W (batched_here) -- contracts Κ_1
//                          at its own node -> Reduction, AccumulateSum
//                          output of {Κ}.
//   Z{a_2;a_3}          -- unrelated Κ-free leaf, root-homed.
//   Top = W * Z          -- Κ-free at its own node (LoopInvariant) -> a
//                          root-level BuildStep that READS W, the {Κ}
//                          block's own output -- the "block precedes its
//                          consumer" direction (water-20's own shape).
// Expected root order: R (and G's/H's/Z's placement is unconstrained by any
// edge) before {Κ}, and {Κ} before Top.
// ===========================================================================

namespace {

sequant::EvalExpr orderedsched_eval_tensor(std::string_view tensor) {
  auto expr = sequant::deserialize<sequant::ExprPtr>(std::string(tensor));
  REQUIRE(static_cast<bool>(expr));
  return sequant::EvalExpr{expr->as<sequant::Tensor>()};
}

sequant::EvalNode<sequant::EvalExpr> orderedsched_leaf(
    std::string_view tensor) {
  return sequant::EvalNode<sequant::EvalExpr>{orderedsched_eval_tensor(tensor)};
}

// An internal node whose own result is `result`'s signature, formed as the
// PRODUCT of `l` and `r` (op_type stamped Product via EvalOpSetter, needed
// for contracted_indices()/is_product() -- test_scope_schedule.cpp's own
// scope_inode never needed this, since build_scope_schedule never consults
// contracted_indices, but build_site_of/analyze_legality do).
sequant::EvalNode<sequant::EvalExpr> orderedsched_inode(
    std::string_view result, sequant::EvalNode<sequant::EvalExpr> l,
    sequant::EvalNode<sequant::EvalExpr> r) {
  sequant::EvalExpr data = orderedsched_eval_tensor(result);
  sequant::EvalOpSetter{}.set(data, sequant::EvalOp::Product);
  return sequant::EvalNode<sequant::EvalExpr>{std::move(data), std::move(l),
                                              std::move(r)};
}

}  // namespace

TEST_CASE(
    "build_ordered_schedule: a root-homed leaf consumed INSIDE the {Κ} "
    "block sorts before it, and the {Κ} block sorts before the root-level "
    "composite that reads its AccumulateSum output -- both directions in "
    "one forest",
    "[ordered-schedule]") {
  auto ctx = sequant::get_default_context().clone();
  auto isr = ctx.mutable_index_space_registry();
  REQUIRE(isr != nullptr);
  sequant::mbpt::add_df_spaces(isr);  // Κ (DF aux)
  auto ctx_resetter = sequant::set_scoped_default_context(std::move(ctx));

  Index const K{L"Κ_1"};

  auto R = orderedsched_leaf("R{a_1;a_2}");
  auto G = orderedsched_leaf("G{Κ_1;a_1}");
  auto V = orderedsched_inode("V{Κ_1;a_2}", G, R);  // contracts a_1; reads R
  auto H = orderedsched_leaf("H{Κ_1;a_3}");
  auto W = orderedsched_inode("W{a_2;a_3}", V, H);  // contracts Κ_1
  W->set_batched_here({{K, sequant::BatchModeType::Contracted}});
  auto Z = orderedsched_leaf("Z{a_2;a_3}");
  auto Top = orderedsched_inode("Top{a_9;a_10}", W, Z);  // Κ-free; reads W

  sequant::BatchPolicy policy;
  policy.is_batchable_contracted_index = [](Index const& ix) {
    return ix.space().base_key() == L"Κ";
  };

  sequant::eval::dryrun::SizeRegime regime;
  regime.space_extent = {{L"a", 8u}, {L"Κ", 6u}};
  auto cm = std::make_shared<sequant::eval::dryrun::CostModel const>(regime);

  std::vector<sequant::EvalNode<sequant::EvalExpr>> forest{Top};
  auto const block_of = [](Index const&) -> std::size_t { return 4; };
  auto rich = sequant::eval::compute_dag_boulevard(forest, *cm, block_of);
  REQUIRE(!rich.cells.empty());

  auto const legality = sequant::eval::analyze_legality(rich, forest, policy);
  REQUIRE(legality.cells.size() == rich.cells.size());

  auto const sched =
      sequant::eval::build_ordered_schedule(rich, legality, policy, {L"Κ"});
  REQUIRE(well_formed(sched));

  auto const value_id_of = [&](std::string_view tensor) -> std::size_t {
    auto const hash = orderedsched_eval_tensor(tensor).hash_value();
    auto const it =
        std::find_if(rich.cells.begin(), rich.cells.end(),
                     [&](auto const& vc) { return vc.hash == hash; });
    REQUIRE(it != rich.cells.end());
    return it->value_id;
  };
  std::size_t const r_id = value_id_of("R{a_1;a_2}");
  std::size_t const w_id = value_id_of("W{a_2;a_3}");
  std::size_t const top_id = value_id_of("Top{a_9;a_10}");

  auto const k_child_idx = orderedsched_index_of_child_block(sched.root, L"Κ");
  REQUIRE(k_child_idx.has_value());
  ScopeBlock const& k_block =
      std::get<ScopeBlock>(sched.root.steps[*k_child_idx].value);
  CHECK(k_block.axis.space().base_key() == L"Κ");

  // W is the {Κ} block's AccumulateSum output; no BuildStep for it anywhere.
  auto const w_out_it =
      std::find_if(k_block.outputs.begin(), k_block.outputs.end(),
                   [&](auto const& p) { return p.first == w_id; });
  REQUIRE(w_out_it != k_block.outputs.end());
  CHECK(w_out_it->second == OutputKind::AccumulateSum);

  // R is a root-level BuildStep, positioned BEFORE the {Κ} block: V (inside
  // {Κ}) directly reads it, so the block must sort after it.
  auto const r_idx = orderedsched_index_of_build_step(sched.root, r_id);
  REQUIRE(r_idx.has_value());
  CHECK(*r_idx < *k_child_idx);

  // Top is a root-level BuildStep, positioned AFTER the {Κ} block: it reads
  // W, the block's own accumulated output.
  auto const top_idx = orderedsched_index_of_build_step(sched.root, top_id);
  REQUIRE(top_idx.has_value());
  CHECK(*top_idx > *k_child_idx);
}

// ===========================================================================
// Task 4: forced loop split. On the synthetic cross-iteration fixture
// B{i_3,i_4} = A{;i_3} * A{;i_4} (occ made batchable in the EXTERNAL role, no
// enclosing occ loop realized), every occ-carrying value is LoopCarried on occ
// (test_legality.cpp's own cross-iteration test pins this). The outer-product
// root B is a strict dependency-ancestor of the two loop-carried leaves it
// reads, so it lands in the CONSUMER pass while the leaves land in the PRODUCER
// pass: build_ordered_schedule must realize the occ loop as TWO ordered sibling
// blocks with distinct ordinals, producer (ordinal 0) before consumer (ordinal
// 1), each escaping its values via AccumulateScatter.
// ===========================================================================
TEST_CASE(
    "build_ordered_schedule: a forced-split occ axis realizes TWO ordered "
    "sibling blocks -- a producer pass (loop-carried operands scattered to "
    "full) before a consumer pass (the cross-iteration read) -- with distinct "
    "ordinals",
    "[ordered-schedule]") {
  using sequant::eval::dryrun::EvalExprDryRun;
  using sequant::eval::dryrun::EvalNodeDryRun;
  using Node = EvalNodeDryRun;

  auto const body =
      orderedsched_witness_slurp(std::string(SEQUANT_UNIT_TESTS_SOURCE_DIR) +
                                 "/data/legality_cross_iteration.txt");
  REQUIRE(!body.empty());
  std::string line = body;
  if (auto nl = line.find('\n'); nl != std::string::npos)
    line = line.substr(0, nl);
  auto expr = sequant::deserialize<sequant::ExprPtr>(line);
  REQUIRE(static_cast<bool>(expr));

  SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
  auto node = sequant::binarize<EvalExprDryRun>(expr);
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
  REQUIRE_FALSE(node.leaf());

  sequant::BatchPolicy policy;
  policy.is_batchable_external_index = [](sequant::Index const& ix) {
    return ix.space().base_key() == L"i";
  };

  sequant::eval::dryrun::SizeRegime regime;
  regime.space_extent = {{L"i", 8u}, {L"a", 16u}};
  auto cm = std::make_shared<sequant::eval::dryrun::CostModel const>(regime);

  std::vector<Node> forest{node};
  auto const block_of = [](sequant::Index const&) -> std::size_t { return 4; };
  auto rich = sequant::eval::compute_dag_boulevard(forest, *cm, block_of);
  REQUIRE(!rich.cells.empty());

  auto const legality = sequant::eval::analyze_legality(rich, forest, policy);

  // No LoopLocal value is read from both passes here (the shared operand is a
  // leaf, materialized regardless), so the demotion source reports nothing:
  // the split is driven purely by the loop-carried structure.
  CHECK(sequant::eval::forced_split_demotions(rich, legality).empty());

  auto const sched =
      sequant::eval::build_ordered_schedule(rich, legality, policy, {L"i"});
  REQUIRE(well_formed(sched));
  CHECK(sched.num_values == rich.cells.size());

  // The root's steps hold exactly TWO occ ("i") child blocks, ordinals 0 & 1.
  std::vector<ScopeBlock const*> occ_blocks;
  for (auto const& step : sched.root.steps)
    if (auto const* c = std::get_if<ScopeBlock>(&step.value))
      if (c->axis.space().base_key() == L"i") occ_blocks.push_back(c);
  REQUIRE(occ_blocks.size() == 2);
  std::vector<int> ordinals{occ_blocks[0]->ordinal, occ_blocks[1]->ordinal};
  std::sort(ordinals.begin(), ordinals.end());
  CHECK(ordinals == std::vector<int>{0, 1});

  ScopeBlock const* producer = nullptr;
  ScopeBlock const* consumer = nullptr;
  for (auto const* blk : occ_blocks)
    (blk->ordinal == 0 ? producer : consumer) = blk;
  REQUIRE(producer != nullptr);
  REQUIRE(consumer != nullptr);

  // The outer-product root B{i_3,i_4} (the node's own value) is the loop-
  // carried CONSUMER: an AccumulateScatter output of the ordinal-1 block, and a
  // BuildStep nowhere.
  auto const b_hash = node->hash_value();
  auto const b_it =
      std::find_if(rich.cells.begin(), rich.cells.end(),
                   [&](auto const& vc) { return vc.hash == b_hash; });
  REQUIRE(b_it != rich.cells.end());
  std::size_t const b_id = b_it->value_id;

  auto const has_scatter = [](ScopeBlock const& blk, std::size_t vid) {
    return std::any_of(
        blk.outputs.begin(), blk.outputs.end(), [&](auto const& p) {
          return p.first == vid && p.second == OutputKind::AccumulateScatter;
        });
  };
  CHECK(has_scatter(*consumer, b_id));
  CHECK_FALSE(has_scatter(*producer, b_id));
  CHECK_FALSE(orderedsched_index_of_build_step(*consumer, b_id).has_value());
  CHECK_FALSE(orderedsched_index_of_build_step(sched.root, b_id).has_value());

  // The producer pass carries the loop-carried leaf operand(s) B reads -- each
  // a scatter output there, none of them B.
  REQUIRE(!producer->outputs.empty());
  for (auto const& [vid, kind] : producer->outputs) {
    CHECK(kind == OutputKind::AccumulateScatter);
    CHECK(vid != b_id);
  }

  // Producer pass is ordered BEFORE the consumer pass in the root's steps.
  std::optional<std::size_t> prod_pos, cons_pos;
  for (std::size_t i = 0; i < sched.root.steps.size(); ++i)
    if (auto const* c = std::get_if<ScopeBlock>(&sched.root.steps[i].value))
      if (c->axis.space().base_key() == L"i")
        (c->ordinal == 0 ? prod_pos : cons_pos) = i;
  REQUIRE(prod_pos.has_value());
  REQUIRE(cons_pos.has_value());
  CHECK(*prod_pos < *cons_pos);
}

// SP2 non-innermost forced split (phase 3 gate): a 2-axis fixture with occ
// OUTER and aux INNER. B{;i_3,i_4} = A{;i_3} * A{;i_4} forces the occ split
// (the outer product reads each A across occ-blocks); each A{;i} is itself
// formed by an aux (Κ) contraction, so it is LoopCarried on occ AND Reduction
// on aux -- the multi-level escape (aux sum, occ scatter) plus the
// non-innermost split. Axes are realized by hand-stamping (set_sliced_modes for
// the external occ, set_batched_here Contracted for aux), the same way the
// aux-only fixtures above stamp Κ -- no optimize() run.
namespace {
// Stamp occ as EXTERNAL and aux (Κ) as Contracted in batched_here, sourcing the
// Index identities from the node's OWN canon/contracted indices. sliced_modes
// is then DERIVED by stamp_lifetime_masks (the cross-occurrence meet), never
// hand-set -- that is how the real optimize()->binarize path realizes an axis.
void orderedsched_stamp_2axis(sequant::EvalNode<sequant::EvalExpr>& n) {
  using sequant::BatchModeType;
  sequant::container::svector<std::pair<Index, BatchModeType>> stamps;
  for (auto const& ix : n->canon_indices())
    if (ix.space().base_key() == L"i")
      stamps.push_back({ix, BatchModeType::External});
  for (auto const& ix : sequant::contracted_indices(n))
    if (ix.space().base_key() == L"Κ")
      stamps.push_back({ix, BatchModeType::Contracted});
  if (!stamps.empty()) n->set_batched_here(stamps);
}

// Post-order walk stamping every node.
void orderedsched_stamp_all(sequant::EvalNode<sequant::EvalExpr>& n) {
  if (!n.leaf()) {
    orderedsched_stamp_all(n.left());
    orderedsched_stamp_all(n.right());
  }
  orderedsched_stamp_2axis(n);
}

sequant::EvalNode<sequant::EvalExpr> orderedsched_2axis_forest_root() {
  auto P3 = orderedsched_leaf("P{Κ_1;i_3}");
  auto Q1 = orderedsched_leaf("Q{;Κ_1}");
  auto A3 =
      orderedsched_inode("A{;i_3}", P3, Q1);  // contracts Κ_1, carries i_3

  auto P4 = orderedsched_leaf("P{Κ_2;i_4}");
  auto Q2 = orderedsched_leaf("Q{;Κ_2}");
  auto A4 =
      orderedsched_inode("A{;i_4}", P4, Q2);  // contracts Κ_2, carries i_4

  auto B = orderedsched_inode("B{;i_3,i_4}", A3, A4);  // outer product on occ
  orderedsched_stamp_all(B);
  return B;
}
}  // namespace

TEST_CASE(
    "build_ordered_schedule: a 2-axis occ-outer/aux-inner term realizes occ as "
    "the outer forced split (phase-3 gate)",
    "[ordered-schedule][sp2-noninner]") {
  auto ctx = sequant::get_default_context().clone();
  auto isr = ctx.mutable_index_space_registry();
  REQUIRE(isr != nullptr);
  sequant::mbpt::add_df_spaces(isr);  // Κ (DF aux)
  auto ctx_resetter = sequant::set_scoped_default_context(std::move(ctx));

  auto B = orderedsched_2axis_forest_root();

  sequant::BatchPolicy policy;
  policy.is_batchable_contracted_index = [](Index const& ix) {
    return ix.space().base_key() == L"Κ";
  };
  policy.is_batchable_external_index = [](Index const& ix) {
    return ix.space().base_key() == L"i";
  };
  policy.batch_spectator_indices = true;
  policy.node_level_placement = true;

  sequant::eval::dryrun::SizeRegime regime;
  regime.space_extent = {{L"i", 8u}, {L"Κ", 6u}};
  auto cm = std::make_shared<sequant::eval::dryrun::CostModel const>(regime);

  std::vector<sequant::EvalNode<sequant::EvalExpr>> forest{B};
  // Derive sliced_modes from the batched_here stamps (the cross-occurrence
  // meet), realizing occ (External) + aux (Contracted) as loop axes.
  sequant::stamp_lifetime_masks(forest);
  auto const block_of = [](Index const&) -> std::size_t { return 4; };
  auto rich = sequant::eval::compute_dag_boulevard(forest, *cm, block_of);
  REQUIRE(!rich.cells.empty());

  auto const legality = sequant::eval::analyze_legality(rich, forest, policy);
  REQUIRE(legality.cells.size() == rich.cells.size());

  // Realization precondition: occ is realized as a LoopCarried axis (which
  // forces a split) and aux as a Reduction axis -- both must be present for the
  // 2-axis occ-outer/aux-inner nest.
  bool occ_carried = false, aux_reduction = false;
  for (auto const& cl : legality.cells)
    for (auto const& ac : cl.per_axis) {
      if (ac.axis.space().base_key() == L"i" &&
          ac.role == sequant::eval::LoopRole::LoopCarried)
        occ_carried = true;
      if (ac.axis.space().base_key() == L"Κ" &&
          ac.role == sequant::eval::LoopRole::Reduction)
        aux_reduction = true;
    }
  CHECK(occ_carried);
  CHECK(aux_reduction);

  // The generalized detection must NOT spuriously assert on an occ-OUTER axis
  // (the old code did); it builds a correct occ-outer / aux-inner NESTED
  // schedule. (This forest's single occ-carried value is the forest root, read
  // by nothing, so no cross-iteration read forces a two-pass split -- the split
  // EMISSION is exercised at C60 aux_occ, where the forcing is real. Here we
  // pin the nesting + the multi-level escape placement.)
  auto const sched =
      sequant::eval::build_ordered_schedule(rich, legality, policy, {});
  REQUIRE(well_formed(sched));
  CHECK(sched.num_values == rich.cells.size());

  // Exactly ONE occ ("i") block at root, enclosing an aux ("Κ") block:
  // occ-outer, aux-inner.
  std::vector<ScopeBlock const*> occ_blocks;
  for (auto const& step : sched.root.steps)
    if (auto const* c = std::get_if<ScopeBlock>(&step.value))
      if (c->axis.space().base_key() == L"i") occ_blocks.push_back(c);
  REQUIRE(occ_blocks.size() == 1);
  ScopeBlock const& occ = *occ_blocks.front();
  bool aux_inside = false;
  for (auto const& s : occ.steps)
    if (auto const* c = std::get_if<ScopeBlock>(&s.value))
      if (c->axis.space().base_key() == L"Κ") aux_inside = true;
  CHECK(aux_inside);

  // The outer-product root B (carried on occ) is an AccumulateScatter output of
  // the occ block; the aux-contracted A is an AccumulateSum output of the aux
  // block (the reduction escapes inner, the carry escapes outer).
  auto const b_hash = B->hash_value();
  auto const b_it =
      std::find_if(rich.cells.begin(), rich.cells.end(),
                   [&](auto const& vc) { return vc.hash == b_hash; });
  REQUIRE(b_it != rich.cells.end());
  std::size_t const b_id = b_it->value_id;
  bool b_scatter_at_occ = false;
  for (auto const& [v, k] : occ.outputs)
    if (v == b_id && k == OutputKind::AccumulateScatter)
      b_scatter_at_occ = true;
  CHECK(b_scatter_at_occ);
  CHECK_FALSE(orderedsched_index_of_build_step(sched.root, b_id).has_value());
}

// Task 4 (DAG-scope runtime slicing plan): populate_cell_mode_to_level.
//
// Extends the 2-axis occ-outer/aux-inner [sp2-noninner] fixture (B{;i_3,i_4},
// same setup as the TEST_CASEs above) with two more forest roots, C and Z,
// each a trivial LOCAL wrapper around an intermediate (M3, Y resp.) that
// itself carries exactly ONE occ mode and has no aux dependence at all --
// M3{;i_3} = M3(La{;i_3}, Lb{;i_3}) uses the SAME physical index (i_3)
// build_ordered_schedule's canonical chain picks as the occ block's own
// representative axis (see the "occ-outer/aux-inner" TEST_CASE above, whose
// dump confirms this); Y{;i_5} = Y(Ya{;i_5}, Yb{;i_5}) uses a DIFFERENT
// physical occ index. Wrapping each in a further root (C{;i_3} = C(M3,
// Lc{;i_3}); Z{;i_5} = Z(Y, Yd{;i_5})) matters: a forest ROOT is always
// forced to materialize in full (an AccumulateScatter escape, never a plain
// BuildStep -- a root without a further consumer, like B in the
// "occ-outer/aux-inner" TEST_CASE, is exactly why THAT fixture alone has no
// BuildStep to test against), so M3 and Y -- each read exactly once, by C/Z
// resp., within the SAME occ iteration -- stay genuine per-iteration
// Transients: plain BuildSteps homed directly inside the occ ScopeBlock.
//
// This reproduces the plan's motivating water-8 DF-leaf structure (a value
// enclosed by an occ-type loop realized via a DIFFERENT physical index than
// the one it carries) at [sp2-noninner]'s own occ+aux setup: TWO cells (M3,
// Y) enclosed by the EXACT SAME ScopeBlock, one whose carried mode IS that
// block's representative Index, one whose carried mode is a different
// physical index of the same axis TYPE.
TEST_CASE(
    "populate_cell_mode_to_level: a BuildStep cell maps its own carried mode "
    "to the enclosing loop's DagScopeLevel; a DIFFERENT physical occ index "
    "enclosed by the SAME loop gets no level (Task 4, DAG-scope runtime "
    "slicing)",
    "[ordered-schedule][sp2-noninner]") {
  auto ctx = sequant::get_default_context().clone();
  auto isr = ctx.mutable_index_space_registry();
  REQUIRE(isr != nullptr);
  sequant::mbpt::add_df_spaces(isr);  // Κ (DF aux)
  auto ctx_resetter = sequant::set_scoped_default_context(std::move(ctx));

  auto B = orderedsched_2axis_forest_root();

  auto La = orderedsched_leaf("La{;i_3}");
  auto Lb = orderedsched_leaf("Lb{;i_3}");
  auto M3 = orderedsched_inode("M3{;i_3}", La, Lb);
  auto Lc = orderedsched_leaf("Lc{;i_3}");
  auto C = orderedsched_inode("C{;i_3}", M3, Lc);
  orderedsched_stamp_all(C);

  auto Ya = orderedsched_leaf("Ya{;i_5}");
  auto Yb = orderedsched_leaf("Yb{;i_5}");
  auto Y = orderedsched_inode("Y{;i_5}", Ya, Yb);
  auto Yd = orderedsched_leaf("Yd{;i_5}");
  auto Z = orderedsched_inode("Z{;i_5}", Y, Yd);
  orderedsched_stamp_all(Z);

  sequant::BatchPolicy policy;
  policy.is_batchable_contracted_index = [](Index const& ix) {
    return ix.space().base_key() == L"Κ";
  };
  policy.is_batchable_external_index = [](Index const& ix) {
    return ix.space().base_key() == L"i";
  };
  policy.batch_spectator_indices = true;
  policy.node_level_placement = true;

  sequant::eval::dryrun::SizeRegime regime;
  regime.space_extent = {{L"i", 8u}, {L"Κ", 6u}};
  auto cm = std::make_shared<sequant::eval::dryrun::CostModel const>(regime);

  std::vector<sequant::EvalNode<sequant::EvalExpr>> forest{B, C, Z};
  sequant::stamp_lifetime_masks(forest);
  auto const block_of = [](Index const&) -> std::size_t { return 4; };
  // NOT const: populate_cell_mode_to_level writes each cell's mode_to_level
  // back into rich, mirroring scope_executor.hpp's own wiring.
  auto rich = sequant::eval::compute_dag_boulevard(forest, *cm, block_of);
  auto const legality = sequant::eval::analyze_legality(rich, forest, policy);
  auto const sched =
      sequant::eval::build_ordered_schedule(rich, legality, policy, {});
  REQUIRE(well_formed(sched));

  // Exactly ONE occ ("i") block at root (same shape as the
  // "occ-outer/aux-inner" TEST_CASE above).
  std::vector<ScopeBlock const*> occ_blocks;
  for (auto const& step : sched.root.steps)
    if (auto const* c = std::get_if<ScopeBlock>(&step.value))
      if (c->axis.space().base_key() == L"i") occ_blocks.push_back(c);
  REQUIRE(occ_blocks.size() == 1);
  ScopeBlock const& occ = *occ_blocks.front();
  auto const occ_level = occ.level;

  auto const cell_of = [&](sequant::EvalNode<sequant::EvalExpr> const& n)
      -> sequant::eval::ValueCell const& {
    auto const h = n->hash_value();
    auto const it = std::find_if(rich.cells.begin(), rich.cells.end(),
                                 [&](auto const& vc) { return vc.hash == h; });
    REQUIRE(it != rich.cells.end());
    return *it;
  };
  sequant::eval::ValueCell const& m3_cell = cell_of(M3);
  sequant::eval::ValueCell const& y_cell = cell_of(Y);

  // Precondition: M3 and Y are each a plain BuildStep directly inside `occ`
  // (not an escape output, not root-level) -- i.e. this test actually
  // exercises the BuildStep-only path populate_cell_mode_to_level consults.
  REQUIRE(orderedsched_index_of_build_step(occ, m3_cell.value_id).has_value());
  REQUIRE(orderedsched_index_of_build_step(occ, y_cell.value_id).has_value());
  REQUIRE(m3_cell.carried == sequant::container::svector<Index>{Index{L"i_3"}});
  REQUIRE(y_cell.carried == sequant::container::svector<Index>{Index{L"i_5"}});

  // RED (pre-implementation): mode_to_level is default-empty on every cell
  // until populate_cell_mode_to_level runs.
  CHECK(m3_cell.mode_to_level.by_mode.empty());
  CHECK(y_cell.mode_to_level.by_mode.empty());

  sequant::eval::populate_cell_mode_to_level(sched, rich);

  // GREEN: M3 carries i_3 -- the SAME physical Index the occ block uses as
  // its own representative axis -- so mode 0 (its only mode) maps to the occ
  // block's DagScopeLevel.
  REQUIRE(m3_cell.mode_to_level.mode_of(occ_level).has_value());
  CHECK(*m3_cell.mode_to_level.mode_of(occ_level) == 0u);

  // Y carries i_5 -- a DIFFERENT physical occ index -- even though Y is a
  // BuildStep directly inside the EXACT SAME occ ScopeBlock as M3, it gets NO
  // level for it: the exact fact whose absence caused the over-slice bug
  // (the plan's water-8 DF-leaf example: carries i_1, but no level for an
  // enclosing i_2 loop it does not carry).
  CHECK_FALSE(y_cell.mode_to_level.mode_of(occ_level).has_value());
}

namespace {
///
/// \brief Task 3 (SP3): recurse through every non-root \c ScopeBlock reachable
/// from \p steps (which start at \p depth, the nesting depth of \p steps'
/// OWN blocks -- 1 for the root's direct children, 2 for their children,
/// etc.), asserting each block's \c level mirrors its \c axis/\c ordinal at
/// the correct nesting depth.
///
void orderedsched_check_levels(
    sequant::container::vector<sequant::eval::Step> const& steps,
    std::size_t depth) {
  for (auto const& step : steps) {
    auto const* block = std::get_if<ScopeBlock>(&step.value);
    if (!block) continue;
    CHECK(block->level.space == block->axis.space().base_key());
    CHECK(block->level.ordinal == block->ordinal);
    CHECK(block->level.depth == depth);
    orderedsched_check_levels(block->steps, depth + 1);
  }
}
}  // namespace

TEST_CASE(
    "build_ordered_schedule: ScopeBlock::level mirrors axis/ordinal at the "
    "correct nesting depth (Task 3, DAG-scope runtime slicing)",
    "[ordered-schedule][sp2-noninner]") {
  auto ctx = sequant::get_default_context().clone();
  auto isr = ctx.mutable_index_space_registry();
  REQUIRE(isr != nullptr);
  sequant::mbpt::add_df_spaces(isr);  // Κ (DF aux)
  auto ctx_resetter = sequant::set_scoped_default_context(std::move(ctx));

  auto B = orderedsched_2axis_forest_root();

  sequant::BatchPolicy policy;
  policy.is_batchable_contracted_index = [](Index const& ix) {
    return ix.space().base_key() == L"Κ";
  };
  policy.is_batchable_external_index = [](Index const& ix) {
    return ix.space().base_key() == L"i";
  };
  policy.batch_spectator_indices = true;
  policy.node_level_placement = true;

  sequant::eval::dryrun::SizeRegime regime;
  regime.space_extent = {{L"i", 8u}, {L"Κ", 6u}};
  auto cm = std::make_shared<sequant::eval::dryrun::CostModel const>(regime);

  std::vector<sequant::EvalNode<sequant::EvalExpr>> forest{B};
  sequant::stamp_lifetime_masks(forest);
  auto const block_of = [](Index const&) -> std::size_t { return 4; };
  auto rich = sequant::eval::compute_dag_boulevard(forest, *cm, block_of);
  REQUIRE(!rich.cells.empty());

  auto const legality = sequant::eval::analyze_legality(rich, forest, policy);
  REQUIRE(legality.cells.size() == rich.cells.size());

  auto const sched =
      sequant::eval::build_ordered_schedule(rich, legality, policy, {});
  REQUIRE(well_formed(sched));

  // The root block itself is the sentinel (default level{}) -- only its
  // reachable non-root descendants (starting at depth 1) are asserted.
  orderedsched_check_levels(sched.root.steps, /*depth=*/1);

  // Sanity: this fixture actually exercises TWO nesting depths (occ outer at
  // depth 1, aux inner at depth 2) -- confirm the recursion isn't vacuous.
  bool saw_depth_1 = false, saw_depth_2 = false;
  std::function<void(sequant::container::vector<sequant::eval::Step> const&)>
      scan = [&](sequant::container::vector<sequant::eval::Step> const& steps) {
        for (auto const& step : steps) {
          auto const* block = std::get_if<ScopeBlock>(&step.value);
          if (!block) continue;
          if (block->level.depth == 1) saw_depth_1 = true;
          if (block->level.depth == 2) saw_depth_2 = true;
          scan(block->steps);
        }
      };
  scan(sched.root.steps);
  CHECK(saw_depth_1);
  CHECK(saw_depth_2);
}

// ===========================================================================
// Task 4 (Fix round 1): the demotion trigger in forced_split_demotions is
// SINGLE-SIDED -- a producer-homed LoopLocal value V (V itself NOT in
// consumer_pass) is demoted the moment ANY direct consumer of V lands in the
// consumer pass, NOT only when V ALSO has a producer-side consumer. A hand-
// built dep graph + legality isolates the four cases (the split-pass BFS is
// one-directional, so a consumer-homed value never leaks back to the producer
// side and needs no demotion -- hence the asymmetric "sole consumer" case is
// the only unsound one, and it is exactly the one the old `&&` trigger missed):
//
//   Lc (id0)  -- LoopCarried leaf on occ  (forces the split)
//   Wc (id5)  -- LoopCarried, reads Lc    => a strict ancestor of Lc, so it is
//                                            in consumer_pass
//   Pu (id4)  -- producer-side root, reads Vboth+Vprod (reads NO carried value,
//                                            so NOT in consumer_pass)
//   Vboth(id1)-- LoopLocal, read by Pu (producer) AND Wc (consumer) => flagged
//   Vsole(id2)-- LoopLocal, read ONLY by Wc (consumer) => flagged (the case the
//                                            old `&&` trigger silently dropped)
//   Vprod(id3)-- LoopLocal, read ONLY by Pu (producer) => NOT flagged
// ===========================================================================
TEST_CASE(
    "forced_split_demotions: single-sided trigger flags a producer-homed "
    "LoopLocal read by any consumer-pass value, including its SOLE consumer",
    "[ordered-schedule]") {
  sequant::Index const i{L"i_1"};

  auto const make_cell =
      [](std::size_t id, std::size_t hash,
         std::vector<std::pair<std::size_t, std::size_t>> const& occs) {
        sequant::eval::ValueCell vc{};
        vc.value_id = id;
        vc.hash = hash;
        vc.first_use = 0;
        vc.last_use = 0;
        for (auto const& [p, cp] : occs) {
          sequant::eval::OccurrenceRec o{};
          o.point = p;
          o.consumer_point = cp;
          vc.occurrences.push_back(std::move(o));
        }
        return vc;
      };

  // consumer_point == point means "forest root" (no structural parent).
  sequant::eval::RichSchedule rich;
  rich.cells.push_back(make_cell(0, 1000, {{0, 50}}));  // Lc  -> Wc
  rich.cells.push_back(
      make_cell(1, 1001, {{10, 40}, {11, 50}}));         // Vboth->Pu,Wc
  rich.cells.push_back(make_cell(2, 1002, {{20, 50}}));  // Vsole-> Wc
  rich.cells.push_back(make_cell(3, 1003, {{30, 40}}));  // Vprod-> Pu
  rich.cells.push_back(make_cell(4, 1004, {{40, 40}}));  // Pu  (root)
  rich.cells.push_back(make_cell(5, 1005, {{50, 50}}));  // Wc  (root)

  auto const make_legality = [&](std::size_t hash,
                                 sequant::eval::LoopRole role) {
    sequant::eval::CellLegality cl;
    cl.hash = hash;
    sequant::eval::AxisClass ac;
    ac.axis = i;
    ac.role = role;
    cl.per_axis.push_back(ac);
    if (role == sequant::eval::LoopRole::LoopCarried)
      cl.forced_split_axes.push_back(i);
    return cl;
  };
  sequant::eval::LegalitySchedule legality;
  legality.cells.push_back(
      make_legality(1000, sequant::eval::LoopRole::LoopCarried));
  legality.cells.push_back(
      make_legality(1001, sequant::eval::LoopRole::LoopLocal));
  legality.cells.push_back(
      make_legality(1002, sequant::eval::LoopRole::LoopLocal));
  legality.cells.push_back(
      make_legality(1003, sequant::eval::LoopRole::LoopLocal));
  legality.cells.push_back(
      make_legality(1004, sequant::eval::LoopRole::LoopLocal));
  legality.cells.push_back(
      make_legality(1005, sequant::eval::LoopRole::LoopCarried));

  auto const dem = sequant::eval::forced_split_demotions(rich, legality);
  auto const flagged = [&](std::size_t hash) {
    return std::any_of(dem.begin(), dem.end(), [&](auto const& p) {
      return p.first == hash && p.second == L"i";
    });
  };

  CHECK(flagged(1001));        // Vboth: producer + consumer reader
  CHECK(flagged(1002));        // Vsole: SOLE consumer is consumer-pass
  CHECK_FALSE(flagged(1003));  // Vprod: only a producer-side consumer
  CHECK_FALSE(flagged(1004));  // Pu: no consumers at all
  CHECK_FALSE(flagged(1000));  // Lc: carried, not a LoopLocal candidate
  CHECK_FALSE(flagged(1005));  // Wc: carried / in consumer_pass
  CHECK(dem.size() == 2);
}

// ===========================================================================
// Task 5: acceptance + executor-shape validation.
//
// Both TEST_CASEs below reuse the SAME two real fixtures the Task 3 and
// Task 4 tests above already build (water-20's aux-only residual and the
// cross-iteration forced-split fixture), factored into two small builder
// functions so the acceptance and executor-shape checks below run against
// literally the same data rather than a re-derived copy. The existing Task
// 3/4 TEST_CASEs above are left untouched (their own inline setup is not
// replaced) to avoid disturbing already-pinned behavior; these builders are
// net-new, consumed only by the two TEST_CASEs that follow them.
// ===========================================================================

namespace {

struct OrderedSchedFixture {
  std::vector<sequant::eval::dryrun::EvalNodeDryRun> forest;
  sequant::eval::RichSchedule rich;
  sequant::eval::LegalitySchedule legality;
  OrderedSchedule sched;
};

// Same recipe as the "water-20 aux-only residual places the Κ-contraction
// result..." TEST_CASE above (lines ~231-333), stopping once the schedule is
// built (no further per-value inspection here -- that is this function's
// caller's job).
OrderedSchedFixture orderedsched_water20_fixture() {
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
      orderedsched_witness_slurp(std::string(SEQUANT_UNIT_TESTS_SOURCE_DIR) +
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

  auto regime = orderedsched_witness_df_regime(kOrderedSchedWater20_pVDZF12);
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
        orderedsched_witness_flatten_product(summands[s]);
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
  auto rich = sequant::eval::compute_dag_boulevard(forest, *cm, block_of);
  REQUIRE(!rich.cells.empty());

  auto legality = sequant::eval::analyze_legality(rich, forest, policy);
  REQUIRE(legality.cells.size() == rich.cells.size());

  auto sched =
      sequant::eval::build_ordered_schedule(rich, legality, policy, {L"Κ"});
  REQUIRE(well_formed(sched));

  return OrderedSchedFixture{std::move(forest), std::move(rich),
                             std::move(legality), std::move(sched)};
}

// Same recipe as the "forced-split occ axis realizes TWO ordered sibling
// blocks..." TEST_CASE above (lines ~608-657), stopping once the schedule is
// built.
OrderedSchedFixture orderedsched_cross_iteration_fixture() {
  using sequant::eval::dryrun::EvalExprDryRun;
  using sequant::eval::dryrun::EvalNodeDryRun;
  using Node = EvalNodeDryRun;

  auto const body =
      orderedsched_witness_slurp(std::string(SEQUANT_UNIT_TESTS_SOURCE_DIR) +
                                 "/data/legality_cross_iteration.txt");
  REQUIRE(!body.empty());
  std::string line = body;
  if (auto nl = line.find('\n'); nl != std::string::npos)
    line = line.substr(0, nl);
  auto expr = sequant::deserialize<sequant::ExprPtr>(line);
  REQUIRE(static_cast<bool>(expr));

  SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
  auto node = sequant::binarize<EvalExprDryRun>(expr);
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
  REQUIRE_FALSE(node.leaf());

  sequant::BatchPolicy policy;
  policy.is_batchable_external_index = [](sequant::Index const& ix) {
    return ix.space().base_key() == L"i";
  };

  sequant::eval::dryrun::SizeRegime regime;
  regime.space_extent = {{L"i", 8u}, {L"a", 16u}};
  auto cm = std::make_shared<sequant::eval::dryrun::CostModel const>(regime);

  std::vector<Node> forest{node};
  auto const block_of = [](sequant::Index const&) -> std::size_t { return 4; };
  auto rich = sequant::eval::compute_dag_boulevard(forest, *cm, block_of);
  REQUIRE(!rich.cells.empty());

  auto legality = sequant::eval::analyze_legality(rich, forest, policy);

  auto sched =
      sequant::eval::build_ordered_schedule(rich, legality, policy, {L"i"});
  REQUIRE(well_formed(sched));

  return OrderedSchedFixture{std::move(forest), std::move(rich),
                             std::move(legality), std::move(sched)};
}

///
/// \brief Collects every value_id PRODUCED anywhere in \p block's subtree,
/// tagged by its production ROLE: \c std::nullopt for a plain \c BuildStep
/// (Transient, per \c build_ordered_schedule's own doc comment -- "Transient
/// is realized as produced-by-BuildStep-and-nothing-else", not as an explicit
/// \c outputs entry), or the stored \c OutputKind for a block \c outputs
/// entry. An independent, test-side mirror of \c
/// detail::collect_production_ids (kept internal to ordered_schedule.hpp),
/// walking only the PUBLIC IR (\c Step / \c BuildStep / \c ScopeBlock /
/// outputs) -- used below to hard-assert COMPLETENESS (every value_id in
/// [0, num_values) is produced at least once), which \c well_formed's own
/// single-producer check does not assert (it only rules out duplicates, not
/// gaps).
///
void orderedsched_collect_productions(
    ScopeBlock const& block,
    std::vector<std::pair<std::size_t, std::optional<OutputKind>>>& out) {
  for (Step const& step : block.steps) {
    if (auto const* build = std::get_if<BuildStep>(&step.value)) {
      out.push_back({build->value_id, std::nullopt});
    } else {
      orderedsched_collect_productions(std::get<ScopeBlock>(step.value), out);
    }
  }
  for (auto const& [vid, kind] : block.outputs) out.push_back({vid, kind});
}

}  // namespace

TEST_CASE(
    "build_ordered_schedule: acceptance -- water-20 + cross-iteration "
    "OrderedSchedules are well_formed with EVERY value_id produced exactly "
    "once, and all four production-site roles (Transient, AccumulateSum, "
    "AccumulateScatter producer pass, AccumulateScatter consumer pass) "
    "reachable across these two inputs are demonstrably exercised",
    "[ordered-schedule]") {
  auto const water20 = orderedsched_water20_fixture();
  auto const cross = orderedsched_cross_iteration_fixture();

  // well_formed + COMPLETENESS: every value_id in [0, num_values) is
  // produced EXACTLY once (well_formed itself only rules out duplicates; the
  // exact-match against [0, num_values) below additionally asserts no gaps).
  for (auto const* fx : {&water20, &cross}) {
    CHECK(well_formed(fx->sched));

    std::vector<std::pair<std::size_t, std::optional<OutputKind>>> prods;
    orderedsched_collect_productions(fx->sched.root, prods);

    std::vector<std::size_t> ids;
    ids.reserve(prods.size());
    for (auto const& [vid, kind] : prods) {
      (void)kind;
      ids.push_back(vid);
    }
    std::sort(ids.begin(), ids.end());
    std::vector<std::size_t> expected(fx->sched.num_values);
    std::iota(expected.begin(), expected.end(), std::size_t{0});
    CHECK(ids == expected);
  }

  // Transient + AccumulateSum (water-20): a plain BuildStep with no outputs
  // entry anywhere (Transient -- individually pinned for the Κ-local
  // intermediate and the I(i,i;a,a) root composite by the water-20 TEST_CASE
  // above; reconfirmed here from the SAME fixture data), and the {Κ} block's
  // own Κ-contraction-result AccumulateSum output (also individually pinned
  // above).
  {
    std::vector<std::pair<std::size_t, std::optional<OutputKind>>> prods;
    orderedsched_collect_productions(water20.sched.root, prods);
    CHECK(std::any_of(prods.begin(), prods.end(),
                      [](auto const& p) { return !p.second.has_value(); }));
    CHECK(std::any_of(prods.begin(), prods.end(), [](auto const& p) {
      return p.second == OutputKind::AccumulateSum;
    }));
  }

  // AccumulateScatter, both the PRODUCER-pass and CONSUMER-pass roles
  // (cross-iteration): the split structure itself (two sibling {i} blocks
  // with distinct ordinals -- individually pinned by the Task 4 TEST_CASE
  // above; reconfirmed here), each escaping its values via AccumulateScatter.
  {
    std::vector<ScopeBlock const*> occ_blocks;
    for (auto const& step : cross.sched.root.steps)
      if (auto const* c = std::get_if<ScopeBlock>(&step.value))
        if (c->axis.space().base_key() == L"i") occ_blocks.push_back(c);
    REQUIRE(occ_blocks.size() == 2);
    std::vector<int> ordinals{occ_blocks[0]->ordinal, occ_blocks[1]->ordinal};
    std::sort(ordinals.begin(), ordinals.end());
    CHECK(ordinals == std::vector<int>{0, 1});

    for (auto const* blk : occ_blocks) {
      REQUIRE(!blk->outputs.empty());
      CHECK(std::all_of(blk->outputs.begin(), blk->outputs.end(),
                        [](auto const& p) {
                          return p.second == OutputKind::AccumulateScatter;
                        }));
    }
  }
}

// ===========================================================================
// Task 5: executor-shape validation -- documents exactly what SP3's executor
// will read off the OrderedSchedule IR. \c scope_executor.hpp's \c walk_scope
// (see its own doc comment, "the value_id -> forest-node bridge") reads a \c
// ScopeNode's \c mode / \c kind / \c homed_values plus EXACTLY ONE child (\c
// node.children.front(), even though \c ScopeNode::children is a vector --
// there is today no consumer of more than one sibling child block at a
// level). \c OrderedSchedule's \c ScopeBlock generalizes this: \c axis / \c
// kind carry the same meaning, \c outputs replaces \c homed_values' implicit
// "built here, never leaves" with an explicit value_id -> \c OutputKind map,
// and -- the piece \c scope_executor.hpp's single-child pattern lacks -- \c
// steps interleaves build steps with an ORDERED LIST of sibling child \c
// ScopeBlocks (Task 4's forced-split producer/consumer passes are two such
// siblings at ONE level). This test asserts that shape is actually present
// and walkable on both real fixtures, and that value_ids resolve through
// \c rich.cells[...].hash exactly as \c scope_executor.hpp's \c
// build_value_node_map bridge already does elsewhere. A structural assertion
// only -- SP3 will consume this; SP2 does not execute anything.
// ===========================================================================
TEST_CASE(
    "build_ordered_schedule: executor-shape -- axis/kind/outputs per block, "
    "an ORDERED list of sibling child blocks (not a single chained child), "
    "and value_ids resolvable through rich.cells[...].hash",
    "[ordered-schedule]") {
  auto const water20 = orderedsched_water20_fixture();
  auto const cross = orderedsched_cross_iteration_fixture();

  // Per-block shape (water-20's {Κ} block): axis (an Index), kind (a
  // BatchModeType), outputs (an svector of (value_id, OutputKind) pairs).
  auto const k_idx =
      orderedsched_index_of_child_block(water20.sched.root, L"Κ");
  REQUIRE(k_idx.has_value());
  ScopeBlock const& k_block =
      std::get<ScopeBlock>(water20.sched.root.steps[*k_idx].value);
  CHECK(k_block.axis.space().base_key() == L"Κ");
  CHECK((k_block.kind == sequant::BatchModeType::Contracted ||
         k_block.kind == sequant::BatchModeType::External));
  for (auto const& [vid, out_kind] : k_block.outputs) {
    CHECK(vid < water20.rich.cells.size());
    CHECK((out_kind == OutputKind::AccumulateSum ||
           out_kind == OutputKind::AccumulateScatter));
  }

  // ORDERED LIST of sibling child blocks, not a single chained child: the
  // cross-iteration root holds TWO {i} ScopeBlock steps side by side in the
  // SAME steps list, each reachable in its own right by walking that list --
  // unlike scope_executor.hpp's ScopeNode::children.front()-only pattern,
  // which would silently see only the first.
  std::vector<ScopeBlock const*> root_children;
  for (auto const& step : cross.sched.root.steps)
    if (auto const* c = std::get_if<ScopeBlock>(&step.value))
      root_children.push_back(c);
  REQUIRE(root_children.size() == 2);  // the split producer/consumer siblings
  CHECK(root_children[0]->ordinal != root_children[1]->ordinal);
  CHECK(root_children[0]->axis.space().base_key() ==
        root_children[1]->axis.space().base_key());  // same axis TYPE ("i")

  // value_ids resolvable through rich.cells[...].hash -- exactly the bridge
  // scope_executor.hpp's build_value_node_map / walk_scope rely on
  // (rich.cells[vid].hash -> vmap lookup, see design integration point 1):
  // every value_id produced anywhere in EITHER schedule indexes validly into
  // its own RichSchedule, that index's stored value_id round-trips, and its
  // hash resolves to an actual forest node.
  for (auto const* fx : {&water20, &cross}) {
    auto const vmap = sequant::eval::build_value_node_map(fx->forest);
    std::vector<std::pair<std::size_t, std::optional<OutputKind>>> prods;
    orderedsched_collect_productions(fx->sched.root, prods);
    REQUIRE(!prods.empty());
    for (auto const& [vid, kind] : prods) {
      (void)kind;
      REQUIRE(vid < fx->rich.cells.size());
      CHECK(fx->rich.cells[vid].value_id == vid);
      CHECK(vmap.find(fx->rich.cells[vid].hash) != vmap.end());
    }
  }
}

// ===========================================================================
// Task 3 of the sliced-value canonical layout / loop-coloring plan (doc/dev/
// specs/2026-08-23-sliced-value-canonical-layout-loop-coloring-design.md):
// prove build_ordered_schedule's realized nest yields a per-(value,
// sliced-mode) -> DAG-scope-loop ASSIGNMENT (compute_sliced_mode_assignment /
// SlicedModeAssignment) -- plain DATA, not a per-cell mode_to_level lookup --
// the coloring input Task 5 feeds to canonicalize_slots.
//
// Reuses [sp2-noninner]'s occ-outer/aux-inner fixture, built inline here so
// this test keeps handles on the DF-leaves P3{Κ_1;i_3} and P4{Κ_2;i_4}: each
// carries BOTH an aux mode and an occ mode, fetched (as an operand of A3/A4
// resp.) under the SAME single realized occ block and the SAME single
// realized aux block nested inside it. This is exactly the shape that must
// map two DIFFERENT physical Index labels (Κ_1 vs Κ_2; i_3 vs i_4) to the
// SAME two loop ids (stability), while Y{;i_5} -- built INSIDE that same occ
// block but carrying a DIFFERENT physical occ index no enclosing loop of ITS
// OWN pushes onto it from outside -- gets NO entry for either loop
// (participation).
// ===========================================================================
TEST_CASE(
    "compute_sliced_mode_assignment: a DF-leaf's aux+occ modes map to the "
    "realized aux/occ loop ids; participation is respected; the SAME "
    "physical loop yields the SAME id across DIFFERENT values (Task 3, "
    "sliced-value canonical layout / loop coloring)",
    "[ordered-schedule][sp2-noninner]") {
  auto ctx = sequant::get_default_context().clone();
  auto isr = ctx.mutable_index_space_registry();
  REQUIRE(isr != nullptr);
  sequant::mbpt::add_df_spaces(isr);  // Κ (DF aux)
  auto ctx_resetter = sequant::set_scoped_default_context(std::move(ctx));

  // Same 2-axis occ-outer/aux-inner shape as orderedsched_2axis_forest_root,
  // built inline so this test keeps a handle on the DF-leaf P3. NOTE: P3/P4
  // (and A3/A4) are the SAME abstract tensor "P"/"A" up to occ-index
  // relabeling, so compute_dag_boulevard CSE-folds them to ONE ValueCell --
  // there is genuinely only ONE DF-leaf VALUE here (P3 and P4 are two
  // handles onto it), which is exactly why the fixture's own doc comment
  // says "B{;i_3,i_4} = A3 x A4 ... reads each A across occ-blocks": ONE
  // shared A/P, fetched from two structural positions.
  auto P3 = orderedsched_leaf("P{Κ_1;i_3}");
  auto Q1 = orderedsched_leaf("Q{;Κ_1}");
  auto A3 =
      orderedsched_inode("A{;i_3}", P3, Q1);  // contracts Κ_1, carries i_3

  auto P4 = orderedsched_leaf("P{Κ_2;i_4}");
  auto Q2 = orderedsched_leaf("Q{;Κ_2}");
  auto A4 =
      orderedsched_inode("A{;i_4}", P4, Q2);  // contracts Κ_2, carries i_4

  auto B = orderedsched_inode("B{;i_3,i_4}", A3, A4);  // outer product on occ
  orderedsched_stamp_all(B);

  // A SECOND, structurally UNRELATED value (M3, tensor labels La/Lb/M3/Lc --
  // none shared with P/Q/A/B) that happens to carry the SAME physical occ
  // Index i_3 as the DF-leaf/A-node above (matching the occ block's own
  // representative axis exactly): the genuine cross-VALUE stability witness
  // -- P and M3 are unrelated cells that must still resolve i_3 to the SAME
  // occ loop id. Wrapped in forest root C (see the [sp2-noninner]
  // populate_cell_mode_to_level TEST_CASE above for why the wrapper matters:
  // a forest root always escapes in full, keeping M3 itself a genuine
  // per-iteration Transient BuildStep inside the occ block).
  auto La = orderedsched_leaf("La{;i_3}");
  auto Lb = orderedsched_leaf("Lb{;i_3}");
  auto M3 = orderedsched_inode("M3{;i_3}", La, Lb);
  auto Lc = orderedsched_leaf("Lc{;i_3}");
  auto C = orderedsched_inode("C{;i_3}", M3, Lc);
  orderedsched_stamp_all(C);

  // A per-iteration Transient (Y) carrying a DIFFERENT physical occ index
  // than any DF-leaf/A-node/M3 above, wrapped in a further forest root Z:
  // the negative "no entry" (participation) fixture.
  auto Ya = orderedsched_leaf("Ya{;i_5}");
  auto Yb = orderedsched_leaf("Yb{;i_5}");
  auto Y = orderedsched_inode("Y{;i_5}", Ya, Yb);
  auto Yd = orderedsched_leaf("Yd{;i_5}");
  auto Z = orderedsched_inode("Z{;i_5}", Y, Yd);
  orderedsched_stamp_all(Z);

  sequant::BatchPolicy policy;
  policy.is_batchable_contracted_index = [](Index const& ix) {
    return ix.space().base_key() == L"Κ";
  };
  policy.is_batchable_external_index = [](Index const& ix) {
    return ix.space().base_key() == L"i";
  };
  policy.batch_spectator_indices = true;
  policy.node_level_placement = true;

  sequant::eval::dryrun::SizeRegime regime;
  regime.space_extent = {{L"i", 8u}, {L"Κ", 6u}};
  auto cm = std::make_shared<sequant::eval::dryrun::CostModel const>(regime);

  std::vector<sequant::EvalNode<sequant::EvalExpr>> forest{B, C, Z};
  sequant::stamp_lifetime_masks(forest);
  auto const block_of = [](Index const&) -> std::size_t { return 4; };
  auto rich = sequant::eval::compute_dag_boulevard(forest, *cm, block_of);
  REQUIRE(!rich.cells.empty());

  auto const legality = sequant::eval::analyze_legality(rich, forest, policy);
  auto const sched =
      sequant::eval::build_ordered_schedule(rich, legality, policy, {});
  REQUIRE(well_formed(sched));

  // Exactly ONE occ ("i") block at root, enclosing exactly ONE aux ("Κ")
  // block -- the [sp2-noninner] shape; no split forced by this forest.
  std::vector<ScopeBlock const*> occ_blocks;
  for (auto const& step : sched.root.steps)
    if (auto const* c = std::get_if<ScopeBlock>(&step.value))
      if (c->axis.space().base_key() == L"i") occ_blocks.push_back(c);
  REQUIRE(occ_blocks.size() == 1);
  ScopeBlock const& occ = *occ_blocks.front();

  std::vector<ScopeBlock const*> aux_blocks;
  for (auto const& s : occ.steps)
    if (auto const* c = std::get_if<ScopeBlock>(&s.value))
      if (c->axis.space().base_key() == L"Κ") aux_blocks.push_back(c);
  REQUIRE(aux_blocks.size() == 1);
  ScopeBlock const& aux = *aux_blocks.front();

  auto const cell_of = [&](sequant::EvalNode<sequant::EvalExpr> const& n)
      -> sequant::eval::ValueCell const& {
    auto const h = n->hash_value();
    auto const it = std::find_if(rich.cells.begin(), rich.cells.end(),
                                 [&](auto const& vc) { return vc.hash == h; });
    REQUIRE(it != rich.cells.end());
    return *it;
  };
  auto const& p_cell = cell_of(P3);  // == cell_of(P4): the ONE shared DF-leaf
  auto const& m3_cell = cell_of(M3);
  auto const& y_cell = cell_of(Y);

  // Preconditions: P3 and P4 really are the SAME folded value (this is what
  // makes the stability check below meaningful rather than vacuous), M3 is a
  // genuinely DIFFERENT value sharing physical index i_3, and Y is a plain
  // BuildStep directly inside `occ` carrying a different physical occ index
  // (mirrors the [sp2-noninner] populate_cell_mode_to_level TEST_CASE above).
  REQUIRE(p_cell.value_id == cell_of(P4).value_id);
  CHECK(p_cell.value_id != m3_cell.value_id);
  REQUIRE(orderedsched_index_of_build_step(occ, m3_cell.value_id).has_value());
  REQUIRE(orderedsched_index_of_build_step(occ, y_cell.value_id).has_value());
  REQUIRE(m3_cell.carried == sequant::container::svector<Index>{Index{L"i_3"}});
  REQUIRE(y_cell.carried == sequant::container::svector<Index>{Index{L"i_5"}});

  auto const assignment =
      sequant::eval::compute_sliced_mode_assignment(sched, rich);

  Index const K1{L"Κ_1"}, i3{L"i_3"}, i5{L"i_5"};

  // GREEN: the DF-leaf's OWN aux mode maps to the aux loop id, and its OWN
  // participating occ mode maps to the occ loop id.
  auto const p_aux = assignment.loop_of(p_cell.value_id, K1);
  auto const p_occ = assignment.loop_of(p_cell.value_id, i3);
  REQUIRE(p_aux.has_value());
  REQUIRE(p_occ.has_value());

  // Correct DagScopeLevel identity, read off the realized nest.
  CHECK(assignment.level_of(*p_aux) == aux.level);
  CHECK(assignment.level_of(*p_occ) == occ.level);

  // STABILITY: the SAME physical occ loop yields the SAME id across TWO
  // DIFFERENT, structurally-unrelated values (the DF-leaf P and M3) that
  // both happen to carry its physical Index i_3.
  auto const m3_occ = assignment.loop_of(m3_cell.value_id, i3);
  REQUIRE(m3_occ.has_value());
  CHECK(*m3_occ == *p_occ);

  // PARTICIPATION: Y carries i_5 (a DIFFERENT physical occ index than the
  // DF-leaf/A-node/M3 above) and is built INSIDE the occ block itself (not
  // fetched across it) -- it gets NO entry for any mode.
  CHECK_FALSE(assignment.loop_of(y_cell.value_id, i5).has_value());
  CHECK_FALSE(assignment.loop_of(y_cell.value_id, K1).has_value());
  CHECK_FALSE(assignment.loop_of(y_cell.value_id, i3).has_value());
}

// ===========================================================================
// Task 7 (sliced-value canonical layout / loop-coloring): the INTEGRATION
// co-pass on a REAL scheduled + sliced forest. Closes the Task-6 gap the
// driver test (test_sliced_canonical_layout.cpp) leaves open: that one runs
// populate_canonical_layouts with an EMPTY assignment (degenerate), so it
// never materializes a NON-empty layout. Here the occ-outer/aux-inner
// [sp2-noninner] fixture has genuine batched_here()/sliced modes, so the
// DF-leaf P (CSE-folded across two occurrences, each carrying an aux + an occ
// mode) gets a real 2-slot canonical_layout and two per-occurrence
// permutations -- and the SAME hash-keyed LoopColoredSliceSeam the ordered
// executor builds resolves each sliced mode to the realized aux/occ loop.
// ===========================================================================
TEST_CASE(
    "populate_canonical_layouts + LoopColoredSliceSeam on the sp2-noninner "
    "occ/aux forest: real per-value layout, per-occurrence permutations, and "
    "loop-resolved slice modes (Task 7)",
    "[ordered-schedule][sp2-noninner][sliced-layout]") {
  // occurrence_key canonicalizes the rank-2 DF-leaf tensors; scope the
  // strict-braket assertion off exactly as the [sliced-layout] probes do.
  auto ctx = sequant::get_default_context().clone();
  auto isr = ctx.mutable_index_space_registry();
  REQUIRE(isr != nullptr);
  sequant::mbpt::add_df_spaces(isr);  // Κ (DF aux)
  ctx.set(sequant::AssertStrictBraKetSymmetry::No);
  auto ctx_resetter = sequant::set_scoped_default_context(std::move(ctx));

  // Same occ-outer/aux-inner shape as the compute_sliced_mode_assignment test
  // above: B{;i_3,i_4} = A3 x A4, each A contracting its own DF aux mode; P3
  // and P4 are two handles onto the SAME CSE-folded DF-leaf value.
  auto P3 = orderedsched_leaf("P{Κ_1;i_3}");
  auto Q1 = orderedsched_leaf("Q{;Κ_1}");
  auto A3 = orderedsched_inode("A{;i_3}", P3, Q1);
  auto P4 = orderedsched_leaf("P{Κ_2;i_4}");
  auto Q2 = orderedsched_leaf("Q{;Κ_2}");
  auto A4 = orderedsched_inode("A{;i_4}", P4, Q2);
  auto B = orderedsched_inode("B{;i_3,i_4}", A3, A4);
  orderedsched_stamp_all(B);

  sequant::BatchPolicy policy;
  policy.is_batchable_contracted_index = [](Index const& ix) {
    return ix.space().base_key() == L"Κ";
  };
  policy.is_batchable_external_index = [](Index const& ix) {
    return ix.space().base_key() == L"i";
  };
  policy.batch_spectator_indices = true;
  policy.node_level_placement = true;

  sequant::eval::dryrun::SizeRegime regime;
  regime.space_extent = {{L"i", 8u}, {L"Κ", 6u}};
  auto cm = std::make_shared<sequant::eval::dryrun::CostModel const>(regime);

  std::vector<sequant::EvalNode<sequant::EvalExpr>> forest{B};
  sequant::stamp_lifetime_masks(forest);
  auto const block_of = [](Index const&) -> std::size_t { return 4; };
  auto rich = sequant::eval::compute_dag_boulevard(forest, *cm, block_of);
  REQUIRE(!rich.cells.empty());

  auto const legality = sequant::eval::analyze_legality(rich, forest, policy);
  auto const sched =
      sequant::eval::build_ordered_schedule(rich, legality, policy, {});
  REQUIRE(well_formed(sched));

  auto const assignment =
      sequant::eval::compute_sliced_mode_assignment(sched, rich);

  // Locate the ONE CSE-folded DF-leaf cell (P3 == P4) -- it must carry two
  // occurrences, exactly the "one value, two use-sites" shape the layout
  // reconciles.
  auto const p_hash = P3->hash_value();
  REQUIRE(p_hash == P4->hash_value());
  auto const p_it =
      std::find_if(rich.cells.begin(), rich.cells.end(),
                   [&](auto const& vc) { return vc.hash == p_hash; });
  REQUIRE(p_it != rich.cells.end());
  REQUIRE(p_it->occurrences.size() == 2);
  auto const p_vid = p_it->value_id;

  // Precondition: every occurrence's layout is empty BEFORE the co-pass.
  for (auto const& c : rich.cells)
    for (auto const& o : c.occurrences) REQUIRE(o.perm_to_canonical.empty());

  // --- THE CO-PASS: materialize canonical_layout + perm_to_canonical. ---
  sequant::eval::populate_canonical_layouts(forest, rich, assignment);

  auto const& p_cell = rich.cells[p_vid];

  // (1) The DF-leaf carries a genuine 2-slot canonical layout (its aux Κ mode
  // and its occ i mode are both sliced), NON-empty -- the real per-value fact
  // the degenerate driver test could not produce.
  REQUIRE(p_cell.canonical_layout.size() == 2);

  // (2) Both occurrences carry a per-occurrence permutation, each aligned with
  // (same slot count as) the value's canonical layout and expressed in THAT
  // occurrence's OWN physical labels (P3: Κ_1/i_3; P4: Κ_2/i_4).
  for (auto const& occ : p_cell.occurrences) {
    REQUIRE(occ.perm_to_canonical.size() == p_cell.canonical_layout.size());
    for (Index const& ix : occ.perm_to_canonical)
      CHECK(std::find(occ.carried.begin(), occ.carried.end(), ix) !=
            occ.carried.end());
  }

  // (3) The two occurrences bind the loops to DIFFERENT physical labels: their
  // permutations are NOT identical (P3's slots name Κ_1/i_3, P4's Κ_2/i_4) --
  // the reconciliation the design's per-occurrence permutation exists for.
  CHECK(p_cell.occurrences[0].perm_to_canonical !=
        p_cell.occurrences[1].perm_to_canonical);

  // --- THE SEAM: build it exactly as run_ordered_schedule_pre_results does,
  // and confirm the ordered executor's runtime resolution reaches the right
  // loop for each of the DF-leaf's own sliced modes. ---
  sequant::LoopColoredSliceSeam seam;
  seam.levels = assignment.levels;
  for (auto const& c : rich.cells) {
    auto const it = assignment.by_value.find(c.value_id);
    if (it != assignment.by_value.end())
      seam.by_hash.emplace(c.hash, it->second);
  }

  Index const K1{L"Κ_1"}, i3{L"i_3"};
  auto const aux_loop = assignment.loop_of(p_vid, K1);
  auto const occ_loop = assignment.loop_of(p_vid, i3);
  REQUIRE(aux_loop.has_value());
  REQUIRE(occ_loop.has_value());

  // loop_of_level inverts the level enumeration, and mode_of recovers P's own
  // sliced Index for that loop from the hash-keyed seam -- the exact two-step
  // slice_to_use performs per fetch (ordered arm).
  auto const aux_level = assignment.level_of(*aux_loop);
  auto const occ_level = assignment.level_of(*occ_loop);
  REQUIRE(seam.loop_of_level(aux_level) == aux_loop);
  REQUIRE(seam.loop_of_level(occ_level) == occ_loop);
  CHECK(seam.mode_of(p_hash, *aux_loop) == std::optional<Index>{K1});
  CHECK(seam.mode_of(p_hash, *occ_loop) == std::optional<Index>{i3});

  // A value hash the seam does not know leaves the fetch unsliced (nullopt).
  CHECK_FALSE(seam.mode_of(p_hash + 987654321u, *aux_loop).has_value());
}
