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
#include <SeQuant/core/eval/cell_table_builder.hpp>
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
#include <catch2/matchers/catch_matchers_string.hpp>

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
#include <unordered_set>
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
  child.latitude_ordinal = 0;
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
  inner.latitude_ordinal = 0;
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
  CHECK(p_loop.latitude_ordinal == 0);
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
  inner_block.latitude_ordinal = 0;
  inner_block.steps.push_back(Step{BuildStep{1}});  // per-iteration partial
  inner_block.outputs.push_back({2, OutputKind::AccumulateSum});

  ScopeBlock outer_block;
  outer_block.axis = outer;
  outer_block.latitude_ordinal = 0;
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
  a.latitude_ordinal = 0;
  a.steps.push_back(Step{BuildStep{1}});
  a.outputs.push_back({2, OutputKind::AccumulateSum});

  ScopeBlock b;
  b.axis = ax;
  b.latitude_ordinal =
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
  child_a.latitude_ordinal = 0;
  child_a.steps.push_back(Step{BuildStep{0}});

  ScopeBlock child_b;
  child_b.axis = i2;
  child_b.latitude_ordinal =
      0;  // duplicate ordinal at the same axis TYPE as child_a
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
    "[.][ordered-schedule][blocked-layers-1-2]") {
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
        auto const it = vmap.find(sequant::eval::value_key_of(vc));
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
//   W{a_2;a_3} = V * H, with Κ realized AT W (node_slice_mask) -- contracts Κ_1
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
    "[.][ordered-schedule][blocked-layers-1-2]") {
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
  W->set_node_slice_mask({{K, sequant::BatchModeType::Contracted}});
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
// task-loopid (2026-09-08 amendment 8): a batched reduction always owns a
// loop identity, even when EVERY operand of the contraction is an INPUT --
// a leaf, or (as here) a value that is never itself home-sliced on the
// reduced mode -- so there is no home-sliced child for compute_dag_boulevard
// to union a loop component through. P{;a_9} = A{i_1;} * B{;i_1} contracts i_1
// at its own node (node_slice_mask stamped Contracted); A and B are LEAVES,
// so neither is ever home-sliced (stamp_lifetime_masks never stamps a leaf's
// sliced_modes -- see lifetime_mask.hpp's "leaves are not stamped"). Before
// the fix, compute_dag_boulevard's union-find seeded a reduction_node ONLY
// from a home-sliced child's edge, so this shape never got a component: P's
// occurrences carried no reduced_slot for i_1, fusion_slot(P, i_1) returned
// -1, and build_ordered_schedule's escape placement guessed loop_slot 0 --
// landing the AccumulateSum escape in whatever nest happened to own slot 0
// rather than the loop A/B are actually read under, so the per-batch read of
// A and B was never sliced on i_1 (every batch contracted the FULL leaf, and
// the sum over n batches overcounted P by n).
// ===========================================================================
TEST_CASE(
    "compute_dag_boulevard seeds a loop identity for a batched reduction "
    "whose operands are ALL inputs (no home-sliced operand to union "
    "through); build_ordered_schedule places its single AccumulateSum "
    "escape by that slot, and the resulting cell-table reads slice both "
    "operands on it",
    "[ordered-schedule][cell_table]") {
  Index const i1{L"i_1"};

  auto A = orderedsched_leaf("A{i_1;}");
  auto B = orderedsched_leaf("B{;i_1}");
  // a_9 is an inert filler free index (EvalExpr's Tensor ctor requires a
  // non-empty index list); it is never batchable, so it plays no role below.
  auto P = orderedsched_inode("P{;a_9}", A, B);  // contracts i_1
  P->set_node_slice_mask({{i1, sequant::BatchModeType::Contracted}});
  P->set_batch_loops_opened_here({{i1, sequant::BatchModeType::Contracted}});

  sequant::BatchPolicy policy;
  policy.is_batchable_contracted_index = [](Index const& ix) {
    return ix.space().base_key() == L"i";
  };

  sequant::eval::dryrun::SizeRegime regime;
  regime.space_extent = {{L"i", 6u}};
  auto cm = std::make_shared<sequant::eval::dryrun::CostModel const>(regime);

  std::vector<sequant::EvalNode<sequant::EvalExpr>> forest{P};
  auto const block_of = [](Index const&) -> std::size_t { return 3; };
  auto rich = sequant::eval::compute_dag_boulevard(forest, *cm, block_of);
  REQUIRE(rich.cells.size() == 3);  // A, B, P

  auto const value_id_of = [&](std::string_view tensor) -> std::size_t {
    auto const hash = orderedsched_eval_tensor(tensor).hash_value();
    auto const it =
        std::find_if(rich.cells.begin(), rich.cells.end(),
                     [&](auto const& vc) { return vc.hash == hash; });
    REQUIRE(it != rich.cells.end());
    return it->value_id;
  };
  std::size_t const a_id = value_id_of("A{i_1;}");
  std::size_t const b_id = value_id_of("B{;i_1}");
  std::size_t const p_id = value_id_of("P{;a_9}");

  // Fix item 1: P's occurrence must carry a reduced_slot entry for i_1 -- a
  // real component, numbered, even though neither A nor B was ever
  // home-sliced to union it through.
  REQUIRE(!rich.cells[p_id].occurrences.empty());
  auto const& p_occ = rich.cells[p_id].occurrences.front();
  auto const rs_it =
      std::find_if(p_occ.reduced_slot.begin(), p_occ.reduced_slot.end(),
                   [&](auto const& e) { return e.first == i1; });
  REQUIRE(rs_it != p_occ.reduced_slot.end());
  CHECK(rs_it->second >= 0);

  auto const legality = sequant::eval::analyze_legality(rich, forest, policy);
  REQUIRE(legality.cells.size() == rich.cells.size());

  // Fix item 2: build_ordered_schedule must NOT throw (fusion_slot resolves
  // now that reduced_slot is stamped) and places P as a single AccumulateSum
  // output of the {i} block, with no BuildStep anywhere (single production
  // site: the escape itself).
  sequant::eval::OrderedSchedule sched;
  REQUIRE_NOTHROW(sched = sequant::eval::build_ordered_schedule(
                      rich, legality, policy, {L"i"}));
  REQUIRE(well_formed(sched));

  auto const i_child_idx = orderedsched_index_of_child_block(sched.root, L"i");
  REQUIRE(i_child_idx.has_value());
  ScopeBlock const& i_block =
      std::get<ScopeBlock>(sched.root.steps[*i_child_idx].value);
  CHECK(i_block.axis.space().base_key() == L"i");

  auto const p_out_it =
      std::find_if(i_block.outputs.begin(), i_block.outputs.end(),
                   [&](auto const& out) { return out.first == p_id; });
  REQUIRE(p_out_it != i_block.outputs.end());
  CHECK(p_out_it->second == OutputKind::AccumulateSum);

  // Exactly one AccumulateSum escape for P anywhere in the schedule (its only
  // production site).
  std::size_t p_sum_outputs = 0;
  std::function<void(ScopeBlock const&)> count_p = [&](ScopeBlock const& b) {
    for (auto const& [ovid, okind] : b.outputs)
      if (ovid == p_id && okind == OutputKind::AccumulateSum) ++p_sum_outputs;
    for (auto const& st : b.steps)
      if (auto const* child = std::get_if<ScopeBlock>(&st.value))
        count_p(*child);
  };
  count_p(sched.root);
  CHECK(p_sum_outputs == 1);

  // Fix item 3 (the actual overcounting bug): the cell table's Read of A and
  // of B, as operands of P's (synthesized) Build cell inside the {i} block,
  // must each carry a NON-EMPTY slice on i_1 -- each batch of the i-loop
  // reads only that batch's block of the leaf. An empty slice there is
  // exactly the pre-fix overcounting shape: the same FULL leaf read every
  // batch, summed n times instead of once.
  auto const sma = sequant::eval::compute_sliced_mode_assignment(sched, rich);
  auto const g = sequant::eval::detail::ordered_schedule_dep_graph(rich);
  sequant::eval::CellTableInputs in;
  in.ordered = &sched;
  in.rich = &rich;
  in.sliced = &sma;
  in.sliced_modes_of = [](std::size_t) {
    return sequant::container::svector<Index>{};
  };
  in.volatile_of = [](std::size_t) { return false; };
  in.n_batches_of = [](sequant::eval::LoopKey const&) -> std::size_t {
    return 2;
  };
  in.operands_of = [&](std::size_t vid) {
    auto const it = g.depends_on.find(vid);
    return it == g.depends_on.end() ? sequant::container::svector<std::size_t>{}
                                    : it->second;
  };
  auto const table = sequant::eval::build_cell_table(in);
  auto const violations =
      sequant::eval::validate_cell_table(table, sched.root, in.n_batches_of);
  CHECK(violations.empty());

  bool a_sliced = false, b_sliced = false;
  for (sequant::eval::Read const& r : table.reads) {
    if (r.operand_value_id == a_id && !r.slice.empty()) a_sliced = true;
    if (r.operand_value_id == b_id && !r.slice.empty()) b_sliced = true;
  }
  CHECK(a_sliced);
  CHECK(b_sliced);
}

// task-loopid fix round 1 (M1): the seeding loop above must be SCOPED to
// modes legality::classify_axis would actually call Reduction for --
// classify_axis reaches Reduction only via its Q1 test (the value carries
// NO index of the SAME SPACE as the contracted mode). A value that
// contracts a batched mode while ALSO carrying another index of that same
// space (e.g. an occ-space contraction inside an occ-carrying result) takes
// the Q2b branch instead (LoopLocal/LoopCarried here, since nothing opens
// an enclosing loop of that space at the root) and must not gain a new
// realized loop depth: it should carry NO reduced_slot and the schedule
// should keep exactly the pre-fix slot-0-fallback shape (an
// AccumulateScatter escape, not a throw -- the throw is Reduction-only).
TEST_CASE(
    "compute_dag_boulevard does NOT seed a reduction loop identity for a "
    "contracted mode the value also carries another index of the SAME "
    "SPACE for; the schedule keeps its previous (slot-0-fallback) shape",
    "[ordered-schedule]") {
  Index const i1{L"i_1"}, i2{L"i_2"}, i3{L"i_3"};

  // R{i_1;i_2} = A{i_1;i_3} * B{i_3;i_2}, contracting i_3 -- but R carries
  // i_1 AND i_2, both the SAME SPACE ("i") as the contracted i_3. A and B
  // are leaves, so neither is ever home-sliced (same as the companion test
  // above).
  auto A = orderedsched_leaf("A{i_1;i_3}");
  auto B = orderedsched_leaf("B{i_3;i_2}");
  auto R = orderedsched_inode("R{i_1;i_2}", A, B);  // contracts i_3
  R->set_node_slice_mask({{i3, sequant::BatchModeType::Contracted}});
  R->set_batch_loops_opened_here({{i3, sequant::BatchModeType::Contracted}});

  sequant::BatchPolicy policy;
  policy.is_batchable_contracted_index = [](Index const& ix) {
    return ix.space().base_key() == L"i";
  };

  sequant::eval::dryrun::SizeRegime regime;
  regime.space_extent = {{L"i", 6u}};
  auto cm = std::make_shared<sequant::eval::dryrun::CostModel const>(regime);

  std::vector<sequant::EvalNode<sequant::EvalExpr>> forest{R};
  auto const block_of = [](Index const&) -> std::size_t { return 3; };
  auto rich = sequant::eval::compute_dag_boulevard(forest, *cm, block_of);
  REQUIRE(rich.cells.size() == 3);  // A, B, R

  auto const value_id_of = [&](std::string_view tensor) -> std::size_t {
    auto const hash = orderedsched_eval_tensor(tensor).hash_value();
    auto const it =
        std::find_if(rich.cells.begin(), rich.cells.end(),
                     [&](auto const& vc) { return vc.hash == hash; });
    REQUIRE(it != rich.cells.end());
    return it->value_id;
  };
  std::size_t const r_id = value_id_of("R{i_1;i_2}");

  // R's occurrence carries a reduced_slot for i_3: the contracted mode owns a
  // loop identity even though R also carries other indices of the same space
  // (roles are per axis; the seeding no longer skips the same-space-carried
  // case).
  REQUIRE(!rich.cells[r_id].occurrences.empty());
  CHECK_FALSE(rich.cells[r_id].occurrences.front().reduced_slot.empty());

  auto const legality = sequant::eval::analyze_legality(rich, forest, policy);
  REQUIRE(legality.cells.size() == rich.cells.size());

  // R's only build-site axis is i_3 itself (the CONTRACTED-at-node test);
  // i_1/i_2 are plain untouched carried "spectator" indices (never in
  // R.sliced_modes(), since nothing opens an "i" loop enclosing the root).
  // Pin the classification directly: i_3 is a Reduction -- decided per axis,
  // the same-space carried i_1/i_2 do not change that.
  auto const r_legality_it = std::find_if(
      legality.cells.begin(), legality.cells.end(),
      [&](auto const& cl) { return cl.hash == rich.cells[r_id].hash; });
  REQUIRE(r_legality_it != legality.cells.end());
  REQUIRE(r_legality_it->per_axis.size() == 1);
  CHECK(r_legality_it->per_axis.front().axis == i3);
  CHECK(r_legality_it->per_axis.front().role ==
        sequant::eval::LoopRole::Reduction);

  // build_ordered_schedule must NOT throw (the Reduction resolves to the
  // stamped reduced_slot) and places R as an AccumulateSum output.
  sequant::eval::OrderedSchedule sched;
  REQUIRE_NOTHROW(sched = sequant::eval::build_ordered_schedule(
                      rich, legality, policy, {L"i"}));
  REQUIRE(well_formed(sched));

  bool found_sum = false;
  std::function<void(ScopeBlock const&)> find_r = [&](ScopeBlock const& b) {
    for (auto const& [ovid, okind] : b.outputs)
      if (ovid == r_id && okind == OutputKind::AccumulateSum) found_sum = true;
    for (auto const& st : b.steps)
      if (auto const* child = std::get_if<ScopeBlock>(&st.value))
        find_r(*child);
  };
  find_r(sched.root);
  CHECK(found_sum);
}

// task-loopid fix round 1 (m2): the escape-placement throw added alongside
// the seeding fix above. Through the real pipeline a Reduction axis always
// gets a reduced_slot now (that is the whole point of the fix), so this
// shape has to be built by hand: a value V reducing i_1 at its own node
// (Reduction role, no carried position for i_1) whose occurrence's
// reduced_slot is left EMPTY -- no loop identity at all, the exact
// precondition the throw guards, in the style of the existing hand-built
// rejection cases above ("outside its nest", "throws instead of silently
// mis-scheduling").
TEST_CASE(
    "build_ordered_schedule: a Reduction axis with no reduced_slot at all "
    "throws instead of guessing loop_slot 0",
    "[ordered-schedule]") {
  using sequant::eval::LoopRole;
  Index const i1{L"i_1"};

  sequant::eval::RichSchedule rich;
  {
    sequant::eval::ValueCell vc{};
    vc.value_id = 0;
    vc.hash = 8100;
    sequant::eval::OccurrenceRec o{};
    o.point = 0;
    o.consumer_point = 0;  // root
    // carried, loop_slot, and reduced_slot are all left default-empty: V
    // reduces i_1 (per the legality entry below) but has no slot recorded
    // for it anywhere.
    vc.occurrences.push_back(o);
    rich.cells.push_back(std::move(vc));
  }

  sequant::eval::LegalitySchedule legality;
  {
    sequant::eval::CellLegality cl;
    cl.hash = 8100;
    cl.per_axis.push_back({i1, LoopRole::Reduction});
    legality.cells.push_back(cl);
  }

  sequant::BatchPolicy policy;
  policy.is_batchable_contracted_index = [](Index const& ix) {
    return ix.space().base_key() == L"i";
  };

  REQUIRE_THROWS_WITH(
      sequant::eval::build_ordered_schedule(rich, legality, policy, {L"i"}),
      Catch::Matchers::ContainsSubstring("no loop identity"));
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
    "[.][ordered-schedule][blocked-layers-1-2]") {
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
  std::vector<int> ordinals{occ_blocks[0]->latitude_ordinal,
                            occ_blocks[1]->latitude_ordinal};
  std::sort(ordinals.begin(), ordinals.end());
  CHECK(ordinals == std::vector<int>{0, 1});

  ScopeBlock const* producer = nullptr;
  ScopeBlock const* consumer = nullptr;
  for (auto const* blk : occ_blocks)
    (blk->latitude_ordinal == 0 ? producer : consumer) = blk;
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
        (c->latitude_ordinal == 0 ? prod_pos : cons_pos) = i;
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
// the external occ, set_node_slice_mask Contracted for aux), the same way the
// aux-only fixtures above stamp Κ -- no optimize() run.
namespace {
// Stamp occ as EXTERNAL and aux (Κ) as Contracted in node_slice_mask, sourcing
// the Index identities from the node's OWN canon/contracted indices.
// sliced_modes is then DERIVED by stamp_lifetime_masks (the cross-occurrence
// meet), never hand-set -- that is how the real optimize()->binarize path
// realizes an axis.
void orderedsched_stamp_2axis(sequant::EvalNode<sequant::EvalExpr>& n) {
  using sequant::BatchModeType;
  sequant::container::svector<std::pair<Index, BatchModeType>> stamps;
  for (auto const& ix : n->canon_indices())
    if (ix.space().base_key() == L"i")
      stamps.push_back({ix, BatchModeType::External});
  for (auto const& ix : sequant::contracted_indices(n))
    if (ix.space().base_key() == L"Κ")
      stamps.push_back({ix, BatchModeType::Contracted});
  if (!stamps.empty()) n->set_node_slice_mask(stamps);
  // Loop-OPENS (peak_profile builds ectx from these): the aux (Κ) Contracted
  // loop opens at ITS contraction node -- this node, where contracted_indices
  // put it. The occ (i) External loop opens at the ROOT only (added below), not
  // at each carrying node, so it is NOT stamped here.
  sequant::container::svector<std::pair<Index, BatchModeType>> opens;
  for (auto const& [ix, kind] : stamps)
    if (kind == BatchModeType::Contracted) opens.push_back({ix, kind});
  if (!opens.empty()) n->set_batch_loops_opened_here(opens);
}

// Post-order walk stamping every node's node_slice_mask + Contracted opens,
// then the occ (External) loop-open at the ROOT only (external mode is on the
// final result, so the root is its outermost carrier).
void orderedsched_stamp_all(sequant::EvalNode<sequant::EvalExpr>& n) {
  using sequant::BatchModeType;
  std::function<void(sequant::EvalNode<sequant::EvalExpr>&)> post =
      [&](sequant::EvalNode<sequant::EvalExpr>& m) {
        if (!m.leaf()) {
          post(m.left());
          post(m.right());
        }
        orderedsched_stamp_2axis(m);
      };
  post(n);
  sequant::container::svector<std::pair<Index, BatchModeType>> opens(
      n->batch_loops_opened_here().begin(), n->batch_loops_opened_here().end());
  for (auto const& ix : n->canon_indices())
    if (ix.space().base_key() == L"i")
      opens.push_back({ix, BatchModeType::External});
  n->set_batch_loops_opened_here(opens);
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
    "[.][ordered-schedule][sp2-noninner][blocked-layers-1-2]") {
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
  // Derive sliced_modes from the node_slice_mask stamps (the cross-occurrence
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
    CHECK(block->level.latitude_ordinal == block->latitude_ordinal);
    CHECK(block->level.depth == depth);
    orderedsched_check_levels(block->steps, depth + 1);
  }
}
}  // namespace

TEST_CASE(
    "build_ordered_schedule: ScopeBlock::level mirrors axis/ordinal at the "
    "correct nesting depth (Task 3, DAG-scope runtime slicing)",
    "[.][ordered-schedule][sp2-noninner][blocked-layers-1-2]") {
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

namespace {

sequant::eval::ValueCell orderedsched_levels_cell(
    std::size_t id, std::size_t hash,
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
}

sequant::eval::CellLegality orderedsched_levels_legality(
    std::size_t hash, sequant::eval::LoopRole role) {
  sequant::Index const i{L"i_1"};
  sequant::eval::CellLegality cl;
  cl.hash = hash;
  sequant::eval::AxisClass ac;
  ac.axis = i;
  ac.role = role;
  cl.per_axis.push_back(ac);
  if (role == sequant::eval::LoopRole::LoopCarried)
    cl.forced_split_axes.push_back(i);
  return cl;
}

// forced_split_levels's `inside` predicate for a fixture with no Reduction
// role at all: no Reduction source exists, so it is never consulted; always
// false is the only well-defined answer.
bool orderedsched_no_reduction_inside(std::size_t, std::size_t) {
  return false;
}

}  // namespace

// Chain: L0 (carried, no carried operand) -> C1 (carried, reads L0's full
// form) -> C2 (carried, reads C1) -> R (non-carried, reads C2). Passes must
// be 0, 1, 2, 3.
TEST_CASE("forced_split_levels: a carried chain gives consecutive passes",
          "[ordered-schedule][levels]") {
  using sequant::eval::LoopRole;
  sequant::eval::RichSchedule rich;
  // occurrence (point, consumer_point); consumer_point == point is a root
  rich.cells.push_back(orderedsched_levels_cell(0, 1000, {{0, 10}}));   // L0
  rich.cells.push_back(orderedsched_levels_cell(1, 1001, {{10, 20}}));  // C1
  rich.cells.push_back(orderedsched_levels_cell(2, 1002, {{20, 30}}));  // C2
  rich.cells.push_back(orderedsched_levels_cell(3, 1003, {{30, 30}}));  // R
  sequant::eval::LegalitySchedule legality;
  legality.cells.push_back(
      orderedsched_levels_legality(1000, LoopRole::LoopCarried));
  legality.cells.push_back(
      orderedsched_levels_legality(1001, LoopRole::LoopCarried));
  legality.cells.push_back(
      orderedsched_levels_legality(1002, LoopRole::LoopCarried));
  legality.cells.push_back(
      orderedsched_levels_legality(1003, LoopRole::LoopLocal));
  auto const g = sequant::eval::detail::ordered_schedule_dep_graph(rich);
  auto const lv = sequant::eval::detail::forced_split_levels(
      rich, legality, g, orderedsched_no_reduction_inside);
  CHECK(lv.carried == std::unordered_set<std::size_t>{0, 1, 2});
  CHECK(lv.pass(0) == 0);
  CHECK(lv.pass(1) == 1);
  CHECK(lv.pass(2) == 2);
  CHECK(lv.pass(3) == 3);
  CHECK(lv.max_pass == 3);
}

// Lift: V (non-carried, base 0) is read only by A and B, both in pass 2
// (they read C1's full form, C1 carried at pass 1). V moves to pass 2.
// Straddle: W (non-carried, base 0) is read by A (pass 2) and by P (pass 0);
// W stays at 0. Edge property: every dependency edge points to an equal or
// earlier pass.
TEST_CASE(
    "forced_split_levels: the lift follows readers that all sit later; "
    "straddling readers keep the base; edges never point later",
    "[ordered-schedule][levels]") {
  using sequant::eval::LoopRole;
  sequant::eval::RichSchedule rich;
  rich.cells.push_back(
      orderedsched_levels_cell(0, 1000, {{0, 10}}));  // L0 carried
  rich.cells.push_back(orderedsched_levels_cell(
      1, 1001, {{10, 40}, {11, 50}}));  // C1 carried, -> A, B
  rich.cells.push_back(
      orderedsched_levels_cell(2, 1002, {{20, 40}, {21, 50}}));  // V -> A, B
  rich.cells.push_back(
      orderedsched_levels_cell(3, 1003, {{30, 40}, {31, 60}}));  // W -> A, P
  rich.cells.push_back(
      orderedsched_levels_cell(4, 1004, {{40, 40}}));  // A root
  rich.cells.push_back(
      orderedsched_levels_cell(5, 1005, {{50, 50}}));  // B root
  rich.cells.push_back(
      orderedsched_levels_cell(6, 1006, {{60, 60}}));  // P root
  sequant::eval::LegalitySchedule legality;
  legality.cells.push_back(
      orderedsched_levels_legality(1000, LoopRole::LoopCarried));
  legality.cells.push_back(
      orderedsched_levels_legality(1001, LoopRole::LoopCarried));
  for (std::size_t h : {1002u, 1003u, 1004u, 1005u, 1006u})
    legality.cells.push_back(
        orderedsched_levels_legality(h, LoopRole::LoopLocal));
  auto const g = sequant::eval::detail::ordered_schedule_dep_graph(rich);
  auto const lv = sequant::eval::detail::forced_split_levels(
      rich, legality, g, orderedsched_no_reduction_inside);
  CHECK(lv.pass(0) == 0);  // L0
  CHECK(lv.pass(1) == 1);  // C1
  CHECK(lv.pass(4) == 2);  // A reads C1
  CHECK(lv.pass(5) == 2);  // B reads C1
  CHECK(lv.pass(6) == 0);  // P reads only W
  CHECK(lv.pass(2) == 2);  // V lifted to its readers' pass
  CHECK(lv.pass(3) == 0);  // W straddles passes 0 and 2: stays
  CHECK(lv.max_pass == 2);
  for (auto const& [v, ops] : g.depends_on)
    for (std::size_t o : ops) CHECK(lv.pass(o) <= lv.pass(v));
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

// Same recipe as the "[w20-auxocc]" TEST_CASE below (the real water-20
// CSV-CCSD doubles residual with AUX+OCC batching -- Kappa
// batchable-contracted, occ ("i") batchable-EXTERNAL, matching MPQC's
// make_csv_batch_policy with occ_target>0), stopping once the schedule is
// built. Unlike orderedsched_water20_fixture()'s aux-ONLY policy (where
// Kappa is always fully contracted and nothing is ever LoopCarried), this
// configuration genuinely carries values on occ -- shared by the
// [w20-auxocc] case and the forced_split_levels equivalence case below so
// the two do not silently drift apart.
OrderedSchedFixture orderedsched_water20_auxocc_fixture() {
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

  // Use the FULL residual (all summands) to match the MPQC w20 run -- a
  // truncated forest can miss the term interaction that produces the malformed
  // schedule. Overridable for bisecting which term first breaks well_formed.
  std::size_t nterms = summands.size();
  if (char const* nt = std::getenv("SEQUANT_UT_DRYRUN_NTERMS"))
    nterms = std::min<std::size_t>(summands.size(), std::atoll(nt));

  auto regime = orderedsched_witness_df_regime(kOrderedSchedWater20_pVDZF12);
  auto cm = std::make_shared<sequant::eval::dryrun::CostModel const>(regime);

  // AUX+OCC: Kappa batchable-contracted (aux), i batchable-EXTERNAL (occ).
  sequant::BatchPolicy policy;
  policy.is_batchable_contracted_index = [](sequant::Index const& ix) {
    return ix.space().base_key() == L"\x39a";
  };
  policy.is_batchable_external_index = [](sequant::Index const& ix) {
    auto const reg = sequant::get_default_context().index_space_registry();
    return reg && ix.space() && reg->is_pure_occupied(ix.space());
  };
  // MATCH MPQC make_csv_batch_policy with occ_target>0: occ batching turns on
  // spectator batching AND node-level placement (which changes the schedule
  // structure -- occ-outer hoisting). Both are gated on occ_target>0 there.
  policy.batch_spectator_indices = true;
  policy.node_level_placement = true;
  policy.batch_target_size = [](sequant::Index const& ix) -> std::size_t {
    return ix.space().base_key() == L"\x39a" ? 256 : 16;
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

  auto const block_of = [](sequant::Index const& ix) -> std::size_t {
    return ix.space().base_key() == L"\x39a" ? 256 : 16;
  };
  auto rich = sequant::eval::compute_dag_boulevard(forest, *cm, block_of);
  REQUIRE(!rich.cells.empty());

  auto legality = sequant::eval::analyze_legality(rich, forest, policy);
  REQUIRE(legality.cells.size() == rich.cells.size());

  // MPQC passes an EMPTY mode_order (cck.ipp: the whole_scope/ordered driver
  // calls evaluate(..., std::initializer_list<std::wstring>{}, ...)), letting
  // build_ordered_schedule DERIVE the forced-split axes from the legality --
  // NOT an explicit {L"Κ", L"i"}. Match that.
  auto sched = sequant::eval::build_ordered_schedule(
      rich, legality, policy, std::initializer_list<std::wstring>{});

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
    "[.][ordered-schedule][blocked-layers-1-2]") {
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
    std::vector<int> ordinals{occ_blocks[0]->latitude_ordinal,
                              occ_blocks[1]->latitude_ordinal};
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
    "[.][ordered-schedule][blocked-layers-1-2]") {
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
  CHECK(root_children[0]->latitude_ordinal !=
        root_children[1]->latitude_ordinal);
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
      CHECK(vmap.find(sequant::eval::value_key_of(fx->rich.cells[vid])) !=
            vmap.end());
    }
  }
}

// Each value's DIRECT operand value_ids are persisted on the schedule, so the
// value-driven ordered executor can fetch each operand by its own cell id.
// The persisted edges are exactly the dep graph the topo-sort used (derived
// from every OccurrenceRec's consumer_point), and only computed (non-leaf)
// values carry operands.
TEST_CASE(
    "build_ordered_schedule persists operand_vids (value/occurrence DAG edges)",
    "[ordered-schedule][value-id]") {
  auto ctx = sequant::get_default_context().clone();
  auto isr = ctx.mutable_index_space_registry();
  REQUIRE(isr != nullptr);
  sequant::mbpt::add_df_spaces(isr);
  auto ctx_resetter = sequant::set_scoped_default_context(std::move(ctx));

  auto B = orderedsched_2axis_forest_root();

  sequant::BatchPolicy policy;
  policy.is_batchable_contracted_index = [](Index const& ix) {
    return ix.space().base_key() == L"\x39a";  // Kappa
  };
  policy.is_batchable_external_index = [](Index const& ix) {
    return ix.space().base_key() == L"i";
  };
  policy.batch_spectator_indices = true;
  policy.node_level_placement = true;

  sequant::eval::dryrun::SizeRegime regime;
  regime.space_extent = {{L"i", 8u}, {L"\x39a", 6u}};
  auto cm = std::make_shared<sequant::eval::dryrun::CostModel const>(regime);

  std::vector<sequant::EvalNode<sequant::EvalExpr>> forest{B};
  sequant::stamp_lifetime_masks(forest);
  auto const block_of = [](Index const&) -> std::size_t { return 4; };
  auto rich = sequant::eval::compute_dag_boulevard(forest, *cm, block_of);
  REQUIRE(!rich.cells.empty());
  auto const legality = sequant::eval::analyze_legality(rich, forest, policy);
  auto const sched =
      sequant::eval::build_ordered_schedule(rich, legality, policy, {L"i"});

  // Persisted edges == the dep graph the topo-sort consumes (its documented
  // source), and there is at least one computed value with operands.
  auto const g = sequant::eval::detail::ordered_schedule_dep_graph(rich);
  REQUIRE(!sched.operand_vids.empty());
  CHECK(sched.operand_vids == g.depends_on);

  for (auto const& [vid, ops] : sched.operand_vids) {
    REQUIRE(vid < rich.cells.size());
    CHECK_FALSE(rich.cells[vid].is_leaf);  // only computed values have operands
    CHECK_FALSE(ops.empty());
    for (auto const op : ops) {
      CHECK(op < rich.cells.size());
      CHECK(op != vid);  // no self-dependency
    }
  }
}

// Task 3 (w20 repro): the SAME real water-20 CSV-CCSD doubles residual, but
// with AUX+OCC batching (Κ contracted + i external) -- the config the MPQC w20
// csv-cck run uses, which trips SEQUANT_ASSERT(well_formed(out)) INSIDE
// build_ordered_schedule. This reproduces that schedule-build failure as a
// local dry-run so the failing well_formed invariant (SEQUANT_DUMP_WF) can be
// diagnosed without the MPI/MPQC run.
TEST_CASE(
    "build_ordered_schedule: water-20 aux+occ residual builds a well-formed "
    "schedule",
    "[ordered-schedule][w20-auxocc]") {
  auto const fx = orderedsched_water20_auxocc_fixture();
  CHECK(well_formed(fx.sched));
  CHECK(fx.sched.num_values == fx.rich.cells.size());
}

namespace {

// A test-local copy of the two-set producer/consumer partition function
// ordered_schedule.hpp used to build before Task 4's per-nest pass-level
// design superseded it (now deleted from production code): the
// LoopCarried-on-axis set and its strict dependency-ancestor closure
// (upward, then the downward LoopLocal-member closure), verbatim. Kept ONLY
// to pin the equivalence this test checks -- that two pass levels (0 and 1)
// reproduce exactly this old two-set partition -- against the CURRENT
// `forced_split_levels`; it has no other caller and is not a claim about
// production behavior.
std::pair<std::unordered_set<std::size_t>, std::unordered_set<std::size_t>>
orderedsched_old_partition(
    std::wstring const& axis_key,
    sequant::eval::LegalitySchedule const& legality,
    sequant::eval::detail::OrderedScheduleDepGraph const& g) {
  std::unordered_set<std::size_t> carried;
  std::unordered_set<std::size_t> consumer_pass;
  for (sequant::eval::CellLegality const& cl : legality.cells) {
    bool const carried_here =
        std::any_of(cl.per_axis.begin(), cl.per_axis.end(),
                    [&](sequant::eval::AxisClass const& ac) {
                      return ac.role == sequant::eval::LoopRole::LoopCarried &&
                             ac.axis.space().base_key() == axis_key;
                    });
    if (!carried_here) continue;
    auto const it = g.value_id_of.find(cl.hash);
    if (it != g.value_id_of.end()) carried.insert(it->second);
  }

  // consumer_pass = strict dependency-ancestors of the carried set: walk UP
  // the consumer edges from each carried value.
  std::vector<std::size_t> stack;
  auto const push = [&](std::size_t v) {
    if (consumer_pass.insert(v).second) stack.push_back(v);
  };
  for (std::size_t c : carried) {
    auto const it = g.consumers_of.find(c);
    if (it != g.consumers_of.end())
      for (std::size_t p : it->second) push(p);
  }
  while (!stack.empty()) {
    std::size_t const v = stack.back();
    stack.pop_back();
    auto const it = g.consumers_of.find(v);
    if (it != g.consumers_of.end())
      for (std::size_t p : it->second) push(p);
  }
  // DOWNWARD closure: a LoopLocal member of this nest whose consumers ALL sit
  // in the consumer pass must move with them.
  for (bool grew = true; grew;) {
    grew = false;
    std::vector<std::size_t> members(consumer_pass.begin(),
                                     consumer_pass.end());
    for (std::size_t v : members) {
      auto const dit = g.depends_on.find(v);
      if (dit == g.depends_on.end()) continue;
      for (std::size_t op : dit->second) {
        if (consumer_pass.count(op) || carried.count(op)) continue;
        auto const cit = g.consumers_of.find(op);
        if (cit == g.consumers_of.end() || cit->second.empty()) continue;
        bool all_in = true;
        for (std::size_t c : cit->second)
          if (!consumer_pass.count(c)) {
            all_in = false;
            break;
          }
        if (all_in && consumer_pass.insert(op).second) grew = true;
      }
    }
  }
  return {std::move(carried), std::move(consumer_pass)};
}

}  // namespace

// Two-level equivalence, on REAL data: orderedsched_water20_auxocc_fixture()
// (the SAME shared fixture the "[w20-auxocc]" TEST_CASE above consumes --
// Kappa batchable-contracted, occ ("i") batchable-EXTERNAL, matching MPQC's
// make_csv_batch_policy with occ_target>0), which genuinely carries values
// on occ (unlike orderedsched_water20_fixture()'s aux-ONLY policy, where
// Kappa is always fully contracted and nothing is ever LoopCarried).
// old.consumer_pass is exactly today's upward-plus-downward pass-1 set; the
// property this pins -- pass(v) >= 1 iff v is in old.consumer_pass -- holds
// for a carried chain of ANY depth, not just the depth-1 case this real
// fixture happens to produce (max_pass here is 1; deeper chains, where it
// would be > 1, are exercised by the two synthetic cases above).
TEST_CASE("forced_split_levels: two levels reproduce the two-set partition",
          "[ordered-schedule][levels]") {
  auto const fx = orderedsched_water20_auxocc_fixture();
  auto const g = sequant::eval::detail::ordered_schedule_dep_graph(fx.rich);
  auto const [old_carried, old_consumer_pass] =
      orderedsched_old_partition(L"i", fx.legality, g);
  auto const lv = sequant::eval::detail::forced_split_levels(
      fx.rich, fx.legality, g, orderedsched_no_reduction_inside);
  REQUIRE(!old_carried.empty());
  CHECK(lv.carried == old_carried);
  CHECK(lv.max_pass >= 1);
  for (std::size_t v = 0; v < fx.rich.cells.size(); ++v)
    CHECK((lv.pass(v) >= 1) == (old_consumer_pass.count(v) != 0));
}

// ===========================================================================
// Fix round 1 (review-task-2.md, Important 2): the one-forced-space
// assertion must count only GENUINE splits (a space whose levels reach
// pass >= 1 somewhere), not every space some cell happens to be
// LoopCarried on. Two spaces, "i" and "a", each carry exactly one value;
// only the "i" value has a reader (so only "i" is genuine -- "a"'s
// max_pass stays 0). build_ordered_schedule must not throw the
// more-than-one-forced-space assertion and must produce a well-formed
// schedule.
// ===========================================================================
TEST_CASE(
    "build_ordered_schedule: two LoopCarried spaces of which only one is "
    "genuine build without the more-than-one-forced-space assertion",
    "[ordered-schedule][levels]") {
  using sequant::eval::LoopRole;
  sequant::Index const i1{L"i_1"};
  sequant::Index const a1{L"a_1"};

  sequant::eval::RichSchedule rich;
  // L_i (carried on "i"): produced at point 0, read by Reader_i at point 1.
  rich.cells.push_back(orderedsched_levels_cell(0, 9000, {{0, 1}}));
  // Reader_i: root (produced at point 1, no further consumer) -- makes
  // space "i" genuine (a real strict-ancestor reader of the carried value).
  rich.cells.push_back(orderedsched_levels_cell(1, 9001, {{1, 1}}));
  // L_a (carried on "a"): itself a root (point == consumer_point) -- no
  // reader at all, so "a" is LoopCarried but NOT genuine (max_pass == 0).
  rich.cells.push_back(orderedsched_levels_cell(2, 9002, {{2, 2}}));

  sequant::eval::LegalitySchedule legality;
  {
    sequant::eval::CellLegality cl;
    cl.hash = 9000;
    cl.per_axis.push_back({i1, LoopRole::LoopCarried});
    cl.forced_split_axes.push_back(i1);
    legality.cells.push_back(cl);
  }
  {
    sequant::eval::CellLegality cl;  // Reader_i: root, no axis roles at all
    cl.hash = 9001;
    legality.cells.push_back(cl);
  }
  {
    sequant::eval::CellLegality cl;
    cl.hash = 9002;
    cl.per_axis.push_back({a1, LoopRole::LoopCarried});
    cl.forced_split_axes.push_back(a1);
    legality.cells.push_back(cl);
  }

  sequant::BatchPolicy policy;
  policy.is_batchable_external_index = [](sequant::Index const& ix) {
    return ix.space().base_key() == L"i" || ix.space().base_key() == L"a";
  };

  sequant::eval::OrderedSchedule sched;
  REQUIRE_NOTHROW(sched = sequant::eval::build_ordered_schedule(
                      rich, legality, policy, {L"i", L"a"}));
  CHECK(sequant::eval::well_formed(sched));
}

// ===========================================================================
// Fix round 1 (review-task-2.md, Important 3): a LoopLocal-only value (no
// role-driven escape) that STRADDLES passes -- read by a same-pass reader
// (here, the carried value C it itself feeds) and by a LATER-pass reader
// homed OUTSIDE its nest (root-homed X) -- is neither materialized (its
// only same-nest reader, C, is not a LATER pass, so `later_same_nest_
// readers` is empty) nor silently accepted: legality promises an
// all-LoopLocal value has no reader outside its own nest, and the schedule
// must throw rather than build an inconsistent schedule.
//
// Nest: V is LoopLocal on BOTH i_1 (fusion slot 0) and i_2 (fusion slot 1),
// which anchors the two loop instances into one nest. V feeds C (LoopCarried
// on i_1) and is ALSO read directly by X, which is LoopLocal on i_3 (fusion
// slot 2) -- a slot that never co-occurs with slots 0/1 in any cell, so it
// forms its OWN, disjoint one-member nest: X's production site resolves
// CONFIDENTLY, and to a nest other than V's. (An earlier version of this
// fixture gave X no per_axis roles at all, so its production site did not
// resolve; production_depth now never guesses a nest for an unresolved
// reader (ruling I1), so an unresolved reader can no longer trip this
// tripwire -- the fixture must give the offending reader a mode that
// resolves, confidently, into a different nest, which is what this version
// does.) Passes (by hand, mirroring the [levels] lift/straddle test): C = 0
// (no carried operand); X = 1 (reads C's completed form, a strict
// ancestor); V's two consumers are C (pass 0) and X (pass 1) -- the
// straddle keeps V at its base pass 0 (no lift), so X (pass 1) is a genuine
// later-pass reader of V, produced in a sibling nest -- outside V's own.
// ===========================================================================
TEST_CASE(
    "build_ordered_schedule: a LoopLocal value straddling passes, read "
    "later by a value outside its nest, throws instead of silently "
    "mis-scheduling",
    "[ordered-schedule][levels]") {
  using sequant::eval::LoopRole;
  sequant::Index const i1{L"i_1"};
  sequant::Index const i2{L"i_2"};
  sequant::Index const i3{L"i_3"};

  sequant::eval::RichSchedule rich;
  {
    sequant::eval::ValueCell vc{};
    vc.value_id = 0;
    vc.hash = 7000;
    sequant::eval::OccurrenceRec o1{};
    o1.point = 11;
    o1.consumer_point = 20;  // feeds C
    o1.carried = {i1, i2};
    o1.loop_slot = {0, 1};
    vc.occurrences.push_back(o1);
    sequant::eval::OccurrenceRec o2{};
    o2.point = 12;
    o2.consumer_point = 30;  // read directly by X
    o2.carried = {i1, i2};
    o2.loop_slot = {0, 1};
    vc.occurrences.push_back(o2);
    rich.cells.push_back(std::move(vc));
  }
  {
    sequant::eval::ValueCell vc{};
    vc.value_id = 1;
    vc.hash = 7001;
    sequant::eval::OccurrenceRec o{};
    o.point = 20;
    o.consumer_point = 30;  // read by X
    o.carried = {i1};
    o.loop_slot = {0};
    vc.occurrences.push_back(o);
    rich.cells.push_back(std::move(vc));
  }
  {
    sequant::eval::ValueCell vc{};  // X: LoopLocal on i_3 (fusion slot 2) --
                                    // a second, disjoint nest from V's
                                    // (slots 0, 1).
    vc.value_id = 2;
    vc.hash = 7002;
    sequant::eval::OccurrenceRec o{};
    o.point = 30;
    o.consumer_point = 30;  // root
    o.carried = {i3};
    o.loop_slot = {2};
    vc.occurrences.push_back(o);
    rich.cells.push_back(std::move(vc));
  }

  sequant::eval::LegalitySchedule legality;
  {
    sequant::eval::CellLegality cl;  // V
    cl.hash = 7000;
    cl.per_axis.push_back({i1, LoopRole::LoopLocal});
    cl.per_axis.push_back({i2, LoopRole::LoopLocal});
    legality.cells.push_back(cl);
  }
  {
    sequant::eval::CellLegality cl;  // C
    cl.hash = 7001;
    cl.per_axis.push_back({i1, LoopRole::LoopCarried});
    cl.forced_split_axes.push_back(i1);
    legality.cells.push_back(cl);
  }
  {
    sequant::eval::CellLegality cl;  // X: LoopLocal on i_3, its own nest
    cl.hash = 7002;
    cl.per_axis.push_back({i3, LoopRole::LoopLocal});
    legality.cells.push_back(cl);
  }

  sequant::BatchPolicy policy;
  policy.is_batchable_external_index = [](sequant::Index const& ix) {
    return ix.space().base_key() == L"i";
  };

  REQUIRE_THROWS_WITH(
      sequant::eval::build_ordered_schedule(rich, legality, policy, {L"i"}),
      Catch::Matchers::ContainsSubstring("outside its nest"));
}

namespace {

// A value with modes `modes` on the batched space, every mode home-sliced,
// loop slots `slots` (parallel to modes), occurrences (point, consumer_point).
sequant::eval::ValueCell orderedsched_nest_cell(
    std::size_t id, std::size_t hash, std::vector<sequant::Index> const& modes,
    std::vector<int> const& slots,
    std::vector<std::pair<std::size_t, std::size_t>> const& occs) {
  sequant::eval::ValueCell vc{};
  vc.value_id = id;
  vc.hash = hash;
  vc.first_use = id;
  vc.last_use = id;
  for (auto const& [p, cp] : occs) {
    sequant::eval::OccurrenceRec o{};
    o.point = p;
    o.consumer_point = cp;
    o.carried.assign(modes.begin(), modes.end());
    o.home.assign(modes.begin(), modes.end());
    o.loop_slot.assign(slots.begin(), slots.end());
    vc.occurrences.push_back(std::move(o));
  }
  return vc;
}

// Legality: one AxisClass per mode with the given role (all the same role).
sequant::eval::CellLegality orderedsched_nest_legality(
    std::size_t hash, std::vector<sequant::Index> const& modes,
    sequant::eval::LoopRole role) {
  sequant::eval::CellLegality cl;
  cl.hash = hash;
  for (auto const& m : modes) {
    sequant::eval::AxisClass ac;
    ac.axis = m;
    ac.role = role;
    cl.per_axis.push_back(ac);
    if (role == sequant::eval::LoopRole::LoopCarried)
      cl.forced_split_axes.push_back(m);
  }
  return cl;
}

// Root-level blocks of the batched space, in schedule order.
std::vector<sequant::eval::ScopeBlock const*> orderedsched_root_blocks(
    sequant::eval::OrderedSchedule const& s) {
  std::vector<sequant::eval::ScopeBlock const*> out;
  for (auto const& st : s.root.steps)
    if (auto const* b = std::get_if<sequant::eval::ScopeBlock>(&st.value))
      out.push_back(b);
  return out;
}

// Derives the cell table of a hand-built per-nest fixture and returns it
// after asserting it validates clean. Pattern follows
// test_ordered_executor.cpp:1552 / :1905, minus the real-forest-derived
// callbacks (no EvalExpr node map exists for a hand-built schedule):
// sliced_modes_of and volatile_of read the fixture's own occurrence
// records / the fixture's own knowledge (nothing here is marked volatile),
// operands_of is the dependency graph already recovered from rich, and
// n_batches_of is a constant stand-in (2 per loop instance; its value plays
// no role in the well-formedness rules checked here).
sequant::eval::CellTable orderedsched_validated_table(
    sequant::eval::RichSchedule const& rich,
    sequant::eval::OrderedSchedule const& sched) {
  auto const g = sequant::eval::detail::ordered_schedule_dep_graph(rich);
  auto const sma = sequant::eval::compute_sliced_mode_assignment(sched, rich);
  sequant::eval::CellTableInputs in;
  in.ordered = &sched;
  in.rich = &rich;
  in.sliced = &sma;
  in.sliced_modes_of = [&](std::size_t vid) {
    return rich.cells[vid].occurrences.front().home;
  };
  in.volatile_of = [](std::size_t) { return false; };
  in.n_batches_of = [](sequant::eval::LoopKey const&) -> std::size_t {
    return 2;
  };
  in.operands_of = [&](std::size_t vid) {
    auto const it = g.depends_on.find(vid);
    return it == g.depends_on.end() ? sequant::container::svector<std::size_t>{}
                                    : it->second;
  };
  auto const table = sequant::eval::build_cell_table(in);
  auto const violations =
      sequant::eval::validate_cell_table(table, sched.root, in.n_batches_of);
  for (auto const& v : violations)
    UNSCOPED_INFO("[" << v.rule << "] " << v.what);
  for (auto const& [cid, pos] : table.unresolved)
    UNSCOPED_INFO("[unresolved] cell#" << cid << " position " << pos
                                       << " (value "
                                       << table.cells[cid].value_id << ")");
  REQUIRE(violations.empty());
  REQUIRE(table.unresolved.empty());
  return table;
}

struct OrderedSchedReductionInLoopFixture {
  sequant::eval::RichSchedule rich;
  sequant::eval::LegalitySchedule legality;
};

// Amendment 8 (design section 9.2): V (id 0) is Reduction on i_1 (reduced
// instance slot 0), read by u (id 1, also Reduction on the same instance --
// produced inside it) and by w (id 2, a root reader with no per_axis at
// all, outside the instance). V and u have empty carried/home (a Reduction
// axis is contracted at the value, no carried position -- fusion_slot falls
// back to reduced_slot, exactly as the escape placement does); reduced_slot
// is patched onto orderedsched_nest_cell's occurrences after the fact since
// that helper only sets carried/home/loop_slot from `modes`.
OrderedSchedReductionInLoopFixture orderedsched_reduction_in_loop_fixture() {
  using sequant::eval::LoopRole;
  sequant::Index const i1{L"i_1"};

  OrderedSchedReductionInLoopFixture fx;
  auto v = orderedsched_nest_cell(0, 9300, {}, {}, {{110, 200}, {111, 300}});
  for (auto& occ : v.occurrences) occ.reduced_slot = {{i1, 0}};
  fx.rich.cells.push_back(std::move(v));

  auto u = orderedsched_nest_cell(1, 9301, {}, {}, {{200, 200}});
  for (auto& occ : u.occurrences) occ.reduced_slot = {{i1, 0}};
  fx.rich.cells.push_back(std::move(u));

  fx.rich.cells.push_back(
      orderedsched_nest_cell(2, 9302, {}, {}, {{300, 300}}));

  fx.legality.cells.push_back(
      orderedsched_nest_legality(9300, {i1}, LoopRole::Reduction));
  fx.legality.cells.push_back(
      orderedsched_nest_legality(9301, {i1}, LoopRole::Reduction));
  fx.legality.cells.push_back(
      orderedsched_nest_legality(9302, {}, LoopRole::LoopLocal));
  return fx;
}

}  // namespace

// V reduced over the instance at depth 0 (Reduction on i_1, reduced_slot 0),
// read by u (also Reduction on i_1: produced inside the loop) and by w (no
// batched mode: root). pass(V)=0, pass(u)=1, pass(w)=0. This raw-levels test
// hand-rolls `inside` (the fixture has exactly one candidate bumping edge,
// V -> u); the per-nest-split test below exercises build_ordered_schedule's
// REAL `inside`, resolved from the loop chain (fusion_slot/depth_of_instance/
// production_depth/type_cluster) against the SAME fixture.
TEST_CASE(
    "forced_split_levels: an in-loop reader of a reduction is bumped, "
    "a root reader is not",
    "[ordered-schedule][levels]") {
  auto const fx = orderedsched_reduction_in_loop_fixture();
  auto const g = sequant::eval::detail::ordered_schedule_dep_graph(fx.rich);
  auto const inside = [](std::size_t reader_vid, std::size_t source_vid) {
    return reader_vid == 1 && source_vid == 0;  // u inside V's instance
  };
  auto const lv = sequant::eval::detail::forced_split_levels(
      fx.rich, fx.legality, g, inside);
  CHECK(lv.pass(0) == 0);  // V
  CHECK(lv.pass(1) == 1);  // u: bumped, produced inside V's reduced instance
  CHECK(lv.pass(2) == 0);  // w: root reader, not bumped
  CHECK(lv.pinned.count(0) == 1);  // V pinned: source of a bumping edge
}

TEST_CASE(
    "per-nest split: a reduction source read inside its own loop by a "
    "later pass gets two pass blocks",
    "[ordered-schedule][per-nest-split]") {
  using sequant::eval::OutputKind;
  auto const fx = orderedsched_reduction_in_loop_fixture();

  sequant::BatchPolicy policy;
  policy.is_batchable_external_index = [](sequant::Index const& ix) {
    return ix.space().base_key() == L"i";
  };
  auto const sched = sequant::eval::build_ordered_schedule(fx.rich, fx.legality,
                                                           policy, {L"i"});
  REQUIRE(well_formed(sched));

  // Two root blocks of the same loop (slot 0), latitudes 0 and 1.
  auto const roots = orderedsched_root_blocks(sched);
  REQUIRE(roots.size() == 2);
  sequant::eval::ScopeBlock const* b0 = nullptr;
  sequant::eval::ScopeBlock const* b1 = nullptr;
  for (auto const* b : roots) {
    CHECK(b->level.loop_slot == 0);
    if (b->latitude_ordinal == 0) b0 = b;
    if (b->latitude_ordinal == 1) b1 = b;
  }
  REQUIRE(b0 != nullptr);
  REQUIRE(b1 != nullptr);

  // Pass 0 lists V's AccumulateSum output; pass 1 holds u's production (also
  // an AccumulateSum output -- u has no BuildStep either, both values escape
  // by their own Reduction role, section 3.2's per-nest realization is
  // latitude-blind to which of build_ids/outputs supplies the pass).
  std::vector<std::pair<std::size_t, std::optional<OutputKind>>> p0, p1;
  orderedsched_collect_productions(*b0, p0);
  orderedsched_collect_productions(*b1, p1);
  auto const has = [](auto const& v, std::size_t id,
                      std::optional<OutputKind> kind) {
    return std::any_of(v.begin(), v.end(), [&](auto const& e) {
      return e.first == id && e.second == kind;
    });
  };
  CHECK(has(p0, 0, OutputKind::AccumulateSum));  // V's sum, pass 0
  CHECK(has(p1, 1, OutputKind::AccumulateSum));  // u's production, pass 1
  CHECK_FALSE(has(p1, 0, OutputKind::AccumulateSum));
  CHECK_FALSE(has(p0, 1, OutputKind::AccumulateSum));

  // Design section 4 / amendment 8's own gate: the derived cell table
  // validates clean -- the 9.1 partial-sum check does not fire (it would
  // have thrown loudly inside build_ordered_schedule above, before this
  // point, had the levels failed to bump u).
  orderedsched_validated_table(fx.rich, sched);
}

TEST_CASE(
    "per-nest split: only the nest holding a later-pass member is "
    "split; a mixed-pass member is scattered to root in its own nest",
    "[ordered-schedule][per-nest-split]") {
  using sequant::eval::LoopRole;
  using sequant::eval::OutputKind;
  sequant::Index const i1{L"i_1"}, i2{L"i_2"}, i3{L"i_3"}, i4{L"i_4"};
  std::vector<sequant::Index> const ab{i1, i2}, cd{i3, i4};

  sequant::eval::RichSchedule rich;
  rich.cells.push_back(orderedsched_nest_cell(0, 2000, ab, {0, 1}, {{0, 10}}));
  rich.cells.push_back(orderedsched_nest_cell(1, 2001, ab, {0, 1}, {{10, 10}}));
  rich.cells.push_back(orderedsched_nest_cell(2, 2002, cd, {2, 3}, {{20, 50}}));
  rich.cells.push_back(
      orderedsched_nest_cell(3, 2003, cd, {2, 3}, {{30, 40}, {31, 50}}));
  rich.cells.push_back(orderedsched_nest_cell(4, 2004, cd, {2, 3}, {{40, 40}}));
  rich.cells.push_back(orderedsched_nest_cell(5, 2005, cd, {2, 3}, {{50, 50}}));

  // ValueCell::carried (distinct from the per-occurrence carried
  // orderedsched_nest_cell already sets) is read directly by
  // cell_table_builder.hpp (slicing_instance) and by
  // compute_sliced_mode_assignment -- but never by build_ordered_schedule
  // itself (checked: that function only reads occ.carried, never
  // ValueCell::carried), so the four fixtures above are unaffected by
  // omitting it. The cell-table derivation below needs it; fill it in from
  // the same (per-value-consistent) modes the occurrences already carry.
  for (auto& vc : rich.cells)
    if (!vc.occurrences.empty()) vc.carried = vc.occurrences.front().carried;

  sequant::eval::LegalitySchedule legality;
  legality.cells.push_back(
      orderedsched_nest_legality(2000, ab, LoopRole::LoopLocal));
  legality.cells.push_back(
      orderedsched_nest_legality(2001, ab, LoopRole::LoopCarried));
  legality.cells.push_back(
      orderedsched_nest_legality(2002, cd, LoopRole::LoopCarried));
  legality.cells.push_back(
      orderedsched_nest_legality(2003, cd, LoopRole::LoopLocal));
  legality.cells.push_back(
      orderedsched_nest_legality(2004, cd, LoopRole::LoopCarried));  // P
  legality.cells.push_back(
      orderedsched_nest_legality(2005, cd, LoopRole::LoopCarried));  // X

  sequant::BatchPolicy policy;
  policy.is_batchable_external_index = [](sequant::Index const& ix) {
    return ix.space().base_key() == L"i";
  };
  auto const sched =
      sequant::eval::build_ordered_schedule(rich, legality, policy, {L"i"});
  REQUIRE(well_formed(sched));

  // Chain: four depths of one space (slots 0..3), two nests {0,1} and {2,3}.
  auto const roots = orderedsched_root_blocks(sched);
  REQUIRE(roots.size() == 3);  // nest A once, nest B twice (passes 0, 1)
  std::vector<int> lats;
  std::vector<int> slots;
  for (auto const* b : roots) {
    lats.push_back(b->latitude_ordinal);
    slots.push_back(b->level.loop_slot);
  }
  // Nest A (slot 0 outermost) is one block at latitude 0.
  auto const a_it = std::find(slots.begin(), slots.end(), 0);
  REQUIRE(a_it != slots.end());
  CHECK(lats[a_it - slots.begin()] == 0);
  CHECK(std::count(slots.begin(), slots.end(), 0) == 1);
  // Nest B (slot 2 outermost) is two blocks, latitudes 0 then 1, in order.
  std::vector<int> b_lats;
  for (std::size_t k = 0; k < roots.size(); ++k)
    if (slots[k] == 2) b_lats.push_back(lats[k]);
  CHECK(b_lats == std::vector<int>{0, 1});

  // V (id 3) is built in nest B's latitude-0 block (at its inner depth) and
  // scattered at both of nest B's levels -- rule 4 fires because X (id 5),
  // though it has no LoopLocal axis of its own (it is a fully LoopCarried
  // forest root, like C), is still PRODUCED inside nest B in pass 1
  // (production_depth, not local_home_depth, decides nest membership for a
  // rule-4 reader). C (id 2) and P (id 4) are each scattered (no BuildStep)
  // in latitude 0; X is scattered (no BuildStep) in latitude 1.
  sequant::eval::ScopeBlock const* b0 = nullptr;
  sequant::eval::ScopeBlock const* b1 = nullptr;
  for (std::size_t k = 0; k < roots.size(); ++k) {
    if (slots[k] == 2 && lats[k] == 0) b0 = roots[k];
    if (slots[k] == 2 && lats[k] == 1) b1 = roots[k];
  }
  REQUIRE(b0 != nullptr);
  REQUIRE(b1 != nullptr);
  std::vector<std::pair<std::size_t, std::optional<OutputKind>>> p0, p1;
  orderedsched_collect_productions(*b0, p0);
  orderedsched_collect_productions(*b1, p1);
  auto const has = [](auto const& v, std::size_t id,
                      std::optional<OutputKind> kind) {
    return std::any_of(v.begin(), v.end(), [&](auto const& e) {
      return e.first == id && e.second == kind;
    });
  };
  CHECK(has(p0, 3, std::nullopt));                   // V built in pass 0
  CHECK(has(p0, 3, OutputKind::AccumulateScatter));  // V scattered
  CHECK(std::count_if(p0.begin(), p0.end(), [](auto const& e) {
          return e.first == 3 && e.second == OutputKind::AccumulateScatter;
        }) == 2);                                    // at both levels
  CHECK(has(p0, 2, OutputKind::AccumulateScatter));  // C scattered
  CHECK_FALSE(has(p0, 2, std::nullopt));
  CHECK(has(p0, 4, OutputKind::AccumulateScatter));  // P scattered, pass 0
  CHECK_FALSE(has(p0, 4, std::nullopt));             // P has no BuildStep
  CHECK(has(p1, 5, OutputKind::AccumulateScatter));  // X scattered, pass 1
  CHECK_FALSE(has(p1, 5, std::nullopt));             // X has no BuildStep
  CHECK_FALSE(has(p1, 3, std::nullopt));             // V not rebuilt

  // Design section 4: the derived cell table validates clean, with zero
  // unresolved positions. R (id 1), P (id 4) and X (id 5) are each true
  // forest roots (no consumer anywhere in this fixture) delivered in full
  // over their own external batched index, so each is LoopCarried like C
  // rather than LoopLocal: a genuine forest root over a batched index does
  // need to escape that index to be delivered in full, which is also what
  // gives each of them a route to the table's root scope, satisfying the
  // life rule's "only at the ROOT scope is a zero-read cell legitimate"
  // (cell_table.hpp).
  auto const table = orderedsched_validated_table(rich, sched);

  // The mixed-pass value (V, id 3) has an Assemble cell at root scope (empty
  // path) -- the form the later pass (X, id 5) reads.
  bool v_root_assemble = false;
  for (auto const& c : table.cells)
    if (c.value_id == 3 &&
        c.production.kind == sequant::eval::ProductionKind::Assemble &&
        c.scope.path.empty())
      v_root_assemble = true;
  CHECK(v_root_assemble);
}

TEST_CASE("per-nest split: a carried chain gives three pass blocks",
          "[ordered-schedule][per-nest-split]") {
  using sequant::eval::LoopRole;
  using sequant::eval::OutputKind;
  sequant::Index const i1{L"i_1"}, i2{L"i_2"};
  std::vector<sequant::Index> const ab{i1, i2};
  sequant::eval::RichSchedule rich;
  rich.cells.push_back(orderedsched_nest_cell(0, 3000, ab, {0, 1}, {{0, 10}}));
  rich.cells.push_back(orderedsched_nest_cell(1, 3001, ab, {0, 1}, {{10, 20}}));
  rich.cells.push_back(orderedsched_nest_cell(2, 3002, ab, {0, 1}, {{20, 20}}));
  sequant::eval::LegalitySchedule legality;
  legality.cells.push_back(
      orderedsched_nest_legality(3000, ab, LoopRole::LoopCarried));
  legality.cells.push_back(
      orderedsched_nest_legality(3001, ab, LoopRole::LoopCarried));
  legality.cells.push_back(
      orderedsched_nest_legality(3002, ab, LoopRole::LoopLocal));
  sequant::BatchPolicy policy;
  policy.is_batchable_external_index = [](sequant::Index const& ix) {
    return ix.space().base_key() == L"i";
  };
  auto const sched =
      sequant::eval::build_ordered_schedule(rich, legality, policy, {L"i"});
  REQUIRE(well_formed(sched));
  auto const roots = orderedsched_root_blocks(sched);
  REQUIRE(roots.size() == 3);
  std::vector<int> lats;
  for (auto const* b : roots) lats.push_back(b->latitude_ordinal);
  CHECK(lats == std::vector<int>{0, 1, 2});
  std::vector<std::pair<std::size_t, std::optional<OutputKind>>> p[3];
  for (int k = 0; k < 3; ++k) orderedsched_collect_productions(*roots[k], p[k]);
  auto const has = [](auto const& v, std::size_t id,
                      std::optional<OutputKind> kind) {
    return std::any_of(v.begin(), v.end(), [&](auto const& e) {
      return e.first == id && e.second == kind;
    });
  };
  CHECK(has(p[0], 0, OutputKind::AccumulateScatter));
  CHECK(has(p[1], 1, OutputKind::AccumulateScatter));
  CHECK(has(p[2], 2, std::nullopt));
  CHECK_FALSE(has(p[0], 1, OutputKind::AccumulateScatter));
}

// NOTE (fixture, not builder): the brief's literal 3-value shape (C, V, X)
// gives V a SINGLE consumer X. forced_split_levels's reverse sweep lifts a
// non-carried value with one consumer exactly to that consumer's own pass
// (pass_of(v) = max(base(v), min over consumers' pass) with a one-element
// min), so V and X always land in the SAME pass and later_same_nest_readers
// (which requires pass_of(reader) > pass_of(v)) never fires. A same-nest,
// SAME-pass second reader P (the same "read by P (same pass) and X (later
// pass)" shape the two-nest test above uses for its V) is added so the
// reverse-sweep min keeps V's pass strictly below X's, giving a genuine
// later-pass reader while V itself (LoopLocal on i_2, slot 1 only) is
// invariant to its nest's outer level (i_1, slot 0): rule 4 (section 7.3)
// skips that level and scatters V only at the inner one it is loop-local on,
// so V's assembled form still reaches root. X is LoopCarried on both axes
// (not LoopLocal): it is a genuine forest root with no consumer of its own
// in this fixture, and the table validator's life rule admits a zero-read
// cell only at root scope, which only a fully escaped (LoopCarried) value
// reaches -- the same pattern the two-nest fixture above uses for its own
// zero-consumer roots.
TEST_CASE(
    "per-nest split: a member loop-local on the inner instance only is "
    "scattered to a root-resident form",
    "[ordered-schedule][per-nest-split]") {
  using sequant::eval::LoopRole;
  sequant::Index const i1{L"i_1"}, i2{L"i_2"};
  std::vector<sequant::Index> const ab{i1, i2}, b{i2};
  sequant::eval::RichSchedule rich;
  rich.cells.push_back(orderedsched_nest_cell(0, 4000, ab, {0, 1}, {{0, 30}}));
  rich.cells.push_back(
      orderedsched_nest_cell(1, 4001, b, {1}, {{10, 20}, {11, 30}}));
  rich.cells.push_back(orderedsched_nest_cell(2, 4002, {}, {}, {{20, 20}}));
  rich.cells.push_back(orderedsched_nest_cell(3, 4003, ab, {0, 1}, {{30, 30}}));
  // ValueCell::carried (distinct from the per-occurrence carried
  // orderedsched_nest_cell already sets) is read by cell_table_builder.hpp
  // (slicing_instance) and compute_sliced_mode_assignment, never by
  // build_ordered_schedule itself; fill it in for orderedsched_validated_
  // table below (see the two-nest fixture above).
  for (auto& vc : rich.cells)
    if (!vc.occurrences.empty()) vc.carried = vc.occurrences.front().carried;
  sequant::eval::LegalitySchedule legality;
  legality.cells.push_back(
      orderedsched_nest_legality(4000, ab, LoopRole::LoopCarried));
  legality.cells.push_back(
      orderedsched_nest_legality(4001, b, LoopRole::LoopLocal));
  legality.cells.push_back(
      orderedsched_nest_legality(4002, {}, LoopRole::LoopLocal));
  legality.cells.push_back(
      orderedsched_nest_legality(4003, ab, LoopRole::LoopCarried));
  sequant::BatchPolicy policy;
  policy.is_batchable_external_index = [](sequant::Index const& ix) {
    return ix.space().base_key() == L"i";
  };
  auto const sched =
      sequant::eval::build_ordered_schedule(rich, legality, policy, {L"i"});
  REQUIRE(well_formed(sched));
  auto const roots = orderedsched_root_blocks(sched);
  REQUIRE(roots.size() == 2);  // one nest, passes 0 and 1
  std::vector<std::pair<std::size_t, std::optional<sequant::eval::OutputKind>>>
      p0;
  orderedsched_collect_productions(*roots[0], p0);
  // V (id 1) built in pass 0 and scattered exactly once (over the inner
  // instance; no escape at the outer one it is invariant to)
  CHECK(std::count_if(p0.begin(), p0.end(), [](auto const& e) {
          return e.first == 1 && e.second == std::nullopt;
        }) == 1);
  CHECK(std::count_if(p0.begin(), p0.end(), [](auto const& e) {
          return e.first == 1 &&
                 e.second == sequant::eval::OutputKind::AccumulateScatter;
        }) == 1);
  auto const table = orderedsched_validated_table(rich, sched);
  // V's Assemble cell resides at root and is formed once per outer batch
  bool found = false;
  for (auto const& c : table.cells)
    if (c.value_id == 1 &&
        c.production.kind == sequant::eval::ProductionKind::Assemble) {
      found = true;
      CHECK(sequant::eval::detail::residency_scope(c).path.empty());
      CHECK(c.produce_if_absent);
      CHECK_FALSE(c.scope.path.empty());  // produced inside the outer loop
    }
  CHECK(found);
}

// NOTE (fixture, not builder): same single-consumer-lift issue as above --
// with only X reading V, V's pass is lifted to exactly X's pass and no
// later-pass reader is ever detected. A same-pass second reader P is added
// for the same reason as in the invariant-outer fixture above. V is
// Reduction on i_2 (the inner instance, slot 1) and LoopLocal on i_1 (the
// outer instance, slot 0): rule 4 (section 7.3) leaves the already-summed
// inner instance as-is and adds a scatter at the outer one, chaining a sum
// then a scatter so V's assembled form reaches root. X is LoopCarried on
// both axes (not LoopLocal), for the same zero-consumer/root-residency
// reason as the invariant-outer fixture above.
TEST_CASE(
    "per-nest split: a member reduced over the inner instance chains a "
    "sum then a scatter to root",
    "[ordered-schedule][per-nest-split]") {
  using sequant::eval::LoopRole;
  sequant::Index const i1{L"i_1"}, i2{L"i_2"};
  std::vector<sequant::Index> const ab{i1, i2};
  sequant::eval::RichSchedule rich;
  rich.cells.push_back(orderedsched_nest_cell(0, 5000, ab, {0, 1}, {{0, 30}}));
  rich.cells.push_back(
      orderedsched_nest_cell(1, 5001, {i1}, {0}, {{10, 20}, {11, 30}}));
  rich.cells.push_back(orderedsched_nest_cell(2, 5002, {}, {}, {{20, 20}}));
  rich.cells.push_back(orderedsched_nest_cell(3, 5003, ab, {0, 1}, {{30, 30}}));
  for (auto& occ : rich.cells[1].occurrences)
    occ.reduced_slot.push_back({i2, 1});
  // ValueCell::carried (distinct from the per-occurrence carried
  // orderedsched_nest_cell already sets) is read by cell_table_builder.hpp
  // (slicing_instance) and compute_sliced_mode_assignment, never by
  // build_ordered_schedule itself; fill it in for orderedsched_validated_
  // table below (see the two-nest fixture above). V's Reduction axis (i_2)
  // is not part of occ.carried (it is a reduced_slot mode, not a sliced
  // one), so this correctly leaves V's own carried set to just i_1.
  for (auto& vc : rich.cells)
    if (!vc.occurrences.empty()) vc.carried = vc.occurrences.front().carried;
  sequant::eval::LegalitySchedule legality;
  legality.cells.push_back(
      orderedsched_nest_legality(5000, ab, LoopRole::LoopCarried));
  {
    sequant::eval::CellLegality cl;
    cl.hash = 5001;
    sequant::eval::AxisClass a1;
    a1.axis = i1;
    a1.role = LoopRole::LoopLocal;
    sequant::eval::AxisClass a2;
    a2.axis = i2;
    a2.role = LoopRole::Reduction;
    cl.per_axis = {a1, a2};
    legality.cells.push_back(cl);
  }
  legality.cells.push_back(
      orderedsched_nest_legality(5002, {}, LoopRole::LoopLocal));
  legality.cells.push_back(
      orderedsched_nest_legality(5003, ab, LoopRole::LoopCarried));
  sequant::BatchPolicy policy;
  policy.is_batchable_external_index = [](sequant::Index const& ix) {
    return ix.space().base_key() == L"i";
  };
  auto const sched =
      sequant::eval::build_ordered_schedule(rich, legality, policy, {L"i"});
  REQUIRE(well_formed(sched));
  auto const roots = orderedsched_root_blocks(sched);
  REQUIRE(roots.size() == 2);  // one nest, passes 0 and 1
  std::vector<std::pair<std::size_t, std::optional<sequant::eval::OutputKind>>>
      p0;
  orderedsched_collect_productions(*roots[0], p0);
  // V (id 1): a Reduction-kind escape (the sum over the inner instance) and
  // an AccumulateScatter escape (over the outer instance), and exactly one
  // BuildStep.
  CHECK(std::count_if(p0.begin(), p0.end(), [](auto const& e) {
          return e.first == 1 &&
                 e.second == sequant::eval::OutputKind::AccumulateSum;
        }) == 1);
  CHECK(std::count_if(p0.begin(), p0.end(), [](auto const& e) {
          return e.first == 1 &&
                 e.second == sequant::eval::OutputKind::AccumulateScatter;
        }) == 1);
  CHECK(std::count_if(p0.begin(), p0.end(), [](auto const& e) {
          return e.first == 1 && e.second == std::nullopt;
        }) == 1);
  auto const table = orderedsched_validated_table(rich, sched);
  // V's OUTERMOST Assemble cell (the end of the sum-then-scatter chain, at
  // the shallowest scope) resides at root; the inner Sum-kind Assemble
  // (still bound to the outer, i_1 instance) legitimately does not.
  bool found = false;
  sequant::eval::TableCell const* outermost = nullptr;
  for (auto const& c : table.cells)
    if (c.value_id == 1 &&
        c.production.kind == sequant::eval::ProductionKind::Assemble) {
      found = true;
      if (!outermost || c.scope.path.size() < outermost->scope.path.size())
        outermost = &c;
    }
  CHECK(found);
  REQUIRE(outermost != nullptr);
  CHECK(sequant::eval::detail::residency_scope(*outermost).path.empty());
}

// Review fix round 1, Important 1: a depth can legitimately carry BOTH a
// LoopLocal mode and a Reduction mode of the SAME value (a reduce+carry-style
// collapse at one depth -- the role loop's own comment at the escape
// emission above already contemplates this for a reduce+carry PAIR). Same
// shape as the inner-escape fixture above, except V's Reduction mode (i_2)
// resolves through reduced_slot to the SAME fusion slot as its own LoopLocal
// mode (i_1, slot 0) instead of a distinct inner slot, so both land at ONE
// depth. Rule 4 must upgrade that depth's escape from AccumulateSum (the
// role loop's own reduction escape) to AccumulateScatter: the loop-local
// mode's batches are disjoint, so summing them would silently combine
// values that must stay separate.
TEST_CASE(
    "per-nest split: a loop-local mode sharing a depth with a reduced mode "
    "upgrades the escape to a scatter",
    "[ordered-schedule][per-nest-split]") {
  using sequant::eval::LoopRole;
  sequant::Index const i1{L"i_1"}, i2{L"i_2"};
  std::vector<sequant::Index> const ab{i1, i2};
  sequant::eval::RichSchedule rich;
  rich.cells.push_back(orderedsched_nest_cell(0, 6000, ab, {0, 1}, {{0, 30}}));
  rich.cells.push_back(
      orderedsched_nest_cell(1, 6001, {i1}, {0}, {{10, 20}, {11, 30}}));
  rich.cells.push_back(orderedsched_nest_cell(2, 6002, {}, {}, {{20, 20}}));
  rich.cells.push_back(orderedsched_nest_cell(3, 6003, ab, {0, 1}, {{30, 30}}));
  for (auto& occ : rich.cells[1].occurrences)
    occ.reduced_slot.push_back({i2, 0});  // SAME slot as i_1: one depth
  for (auto& vc : rich.cells)
    if (!vc.occurrences.empty()) vc.carried = vc.occurrences.front().carried;
  sequant::eval::LegalitySchedule legality;
  legality.cells.push_back(
      orderedsched_nest_legality(6000, ab, LoopRole::LoopCarried));
  {
    sequant::eval::CellLegality cl;
    cl.hash = 6001;
    sequant::eval::AxisClass a1;
    a1.axis = i1;
    a1.role = LoopRole::LoopLocal;
    sequant::eval::AxisClass a2;
    a2.axis = i2;
    a2.role = LoopRole::Reduction;
    cl.per_axis = {a1, a2};
    legality.cells.push_back(cl);
  }
  legality.cells.push_back(
      orderedsched_nest_legality(6002, {}, LoopRole::LoopLocal));
  legality.cells.push_back(
      orderedsched_nest_legality(6003, ab, LoopRole::LoopCarried));
  sequant::BatchPolicy policy;
  policy.is_batchable_external_index = [](sequant::Index const& ix) {
    return ix.space().base_key() == L"i";
  };
  auto const sched =
      sequant::eval::build_ordered_schedule(rich, legality, policy, {L"i"});
  REQUIRE(well_formed(sched));
  auto const roots = orderedsched_root_blocks(sched);
  REQUIRE(roots.size() == 2);
  std::vector<std::pair<std::size_t, std::optional<sequant::eval::OutputKind>>>
      p0;
  orderedsched_collect_productions(*roots[0], p0);
  // V (id 1): exactly one escape at the shared depth, and it is a Scatter,
  // not a Sum -- the dominance upgrade fired.
  CHECK(std::count_if(p0.begin(), p0.end(), [](auto const& e) {
          return e.first == 1 &&
                 e.second == sequant::eval::OutputKind::AccumulateScatter;
        }) == 1);
  CHECK(std::count_if(p0.begin(), p0.end(), [](auto const& e) {
          return e.first == 1 &&
                 e.second == sequant::eval::OutputKind::AccumulateSum;
        }) == 0);
}

// Review fix round 1, Important 2: the outside-nest tripwire fires only when
// the value has NO role-driven escape (checked BEFORE rule 4 runs); a value
// that DOES have one is exempt, since that escape already assembles a full
// form with root residency that any later-pass reader, in any nest, can see.
// This pins the POSITIVE case (mirrors the existing [levels] "straddling
// passes" fixture, but checks the exact throw message under this task's own
// tag): V is LoopLocal only (no role escape) in nest A (i_4, i_5), feeds C
// (LoopCarried on i_4, same pass -- makes "i" genuine) and is ALSO read
// directly by X, LoopLocal on the disjoint i_6 (its own nest B), at a later
// pass -- V's own nest never opens for X, so legality and the schedule
// disagree.
TEST_CASE(
    "per-nest split: a loop-local value with no role escape, read later by "
    "a value outside its nest, throws",
    "[ordered-schedule][per-nest-split]") {
  using sequant::eval::LoopRole;
  sequant::Index const i4{L"i_4"}, i5{L"i_5"}, i6{L"i_6"};

  sequant::eval::RichSchedule rich;
  {
    sequant::eval::ValueCell vc{};
    vc.value_id = 0;
    vc.hash = 7100;
    sequant::eval::OccurrenceRec o1{};
    o1.point = 11;
    o1.consumer_point = 20;  // feeds C
    o1.carried = {i4, i5};
    o1.loop_slot = {0, 1};
    vc.occurrences.push_back(o1);
    sequant::eval::OccurrenceRec o2{};
    o2.point = 12;
    o2.consumer_point = 30;  // read directly by X
    o2.carried = {i4, i5};
    o2.loop_slot = {0, 1};
    vc.occurrences.push_back(o2);
    rich.cells.push_back(std::move(vc));
  }
  {
    sequant::eval::ValueCell vc{};
    vc.value_id = 1;
    vc.hash = 7101;
    sequant::eval::OccurrenceRec o{};
    o.point = 20;
    o.consumer_point = 30;  // read by X
    o.carried = {i4};
    o.loop_slot = {0};
    vc.occurrences.push_back(o);
    rich.cells.push_back(std::move(vc));
  }
  {
    sequant::eval::ValueCell vc{};  // X: LoopLocal on i_6, a disjoint nest
                                    // from V's own (slots 0, 1)
    vc.value_id = 2;
    vc.hash = 7102;
    sequant::eval::OccurrenceRec o{};
    o.point = 30;
    o.consumer_point = 30;  // root
    o.carried = {i6};
    o.loop_slot = {2};
    vc.occurrences.push_back(o);
    rich.cells.push_back(std::move(vc));
  }

  sequant::eval::LegalitySchedule legality;
  {
    sequant::eval::CellLegality cl;  // V
    cl.hash = 7100;
    cl.per_axis.push_back({i4, LoopRole::LoopLocal});
    cl.per_axis.push_back({i5, LoopRole::LoopLocal});
    legality.cells.push_back(cl);
  }
  {
    sequant::eval::CellLegality cl;  // C
    cl.hash = 7101;
    cl.per_axis.push_back({i4, LoopRole::LoopCarried});
    cl.forced_split_axes.push_back(i4);
    legality.cells.push_back(cl);
  }
  {
    sequant::eval::CellLegality cl;  // X: LoopLocal on i_6, its own nest
    cl.hash = 7102;
    cl.per_axis.push_back({i6, LoopRole::LoopLocal});
    legality.cells.push_back(cl);
  }

  sequant::BatchPolicy policy;
  policy.is_batchable_external_index = [](sequant::Index const& ix) {
    return ix.space().base_key() == L"i";
  };

  REQUIRE_THROWS_WITH(
      sequant::eval::build_ordered_schedule(rich, legality, policy, {L"i"}),
      Catch::Matchers::ContainsSubstring("legality and the schedule disagree"));
}

// Re-review fix round 2 (residual gap in fix round 1's Important 2): the
// tripwire's precondition must be scoped to the SPECIFIC LoopLocal instance
// that would otherwise silently vanish, not to "any role escape anywhere in
// the nest" -- a value can have TWO instances of its own nest at DIFFERENT
// depths, one role-escaped and one not, and the unescaped one is exactly as
// invisible to an outside-nest later-pass reader as it would be if the value
// had no role escape at all. V has i_1 (LoopLocal, its own home) and i_2
// (Reduction, escaped via reduced_slot) at TWO DIFFERENT depths of one nest;
// X reads V from a disjoint nest at a later pass, with no same-nest later-
// pass reader to trigger rule 4 for i_1's own instance -- so i_1's per-batch
// form is never delivered to root, and legality and the schedule disagree,
// even though i_2 IS role-escaped.
TEST_CASE(
    "per-nest split: a value with one escaped and one unescaped loop-local "
    "instance in one nest, read later from a different nest, throws",
    "[ordered-schedule][per-nest-split]") {
  using sequant::eval::LoopRole;
  sequant::Index const i1{L"i_1"}, i2{L"i_2"}, i3{L"i_3"};
  std::vector<sequant::Index> const ab{i1, i2};
  sequant::eval::RichSchedule rich;
  rich.cells.push_back(orderedsched_nest_cell(0, 8000, ab, {0, 1}, {{0, 30}}));
  rich.cells.push_back(
      orderedsched_nest_cell(1, 8001, {i1}, {0}, {{10, 20}, {11, 30}}));
  rich.cells.push_back(orderedsched_nest_cell(2, 8002, {}, {}, {{20, 20}}));
  rich.cells.push_back(orderedsched_nest_cell(3, 8003, {i3}, {2}, {{30, 30}}));
  for (auto& occ : rich.cells[1].occurrences)
    occ.reduced_slot.push_back({i2, 1});
  sequant::eval::LegalitySchedule legality;
  legality.cells.push_back(
      orderedsched_nest_legality(8000, ab, LoopRole::LoopCarried));
  {
    sequant::eval::CellLegality cl;  // V: i_1 LoopLocal (unescaped, its own
                                     // home), i_2 Reduction (escaped)
    cl.hash = 8001;
    sequant::eval::AxisClass a1;
    a1.axis = i1;
    a1.role = LoopRole::LoopLocal;
    sequant::eval::AxisClass a2;
    a2.axis = i2;
    a2.role = LoopRole::Reduction;
    cl.per_axis = {a1, a2};
    legality.cells.push_back(cl);
  }
  legality.cells.push_back(
      orderedsched_nest_legality(8002, {}, LoopRole::LoopLocal));
  legality.cells.push_back(  // X: LoopLocal on i_3, a disjoint nest from
                             // V's own (slots 0, 1)
      orderedsched_nest_legality(8003, {i3}, LoopRole::LoopLocal));
  sequant::BatchPolicy policy;
  policy.is_batchable_external_index = [](sequant::Index const& ix) {
    return ix.space().base_key() == L"i";
  };
  REQUIRE_THROWS_WITH(
      sequant::eval::build_ordered_schedule(rich, legality, policy, {L"i"}),
      Catch::Matchers::ContainsSubstring("legality and the schedule disagree"));
}
