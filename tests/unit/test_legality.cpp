#include <SeQuant/core/basis.hpp>
#include <SeQuant/core/batch_policy.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/eval/backends/dryrun/cost_model_object.hpp>
#include <SeQuant/core/eval/backends/dryrun/eval_expr.hpp>
#include <SeQuant/core/eval/backends/dryrun/size_regime.hpp>

#include <SeQuant/core/eval/eval.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/eval/legality.hpp>
#include <SeQuant/core/eval/ordered_executor.hpp>
#include <SeQuant/core/eval/peak_profile.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/optimize/optimize.hpp>
#include <SeQuant/core/optimize/options.hpp>
#include <SeQuant/core/runtime.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/space_qns.hpp>
#include "csv_size_fixtures.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <optional>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

TEST_CASE(
    "build_site_of collects the batchable contracted axis of a binary "
    "contraction and excludes a non-batchable one",
    "[legality]") {
  using sequant::eval::dryrun::EvalExprDryRun;

  // g{i_1;a_3} * t{a_3;i_1}: contracts BOTH i_1 (occupied) and a_3
  // (virtual), fully -- so canon_indices() on the product node is empty and
  // contracted_indices() reports both i_1 and a_3. Mark only a_3 batchable
  // so the assertion distinguishes "in the build-site" from "in the
  // subtree at all".
  auto expr = sequant::deserialize<sequant::ExprPtr>("g{i_1;a_3} * t{a_3;i_1}");
  REQUIRE(static_cast<bool>(expr));

  SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
  auto node = sequant::binarize<EvalExprDryRun>(expr);
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
  REQUIRE_FALSE(node.leaf());

  sequant::BatchPolicy policy;
  policy.is_batchable_contracted_index = [](sequant::Index const& ix) {
    return ix.basis().base_key() == L"a";
  };

  // build_site_of now reads the DP's ACTUAL per-node decision (a Contracted
  // node_slice_mask stamp == the node's chosen aprime), NOT policy.is_batchable
  // (legality.hpp) -- so stamp a_3 as the DP would when it batches this
  // contraction here. i_1 is left unstamped (not batched), so it stays out.
  node->set_node_slice_mask(
      {{sequant::Index{L"a_3"}, sequant::BatchModeType::Contracted}});

  auto const site = sequant::eval::build_site_of(node, policy);

  auto const contains = [&](std::wstring const& key) {
    return std::any_of(site.begin(), site.end(), [&](sequant::Index const& ix) {
      return ix.basis().base_key() == key;
    });
  };
  CHECK(contains(L"a"));
  CHECK_FALSE(contains(L"i"));
}

TEST_CASE(
    "classify_axis decides per axis: a value carrying one occupied index and "
    "contracting another in batches is loop-local on the first and a "
    "reduction on the second",
    "[legality]") {
  using sequant::Index;
  using sequant::eval::LoopRole;
  using sequant::eval::OccurrenceRec;
  // value V{i_1;a_1} = sum_{i_2} X{i_1,i_2;a_1} * Y{i_2}: carries i_1 (sliced,
  // lockstep with an enclosing i_1 loop) and contracts i_2 at its own node
  // (sliced there by the DP). Both indices are of the occupied space.
  Index const i1{L"i_1"}, i2{L"i_2"}, a1{L"a_1"};
  sequant::container::svector<Index> const carried{i1, a1};
  sequant::container::svector<Index> const contracted_below{i2};
  OccurrenceRec occ;
  occ.carried = carried;
  occ.home = {i1};  // this occurrence's own home (positional, same frame)
  occ.ectx = {{i1, std::pair<std::size_t, std::size_t>{0, 16}}};
  occ.ectx_opener_point = {1};
  occ.loop_slot = {0, -1};
  sequant::container::svector<OccurrenceRec> const occurrences{occ};
  // the enclosing i_1 loop is the instance (i, 0) that slices position 0
  auto const enclosing = [](OccurrenceRec const&) {
    return sequant::container::svector<sequant::eval::LoopInstance>{{L"i", 0}};
  };
  CHECK(sequant::eval::classify_axis(carried, contracted_below, i1, occurrences,
                                     enclosing) == LoopRole::LoopLocal);
  CHECK(sequant::eval::classify_axis(carried, contracted_below, i2, occurrences,
                                     enclosing) == LoopRole::Reduction);
  // and an axis the value neither carries nor contracts is invariant
  CHECK(sequant::eval::classify_axis(carried, contracted_below, Index{L"i_3"},
                                     occurrences,
                                     enclosing) == LoopRole::LoopInvariant);
}

// ===========================================================================
// classify_axis / analyze_legality: four-way per-axis classification,
// validated on the real water-20 CSV-CCSD doubles residual (DF/aux-only
// batching). No shared test header exists for these DryRun fixtures, so the
// loading recipe is duplicated here under a `legality_` prefix to avoid any
// symbol collision under CMake UNITY_BUILD grouping with other test files'
// identically-shaped (but differently-named) anonymous-namespace helpers.
// ===========================================================================

namespace {

std::string legality_witness_slurp(std::string const& path) {
  std::ifstream in(path);
  std::stringstream ss;
  ss << in.rdbuf();
  return ss.str();
}

// Problem size, provenance and occupancy adjustment: csv_size_fixtures.hpp.
using LegalityWater20ProblemSize = sequant::tests::csv_sizes::CsvProblemSize;
inline constexpr auto const& kLegalityWater20_pVDZF12 =
    sequant::tests::csv_sizes::kWater20_pVDZF12;

sequant::eval::dryrun::SizeRegime legality_witness_df_regime(
    LegalityWater20ProblemSize const& p) {
  return sequant::tests::csv_sizes::df_size_regime(p);
}

sequant::ExprPtr legality_witness_flatten_product(sequant::ExprPtr const& e) {
  if (!e->is<sequant::Product>()) return e;
  auto const& p = e->as<sequant::Product>();
  return sequant::ex<sequant::Product>(p.scalar(), p.factors(),
                                       sequant::Product::Flatten::Yes);
}

}  // namespace

TEST_CASE(
    "classify_axis / analyze_legality: four-way axis classification on the "
    "water-20 aux-only residual",
    "[.][legality][blocked-stamp-based-legality]") {
  using sequant::eval::dryrun::EvalExprDryRun;
  using sequant::eval::dryrun::EvalNodeDryRun;
  using Node = EvalNodeDryRun;

  auto ctx = sequant::get_default_context_snapshot();
  ctx.set_first_dummy_index_ordinal(1000000);
  REQUIRE(ctx.index_basis_registry() != nullptr);
  auto ibr = std::make_shared<sequant::IndexBasisRegistry>(
      *ctx.index_basis_registry());
  sequant::mbpt::add_pao_basis(
      ibr, sequant::IndexSpace::QuantumNumbers{sequant::mbpt::Spin::any});
  sequant::mbpt::add_df_spaces(ibr);
  ctx.set(ibr);
  auto ctx_resetter = sequant::set_scoped_default_context(std::move(ctx));

  auto const body =
      legality_witness_slurp(std::string(SEQUANT_UNIT_TESTS_SOURCE_DIR) +
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

  // Same term cap as the water-20 witness: enough terms for aux
  // composites to be shared cross-tree, small enough to stay quick.
  std::size_t nterms = std::min<std::size_t>(summands.size(), 40);
  if (char const* nt = std::getenv("SEQUANT_UT_DRYRUN_NTERMS"))
    nterms = std::min<std::size_t>(summands.size(), std::atoll(nt));

  auto regime = legality_witness_df_regime(kLegalityWater20_pVDZF12);
  auto cm = std::make_shared<sequant::eval::dryrun::CostModel const>(regime);

  // EXACT MPQC aux-only config (make_csv_batch_policy, aux_target=256): Κ is
  // the only batchable mode, contracted role.
  sequant::BatchPolicy policy;
  policy.is_batchable_contracted_index = [](sequant::Index const& ix) {
    return ix.basis().base_key() == L"Κ";
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
    sequant::ExprPtr const term = legality_witness_flatten_product(summands[s]);
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

  // First real LegalitySchedule -- Step 3's wiring requirement. Every cell
  // must get one, with a per_axis entry for every axis in its OWN (at-node)
  // build site.
  auto const legality = sequant::eval::analyze_legality(rich, forest, policy);
  REQUIRE(legality.cells.size() == rich.cells.size());

  auto const is_K = [](sequant::Index const& ix) {
    return ix.basis().base_key() == L"Κ";
  };
  auto const carries_type =
      [&](sequant::container::svector<sequant::Index> const& v,
          auto const& pred) { return std::any_of(v.begin(), v.end(), pred); };

  auto const vmap = sequant::eval::build_value_node_map(forest);

  // A concrete Κ Index instance to pass as `axis` to classify_axis --
  // classify_axis only ever compares TYPE (space().base_key()), never
  // identity, so any Κ-typed Index found in the forest is a valid probe.
  std::optional<sequant::Index> k_axis;

  // The Κ-CONTRACTION RESULT -- a non-leaf value that does NOT
  // itself carry Κ (Q1 false) but DOES contract Κ at its own node (Q2a
  // true) => Reduction. Also record its hash so target 2 can find its
  // consumer.
  std::optional<std::size_t> mu_mu_hash;
  {
    // Two passes over the same predicate (Q1 false, Q2a true): the first
    // pass only accepts an EXACT μ̃,μ̃ carried pair (the canonical
    // example); if the fixture has none, the second pass falls back to the
    // first Κ-contraction-result candidate of any shape (the
    // classify_axis == Reduction assertion itself is the load-bearing part
    // either way, not the exact carried signature).
    auto const is_mu_mu_pair = [](auto const& carried) {
      return carried.size() == 2 &&
             std::all_of(carried.begin(), carried.end(),
                         [](sequant::Index const& ix) {
                           return ix.basis().base_key() == L"μ̃";
                         });
    };

    for (bool const require_mu_mu : {true, false}) {
      if (mu_mu_hash) break;
      for (auto const& vc : rich.cells) {
        auto const it = vmap.find(sequant::eval::value_key_of(vc));
        if (it == vmap.end() || it->second->leaf()) continue;
        if (carries_type(vc.carried, is_K)) continue;  // Q1 would be true
        auto const contracted = sequant::contracted_indices(*it->second);
        auto const k_it =
            std::find_if(contracted.begin(), contracted.end(), is_K);
        if (k_it == contracted.end()) continue;  // does not reduce Κ here
        if (require_mu_mu && !is_mu_mu_pair(vc.carried)) continue;
        if (!k_axis) k_axis = *k_it;
        mu_mu_hash = vc.hash;
        if (!require_mu_mu)
          WARN(
              "no exact μ̃,μ̃ Κ-contraction result in this fixture; using "
              "carried.size()="
              << vc.carried.size() << " instead");

        sequant::container::svector<sequant::Index> contracted_below(
            contracted.begin(), contracted.end());
        auto const no_enclosing = [](sequant::eval::OccurrenceRec const&) {
          return sequant::container::svector<sequant::eval::LoopInstance>{};
        };
        CHECK(sequant::eval::classify_axis(
                  vc.carried, contracted_below, *k_axis, vc.occurrences,
                  no_enclosing) == sequant::eval::LoopRole::Reduction);

        // End-to-end wiring check: analyze_legality's own per_axis for this
        // cell (Κ IS in its at-node build site: contracted-here) must
        // agree.
        auto const cl_it =
            std::find_if(legality.cells.begin(), legality.cells.end(),
                         [&](auto const& cl) { return cl.hash == vc.hash; });
        REQUIRE(cl_it != legality.cells.end());
        auto const ax_it =
            std::find_if(cl_it->per_axis.begin(), cl_it->per_axis.end(),
                         [&](auto const& ac) { return is_K(ac.axis); });
        REQUIRE(ax_it != cl_it->per_axis.end());
        CHECK(ax_it->role == sequant::eval::LoopRole::Reduction);

        // I(μ̃,μ̃)'s home_floor must EXCLUDE Κ -- a Reduction axis is
        // lifted out, never a home-floor member. Cross-check against the
        // existing meet-residency home (ValueCell::home_modes): both
        // "homes" should agree that Κ is not in the floor.
        CHECK_FALSE(carries_type(cl_it->home_floor, is_K));
        CHECK(carries_type(cl_it->home_floor, is_K) ==
              carries_type(vc.home_modes, is_K));
        break;
      }
    }
    REQUIRE(mu_mu_hash.has_value());
  }
  REQUIRE(k_axis.has_value());

  // The Κ-FREE COMPOSITE ABOVE THE CONTRACTION -- the immediate
  // consumer of the Κ-contraction result found above. Its OWN node does not
  // carry Κ (it was already summed away below) and does not contract Κ
  // AGAIN at its own node, so classify_axis(..., Κ, ...) must be
  // LoopInvariant -- this is exactly the case an uncorrected
  // (whole-subtree) build_site_of would have gotten wrong (it would have
  // seen Κ in the subtree and reported Reduction).
  {
    // A POINTER into `forest` (a local declared above, which outlives it), not
    // a subtree copy, found by an iterative walk in the SAME order the
    // recursion this replaced used: the trees are pushed in REVERSE so tree 0
    // pops first, and each node's RIGHT child is pushed before its LEFT, so the
    // pop order is forest-forward, node-then-left-subtree-then-right-subtree
    // pre-order. The forest reaches us as a single in-place Sum tree whose left
    // spine is as deep as the term count, so neither the subtree copy nor a
    // recursive search belongs here.
    Node const* parent = nullptr;
    {
      sequant::container::svector<Node const*> stack;
      for (auto it = forest.rbegin(); it != forest.rend(); ++it)
        stack.push_back(&*it);
      while (!parent && !stack.empty()) {
        Node const& n = *stack.back();
        stack.pop_back();
        if (n.leaf()) continue;
        if (n.left()->hash_value() == *mu_mu_hash ||
            n.right()->hash_value() == *mu_mu_hash) {
          parent = &n;
          break;
        }
        stack.push_back(&n.right());
        stack.push_back(&n.left());
      }
    }
    REQUIRE(parent != nullptr);

    auto const parent_hash = (*parent)->hash_value();
    auto const cell_it =
        std::find_if(rich.cells.begin(), rich.cells.end(),
                     [&](auto const& vc) { return vc.hash == parent_hash; });
    REQUIRE(cell_it != rich.cells.end());

    // Sanity: the composite must not itself carry Κ (else Q1 would route it
    // to LoopLocal/LoopCarried, not the LoopInvariant path being tested).
    REQUIRE_FALSE(carries_type(cell_it->carried, is_K));

    auto const contracted = sequant::contracted_indices(*parent);
    sequant::container::svector<sequant::Index> contracted_below(
        contracted.begin(), contracted.end());
    // Also not reduced AGAIN at this node.
    REQUIRE_FALSE(carries_type(contracted_below, is_K));

    CHECK(sequant::eval::classify_axis(
              cell_it->carried, contracted_below, *k_axis, cell_it->occurrences,
              sequant::eval::enclosing_instances_of(rich)) ==
          sequant::eval::LoopRole::LoopInvariant);

    // Illustrative signature check ([i i a a]); informational
    // -- the LoopInvariant assertion above is the load-bearing one.
    INFO("I(i,i;a,a) candidate carried.size()=" << cell_it->carried.size());

    // Κ must NOT appear in
    // this cell's own (at-node) build_site / per_axis at all (the
    // uncorrected whole-subtree build_site_of would have wrongly included
    // it as Reduction).
    auto const cl_it =
        std::find_if(legality.cells.begin(), legality.cells.end(),
                     [&](auto const& cl) { return cl.hash == parent_hash; });
    REQUIRE(cl_it != legality.cells.end());
    CHECK_FALSE(carries_type(cl_it->build_site, is_K));
    CHECK(std::none_of(cl_it->per_axis.begin(), cl_it->per_axis.end(),
                       [&](auto const& ac) { return is_K(ac.axis); }));

    // I(i,i;a,a)'s home_floor must EXCLUDE Κ too -- Κ is
    // LoopInvariant here (absent from build_site altogether, the IMPLICIT
    // case), so it can never surface in home_floor either. Cross-check
    // against the meet-residency home (ValueCell::home_modes).
    CHECK_FALSE(carries_type(cl_it->home_floor, is_K));
    CHECK(carries_type(cl_it->home_floor, is_K) ==
          carries_type(cell_it->home_modes, is_K));
  }

  // A Κ-CARRYING leaf/intermediate used AT the loop variable
  // (Q1 true, Q2b lockstep) -> LoopLocal. Scan every value that itself
  // carries a Κ-typed index (leaves included -- e.g. the DF integral
  // g{Κ;μ̃,ν̃}) and classify it; at least one LoopLocal
  // candidate must be found in this fixture (REQUIRE'd below, same rigor as
  // the Κ-contraction-result and Κ-free-composite checks). The loop-CARRIED
  // counterpart is covered by the synthetic cross-iteration fixture below, so
  // found_loop_carried is reported but NOT asserted either way.
  {
    bool found_loop_local = false;
    bool found_loop_carried = false;
    for (auto const& vc : rich.cells) {
      if (!carries_type(vc.carried, is_K)) continue;
      auto const it = vmap.find(sequant::eval::value_key_of(vc));
      if (it == vmap.end()) continue;
      auto const contracted = sequant::contracted_indices(*it->second);
      sequant::container::svector<sequant::Index> contracted_below(
          contracted.begin(), contracted.end());
      auto const role = sequant::eval::classify_axis(
          vc.carried, contracted_below, *k_axis, vc.occurrences,
          sequant::eval::enclosing_instances_of(rich));
      if (role == sequant::eval::LoopRole::LoopLocal) {
        found_loop_local = true;
        // End-to-end wiring check for this cell too (Κ IS in its build
        // site: carried).
        auto const cl_it =
            std::find_if(legality.cells.begin(), legality.cells.end(),
                         [&](auto const& cl) { return cl.hash == vc.hash; });
        REQUIRE(cl_it != legality.cells.end());
        auto const ax_it =
            std::find_if(cl_it->per_axis.begin(), cl_it->per_axis.end(),
                         [&](auto const& ac) { return is_K(ac.axis); });
        REQUIRE(ax_it != cl_it->per_axis.end());
        CHECK(ax_it->role == sequant::eval::LoopRole::LoopLocal);

        // A Κ-local value's home_floor must INCLUDE Κ -- it stays
        // sliced on (homed inside) the axis it is locked to lockstep with
        // its enclosing loop. Cross-check against the meet-residency home
        // (ValueCell::home_modes); INFO records any divergence for the
        // report rather than forcing agreement.
        CHECK(carries_type(cl_it->home_floor, is_K));
        WARN("Κ-local value hash="
             << vc.hash << " leaf=" << it->second->leaf()
             << " home_floor has Κ=" << carries_type(cl_it->home_floor, is_K)
             << " home_modes has Κ=" << carries_type(vc.home_modes, is_K));
      } else if (role == sequant::eval::LoopRole::LoopCarried) {
        found_loop_carried = true;
      }
    }
    // Every Κ-carrying candidate that WAS classified LoopLocal was already
    // checked above (per_axis wiring agreement); this REQUIRE guarantees the
    // LoopLocal path is actually exercised at least once by this fixture
    // (matching the Κ-contraction-result and Κ-free-composite checks, which
    // REQUIRE their own candidate was found) --
    // without it, a future shift in the fixture/optimizer output could
    // silently drop LoopLocal coverage entirely.
    INFO("water-20 Κ-carrying values: LoopLocal seen="
         << found_loop_local << " LoopCarried seen=" << found_loop_carried);
    REQUIRE(found_loop_local);
    // The loop-carried case is covered by the synthetic cross-iteration
    // fixture -- found_loop_carried is recorded (via INFO above) but not
    // asserted either way.
  }
}

// ===========================================================================
// forced_split_axes + the monotone fixpoint, on a SYNTHETIC
// cross-iteration fixture. Water-20's aux-only residual has no LoopCarried
// value (every Κ is contracted immediately), so the loop-carried case must be
// authored here. The fixture is B{i_3,i_4} = A{;i_3} * A{;i_4}: an outer
// product that carries TWO distinct occupied indices into its own result, so
// with the occupied space made batchable in the EXTERNAL role and no enclosing
// occ batch loop realized, the occ axis survives into the value's result with
// nothing to lock it to a loop iteration -> classify_axis == LoopCarried.
// Saved to data/legality_cross_iteration.txt for reuse below.
// ===========================================================================
TEST_CASE(
    "analyze_legality: a loop-carried value records its axis in "
    "forced_split_axes; the fixpoint terminates (synthetic cross-iteration "
    "fixture)",
    "[.][legality][blocked-stamp-based-legality]") {
  using sequant::eval::dryrun::EvalExprDryRun;
  using sequant::eval::dryrun::EvalNodeDryRun;
  using Node = EvalNodeDryRun;

  auto const body =
      legality_witness_slurp(std::string(SEQUANT_UNIT_TESTS_SOURCE_DIR) +
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

  // occ ("i") batchable in the EXTERNAL role -> i_3/i_4 are batch axes.
  sequant::BatchPolicy policy;
  policy.is_batchable_external_index = [](sequant::Index const& ix) {
    return ix.basis().base_key() == L"i";
  };

  // Minimal DryRun cost model: compute_dag_boulevard prices no footprints, and
  // block_of only sizes enclosing-loop entries (there are none here).
  sequant::eval::dryrun::SizeRegime regime;
  regime.space_extent = {{L"i", 8u}, {L"a", 16u}};
  auto cm = std::make_shared<sequant::eval::dryrun::CostModel const>(regime);

  std::vector<Node> forest{node};
  auto const block_of = [](sequant::Index const&) -> std::size_t { return 4; };
  auto rich = sequant::eval::compute_dag_boulevard(forest, *cm, block_of);
  REQUIRE(!rich.cells.empty());

  // The fixpoint must terminate (a run past the cap SEQUANT_ASSERTs inside).
  auto const legality = sequant::eval::analyze_legality(rich, forest, policy);
  REQUIRE(legality.cells.size() == rich.cells.size());

  auto const is_i = [](sequant::Index const& ix) {
    return ix.basis().base_key() == L"i";
  };
  auto const has_i = [&](auto const& v) {
    return std::any_of(v.begin(), v.end(), is_i);
  };

  // Every occ-carrying value here is LoopCarried on occ (survives the axis,
  // no enclosing occ loop): its occ axis must land in forced_split_axes and
  // be EXCLUDED from home_floor (a LoopCarried axis is lifted out).
  bool found_loop_carried = false;
  for (auto const& cl : legality.cells) {
    if (!has_i(cl.build_site)) continue;
    for (auto const& ac : cl.per_axis)
      if (is_i(ac.axis)) CHECK(ac.role == sequant::eval::LoopRole::LoopCarried);
    if (has_i(cl.forced_split_axes)) {
      found_loop_carried = true;
      CHECK_FALSE(has_i(cl.home_floor));
    }
  }
  // The B{i_3,i_4} outer-product root itself carries both occ indices; its
  // forced_split_axes must name occ.
  REQUIRE(found_loop_carried);

  // A value that is NOT loop-carried on occ has no occ forced split (sanity:
  // forced_split_axes is not spuriously populated). No such value exists in
  // this all-LoopCarried fixture, so we only assert the record never contains
  // a non-build-site axis.
  for (auto const& cl : legality.cells)
    for (auto const& fx : cl.forced_split_axes)
      CHECK(std::find(cl.build_site.begin(), cl.build_site.end(), fx) !=
            cl.build_site.end());
}

// ===========================================================================
// classify_axis: a value carrying TWO same-space indices sliced by two loop
// instances, only one of which encloses the occurrence. The per-occurrence
// verdict is per axis and per recorded instance: the position whose instance
// encloses the read is LoopLocal, the other is LoopCarried. Direct unit test
// on classify_axis with a HAND-CONSTRUCTED OccurrenceRec -- no
// forest/optimize pipeline needed.
// ===========================================================================
TEST_CASE(
    "classify_axis: a value carrying TWO same-space indices, one sliced by an "
    "enclosing loop instance and one by a non-enclosing one, classifies "
    "LoopLocal / LoopCarried per axis",
    "[legality]") {
  using sequant::Index;
  using sequant::container::svector;
  using sequant::eval::classify_axis;
  using sequant::eval::LoopInstance;
  using sequant::eval::LoopRole;
  using sequant::eval::OccurrenceAxisRole;
  using sequant::eval::OccurrenceRec;
  using sequant::eval::OccurrenceRole;

  auto const i_1 = Index{L"i_1"};
  auto const i_2 = Index{L"i_2"};

  // B{i_1,i_2} = A{;i_1} * A{;i_2}: an outer product carrying two DISTINCT
  // occ indices into its own result, each sliced by its own loop instance
  // (slots 0 and 1), read under a realized loop that is instance 0 only.
  svector<Index> const carried{i_1, i_2};
  svector<Index> const contracted_below{};

  OccurrenceRec occ;
  occ.point = 7;
  occ.consumer_point = 8;
  occ.carried = carried;
  occ.home = {i_1, i_2};
  occ.loop_slot = {0, 1};
  occ.ectx.push_back({i_1, {0u, 4u}});
  occ.ectx_opener_point.push_back(9);
  svector<OccurrenceRec> const occurrences{occ};

  // The enclosing instances, as analyze_legality resolves them through the
  // openers: only (i, 0) is running at this occurrence.
  auto const enclosing = [](OccurrenceRec const&) {
    return svector<LoopInstance>{{L"i", 0}};
  };

  svector<OccurrenceAxisRole> per_occ;
  CHECK(classify_axis(carried, contracted_below, i_1, occurrences, enclosing,
                      &per_occ) == LoopRole::LoopLocal);
  REQUIRE(per_occ.size() == 1);
  CHECK(per_occ.front().point == 7);
  CHECK(per_occ.front().role == OccurrenceRole::Lockstep);

  per_occ.clear();
  CHECK(classify_axis(carried, contracted_below, i_2, occurrences, enclosing,
                      &per_occ) == LoopRole::LoopCarried);
  REQUIRE(per_occ.size() == 1);
  CHECK(per_occ.front().role == OccurrenceRole::Assembled);

  // A batched position with no recorded loop identity is a defect of the
  // record, not a case to guess at.
  OccurrenceRec unstamped = occ;
  unstamped.loop_slot = {0, -1};
  REQUIRE_THROWS_AS(classify_axis(carried, contracted_below, i_2,
                                  svector<OccurrenceRec>{unstamped}, enclosing),
                    sequant::Exception);
}

// ===========================================================================
// classify_axis: a value with two occurrences, one read inside the loop
// instance that slices it (lockstep) and one read where that instance is not
// running (the value opens the loop itself; its consumer reads it whole).
// The value is LoopLocal -- it has a lockstep reader -- and the assembled
// read is recorded per occurrence, so the scheduler also escapes it.
// (Issue #720: the outside occurrence used to be skipped, and the value was
// built inside the loop only.)
// ===========================================================================
TEST_CASE(
    "classify_axis: one lockstep and one assembled occurrence of a value give "
    "LoopLocal with a per-occurrence Assembled verdict",
    "[legality]") {
  using sequant::Index;
  using sequant::container::svector;
  using sequant::eval::classify_axis;
  using sequant::eval::LoopInstance;
  using sequant::eval::LoopRole;
  using sequant::eval::OccurrenceAxisRole;
  using sequant::eval::OccurrenceRec;
  using sequant::eval::OccurrenceRole;

  auto const a_2 = Index{L"a_2"};
  svector<Index> const carried{Index{L"i_3"}, Index{L"i_1"}, a_2,
                               Index{L"a_3"}};
  svector<Index> const contracted_below{};

  OccurrenceRec in_loop;
  in_loop.point = 3;
  in_loop.consumer_point = 5;
  in_loop.carried = carried;
  in_loop.home = {a_2};
  in_loop.loop_slot = {-1, -1, 0, -1};
  in_loop.ectx.push_back({a_2, {0u, 4u}});
  in_loop.ectx_opener_point.push_back(6);

  OccurrenceRec at_root = in_loop;
  at_root.point = 11;
  at_root.consumer_point = 13;
  at_root.ectx.clear();
  at_root.ectx_opener_point.clear();
  at_root.opens.push_back({a_2, sequant::BatchModeType::External});

  svector<OccurrenceRec> const occurrences{in_loop, at_root};
  auto const enclosing = [](OccurrenceRec const& occ) {
    return occ.ectx.empty() ? svector<LoopInstance>{}
                            : svector<LoopInstance>{{L"a", 0}};
  };

  svector<OccurrenceAxisRole> per_occ;
  CHECK(classify_axis(carried, contracted_below, a_2, occurrences, enclosing,
                      &per_occ) == LoopRole::LoopLocal);
  REQUIRE(per_occ.size() == 2);
  CHECK(per_occ[0].point == 3);
  CHECK(per_occ[0].role == OccurrenceRole::Lockstep);
  CHECK(per_occ[1].point == 11);
  CHECK(per_occ[1].role == OccurrenceRole::Assembled);

  // Only the outside occurrence: the value is LoopCarried as a whole.
  per_occ.clear();
  CHECK(classify_axis(carried, contracted_below, a_2,
                      svector<OccurrenceRec>{at_root}, enclosing,
                      &per_occ) == LoopRole::LoopCarried);
  REQUIRE(per_occ.size() == 1);
  CHECK(per_occ.front().role == OccurrenceRole::Assembled);
}

// ===========================================================================
// enclosing_instances resolves each enclosing loop through the occurrence
// that opened it and that opener's stamps, and throws on every way the
// record can fail to name the instance -- never falling back to a label.
// ===========================================================================
TEST_CASE(
    "enclosing_instances: resolved through the opener's stamps; an opener "
    "the record does not name, cannot find, or that carries no stamp throws",
    "[legality]") {
  using sequant::Index;
  using sequant::container::svector;
  using sequant::eval::enclosing_instances;
  using sequant::eval::LoopInstance;
  using sequant::eval::OccurrenceRec;
  using sequant::eval::PointOccurrence;

  auto const i_1 = Index{L"i_1"}, k_1 = Index{L"a_1"};

  // The opener (point 1): carries i_1 sliced by (i, 2) and opens that loop as
  // External; also opens a Contracted loop over a_1 stamped (a, 0).
  OccurrenceRec opener;
  opener.point = 1;
  opener.consumer_point = 1;
  opener.carried = {i_1};
  opener.home = {i_1};
  opener.loop_slot = {2};
  opener.opens = {{i_1, sequant::BatchModeType::External},
                  {k_1, sequant::BatchModeType::Contracted}};
  opener.reduced_slot = {{k_1, 0}};
  std::unordered_map<std::size_t, PointOccurrence> points;
  points[1] = {0, &opener};

  // An occurrence under both loops names the opener for each.
  OccurrenceRec occ;
  occ.point = 5;
  occ.consumer_point = 1;
  occ.ectx = {{i_1, {0u, 4u}}, {k_1, {0u, 8u}}};
  occ.ectx_opener_point = {1, 1};
  auto const inst = enclosing_instances(occ, points);
  REQUIRE(inst.size() == 2);
  CHECK(inst[0] == LoopInstance{L"i", 2});
  CHECK(inst[1] == LoopInstance{L"a", 0});

  // Fewer opener points than enclosing loops: the record is inconsistent.
  OccurrenceRec short_rec = occ;
  short_rec.ectx_opener_point = {1};
  REQUIRE_THROWS_AS(enclosing_instances(short_rec, points), sequant::Exception);

  // An opener point that is no occurrence.
  OccurrenceRec lost = occ;
  lost.ectx_opener_point = {1, 7};
  REQUIRE_THROWS_AS(enclosing_instances(lost, points), sequant::Exception);

  // An opener that does not open the loop the occurrence says it does.
  OccurrenceRec mislabeled = occ;
  mislabeled.ectx = {{Index{L"i_2"}, {0u, 4u}}};
  mislabeled.ectx_opener_point = {1};
  REQUIRE_THROWS_AS(enclosing_instances(mislabeled, points),
                    sequant::Exception);

  // An opener that opens the loop but carries no stamp for it.
  OccurrenceRec unstamped_opener = opener;
  unstamped_opener.loop_slot = {-1};
  std::unordered_map<std::size_t, PointOccurrence> points2;
  points2[1] = {0, &unstamped_opener};
  OccurrenceRec under_i = occ;
  under_i.ectx = {{i_1, {0u, 4u}}};
  under_i.ectx_opener_point = {1};
  REQUIRE_THROWS_AS(enclosing_instances(under_i, points2), sequant::Exception);
}

// ===========================================================================
// classify_axis (spectator): a value carrying a batched, lockstep-bound index
// (i_1) BESIDE a same-space FULL spectator dimension (i_4 -- NOT batched, no
// enclosing loop) must classify the batched axis LoopLocal. The spectator is
// a retained full dimension, not a loop, so it must not drag the axis to
// LoopCarried. (Regression: w8 aux+occ value [i_4 mu mu i_1], sliced only on
// i_1, was wrongly escaped and evicted before an in-loop consumer read it.)
// ===========================================================================
TEST_CASE(
    "classify_axis (spectator): a same-space FULL non-batched carried index "
    "does not force LoopCarried on a lockstep-bound batched axis",
    "[legality]") {
  using sequant::Index;
  using sequant::container::svector;
  using sequant::eval::classify_axis;
  using sequant::eval::LoopRole;
  using sequant::eval::OccurrenceRec;

  auto const i_1 = Index{L"i_1"};
  auto const i_4 = Index{L"i_4"};

  // Value carries i_4 (a full, un-batched occ dimension) and i_1 (batched).
  svector<Index> const carried{i_4, i_1};
  svector<Index> const contracted_below{};

  OccurrenceRec occ;
  occ.point = 0;
  occ.consumer_point = 1;
  occ.carried = carried;
  occ.home = {i_1};
  occ.loop_slot = {-1, 0};  // i_4 has no loop identity; i_1 is instance (i, 0)
  occ.ectx.push_back({i_1, {0u, 4u}});  // only i_1 has a realized loop
  occ.ectx_opener_point.push_back(2);

  svector<OccurrenceRec> const occurrences{occ};

  // Only i_1 is a batched (sliced) mode, by loop instance (i, 0); i_4 is a
  // full spectator dimension with no loop identity.
  auto const enclosing = [](OccurrenceRec const&) {
    return svector<sequant::eval::LoopInstance>{{L"i", 0}};
  };

  // i_1 is lockstep with its enclosing loop instance, so the axis is
  // LoopLocal -- the value is homed inside the i_1 loop, carrying full i_4
  // through. (i_4 is not a build-site axis and is never classified.)
  CHECK(classify_axis(carried, contracted_below, i_1, occurrences, enclosing) ==
        LoopRole::LoopLocal);
}

// ===========================================================================
// forced_split_types groups CellLegality::forced_split_axes by
// axis TYPE (IndexBasis::base_key()), collapsing multiple same-space Index
// INSTANCES (e.g. the cross-iteration fixture's i_3/i_4 outer product, both
// LoopCarried on occ) into the single loop (axis) they jointly force to
// split.
// ===========================================================================
TEST_CASE(
    "forced_split_types collapses two distinct same-space LoopCarried "
    "indices into ONE occ-space entry",
    "[legality]") {
  using sequant::Index;
  using sequant::eval::CellLegality;
  using sequant::eval::forced_split_types;

  CellLegality cl;
  cl.hash = 0;
  // Two DISTINCT occ Index instances, both LoopCarried on the same space --
  // the per-instance record the synthetic B{i_3,i_4} = A{;i_3} * A{;i_4}
  // outer product produces.
  cl.forced_split_axes = {Index{L"i_3"}, Index{L"i_4"}};

  auto const types = forced_split_types(cl);
  REQUIRE(types.size() == 1);
  CHECK(types.front().basis().base_key() == L"i");

  // A cell with no forced splits collapses to an empty list.
  CellLegality const empty_cl;
  CHECK(forced_split_types(empty_cl).empty());

  // Two DIFFERENT-space LoopCarried axes (occ + virtual) must NOT collapse --
  // forced_split_types groups WITHIN a space type only. ("i"/"a" are present
  // in the default context without any scoped registry setup, unlike the DF
  // aux space "Κ", which needs a scoped registry.)
  CellLegality mixed_cl;
  mixed_cl.forced_split_axes = {Index{L"i_3"}, Index{L"i_4"}, Index{L"a_1"}};
  auto const mixed_types = forced_split_types(mixed_cl);
  REQUIRE(mixed_types.size() == 2);
}

TEST_CASE("forced_split_types keeps a named basis apart from its space",
          "[legality][basis]") {
  using namespace sequant;
  auto ibr = mbpt::make_min_sr_spaces();
  const IndexSpace a = ibr->retrieve_space(L"a"), i = ibr->retrieve_space(L"i");
  ibr->add(L"μ̃", a, std::numeric_limits<IndexBasis::instance_type>::max());
  auto ctx = set_scoped_default_context(
      Context({.index_basis_registry_shared_ptr = std::move(ibr),
               .vacuum = Vacuum::SingleProduct}));
  const Index pao(get_default_context().index_basis_registry()->retrieve(L"μ̃"),
                  1);
  const Index composite(a, 2, {Index(i, 1), Index(i, 2)});
  eval::CellLegality cell;
  cell.forced_split_axes = {composite, pao,
                            Index(a, 3, {Index(i, 1), Index(i, 3)})};
  const auto types = eval::forced_split_types(cell);
  REQUIRE(types.size() == 2);  // a (the composites) and μ̃
  CHECK(types[0] == composite);
  CHECK(types[1] == pao);
}

// ===========================================================================
// Legality acceptance: a single, explicitly-named test documenting that
// all FOUR LoopRole values are demonstrably produced by analyze_legality /
// classify_axis across the two canonical inputs. This test does NOT
// re-derive the detailed per-role assertions already proven end-to-end
// elsewhere in this file -- it only re-tallies which LoopRole each already-
// exercised input actually produces, as one consolidated acceptance point:
//   - "classify_axis / analyze_legality: four-way axis classification on the
//     water-20 aux-only residual" already REQUIREs/CHECKs, against
//     analyze_legality's own LegalitySchedule: Reduction (the mu~,mu~
//     Kappa-contraction result), LoopInvariant (the I(i,i;a,a) composite
//     above it -- the whole-subtree-vs-at-node case;
//     home_floor excludes Kappa in both), and LoopLocal (at least one
//     Kappa-carrying value, home_floor includes Kappa).
//   - "analyze_legality: a loop-carried value records its axis in
//     forced_split_axes; ..." already REQUIREs LoopCarried on the synthetic
//     cross-iteration fixture, plus forced_split_axes is
//     non-empty and home_floor EXCLUDES the split axis for that value.
// Both fixtures are re-run here (same recipe, same shared helpers) purely to
// tally LoopRole coverage; LoopInvariant is IMPLICIT in analyze_legality's
// own output (an axis classified LoopInvariant is, by construction, never a
// build_site member -- see CellLegality::home_floor's doc comment), so it is
// witnessed via one direct classify_axis call on a cell with no Kappa
// dependence at its own node, exactly as the water-20 test's Κ-free-composite
// check does.
// ===========================================================================
TEST_CASE(
    "legality acceptance: analyze_legality / classify_axis demonstrably "
    "produce all four LoopRole values (LoopLocal, Reduction, LoopCarried, "
    "LoopInvariant) across the water-20 and cross-iteration canonical "
    "inputs",
    "[.][legality][blocked-stamp-based-legality]") {
  using sequant::eval::LoopRole;

  bool saw_local = false;
  bool saw_reduction = false;
  bool saw_invariant = false;

  // --- Water-20 aux-only residual: LoopLocal, Reduction, LoopInvariant ---
  {
    using sequant::eval::dryrun::EvalExprDryRun;
    using sequant::eval::dryrun::EvalNodeDryRun;
    using Node = EvalNodeDryRun;

    auto ctx = sequant::get_default_context_snapshot();
    ctx.set_first_dummy_index_ordinal(1000000);
    REQUIRE(ctx.index_basis_registry() != nullptr);
    auto ibr = std::make_shared<sequant::IndexBasisRegistry>(
        *ctx.index_basis_registry());
    sequant::mbpt::add_pao_basis(
        ibr, sequant::IndexSpace::QuantumNumbers{sequant::mbpt::Spin::any});
    sequant::mbpt::add_df_spaces(ibr);
    ctx.set(ibr);
    auto ctx_resetter = sequant::set_scoped_default_context(std::move(ctx));

    auto const body =
        legality_witness_slurp(std::string(SEQUANT_UNIT_TESTS_SOURCE_DIR) +
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

    auto regime = legality_witness_df_regime(kLegalityWater20_pVDZF12);
    auto cm = std::make_shared<sequant::eval::dryrun::CostModel const>(regime);

    sequant::BatchPolicy policy;
    policy.is_batchable_contracted_index = [](sequant::Index const& ix) {
      return ix.basis().base_key() == L"Κ";
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
          legality_witness_flatten_product(summands[s]);
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

    for (auto const& cl : legality.cells) {
      for (auto const& ac : cl.per_axis) {
        if (ac.role == LoopRole::LoopLocal) saw_local = true;
        if (ac.role == LoopRole::Reduction) saw_reduction = true;
      }
    }

    auto const is_K = [](sequant::Index const& ix) {
      return ix.basis().base_key() == L"Κ";
    };
    auto const carries_K =
        [&](sequant::container::svector<sequant::Index> const& v) {
          return std::any_of(v.begin(), v.end(), is_K);
        };

    // A concrete Κ Index instance to probe with -- classify_axis only ever
    // compares TYPE, never identity.
    sequant::Index k_axis;
    bool found_k = false;
    for (auto const& vc : rich.cells) {
      auto const it = std::find_if(vc.carried.begin(), vc.carried.end(), is_K);
      if (it != vc.carried.end()) {
        k_axis = *it;
        found_k = true;
        break;
      }
    }
    REQUIRE(found_k);

    // LoopInvariant, witnessed directly (implicit case, never in per_axis):
    // a non-leaf cell whose OWN node neither carries Κ nor contracts it.
    auto const vmap = sequant::eval::build_value_node_map(forest);
    for (auto const& vc : rich.cells) {
      if (carries_K(vc.carried)) continue;
      auto const it = vmap.find(sequant::eval::value_key_of(vc));
      if (it == vmap.end() || it->second->leaf()) continue;
      auto const contracted = sequant::contracted_indices(*it->second);
      if (std::any_of(contracted.begin(), contracted.end(), is_K)) continue;
      sequant::container::svector<sequant::Index> contracted_below(
          contracted.begin(), contracted.end());
      if (sequant::eval::classify_axis(
              vc.carried, contracted_below, k_axis, vc.occurrences,
              sequant::eval::enclosing_instances_of(rich)) ==
          LoopRole::LoopInvariant) {
        saw_invariant = true;
        break;
      }
    }
  }

  CHECK(saw_local);
  CHECK(saw_reduction);
  CHECK(saw_invariant);

  // --- Synthetic cross-iteration fixture: LoopCarried ---
  bool saw_carried = false;
  {
    using sequant::eval::dryrun::EvalExprDryRun;
    using sequant::eval::dryrun::EvalNodeDryRun;
    using Node = EvalNodeDryRun;

    auto const body =
        legality_witness_slurp(std::string(SEQUANT_UNIT_TESTS_SOURCE_DIR) +
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
      return ix.basis().base_key() == L"i";
    };

    sequant::eval::dryrun::SizeRegime regime;
    regime.space_extent = {{L"i", 8u}, {L"a", 16u}};
    auto cm = std::make_shared<sequant::eval::dryrun::CostModel const>(regime);

    std::vector<Node> forest{node};
    auto const block_of = [](sequant::Index const&) -> std::size_t {
      return 4;
    };
    auto rich = sequant::eval::compute_dag_boulevard(forest, *cm, block_of);
    REQUIRE(!rich.cells.empty());

    auto const legality = sequant::eval::analyze_legality(rich, forest, policy);
    REQUIRE(legality.cells.size() == rich.cells.size());

    for (auto const& cl : legality.cells)
      for (auto const& ac : cl.per_axis)
        if (ac.role == LoopRole::LoopCarried) saw_carried = true;
  }

  CHECK(saw_carried);

  // Consolidated acceptance: all four roles demonstrably reachable.
  CHECK((saw_local && saw_reduction && saw_invariant && saw_carried));
}
