#include <SeQuant/core/batch_policy.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/eval/backends/dryrun/cost_model_object.hpp>
#include <SeQuant/core/eval/backends/dryrun/eval_expr.hpp>
#include <SeQuant/core/eval/backends/dryrun/size_regime.hpp>
#include <SeQuant/core/eval/eval.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/eval/legality.hpp>
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

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <optional>
#include <sstream>
#include <string>
#include <unordered_map>
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
    return ix.space().base_key() == L"a";
  };

  auto const site = sequant::eval::build_site_of(node, policy);

  auto const contains = [&](std::wstring const& key) {
    return std::any_of(site.begin(), site.end(), [&](sequant::Index const& ix) {
      return ix.space().base_key() == key;
    });
  };
  CHECK(contains(L"a"));
  CHECK_FALSE(contains(L"i"));
}

// ===========================================================================
// classify_axis / analyze_legality: four-way per-axis classification,
// validated on the real water-20 CSV-CCSD doubles residual (DF/aux-only
// batching, exactly the config exercised by test_scope_executor.cpp's
// "[.][scope-executor-witness-water20]" witness). No shared test header
// exists for these DryRun fixtures (see that file's own doc comments), so
// the loading recipe is duplicated here under a `legality_` prefix to avoid
// any symbol collision under CMake UNITY_BUILD grouping with
// test_scope_executor.cpp's identically-shaped (but differently-named)
// anonymous-namespace helpers.
// ===========================================================================

namespace {

std::string legality_witness_slurp(std::string const& path) {
  std::ifstream in(path);
  std::stringstream ss;
  ss << in.rdbuf();
  return ss.str();
}

// water-20 pVDZ-F12 problem size (duplicated from test_scope_executor.cpp's
// kWitnessWater20_pVDZF12 / test_eval_dryrun.cpp's kWater20_pVDZF12).
struct LegalityWater20ProblemSize {
  std::size_t mu_tilde;
  std::size_t aux;
  std::size_t i_occ;
  std::array<double, 5> pno_M;
  std::array<double, 5> osv_M;
};

inline constexpr LegalityWater20ProblemSize kLegalityWater20_pVDZF12{
    /*mu_tilde=*/896u,
    /*aux=*/1682u,
    /*i_occ=*/80u,
    /*pno_M=*/
    {1.0, 23.175775480059084, 25.865548281212597, 28.171416142614103,
     30.03848680550367},
    /*osv_M=*/
    {1.0, 58.987499999999997, 59.289227520688783, 59.584437469011633,
     59.872014818179686}};

sequant::eval::dryrun::SizeRegime legality_witness_df_regime(
    LegalityWater20ProblemSize const& p) {
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
    "[legality]") {
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

  // Same term cap as the whole-scope witness: enough terms for aux
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
    return ix.space().base_key() == L"Κ";
  };
  auto const carries_type =
      [&](sequant::container::svector<sequant::Index> const& v,
          auto const& pred) { return std::any_of(v.begin(), v.end(), pred); };

  auto const vmap = sequant::eval::build_value_node_map(forest);

  // A concrete Κ Index instance to pass as `axis` to classify_axis --
  // classify_axis only ever compares TYPE (space().base_key()), never
  // identity, so any Κ-typed Index found in the forest is a valid probe.
  std::optional<sequant::Index> k_axis;

  // Target 1: the Κ-CONTRACTION RESULT -- a non-leaf value that does NOT
  // itself carry Κ (Q1 false) but DOES contract Κ at its own node (Q2a
  // true) => Reduction. Also record its hash so target 2 can find its
  // consumer.
  std::optional<std::size_t> mu_mu_hash;
  {
    // Two passes over the same predicate (Q1 false, Q2a true): the first
    // pass only accepts an EXACT μ̃,μ̃ carried pair (the brief's canonical
    // example); if the fixture has none, the second pass falls back to the
    // first Κ-contraction-result candidate of any shape (the
    // classify_axis == Reduction assertion itself is the load-bearing part
    // either way, not the exact carried signature).
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
        if (carries_type(vc.carried, is_K)) continue;  // Q1 would be true
        auto const contracted = sequant::contracted_indices(it->second);
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
        CHECK(sequant::eval::classify_axis(vc.carried, contracted_below,
                                           *k_axis, vc.occurrences) ==
              sequant::eval::LoopRole::Reduction);

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

        // Task 3: I(μ̃,μ̃)'s home_floor must EXCLUDE Κ -- a Reduction axis is
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

  // Target 2: the Κ-FREE COMPOSITE ABOVE THE CONTRACTION -- the immediate
  // consumer of the Κ-contraction result found above. Its OWN node does not
  // carry Κ (it was already summed away below) and does not contract Κ
  // AGAIN at its own node, so classify_axis(..., Κ, ...) must be
  // LoopInvariant -- this is exactly the case Task 1's uncorrected
  // (whole-subtree) build_site_of would have gotten wrong (it would have
  // seen Κ in the subtree and reported Reduction).
  {
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

    // Sanity: the composite must not itself carry Κ (else Q1 would route it
    // to LoopLocal/LoopCarried, not the LoopInvariant path being tested).
    REQUIRE_FALSE(carries_type(cell_it->carried, is_K));

    auto const contracted = sequant::contracted_indices(*parent);
    sequant::container::svector<sequant::Index> contracted_below(
        contracted.begin(), contracted.end());
    // Also not reduced AGAIN at this node.
    REQUIRE_FALSE(carries_type(contracted_below, is_K));

    CHECK(sequant::eval::classify_axis(cell_it->carried, contracted_below,
                                       *k_axis, cell_it->occurrences) ==
          sequant::eval::LoopRole::LoopInvariant);

    // Illustrative signature check from the brief ([i i a a]); informational
    // -- the LoopInvariant assertion above is the load-bearing one.
    INFO("I(i,i;a,a) candidate carried.size()=" << cell_it->carried.size());

    // Confirm the Task-1 correction actually landed: Κ must NOT appear in
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

    // Task 3: I(i,i;a,a)'s home_floor must EXCLUDE Κ too -- Κ is
    // LoopInvariant here (absent from build_site altogether, the IMPLICIT
    // case), so it can never surface in home_floor either. Cross-check
    // against the meet-residency home (ValueCell::home_modes).
    CHECK_FALSE(carries_type(cl_it->home_floor, is_K));
    CHECK(carries_type(cl_it->home_floor, is_K) ==
          carries_type(cell_it->home_modes, is_K));
  }

  // Target 3: a Κ-CARRYING leaf/intermediate used AT the loop variable
  // (Q1 true, Q2b lockstep) -> LoopLocal. Scan every value that itself
  // carries a Κ-typed index (leaves included -- e.g. the DF integral
  // g{Κ;μ̃,ν̃}) and classify it; the brief requires at least one LoopLocal
  // candidate be found in this fixture (REQUIRE'd below, same rigor as
  // Targets 1/2). The loop-CARRIED counterpart is deferred to Task 5's
  // synthetic fixture per the brief, so found_loop_carried is reported but
  // NOT asserted either way.
  {
    bool found_loop_local = false;
    bool found_loop_carried = false;
    for (auto const& vc : rich.cells) {
      if (!carries_type(vc.carried, is_K)) continue;
      auto const it = vmap.find(vc.hash);
      if (it == vmap.end()) continue;
      auto const contracted = sequant::contracted_indices(it->second);
      sequant::container::svector<sequant::Index> contracted_below(
          contracted.begin(), contracted.end());
      auto const role = sequant::eval::classify_axis(
          vc.carried, contracted_below, *k_axis, vc.occurrences);
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

        // Task 3: a Κ-local value's home_floor must INCLUDE Κ -- it stays
        // sliced on (homed inside) the axis it is locked to lockstep with
        // its enclosing loop. Cross-check against the meet-residency home
        // (ValueCell::home_modes); INFO records any divergence for the
        // report rather than forcing agreement.
        CHECK(carries_type(cl_it->home_floor, is_K));
        WARN("Κ-local value hash="
             << vc.hash << " leaf=" << it->second.leaf()
             << " home_floor has Κ=" << carries_type(cl_it->home_floor, is_K)
             << " home_modes has Κ=" << carries_type(vc.home_modes, is_K));
      } else if (role == sequant::eval::LoopRole::LoopCarried) {
        found_loop_carried = true;
      }
    }
    // Every Κ-carrying candidate that WAS classified LoopLocal was already
    // checked above (per_axis wiring agreement); this REQUIRE guarantees the
    // LoopLocal path is actually exercised at least once by this fixture
    // (matching Targets 1/2, which REQUIRE their own candidate was found) --
    // without it, a future shift in the fixture/optimizer output could
    // silently drop LoopLocal coverage entirely.
    INFO("water-20 Κ-carrying values: LoopLocal seen="
         << found_loop_local << " LoopCarried seen=" << found_loop_carried);
    REQUIRE(found_loop_local);
    // Loop-carried case remains deferred to Task 5's synthetic fixture per
    // the brief -- found_loop_carried is recorded (via INFO above) but not
    // asserted either way.
  }
}

// ===========================================================================
// Task 4: forced_split_axes + the monotone fixpoint, on a SYNTHETIC
// cross-iteration fixture. Water-20's aux-only residual has no LoopCarried
// value (every Κ is contracted immediately), so the loop-carried case must be
// authored here. The fixture is B{i_3,i_4} = A{;i_3} * A{;i_4}: an outer
// product that carries TWO distinct occupied indices into its own result, so
// with the occupied space made batchable in the EXTERNAL role and no enclosing
// occ batch loop realized, the occ axis survives into the value's result with
// nothing to lock it to a loop iteration -> classify_axis == LoopCarried.
// Saved to data/legality_cross_iteration.txt for Task 5 to reuse.
// ===========================================================================
TEST_CASE(
    "analyze_legality: a loop-carried value records its axis in "
    "forced_split_axes; the fixpoint terminates (synthetic cross-iteration "
    "fixture)",
    "[legality]") {
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
    return ix.space().base_key() == L"i";
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
    return ix.space().base_key() == L"i";
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
