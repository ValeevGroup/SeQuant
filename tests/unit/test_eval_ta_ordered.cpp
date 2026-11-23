#include "test_eval_ta_fixtures.hpp"

using namespace sequant::tests::eval_ta;

// Ordered executor (see ordered_schedule.hpp's own doc comments for the
// OrderedSchedule IR and
// ordered_executor.hpp's detail::run_ordered_contracted_block for the
// executor this test drives): EQUIVALENCE for a realized Contracted loop
// block (AccumulateSum reduction). A two-root forest sharing an aux-carrying
// composite S = g*h, both roots contracting the aux batch mode x_1 AT THEIR
// OWN ROOT, so each root itself is the block's AccumulateSum output -- "a
// reduction value consumed by a value built after the loop",
// realized here as the shared per-root combine step that runs after
// evaluate_ordered_schedule's loop-block walk closes -- driven through the
// pipeline (analyze_legality -> build_ordered_schedule ->
// evaluate_ordered_schedule). Real TA data at a small tractable size (the
// DryRun backend is zero-data and cannot witness a dropped/mis-accumulated
// batch).
// Numeric regression for the w8 aux-c/occ-e/occ-c loss (2026-09-11): the DP
// opens two External occ loops AND the Contracted occ loop at ONE node (the
// i_3-reducing node of g{i_3;i_1} * t{a_1,a_2;i_3,i_2}); the builder nests the
// contracted loop between the externals and the node's value is built inside
// all three, partial over the contracted loop, escaping scatter -> sum ->
// scatter. Ground truth = the unbatched evaluation.
TEST_CASE("ordered-mixed-open", "[eval][ordered-executor][mixed-open]") {
  using sequant::evaluate;
  using sequant::eval::analyze_legality;
  using sequant::eval::build_ordered_schedule;
  using sequant::eval::compute_dag_boulevard;
  using sequant::eval::evaluate_ordered_schedule;
  using sequant::eval::OrderedSchedule;
  using sequant::eval::RichSchedule;
  using TArrayD = TA::TSpArrayD;
  using node_t = sequant::FullBinaryNode<sequant::EvalExprTA>;

  auto& world = TA::get_default_world();
  // nocc 4 tiled by 2 -> two occ batches of one tile each.
  rand_tensor_yield<double, TA::SparsePolicy> yield_{world, 4, 6, 12};
  yield_.set_max_tile(2);
  auto const occ =
      sequant::get_default_context().index_basis_registry()->retrieve_space(
          L"i");

  auto const t1 = sequant::deserialize<sequant::ExprPtr>(
      L"g{i_3;i_1} * t{a_1,a_2;i_3,i_2}");
  std::vector<node_t> forest{eval_node(t1)};
  std::string target;
  for (auto const& ix : forest[0]->canon_indices()) {
    if (!target.empty()) target += ",";
    target += sequant::toUtf8(ix.full_label());
  }

  auto const ref = evaluate(forest, target, yield_)->get<TArrayD>();

  // The i_3-reducing root opens all three occ loops; its leaves carry the
  // externals they hold (the DP stamps External on every carrying node).
  auto const find = [&](std::wstring const& lbl) {
    for (auto const& ix : forest[0]->canon_indices())
      if (ix.full_label() == lbl) return ix;
    for (auto const& ix : forest[0].left()->canon_indices())
      if (ix.full_label() == lbl) return ix;
    for (auto const& ix : forest[0].right()->canon_indices())
      if (ix.full_label() == lbl) return ix;
    throw sequant::Exception{"index not found"};
  };
  auto const i1 = find(L"i_1"), i2 = find(L"i_2"), i3 = find(L"i_3");
  using sequant::BatchModeType;
  forest[0]->set_node_slice_mask({{i1, BatchModeType::External},
                                  {i2, BatchModeType::External},
                                  {i3, BatchModeType::Contracted}});
  forest[0]->set_batch_loops_opened_here({{i1, BatchModeType::External},
                                          {i2, BatchModeType::External},
                                          {i3, BatchModeType::Contracted}});
  for (auto* leaf : {&forest[0].left(), &forest[0].right()}) {
    sequant::container::svector<std::pair<sequant::Index, BatchModeType>> m;
    for (auto const& ix : (*leaf)->canon_indices())
      if (ix == i1 || ix == i2) m.push_back({ix, BatchModeType::External});
    if (!m.empty()) (*leaf)->set_node_slice_mask(m);
  }
  sequant::eval::stamp_lifetime_masks(forest);

  sequant::BatchPolicy policy;
  policy.is_batchable_contracted_index = [occ](sequant::Index const& ix) {
    return ix.space() == occ;
  };
  policy.is_batchable_external_index = [occ](sequant::Index const& ix) {
    return ix.space() == occ;
  };

  sequant::eval::dryrun::SizeRegime const regime;
  sequant::eval::dryrun::CostModel const cm{regime};
  auto const block_of = [](sequant::Index const&) -> std::size_t { return 2; };
  RichSchedule rich = compute_dag_boulevard(forest, cm, block_of);
  // Force the culprit's nesting: External > Contracted > External. Same-node
  // opens carry no order of their own (the builder tie-breaks by slot), so
  // impose it through loop_order exactly as cross-term constraints did in
  // the w8 schedule.
  {
    std::vector<int> ext, con;
    for (auto const& [key, kind] : rich.loop_kind)
      if (key.first == L"i")
        (kind == BatchModeType::External ? ext : con).push_back(key.second);
    REQUIRE(ext.size() == 2);
    REQUIRE(con.size() == 1);
    std::pair<std::wstring, int> const a{L"i", ext[0]}, c{L"i", con[0]},
        b{L"i", ext[1]};
    rich.loop_order.emplace(std::make_pair(a, c), std::size_t{0});
    rich.loop_order.emplace(std::make_pair(c, b), std::size_t{0});
  }
  auto const legality = analyze_legality(rich, forest, policy);
  OrderedSchedule const ordered =
      build_ordered_schedule(rich, legality, policy, {L"i"});
  if (std::getenv("SEQUANT_UT_SCHED_TREE"))
    sequant::eval::detail::dump_schedule_tree(ordered.root, 0);

  auto cache = sequant::CacheManager<node_t>::empty();
  auto aops = yield_.array_ops();
  cache.set_array_ops(&aops);
  std::function<std::size_t(sequant::Index const&)> const target_batch =
      [](sequant::Index const&) -> std::size_t { return 2; };
  auto const got = evaluate_ordered_schedule(forest, ordered, rich, target,
                                             yield_, cache, target_batch)
                       ->get<TArrayD>();

  TArrayD diff;
  diff(target) = got(target) - ref(target);
  double const rel = TA::norm2(diff) / TA::norm2(ref);
  INFO("relative L2 diff = " << rel);
  CHECK(rel < 1e-10);
}

TEST_CASE("ordered-mixed-open-aux-sparse",
          "[eval][ordered-executor][mixed-open][mixed-open-pia-sparse]") {
  using sequant::evaluate;
  using sequant::eval::analyze_legality;
  using sequant::eval::build_ordered_schedule;
  using sequant::eval::compute_dag_boulevard;
  using sequant::eval::evaluate_ordered_schedule;
  using sequant::eval::OrderedSchedule;
  using sequant::eval::RichSchedule;
  using TArrayD = TA::DistArray<TA::Tensor<double>, TA::SparsePolicy>;
  using node_t = sequant::FullBinaryNode<sequant::EvalExprTA>;

  // The fixture's a_3 sits in the BRA of both h and t, which trips the
  // tensor-network canonicalizer's strict bra<->ket policy (an internal index
  // may appear at most once per side under Conjugate). This case is about the
  // array operations of a mixed same-node open, not about the canonicalizer's
  // covariance assumptions, so relax the policy for its scope -- the same
  // treatment the other multi-bra fixtures in this file already use.
  auto ctx_relaxed = sequant::get_default_context_snapshot();
  ctx_relaxed.set(sequant::AssertStrictBraKetSymmetry::No);
  auto const ctx_resetter =
      sequant::set_scoped_default_context(std::move(ctx_relaxed));
  auto& world = TA::get_default_world();
  // nocc 4 tiled by 2 -> two occ batches of one tile each.
  auto const envi = [](char const* k, std::size_t d) {
    char const* e = std::getenv(k);
    return e ? static_cast<std::size_t>(std::atol(e)) : d;
  };
  std::size_t const nocc = envi("SEQUANT_UT_FIX_NOCC", 4),
                    tile = envi("SEQUANT_UT_FIX_TILE", 2),
                    occblk = envi("SEQUANT_UT_FIX_OCCBLK", 2),
                    naux = envi("SEQUANT_UT_FIX_NAUX", 12),
                    auxblk = envi("SEQUANT_UT_FIX_AUXBLK", 4);
  rand_tensor_yield<double, TA::SparsePolicy> yield_{world, nocc, 6, naux};
  yield_.set_max_tile(tile);
  auto const occ =
      sequant::get_default_context().index_basis_registry()->retrieve_space(
          L"i");

  auto const t1 = sequant::deserialize<sequant::ExprPtr>(
      L"(g{i_3;i_1;x_1} * h{a_3;a_1;x_1}) * t{a_3,a_2;i_3,i_2}");
  std::vector<node_t> forest{eval_node(t1)};
  std::string target;
  for (auto const& ix : forest[0]->canon_indices()) {
    if (!target.empty()) target += ",";
    target += sequant::toUtf8(ix.full_label());
  }

  auto const ref = evaluate(forest, target, yield_)->get<TArrayD>();

  // The i_3-reducing root opens all three occ loops; its leaves carry the
  // externals they hold (the DP stamps External on every carrying node).
  auto const find = [&](std::wstring const& lbl) {
    std::optional<sequant::Index> out;
    forest[0].visit([&](node_t const& n) {
      for (auto const& ix : n->canon_indices())
        if (ix.full_label() == lbl) out = ix;
    });
    if (!out) throw sequant::Exception{"index not found"};
    return *out;
  };
  auto const i1 = find(L"i_1"), i2 = find(L"i_2"), i3 = find(L"i_3");
  using sequant::BatchModeType;
  forest[0]->set_node_slice_mask({{i1, BatchModeType::External},
                                  {i2, BatchModeType::External},
                                  {i3, BatchModeType::Contracted}});
  forest[0]->set_batch_loops_opened_here({{i1, BatchModeType::External},
                                          {i2, BatchModeType::External},
                                          {i3, BatchModeType::Contracted}});
  auto const aux =
      sequant::get_default_context().index_basis_registry()->retrieve_space(
          L"x");
  forest[0].visit([&](node_t const& cn) {
    if (&cn == &forest[0]) return;
    auto& n = const_cast<node_t&>(cn);
    sequant::container::svector<std::pair<sequant::Index, BatchModeType>> m;
    for (auto const& ix : n->canon_indices())
      if (ix == i1 || ix == i2) m.push_back({ix, BatchModeType::External});
    if (!n.leaf()) {
      // the x_1-contracting node opens the aux loop
      for (auto const& ix : n.left()->canon_indices())
        if (ix.space() == aux) {
          m.push_back({ix, BatchModeType::Contracted});
          n->set_batch_loops_opened_here({{ix, BatchModeType::Contracted}});
        }
    }
    if (!m.empty()) n->set_node_slice_mask(m);
  });
  sequant::eval::stamp_lifetime_masks(forest);

  sequant::BatchPolicy policy;
  policy.is_batchable_contracted_index = [occ, aux](sequant::Index const& ix) {
    return ix.space() == occ || ix.space() == aux;
  };
  policy.is_batchable_external_index = [occ](sequant::Index const& ix) {
    return ix.space() == occ;
  };

  sequant::eval::dryrun::SizeRegime const regime;
  sequant::eval::dryrun::CostModel const cm{regime};
  auto const block_of = [occ, occblk,
                         auxblk](sequant::Index const& ix) -> std::size_t {
    return ix.space() == occ ? occblk : auxblk;
  };
  RichSchedule rich = compute_dag_boulevard(forest, cm, block_of);
  // Force the culprit's nesting: External > Contracted > External. Same-node
  // opens carry no order of their own (the builder tie-breaks by slot), so
  // impose it through loop_order exactly as cross-term constraints did in
  // the w8 schedule.
  {
    std::vector<int> ext, con;
    for (auto const& [key, kind] : rich.loop_kind)
      if (key.first == L"i")
        (kind == BatchModeType::External ? ext : con).push_back(key.second);
    REQUIRE(ext.size() == 2);
    REQUIRE(con.size() == 1);
    std::pair<std::wstring, int> const a{L"i", ext[0]}, c{L"i", con[0]},
        b{L"i", ext[1]};
    rich.loop_order.emplace(std::make_pair(a, c), std::size_t{0});
    rich.loop_order.emplace(std::make_pair(c, b), std::size_t{0});
  }
  auto const legality = analyze_legality(rich, forest, policy);
  OrderedSchedule const ordered =
      build_ordered_schedule(rich, legality, policy, {L"i", L"x"});
  if (std::getenv("SEQUANT_UT_SCHED_TREE"))
    sequant::eval::detail::dump_schedule_tree(ordered.root, 0);

  auto cache = sequant::CacheManager<node_t>::empty();
  auto aops = yield_.array_ops();
  cache.set_array_ops(&aops);
  std::function<std::size_t(sequant::Index const&)> const target_batch =
      [occ, occblk, auxblk](sequant::Index const& ix) -> std::size_t {
    return ix.space() == occ ? occblk : auxblk;
  };
  auto const got = evaluate_ordered_schedule(forest, ordered, rich, target,
                                             yield_, cache, target_batch)
                       ->get<TArrayD>();

  TArrayD diff;
  diff(target) = got(target) - ref(target);
  double const rel = TA::norm2(diff) / TA::norm2(ref);
  INFO("relative L2 diff = " << rel);
  CHECK(rel < 1e-10);
}

// Nested-array (ToT, pair composites) replica of the w8 aux-c/occ-e/occ-c
// culprit: the occupied-Fock-like CSV term with the amplitude t{a<i_2,i_3>,
// a<i_2,i_3>; i_3, i_2} read sliced on BOTH outer modes (external i_2,
// contracted i_3) inside external i_1 > contracted i_3 > external i_2, the
// reduction node N1 opening all three loops, the aux-contracted integral
// product invariant on the inner external loop. Ground truth = the unbatched
// evaluation.
TEST_CASE("ordered-tot-culprit",
          "[eval][ordered-executor][mixed-open][mixed-open-tot]") {
  using sequant::evaluate;
  using sequant::eval::analyze_legality;
  using sequant::eval::build_ordered_schedule;
  using sequant::eval::compute_dag_boulevard;
  using sequant::eval::evaluate_ordered_schedule;
  using sequant::eval::OrderedSchedule;
  using sequant::eval::RichSchedule;
  using ToTArray = TA::DistArray<tot_tile_t<double>, TA::SparsePolicy>;
  using node_t = sequant::FullBinaryNode<sequant::EvalExprTA>;
  using sequant::BatchModeType;

  // the rank-(2,2) non-symmetric amplitude trips the strict-braket assertion
  // of the tensor-network canonicalizer in Debug; the fixture is about the
  // array operations, so relax it here.
  auto ctx_relaxed = sequant::get_default_context_snapshot();
  ctx_relaxed.set(sequant::AssertStrictBraKetSymmetry::No);
  auto const ctx_resetter =
      sequant::set_scoped_default_context(std::move(ctx_relaxed));
  auto& world = TA::get_default_world();
  // nocc 4 tiled by 2 -> two occ batches; aux 8 tiled by 4 -> two aux batches.
  rand_tensor_yield<double, TA::SparsePolicy> yield_{world, 4, 4, 8};
  yield_.set_max_tile(2);
  auto const ibr = sequant::get_default_context().index_basis_registry();
  auto const occ = ibr->retrieve_space(L"i");
  auto const aux = ibr->retrieve_space(L"x");

  // P = g*h (K node) ; Q = (P*C)*C{a_1;m_1} (carries a_1<i_1,i_2>, so i_2
  // survives the i_3 reduction as in the MPQC tree) ; U = t*C ; N1 = Q*U
  // (reduces i_3, a_3) ; R = (N1*s)*C{m_4;a_2}.
  auto const t1 = sequant::deserialize<sequant::ExprPtr>(
      L"(((((g{i_3;i_1;x_1} * h{m_1;m_2;x_1}) * C{m_2;a_3<i_2,i_3>}) * "
      L"C{a_1<i_1,i_2>;m_1}) * "
      L"(t{a_3<i_2,i_3>,a_4<i_2,i_3>;i_3,i_2} * C{a_4<i_2,i_3>;m_3})) * "
      L"s{m_3;m_4}) * C{m_4;a_2<i_1,i_2>}");
  std::vector<node_t> forest{eval_node(t1)};
  std::string const target = forest.front()->annot();

  auto const ref = evaluate(forest, target, yield_)->get<ToTArray>();
  REQUIRE(TA::norm2(ref) > 0.0);

  auto const find = [&](std::wstring const& lbl) {
    std::optional<sequant::Index> out;
    forest[0].visit([&](node_t const& n) {
      for (auto const& ix : n->canon_indices())
        if (ix.full_label() == lbl) out = ix;
    });
    if (!out) throw sequant::Exception{"index not found"};
    return *out;
  };
  auto const i1 = find(L"i_1"), i2 = find(L"i_2"), i3 = find(L"i_3");
  // N1 = the node whose result carries no i_3 but whose children do.
  node_t* n1 = nullptr;
  forest[0].visit([&](node_t const& cn) {
    auto& n = const_cast<node_t&>(cn);
    if (n.leaf()) return;
    auto const has = [](node_t const& x, sequant::Index const& ix) {
      for (auto const& c : x->canon_indices())
        if (c == ix) return true;
      return false;
    };
    if (!has(n, i3) && has(n.left(), i3) && has(n.right(), i3)) n1 = &n;
  });
  REQUIRE(n1 != nullptr);
  (*n1)->set_node_slice_mask({{i1, BatchModeType::External},
                              {i2, BatchModeType::External},
                              {i3, BatchModeType::Contracted}});
  (*n1)->set_batch_loops_opened_here({{i1, BatchModeType::External},
                                      {i2, BatchModeType::External},
                                      {i3, BatchModeType::Contracted}});
  forest[0].visit([&](node_t const& cn) {
    auto& n = const_cast<node_t&>(cn);
    if (&n == n1) return;
    sequant::container::svector<std::pair<sequant::Index, BatchModeType>> m;
    for (auto const& ix : n->canon_indices())
      if (ix == i1 || ix == i2) m.push_back({ix, BatchModeType::External});
    if (!n.leaf()) {
      for (auto const& ix : n.left()->canon_indices())
        if (ix.space() == aux) {
          bool carried = false;
          for (auto const& c : n->canon_indices()) carried |= (c == ix);
          if (!carried) {
            m.push_back({ix, BatchModeType::Contracted});
            n->set_batch_loops_opened_here({{ix, BatchModeType::Contracted}});
          }
        }
    }
    if (!m.empty()) n->set_node_slice_mask(m);
  });
  sequant::eval::stamp_lifetime_masks(forest);

  sequant::BatchPolicy policy;
  policy.is_batchable_contracted_index = [occ, aux](sequant::Index const& ix) {
    return !ix.has_proto_indices() && (ix.space() == occ || ix.space() == aux);
  };
  policy.is_batchable_external_index = [occ](sequant::Index const& ix) {
    return !ix.has_proto_indices() && ix.space() == occ;
  };

  sequant::eval::dryrun::SizeRegime const regime;
  sequant::eval::dryrun::CostModel const cm{regime};
  auto const block_of = [occ](sequant::Index const& ix) -> std::size_t {
    return ix.space() == occ ? 2 : 4;
  };
  RichSchedule rich = compute_dag_boulevard(forest, cm, block_of);
  {
    std::vector<int> ext, con;
    for (auto const& [key, kind] : rich.loop_kind)
      if (key.first == L"i")
        (kind == BatchModeType::External ? ext : con).push_back(key.second);
    REQUIRE(ext.size() == 2);
    REQUIRE(con.size() == 1);
    std::pair<std::wstring, int> const a{L"i", ext[0]}, c{L"i", con[0]},
        b{L"i", ext[1]};
    rich.loop_order.emplace(std::make_pair(a, c), std::size_t{0});
    rich.loop_order.emplace(std::make_pair(c, b), std::size_t{0});
  }
  auto const legality = analyze_legality(rich, forest, policy);
  OrderedSchedule const ordered =
      build_ordered_schedule(rich, legality, policy, {L"i", L"x"});
  if (std::getenv("SEQUANT_UT_SCHED_TREE"))
    sequant::eval::detail::dump_schedule_tree(ordered.root, 0);

  auto cache = sequant::CacheManager<node_t>::empty();
  auto aops = yield_.template array_ops<ToTArray>();
  cache.set_array_ops(&aops);
  std::function<std::size_t(sequant::Index const&)> const target_batch =
      [occ](sequant::Index const& ix) -> std::size_t {
    return ix.space() == occ ? 2 : 4;
  };
  auto const got = evaluate_ordered_schedule(forest, ordered, rich, target,
                                             yield_, cache, target_batch)
                       ->get<ToTArray>();

  ToTArray diff;
  diff(target) = got(target) - ref(target);
  double const rel =
      std::sqrt(diff(target).dot(diff(target)) / ref(target).dot(ref(target)));
  INFO("relative L2 diff = " << rel);
  if (std::getenv("SEQUANT_UT_REPLICA_DIAG"))
    std::cerr << "[replica] batched ordered vs forest: rel=" << rel
              << std::endl;

  // Control: the same DAG through the ordered executor with NO batching
  // (policy declines every index) must match the forest to FP noise; if it
  // does not, the discrepancy is in the ordered path itself, not in slicing.
  {
    sequant::BatchPolicy none;
    none.is_batchable_contracted_index = [](sequant::Index const&) {
      return false;
    };
    none.is_batchable_external_index = [](sequant::Index const&) {
      return false;
    };
    std::vector<node_t> forest2{eval_node(t1)};
    RichSchedule const rich2 = compute_dag_boulevard(forest2, cm, block_of);
    auto const legality2 = analyze_legality(rich2, forest2, none);
    OrderedSchedule const ordered2 =
        build_ordered_schedule(rich2, legality2, none, {});
    auto cache2 = sequant::CacheManager<node_t>::empty();
    cache2.set_array_ops(&aops);
    auto const got2 =
        evaluate_ordered_schedule(forest2, ordered2, rich2, target, yield_,
                                  cache2, target_batch)
            ->get<ToTArray>();
    ToTArray diff2;
    diff2(target) = got2(target) - ref(target);
    double const rel2 = std::sqrt(diff2(target).dot(diff2(target)) /
                                  ref(target).dot(ref(target)));
    if (std::getenv("SEQUANT_UT_REPLICA_DIAG"))
      std::cerr << "[replica] UNBATCHED ordered vs forest: rel=" << rel2
                << std::endl;
    CHECK(rel2 < 1e-10);
  }
  CHECK(rel < 1e-10);
}

TEST_CASE("ordered-vs-forest-contracted-block", "[eval][ordered-executor]") {
  using sequant::evaluate;
  using sequant::eval::analyze_legality;
  using sequant::eval::build_ordered_schedule;
  using sequant::eval::compute_dag_boulevard;
  using sequant::eval::evaluate_ordered_schedule;
  using sequant::eval::OrderedSchedule;
  using sequant::eval::RichSchedule;
  using sequant::eval::ScopeBlock;
  using TArrayD = TA::TSpArrayD;
  using node_t = sequant::FullBinaryNode<sequant::EvalExprTA>;

  auto& world = TA::get_default_world();
  // aux multi-tiled (naux 12 in tiles of 4 -> 3 tiles) so x_1 batches.
  rand_tensor_yield<double, TA::SparsePolicy> yield_{world, 4, 6, 12};
  yield_.set_max_tile(4);

  auto const aux =
      sequant::get_default_context().index_basis_registry()->retrieve_space(
          L"x");

  // Two summable roots sharing S = g*h (carries aux x_1, free at the child,
  // contracted at each root -- an AccumulateSum output of the x_1 block):
  //   F1 = (g{a2;i1;x1} * h{i3;a2}) * (p{a3;i2;x1} * q{i4;a3})
  //   F2 = (g{a2;i1;x1} * h{i3;a2}) * (r{a3;i2;x1} * w{i4;a3})
  // both contract x_1 (aux-aux) -> results over {i1,i2,i3,i4}.
  auto const t1 = sequant::deserialize<sequant::ExprPtr>(
      L"(g{a_2;i_1;x_1} * h{i_3;a_2}) * (p{a_3;i_2;x_1} * q{i_4;a_3})");
  auto const t2 = sequant::deserialize<sequant::ExprPtr>(
      L"(g{a_2;i_1;x_1} * h{i_3;a_2}) * (r{a_3;i_2;x_1} * w{i_4;a_3})");
  std::string const target = "i_1,i_2,i_3,i_4";
  std::vector<node_t> forest{eval_node(t1), eval_node(t2)};

  // Reference: unbatched forest descent (ground truth), taken BEFORE
  // node_slice_mask() is stamped below -- forest descent ignores it regardless
  // (no custom evaluator installed), but this keeps the reference call
  // manifestly independent of the annotation.
  auto const ref = evaluate(forest, target, yield_)->get<TArrayD>();

  // compute_dag_boulevard (peak_profile.hpp) builds each occurrence's
  // enclosing-loop context (OccurrenceRec::ectx) off node_slice_mask() during
  // its descent -- the SAME annotation the runtime batched evaluator consults,
  // and the one build_ordered_schedule's own fixture
  // (test_ordered_schedule.cpp's W/{Κ} forest) stamps on the node that REALIZES
  // the loop. Without it no occurrence has an enclosing x_1 loop, so
  // classify_axis's LoopLocal check (bound lockstep to an enclosing loop of the
  // SAME index) can never find one and every x_1-carrying value falls to
  // LoopCarried (AccumulateScatter) instead of Reduction/LoopLocal -- exactly
  // the wrong classification this test's own precondition check caught. Each
  // root itself is the node that contracts x_1 (batch_axis finds it as the
  // shared index between the root's two factors), so node_slice_mask() goes
  // there.
  auto accept_aux = [aux](sequant::Index const& ix) {
    return ix.space() == aux;
  };
  sequant::Index x1;
  for (auto& nd : forest) {
    auto const ax = sequant::batch_axis(nd, accept_aux);
    REQUIRE(ax.has_value());
    x1 = *ax;
    nd->set_node_slice_mask({{*ax, sequant::BatchModeType::Contracted}});
    // This node realizes (opens) the loop, so mirror the DP: peak_profile's
    // ectx is now built from batch_loops_opened_here, not node_slice_mask.
    nd->set_batch_loops_opened_here(
        {{*ax, sequant::BatchModeType::Contracted}});
  }
  REQUIRE(x1.space() == aux);

  // x_1 (the aux space) is the only batchable Contracted axis; no
  // External axis at all -- the natural fixture for a Contracted-only block.
  sequant::BatchPolicy policy;
  policy.is_batchable_contracted_index = [aux](sequant::Index const& ix) {
    return ix.space() == aux;
  };
  policy.is_batchable_external_index = [](sequant::Index const&) {
    return false;
  };

  sequant::eval::dryrun::SizeRegime const regime;
  sequant::eval::dryrun::CostModel const cm{regime};
  auto const block_of = [](sequant::Index const&) -> std::size_t { return 4; };
  RichSchedule const rich = compute_dag_boulevard(forest, cm, block_of);
  auto const legality = analyze_legality(rich, forest, policy);
  OrderedSchedule const ordered =
      build_ordered_schedule(rich, legality, policy, {L"x"});

  // Precondition this test exists to exercise: a realized Contracted loop
  // block at the root, with at least one AccumulateSum output (the shared
  // reduction each root's own build resolves to) -- the shape the unbatched
  // tests do NOT produce.
  ScopeBlock const* x_block = nullptr;
  for (auto const& step : ordered.root.steps)
    if (auto const* b = std::get_if<ScopeBlock>(&step.value)) x_block = b;
  REQUIRE(x_block != nullptr);
  REQUIRE(x_block->kind == sequant::BatchModeType::Contracted);
  REQUIRE(!x_block->outputs.empty());
  for (auto const& [vid, kind] : x_block->outputs)
    REQUIRE(kind == sequant::eval::OutputKind::AccumulateSum);

  // The x_1 block must home the shared composite S as a plain BuildStep
  // (LoopLocal).
  std::size_t n_build_ids = 0;
  for (auto const& step : x_block->steps)
    if (std::holds_alternative<sequant::eval::BuildStep>(step.value))
      ++n_build_ids;
  REQUIRE(n_build_ids > 0);

  // evaluate_ordered_schedule over the single x_1 loop block.
  auto cache = sequant::CacheManager<node_t>::empty();
  auto aops = yield_.array_ops();
  cache.set_array_ops(&aops);
  std::function<std::size_t(sequant::Index const&)> const target_batch =
      [](sequant::Index const&) -> std::size_t { return 4; };
  auto const got = evaluate_ordered_schedule(forest, ordered, rich, target,
                                             yield_, cache, target_batch)
                       ->get<TArrayD>();

  // Agreement to within FP noise (the batched accumulation reorders the
  // x_1 contraction vs the unbatched reference).
  TArrayD diff;
  diff(target) = got(target) - ref(target);
  double const rel = TA::norm2(diff) / TA::norm2(ref);
  INFO("relative L2 diff = " << rel);
  CHECK(rel < 1e-10);
}

// Numerical equivalence of evaluate_ordered_schedule vs forest descent for an
// UNBATCHED flat forest: NO batch axis is stamped, so build_ordered_schedule
// emits only root-level BuildSteps and no ScopeBlock. This exercises the
// ordered executor's pure root-walk + combine_forest_roots + homing path with
// real TA data -- the path the water-8 CSV-CCk run diverges on even with
// batching OFF, and which the batched flat tests above never cover. If the
// ordered result disagrees with forest descent here, the bug is in that
// unbatched core, not in the batched loop path or in ToT tiles.
TEST_CASE("ordered-vs-forest-flat", "[eval][ordered-executor]") {
  using sequant::evaluate;
  using sequant::eval::analyze_legality;
  using sequant::eval::build_ordered_schedule;
  using sequant::eval::compute_dag_boulevard;
  using sequant::eval::evaluate_ordered_schedule;
  using sequant::eval::OrderedSchedule;
  using sequant::eval::RichSchedule;
  using TArrayD = TA::TSpArrayD;
  using node_t = sequant::FullBinaryNode<sequant::EvalExprTA>;

  auto& world = TA::get_default_world();
  rand_tensor_yield<double, TA::SparsePolicy> yield_{world, 4, 6, 12};
  yield_.set_max_tile(4);

  // Two summable roots sharing S = g*h, both contracting the aux x_1 (an
  // aux-aux edge), results over {i1,i2,i3,i4} -- same forest as the batched
  // Contracted test, but here NOT stamped node_slice_mask and NOT made
  // batchable, so the schedule is a flat sequence of BuildSteps (no loop
  // block).
  auto const t1 = sequant::deserialize<sequant::ExprPtr>(
      L"(g{a_2;i_1;x_1} * h{i_3;a_2}) * (p{a_3;i_2;x_1} * q{i_4;a_3})");
  auto const t2 = sequant::deserialize<sequant::ExprPtr>(
      L"(g{a_2;i_1;x_1} * h{i_3;a_2}) * (r{a_3;i_2;x_1} * w{i_4;a_3})");
  std::string const target = "i_1,i_2,i_3,i_4";
  std::vector<node_t> forest{eval_node(t1), eval_node(t2)};

  auto const ref = evaluate(forest, target, yield_)->get<TArrayD>();

  // Empty policy: nothing is batchable, so no node_slice_mask is consulted and
  // the schedule has no ScopeBlock.
  sequant::BatchPolicy policy;
  policy.is_batchable_contracted_index = [](sequant::Index const&) {
    return false;
  };
  policy.is_batchable_external_index = [](sequant::Index const&) {
    return false;
  };

  sequant::eval::dryrun::SizeRegime const regime;
  sequant::eval::dryrun::CostModel const cm{regime};
  auto const block_of = [](sequant::Index const&) -> std::size_t { return 4; };
  RichSchedule const rich = compute_dag_boulevard(forest, cm, block_of);
  auto const legality = analyze_legality(rich, forest, policy);
  OrderedSchedule const ordered =
      build_ordered_schedule(rich, legality, policy, {});

  // Precondition: no ScopeBlock -- the schedule is flat root BuildSteps.
  for (auto const& step : ordered.root.steps)
    REQUIRE(!std::holds_alternative<sequant::eval::ScopeBlock>(step.value));

  auto cache = sequant::CacheManager<node_t>::empty();
  std::function<std::size_t(sequant::Index const&)> const target_batch =
      [](sequant::Index const&) -> std::size_t { return 4; };
  auto const got = evaluate_ordered_schedule(forest, ordered, rich, target,
                                             yield_, cache, target_batch)
                       ->get<TArrayD>();

  TArrayD diff;
  diff(target) = got(target) - ref(target);
  double const rel = TA::norm2(diff) / TA::norm2(ref);
  INFO("relative L2 diff (unbatched ordered vs forest descent) = " << rel);
  CHECK(rel < 1e-12);
}

// Numerical equivalence of evaluate_ordered_schedule vs forest descent for a
// ToT (Tensor-of-Tensor / composite-proto-index) forest, UNBATCHED. This is the
// CSV-CCk tile type: composite legs C{m;a<i,j>} make the result a nested ToT.
// The flat-tensor ordered tests (batched and unbatched) all pass, and the
// water-8 CSV-CCk run diverges even with batching OFF, so this reproduces that
// divergence in the fast unit suite: if the ordered ToT result disagrees with
// forest descent, the bug is in the ordered executor's ToT path (root-walk /
// combine_forest_roots / homing), independent of batching.
TEST_CASE("ordered-vs-forest-tot", "[eval][ordered-executor]") {
  using sequant::evaluate;
  using sequant::eval::analyze_legality;
  using sequant::eval::build_ordered_schedule;
  using sequant::eval::compute_dag_boulevard;
  using sequant::eval::evaluate_ordered_schedule;
  using sequant::eval::OrderedSchedule;
  using sequant::eval::RichSchedule;
  using node_t = sequant::FullBinaryNode<sequant::EvalExprTA>;
  using ToTArray = TA::DistArray<tot_tile_t<double>, TA::SparsePolicy>;

  auto& world = TA::get_default_world();
  rand_tensor_yield<double, TA::SparsePolicy> yield_{world, /*nocc=*/12,
                                                     /*nvirt=*/4, /*naux=*/8};
  yield_.set_max_tile(4);

  // Two summable ToT roots (composite proto occ a<i_1,i_2>), sharing the g leg;
  // both results carry the same outer occ target, so the root-combine sums
  // them -- exercising combine_forest_roots on ToT tiles too.
  auto const t1 = sequant::deserialize<sequant::ExprPtr>(
      L"(g{m_1;m_2} * C{m_2;a1<i_1,i_2>}) * C{a2<i_1,i_2>;m_1}");
  auto const t2 = sequant::deserialize<sequant::ExprPtr>(
      L"(g{m_1;m_2} * D{m_2;a1<i_1,i_2>}) * D{a2<i_1,i_2>;m_1}");
  std::vector<node_t> forest{eval_node(t1), eval_node(t2)};
  std::string const target = forest.front()->annot();

  auto const ref = evaluate(forest, target, yield_)->get<ToTArray>();
  REQUIRE(TA::norm2(ref) > 0.0);

  sequant::BatchPolicy policy;
  policy.is_batchable_contracted_index = [](sequant::Index const&) {
    return false;
  };
  policy.is_batchable_external_index = [](sequant::Index const&) {
    return false;
  };

  sequant::eval::dryrun::SizeRegime const regime;
  sequant::eval::dryrun::CostModel const cm{regime};
  auto const block_of = [](sequant::Index const&) -> std::size_t { return 4; };
  RichSchedule const rich = compute_dag_boulevard(forest, cm, block_of);
  auto const legality = analyze_legality(rich, forest, policy);
  OrderedSchedule const ordered =
      build_ordered_schedule(rich, legality, policy, {});
  for (auto const& step : ordered.root.steps)
    REQUIRE(!std::holds_alternative<sequant::eval::ScopeBlock>(step.value));

  auto cache = sequant::CacheManager<node_t>::empty();
  std::function<std::size_t(sequant::Index const&)> const target_batch =
      [](sequant::Index const&) -> std::size_t { return 4; };
  auto const got = evaluate_ordered_schedule(forest, ordered, rich, target,
                                             yield_, cache, target_batch)
                       ->get<ToTArray>();

  ToTArray diff;
  diff(target) = got(target) - ref(target);
  double const rel =
      std::sqrt(diff(target).dot(diff(target)) / ref(target).dot(ref(target)));
  INFO("relative L2 diff (unbatched ToT ordered vs forest descent) = " << rel);
  CHECK(rel < 1e-12);
}

// CROSS-ITERATION equivalence of evaluate_ordered_schedule vs forest descent:
// the MPQC CC pattern -- ONE cache reused across iterations, cache.reset()
// between them, and a "volatile" amplitude leaf ('t') that changes each
// iteration while every other leaf (the integrals) stays fixed. The ordered
// executor's volatile-vs-persistent homing MUST drop 't'-dependent values on
// reset() (recompute them with the new 't') and keep the 't'-independent ones.
// The water-8 CSV-CCk run's ordered path is correct on iteration 1 but wrong on
// iteration 2 (energy -1.6028 == forest descent at iter 1, then diverges), so
// this reproduces that cross-iteration divergence in the fast unit suite. If
// the second ordered result disagrees with forest descent, the ordered homing
// leaks stale 't'-dependent state (or mutates a persistent buffer) across
// reset().
TEST_CASE("ordered-vs-forest-reset", "[eval][ordered-executor]") {
  using sequant::evaluate;
  using sequant::eval::analyze_legality;
  using sequant::eval::build_ordered_schedule;
  using sequant::eval::compute_dag_boulevard;
  using sequant::eval::evaluate_ordered_schedule;
  using sequant::eval::OrderedSchedule;
  using sequant::eval::RichSchedule;
  using TArrayD = TA::TSpArrayD;
  using node_t = sequant::FullBinaryNode<sequant::EvalExprTA>;

  auto& world = TA::get_default_world();
  rand_tensor_yield<double, TA::SparsePolicy> yield_{world, 8, 4, 8};
  yield_.set_max_tile(4);

  // Shared composite S = g*h (t-INDEPENDENT -> persistent, survives reset).
  // Root F1 = S * t   (contains 't' -> volatile, dropped+recomputed on reset).
  // Root F2 = S * u   (no 't'       -> persistent, kept across reset).
  // Both results over {i_1,i_2}; the root-combine sums them.
  auto const t1 = sequant::deserialize<sequant::ExprPtr>(
      L"(g{i_1;a_1} * h{a_1;i_3}) * t{i_3;i_2}");
  auto const t2 = sequant::deserialize<sequant::ExprPtr>(
      L"(g{i_1;a_1} * h{a_1;i_3}) * u{i_3;i_2}");
  std::string const target = "i_1,i_2";
  std::vector<node_t> forest{eval_node(t1), eval_node(t2)};

  // node-level volatility: leaf tensor labeled 't' (as MPQC lifts
  // policy.is_volatile_leaf).
  std::function<bool(node_t const&)> const is_vol = [](node_t const& n) {
    return n.leaf() && n->is_tensor() && n->as_tensor().label() == L"t";
  };
  auto const is_vol_leaf = [](node_t const& n) {
    return n.leaf() && n->is_tensor() && n->as_tensor().label() == L"t";
  };

  sequant::BatchPolicy policy;
  policy.is_batchable_contracted_index = [](sequant::Index const&) {
    return false;
  };
  policy.is_batchable_external_index = [](sequant::Index const&) {
    return false;
  };
  sequant::eval::dryrun::SizeRegime const regime;
  sequant::eval::dryrun::CostModel const cm{regime};
  auto const block_of = [](sequant::Index const&) -> std::size_t { return 4; };
  RichSchedule const rich = compute_dag_boulevard(forest, cm, block_of);
  auto const legality = analyze_legality(rich, forest, policy);
  OrderedSchedule const ordered =
      build_ordered_schedule(rich, legality, policy, {});

  std::function<std::size_t(sequant::Index const&)> const target_batch =
      [](sequant::Index const&) -> std::size_t { return 4; };

  // ONE cache reused across both "iterations", persistence-classified by the
  // SAME volatile predicate MPQC uses (shared_cache_for wraps this).
  auto cache = sequant::cache_manager(forest, is_vol_leaf);

  // ---- iteration 1 (t = A) ----
  auto const ref1 = evaluate(forest, target, yield_)->get<TArrayD>();
  auto const got1 =
      evaluate_ordered_schedule(forest, ordered, rich, target, yield_, cache,
                                target_batch, {}, is_vol)
          ->get<TArrayD>();
  {
    TArrayD d;
    d(target) = got1(target) - ref1(target);
    double const rel1 = TA::norm2(d) / TA::norm2(ref1);
    INFO("iter-1 relative L2 diff = " << rel1);
    CHECK(rel1 < 1e-12);
  }

  // ---- between iterations: reset the cache and change ONLY 't' ----
  cache.reset();
  yield_.reset_label(L"t");

  // ---- iteration 2 (t = B, integrals unchanged) ----
  auto const ref2 = evaluate(forest, target, yield_)->get<TArrayD>();
  auto const got2 =
      evaluate_ordered_schedule(forest, ordered, rich, target, yield_, cache,
                                target_batch, {}, is_vol)
          ->get<TArrayD>();
  TArrayD d2;
  d2(target) = got2(target) - ref2(target);
  double const rel2 = TA::norm2(d2) / TA::norm2(ref2);
  INFO("iter-2 relative L2 diff (after reset + changed 't') = " << rel2);
  CHECK(rel2 < 1e-12);
}

// The water-8 intersection: a BATCHED (aux/Contracted) loop AND cross-iteration
// cache reuse with a changed volatile 't'. All prior ordered tests pass; the
// batched ones are single-eval and the cross-iteration one is unbatched, so
// this combines them -- the untested overlap where the batched homing / scratch
// might leak state across reset(). A shared composite S=g*h is persistent
// (survives reset); root F1 = S*(t*q) contracts aux x_1 AND contains 't' (so
// volatile, recomputed with the new 't'); root F2 = S*(p*q) has no 't'
// (persistent). If ordered disagrees with forest descent on iteration 2, the
// batched path leaks stale 't'-dependent state across reset.
//
// It runs once per cache reuse threshold (min_repeats): 1 registers every
// eligible internal node in the cache's node-keyed map, 2 only the repeated
// ones plus the persistent frontier, and 1000000 the frontier alone. The
// ordered executor takes values, lives and persistence from its cell table
// and the persistent store, not from that map, so its result must not depend
// on the threshold -- which is what lets MPQC hand it the same
// eval:cache:min_repeats (default 2) as forest descent.
// The ordered executor with values kept on disk (eval::DiskSpill) on the
// TiledArray backend: same numbers as forest descent, with the qualifying
// values really written to and read back from scratch files.
TEST_CASE("ordered-disk-spill-ta", "[eval][ordered-executor][disk]") {
  using sequant::evaluate;
  using sequant::eval::analyze_legality;
  using sequant::eval::build_ordered_schedule;
  using sequant::eval::compute_dag_boulevard;
  using sequant::eval::evaluate_ordered_schedule;
  using sequant::eval::OrderedSchedule;
  using sequant::eval::RichSchedule;
  using TArrayD = TA::TSpArrayD;
  using node_t = sequant::FullBinaryNode<sequant::EvalExprTA>;

  auto& world = TA::get_default_world();
  rand_tensor_yield<double, TA::SparsePolicy> yield_{world, 4, 6, 12};
  yield_.set_max_tile(4);
  auto const aux =
      sequant::get_default_context().index_basis_registry()->retrieve_space(
          L"x");

  // two roots sharing the aux-carrying intermediate g*h; x_1 is batched
  auto const t1 = sequant::deserialize<sequant::ExprPtr>(
      L"(g{a_2;i_1;x_1} * h{i_3;a_2}) * (t{a_3;i_2;x_1} * q{i_4;a_3})");
  auto const t2 = sequant::deserialize<sequant::ExprPtr>(
      L"(g{a_2;i_1;x_1} * h{i_3;a_2}) * (p{a_3;i_2;x_1} * q{i_4;a_3})");
  std::string const target = "i_1,i_2,i_3,i_4";
  std::vector<node_t> forest{eval_node(t1), eval_node(t2)};
  auto accept_aux = [aux](sequant::Index const& ix) {
    return ix.space() == aux;
  };
  for (auto& nd : forest) {
    auto const ax = sequant::batch_axis(nd, accept_aux);
    REQUIRE(ax.has_value());
    nd->set_node_slice_mask({{*ax, sequant::BatchModeType::Contracted}});
    nd->set_batch_loops_opened_here(
        {{*ax, sequant::BatchModeType::Contracted}});
  }
  std::function<bool(node_t const&)> const is_vol = [](node_t const& n) {
    return n.leaf() && n->is_tensor() && n->as_tensor().label() == L"t";
  };
  sequant::BatchPolicy policy;
  policy.is_batchable_contracted_index = accept_aux;
  policy.is_batchable_external_index = [](sequant::Index const&) {
    return false;
  };
  sequant::eval::dryrun::SizeRegime const regime;
  sequant::eval::dryrun::CostModel const cm{regime};
  auto const block_of = [](sequant::Index const&) -> std::size_t { return 4; };
  RichSchedule const rich = compute_dag_boulevard(forest, cm, block_of);
  auto const legality = analyze_legality(rich, forest, policy);
  OrderedSchedule const ordered =
      build_ordered_schedule(rich, legality, policy, {L"x"});
  std::function<std::size_t(sequant::Index const&)> const target_batch =
      [](sequant::Index const&) -> std::size_t { return 4; };

  auto const ref = evaluate(forest, target, yield_)->get<TArrayD>();
  auto const rel_error = [&ref](TArrayD const& got) {
    TArrayD d;
    d(std::string("0,1,2,3")) =
        got(std::string("0,1,2,3")) - ref(std::string("0,1,2,3"));
    return TA::norm2(d) / TA::norm2(ref);
  };
  auto const disk_stats = [&world] {
    world.gop.fence();
    auto const s = TA::disk_io_stats();
    std::array<std::size_t, 2> n{s.tiles_written, s.tiles_read};
    world.gop.sum(n.data(), 2);
    return n;
  };

  auto cache = sequant::cache_manager(forest, is_vol, std::size_t{2});
  auto aops = yield_.array_ops();
  cache.set_array_ops(&aops);
  sequant::eval::DiskSpill spill;
  spill.min_bytes = 1;  // everything that qualifies
  cache.set_disk_spill(&spill);

  TA::reset_disk_io_stats();
  world.gop.fence();
  for (int iteration = 0; iteration < 2; ++iteration) {
    INFO("iteration " << iteration);
    auto const got =
        evaluate_ordered_schedule(forest, ordered, rich, target, yield_, cache,
                                  target_batch, {}, is_vol)
            ->get<TArrayD>();
    CHECK(rel_error(got) < 1e-12);
    cache.reset();
  }
  // values went to disk and came back
  CHECK(spill.n_writes > 0);
  CHECK(spill.n_reads > 0);
  auto const io = disk_stats();
  CHECK(io[0] > 0);  // tiles written
  CHECK(io[1] > 0);  // tiles read
}

// A value assembled by scattering the batches of a loop into disjoint blocks,
// and read after the loop, is assembled on disk block by block once its size
// is known from an earlier evaluation (DiskSpill::measured_bytes): it is
// never whole in memory.
TEST_CASE("ordered-disk-scatter-ta", "[eval][ordered-executor][disk]") {
  using sequant::evaluate;
  using sequant::eval::analyze_legality;
  using sequant::eval::build_ordered_schedule;
  using sequant::eval::compute_dag_boulevard;
  using sequant::eval::evaluate_ordered_schedule;
  using sequant::eval::OrderedSchedule;
  using sequant::eval::OutputKind;
  using sequant::eval::RichSchedule;
  using sequant::eval::ScopeBlock;
  using TArrayD = TA::TSpArrayD;
  using node_t = sequant::FullBinaryNode<sequant::EvalExprTA>;

  auto& world = TA::get_default_world();
  rand_tensor_yield<double, TA::SparsePolicy> yield_{world, 8, 6, 12};
  yield_.set_max_tile(4);
  auto const occ =
      sequant::get_default_context().index_basis_registry()->retrieve_space(
          L"i");

  // S = g*h carries i_1 and is built inside a loop over i_1 (opened at S);
  // the root reads S whole, after the loop. g is volatile, so S is rebuilt
  // by every evaluation rather than carried over from the previous one.
  auto const t = sequant::deserialize<sequant::ExprPtr>(
      L"(g{a_2;i_1} * h{a_3;a_2}) * q{i_1;a_5}");
  std::string const target = "a_3,a_5";
  std::vector<node_t> forest{eval_node(t)};
  auto const ref = evaluate(forest, target, yield_)->get<TArrayD>();

  sequant::Index const i1{L"i_1"};
  auto& S = const_cast<node_t&>(forest.front().left());
  REQUIRE(!S.leaf());
  S->set_node_slice_mask({{i1, sequant::BatchModeType::External}});
  S->set_batch_loops_opened_here({{i1, sequant::BatchModeType::External}});

  sequant::BatchPolicy policy;
  policy.is_batchable_contracted_index = [](sequant::Index const&) {
    return false;
  };
  policy.is_batchable_external_index = [occ](sequant::Index const& ix) {
    return ix.space() == occ;
  };
  sequant::eval::dryrun::SizeRegime const regime;
  sequant::eval::dryrun::CostModel const cm{regime};
  auto const block_of = [](sequant::Index const&) -> std::size_t { return 4; };
  RichSchedule const rich = compute_dag_boulevard(forest, cm, block_of);
  auto const legality = analyze_legality(rich, forest, policy);
  OrderedSchedule const ordered =
      build_ordered_schedule(rich, legality, policy, {L"i"});
  ScopeBlock const* i_block = nullptr;
  for (auto const& step : ordered.root.steps)
    if (auto const* b = std::get_if<ScopeBlock>(&step.value)) i_block = b;
  REQUIRE(i_block != nullptr);
  REQUIRE(i_block->outputs.size() == 1);
  REQUIRE(i_block->outputs.front().second == OutputKind::AccumulateScatter);

  std::function<bool(node_t const&)> const is_vol = [](node_t const& n) {
    return n.leaf() && n->is_tensor() && n->as_tensor().label() == L"g";
  };
  auto cache = sequant::CacheManager<node_t>::empty();
  auto aops = yield_.array_ops();
  REQUIRE(static_cast<bool>(aops.make_zeros_on_disk));
  cache.set_array_ops(&aops);
  std::function<std::size_t(sequant::Index const&)> const target_batch =
      [](sequant::Index const&) -> std::size_t { return 4; };
  sequant::eval::DiskSpill spill;
  spill.min_bytes = 1;
  cache.set_disk_spill(&spill);

  auto const rel_error = [&](TArrayD const& got) {
    TArrayD d;
    d(target) = got(target) - ref(target);
    return TA::norm2(d) / TA::norm2(ref);
  };
  for (int iteration = 0; iteration < 3; ++iteration) {
    INFO("iteration " << iteration);
    auto const got =
        evaluate_ordered_schedule(forest, ordered, rich, target, yield_, cache,
                                  target_batch, {}, is_vol)
            ->get<TArrayD>();
    CHECK(rel_error(got) < 1e-12);
    // the first evaluation learns the value's size while assembling it in
    // memory; later ones assemble it on disk
    if (iteration == 0)
      CHECK(spill.n_assembled_on_disk == 0);
    else
      CHECK(spill.n_assembled_on_disk == static_cast<std::size_t>(iteration));
    cache.reset();
  }
  CHECK(spill.n_writes > 0);
  CHECK(spill.n_reads > 0);
}

TEST_CASE("ordered-vs-forest-reset-batched",
          "[eval][ordered-executor][ordered-crossiter-bug][cache_manager]") {
  using sequant::evaluate;
  using sequant::eval::analyze_legality;
  using sequant::eval::build_ordered_schedule;
  using sequant::eval::compute_dag_boulevard;
  using sequant::eval::evaluate_ordered_schedule;
  using sequant::eval::OrderedSchedule;
  using sequant::eval::RichSchedule;
  using TArrayD = TA::TSpArrayD;
  using node_t = sequant::FullBinaryNode<sequant::EvalExprTA>;

  auto& world = TA::get_default_world();
  rand_tensor_yield<double, TA::SparsePolicy> yield_{world, 4, 6, 12};
  yield_.set_max_tile(4);
  auto const aux =
      sequant::get_default_context().index_basis_registry()->retrieve_space(
          L"x");

  auto const t1 = sequant::deserialize<sequant::ExprPtr>(
      L"(g{a_2;i_1;x_1} * h{i_3;a_2}) * (t{a_3;i_2;x_1} * q{i_4;a_3})");
  auto const t2 = sequant::deserialize<sequant::ExprPtr>(
      L"(g{a_2;i_1;x_1} * h{i_3;a_2}) * (p{a_3;i_2;x_1} * q{i_4;a_3})");
  std::string const target = "i_1,i_2,i_3,i_4";
  std::vector<node_t> forest{eval_node(t1), eval_node(t2)};

  // stamp aux Contracted on each root (the node that contracts x_1).
  auto accept_aux = [aux](sequant::Index const& ix) {
    return ix.space() == aux;
  };
  for (auto& nd : forest) {
    auto const ax = sequant::batch_axis(nd, accept_aux);
    REQUIRE(ax.has_value());
    nd->set_node_slice_mask({{*ax, sequant::BatchModeType::Contracted}});
    // realizer node opens the loop; ectx is built from opens (peak_profile).
    nd->set_batch_loops_opened_here(
        {{*ax, sequant::BatchModeType::Contracted}});
  }

  std::function<bool(node_t const&)> const is_vol = [](node_t const& n) {
    return n.leaf() && n->is_tensor() && n->as_tensor().label() == L"t";
  };
  auto const is_vol_leaf = is_vol;

  sequant::BatchPolicy policy;
  policy.is_batchable_contracted_index = [aux](sequant::Index const& ix) {
    return ix.space() == aux;
  };
  policy.is_batchable_external_index = [](sequant::Index const&) {
    return false;
  };
  sequant::eval::dryrun::SizeRegime const regime;
  sequant::eval::dryrun::CostModel const cm{regime};
  auto const block_of = [](sequant::Index const&) -> std::size_t { return 4; };
  RichSchedule const rich = compute_dag_boulevard(forest, cm, block_of);
  auto const legality = analyze_legality(rich, forest, policy);
  OrderedSchedule const ordered =
      build_ordered_schedule(rich, legality, policy, {L"x"});
  bool has_block = false;
  for (auto const& step : ordered.root.steps)
    if (std::holds_alternative<sequant::eval::ScopeBlock>(step.value))
      has_block = true;
  REQUIRE(has_block);

  std::function<std::size_t(sequant::Index const&)> const target_batch =
      [](sequant::Index const&) -> std::size_t { return 4; };
  auto const min_repeats =
      GENERATE(std::size_t{1}, std::size_t{2}, std::size_t{1000000});
  INFO("min_repeats = " << min_repeats);
  auto cache = sequant::cache_manager(forest, is_vol_leaf, min_repeats);
  auto aops = yield_.array_ops();
  cache.set_array_ops(&aops);

  // The threshold really changes what the map holds. At 1 it registers the two
  // roots; S, t*q and p*q carry the batched x_1 and are refused run-scope
  // residence whatever the threshold. Above 1 nothing qualifies, so the map is
  // empty and the ordered executor must manage on its cell table and
  // persistent store alone.
  auto const n_keys = [](auto const& c) {
    std::size_t n = 0;
    c.for_each_key([&n](node_t const&) { ++n; });
    return n;
  };
  auto const n_keys_at_1 =
      n_keys(sequant::cache_manager(forest, is_vol_leaf, std::size_t{1}));
  REQUIRE(n_keys_at_1 > 0);
  if (min_repeats > 1) CHECK(n_keys(cache) < n_keys_at_1);

  auto const ref1 = evaluate(forest, target, yield_)->get<TArrayD>();
  auto const got1 =
      evaluate_ordered_schedule(forest, ordered, rich, target, yield_, cache,
                                target_batch, {}, is_vol)
          ->get<TArrayD>();
  {
    TArrayD d;
    d(target) = got1(target) - ref1(target);
    double const r = TA::norm2(d) / TA::norm2(ref1);
    INFO("iter-1 rel diff = " << r);
    CHECK(r < 1e-10);
  }

  cache.reset();
  yield_.reset_label(L"t");

  auto const ref2 = evaluate(forest, target, yield_)->get<TArrayD>();
  auto const got2 =
      evaluate_ordered_schedule(forest, ordered, rich, target, yield_, cache,
                                target_batch, {}, is_vol)
          ->get<TArrayD>();
  TArrayD d2;
  d2(target) = got2(target) - ref2(target);
  double const r2 = TA::norm2(d2) / TA::norm2(ref2);
  INFO("iter-2 rel diff (batched, after reset + changed 't') = " << r2);
  CHECK(r2 < 1e-10);
}

// Ordered executor, the External counterpart of the Contracted-block test
// (see ordered_schedule.hpp's own doc comments for the
// OrderedSchedule IR and ordered_executor.hpp's
// detail::run_ordered_contracted_block for the executor this test drives):
// EQUIVALENCE for a realized External loop block (AccumulateScatter,
// loop-carried output). Mirrors the Contracted-block test's harness almost
// exactly (same shared-composite structure, same pipeline), but
// with the batch axis a SPECTATOR that survives into each root's own result
// (occ i_1) instead of one both roots reduce away: F1/F2 = S ⊗ T_k, an OUTER
// PRODUCT between the i_1-carrying shared composite S and an all-virtual
// factor T_k that shares no index with S at all, so i_1 is never contracted
// anywhere -- it is free on F1/F2's own result, which is exactly what makes
// build_ordered_schedule classify it LoopCarried (escaping as
// AccumulateScatter) rather than Reduction (AccumulateSum, the Contracted-block
// test's shape). Every i-space index anywhere in the forest is the SAME literal
// i_1 (every other free index is virtual, a different space entirely,
// deliberately kept out of the "i" axis TYPE's bucket) so the
// single-physical-label simplification run_ordered_contracted_block's own \note
// documents holds throughout -- see that \note for why a forest naming the occ
// TYPE under more than one physical Index in one block is out of scope here.
// Real TA data at a small tractable size (the DryRun backend is zero-data and
// cannot witness a dropped/double-written scatter slice).
TEST_CASE("ordered-vs-forest-external-block", "[eval][ordered-executor]") {
  using sequant::evaluate;
  using sequant::eval::analyze_legality;
  using sequant::eval::build_ordered_schedule;
  using sequant::eval::compute_dag_boulevard;
  using sequant::eval::evaluate_ordered_schedule;
  using sequant::eval::OrderedSchedule;
  using sequant::eval::OutputKind;
  using sequant::eval::RichSchedule;
  using sequant::eval::ScopeBlock;
  using TArrayD = TA::TSpArrayD;
  using node_t = sequant::FullBinaryNode<sequant::EvalExprTA>;

  auto& world = TA::get_default_world();
  // occ multi-tiled (nocc 8 in tiles of 4 -> 2 tiles) so i_1 batches into
  // more than one disjoint scatter slice.
  rand_tensor_yield<double, TA::SparsePolicy> yield_{world, 8, 6, 12};
  yield_.set_max_tile(4);

  auto const occ =
      sequant::get_default_context().index_basis_registry()->retrieve_space(
          L"i");

  // S = g*h carries i_1 (free, LoopLocal under the loop F1/F2 realize) and a
  // virtual spectator a_3 (contracts a_2). T1/T2 are ALL-virtual (no occ
  // index anywhere), so F_k = S ⊗ T_k is an OUTER PRODUCT: i_1 shares nothing
  // to contract against and survives whole into F_k's own result --
  // LoopCarried, not Reduction.
  //   S  = g{a_2;i_1} * h{a_3;a_2}
  //   T1 = p{a_4;a_5} * q{a_6;a_4}
  //   T2 = r{a_4;a_5} * w{a_6;a_4}
  //   F1 = S * T1, F2 = S * T2  ->  results over {i_1, a_3, a_5, a_6}
  auto const t1 = sequant::deserialize<sequant::ExprPtr>(
      L"(g{a_2;i_1} * h{a_3;a_2}) * (p{a_4;a_5} * q{a_6;a_4})");
  auto const t2 = sequant::deserialize<sequant::ExprPtr>(
      L"(g{a_2;i_1} * h{a_3;a_2}) * (r{a_4;a_5} * w{a_6;a_4})");
  std::string const target = "i_1,a_3,a_5,a_6";
  std::vector<node_t> forest{eval_node(t1), eval_node(t2)};

  // Reference: unbatched forest descent (ground truth), taken BEFORE
  // node_slice_mask() is stamped below -- forest descent ignores it regardless
  // (no custom evaluator installed), but this keeps the reference call
  // manifestly independent of the annotation.
  auto const ref = evaluate(forest, target, yield_)->get<TArrayD>();

  // compute_dag_boulevard (peak_profile.hpp) builds each occurrence's
  // enclosing-loop context (OccurrenceRec::ectx) off node_slice_mask() during
  // its descent (see the Contracted-block test's own comment for the full
  // mechanics). Both roots themselves realize the i_1 loop (mirroring the
  // Contracted-block test's own precedent of stamping node_slice_mask directly
  // on each forest root): descendants below them (S, g, h) see an enclosing i_1
  // loop and lock-step to it (LoopLocal); the roots' OWN occurrences have no
  // enclosing loop at all (nothing wraps a forest root) and still carry i_1 on
  // their own result, which is exactly the LoopCarried (AccumulateScatter)
  // shape this test exists to exercise.
  sequant::Index const i1{L"i_1"};
  for (auto& nd : forest) {
    nd->set_node_slice_mask({{i1, sequant::BatchModeType::External}});
    // each root realizes (opens) the i_1 loop; ectx is built from opens.
    nd->set_batch_loops_opened_here({{i1, sequant::BatchModeType::External}});
  }

  // i_1 (the occ space) is the only batchable index, marked
  // EXTERNAL (spectator) rather than contracted -- the natural fixture for
  // an External-only block.
  sequant::BatchPolicy policy;
  policy.is_batchable_contracted_index = [](sequant::Index const&) {
    return false;
  };
  policy.is_batchable_external_index = [occ](sequant::Index const& ix) {
    return ix.space() == occ;
  };

  sequant::eval::dryrun::SizeRegime const regime;
  sequant::eval::dryrun::CostModel const cm{regime};
  auto const block_of = [](sequant::Index const&) -> std::size_t { return 4; };
  RichSchedule const rich = compute_dag_boulevard(forest, cm, block_of);
  auto const legality = analyze_legality(rich, forest, policy);
  OrderedSchedule const ordered =
      build_ordered_schedule(rich, legality, policy, {L"i"});

  // Precondition this test exists to exercise: a realized External loop
  // block at the root, with an AccumulateScatter output for EACH root (F1
  // and F2 each escape the loop on their own).
  ScopeBlock const* i_block = nullptr;
  for (auto const& step : ordered.root.steps)
    if (auto const* b = std::get_if<ScopeBlock>(&step.value)) i_block = b;
  REQUIRE(i_block != nullptr);
  REQUIRE(i_block->kind == sequant::BatchModeType::External);
  REQUIRE(i_block->outputs.size() == 2);
  for (auto const& [vid, kind] : i_block->outputs) {
    (void)vid;
    REQUIRE(kind == OutputKind::AccumulateScatter);
  }

  // The i_1 block must home the shared composite S as a plain BuildStep
  // (LoopLocal), exactly like the Contracted-block test's n_build_ids > 0
  // check.
  std::size_t n_build_ids = 0;
  for (auto const& step : i_block->steps)
    if (std::holds_alternative<sequant::eval::BuildStep>(step.value))
      ++n_build_ids;
  REQUIRE(n_build_ids > 0);

  // evaluate_ordered_schedule over the single i_1 loop block.
  auto cache = sequant::CacheManager<node_t>::empty();
  auto aops = yield_.array_ops();
  cache.set_array_ops(&aops);
  std::function<std::size_t(sequant::Index const&)> const target_batch =
      [](sequant::Index const&) -> std::size_t { return 4; };
  auto const got = evaluate_ordered_schedule(forest, ordered, rich, target,
                                             yield_, cache, target_batch)
                       ->get<TArrayD>();

  // Agreement to within FP noise (the batched scatter reassembles i_1 from
  // disjoint slices in a different order than the unbatched reference).
  TArrayD diff;
  diff(target) = got(target) - ref(target);
  double const rel = TA::norm2(diff) / TA::norm2(ref);
  INFO("relative L2 diff = " << rel);
  CHECK(rel < 1e-10);
}

TEST_CASE("ta_tot_kconjugation_end_to_end", "[eval]") {
  // The two core conjugation states on a tensor-of-tensors leaf, end to end.
  // Neither is a node shape: binarize serves a marked spelling as a plain leaf
  // whose CanonTransform carries the conjugation, and retrieval applies that
  // transform exactly once, so the engine hands back the elementwise conjugate
  // of what the yielder serves for the bare leaf.
  //
  // Both states land on that one reference. A '꙳' over a real basis with an
  // indefinite hermiticity is kept symbolically and decodes to a pure {conj}
  // with the slots in place; a '⁺' on a non-Hermitian tensor decodes to
  // {conj, braket_swap}, and a bra<->ket exchange is layout-invariant for a
  // ToT array because the outer/inner split is by proto indices, not by
  // bundle.
  using namespace sequant;
  auto& world = TA::get_default_world();
  size_t const nocc = 2, nvirt = 3;
  rand_tensor_yield<std::complex<double>, TA::DensePolicy> yield{world, nocc,
                                                                 nvirt};
  using ArrayToT = typename decltype(yield)::array_tot_type;

  // every element of a ToT array, in outer-then-inner traversal order
  auto const flatten = [](ArrayToT const& arr) {
    std::vector<std::complex<double>> out;
    for (auto it = arr.begin(); it != arr.end(); ++it) {
      auto const& outer = it->get();
      for (std::size_t o = 0; o < outer.size(); ++o) {
        auto const& inner = outer[o];
        for (std::size_t k = 0; k < inner.size(); ++k) out.push_back(inner[k]);
      }
    }
    return out;
  };

  // braket symmetry pinned explicitly (:C): the two orientations below are a
  // Conjugate (Hermitian) ToT leaf, independent of the ambient deserializer
  // defaults (which become conservative NonHermitian with the
  // default-tensor-symmetry rework, PR #596)
  auto const swapped =
      deserialize<sequant::ExprPtr>(L"t{i2,i3;a3<i2,i3>,a4<i2,i3>}:N-C-S");
  auto const canonical =
      deserialize<sequant::ExprPtr>(L"t{a3<i2,i3>,a4<i2,i3>;i2,i3}:N-C-S");

  // eval-boundary precondition: leaves keep their orientation (a Conjugate
  // tensor's bundles are never exchanged, no state) -- the two orientations
  // are distinct nodes
  EvalExpr const swapped_leaf{swapped->as<Tensor>()};
  EvalExpr const canon_leaf{canonical->as<Tensor>()};
  auto const is_conj = [](EvalExpr const& leaf) {
    return leaf.expr()->as<Tensor>().kconjugated();
  };
  REQUIRE_FALSE(is_conj(swapped_leaf));
  REQUIRE_FALSE(is_conj(canon_leaf));
  REQUIRE(swapped_leaf.hash_value() != canon_leaf.hash_value());

  SECTION("'꙳' over a real basis: a pure conjugation, slots in place") {
    // the same slot layout over a real basis, with parity None and an
    // indefinite hermiticity: there the '꙳' is kept and is a pure conjugation
    auto rsp = [](std::wstring_view label) {
      IndexSpace sp = Index(label).space();
      sp.field(Field::Real);
      return sp;
    };
    const Index i2(rsp(L"i_2"), 2);
    const Index i3(rsp(L"i_3"), 3);
    const Index a3(rsp(L"a_3"), 3, container::vector<Index>{i2, i3});
    const Index a4(rsp(L"a_4"), 4, container::vector<Index>{i2, i3});
    auto const real_canonical = ex<Tensor>(
        L"t", bra{a3, a4}, ket{i2, i3},
        TensorSymmetries{.perm = Symmetry::Nonsymm,
                         .hermiticity = Hermiticity::NonHermitian,
                         .conjugation_parity = ConjugationParity::None,
                         .column = ColumnSymmetry::Symm});
    REQUIRE(real_canonical->as<Tensor>().base_field() == Field::Real);

    auto conj_side = real_canonical->clone();
    REQUIRE(conj_side->as<Tensor>().kconjugate() == 1);
    REQUIRE(conj_side->as<Tensor>().kconjugated());
    // the slots stay in place: the conjugation is not the bundle exchange
    REQUIRE(conj_side->as<Tensor>().bra()[0].label() == L"a_3");
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto const node = binarize<EvalExprTA>(conj_side);
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    REQUIRE(node.leaf());
    REQUIRE(node->canon_transform() == CanonTransform{.conj = true});
    REQUIRE_FALSE(node->expr()->as<Tensor>().kconjugated());
    // the identity layout: the stored array's own modes
    REQUIRE(node->canon_indices() ==
            EvalExprTA{real_canonical->as<Tensor>()}.canon_indices());

    // the array served for the bare leaf, read before the evaluation so that
    // a kernel writing into it could not fake the comparison
    auto const served =
        flatten(yield(node->expr()->as<Tensor>())->get<ArrayToT>());
    REQUIRE_FALSE(served.empty());

    auto cache = CacheManager<FullBinaryNode<EvalExprTA>>::empty();
    auto const res = evaluate(node, node->annot(), yield, cache);
    auto const got = flatten(res->get<ArrayToT>());

    REQUIRE(got.size() == served.size());
    for (std::size_t n = 0; n < served.size(); ++n) {
      auto const expected = std::conj(served[n]);
      CHECK(got[n].real() == Catch::Approx(expected.real()));
      CHECK(got[n].imag() == Catch::Approx(expected.imag()));
    }
  }

  SECTION("'⁺' on a non-Hermitian leaf: {conj, braket_swap}") {
    // an indefinite hermiticity keeps the '⁺', which exchanges the bundles;
    // for a ToT array that exchange is the identity layout, so the served
    // values come back conjugated and in place
    auto const bare =
        deserialize<sequant::ExprPtr>(L"t{a3<i2,i3>,a4<i2,i3>;i2,i3}:N-N-S");
    auto adj_side = bare->clone();
    REQUIRE(adj_side->as<Tensor>().adjoint() == 1);
    REQUIRE(adj_side->as<Tensor>().adjointed());
    REQUIRE(adj_side->as<Tensor>().bra()[0].label() == L"i_2");
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto const node = binarize<EvalExprTA>(adj_side);
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    REQUIRE(node.leaf());
    REQUIRE(node->canon_transform() ==
            CanonTransform{.conj = true, .braket_swap = true});
    REQUIRE_FALSE(node->expr()->as<Tensor>().adjointed());
    // the exchange moves no ToT mode: the layout is the bare leaf's
    REQUIRE(node->canon_indices() ==
            EvalExprTA{bare->as<Tensor>()}.canon_indices());

    auto const served =
        flatten(yield(node->expr()->as<Tensor>())->get<ArrayToT>());
    REQUIRE_FALSE(served.empty());

    auto cache = CacheManager<FullBinaryNode<EvalExprTA>>::empty();
    auto const res = evaluate(node, node->annot(), yield, cache);
    auto const got = flatten(res->get<ArrayToT>());

    REQUIRE(got.size() == served.size());
    for (std::size_t n = 0; n < served.size(); ++n) {
      CHECK(got[n].real() == Catch::Approx(served[n].real()));
      CHECK(got[n].imag() == Catch::Approx(-served[n].imag()));
    }
  }
}
// ---------------------------------------------------------------------------
// External-placement correctness reproducer.
//
// Drives the FULL optimizer -> binarize -> stamp_lifetime_masks -> batched TA
// evaluate pipeline (unlike the eval_batched_* cases, which stamp batch modes
// by hand), so the External BatchModeType stamps come from the DP's own
// per-node external opens, not a fixture. This is the only case in which the
// External stamps are produced by the optimizer rather than hand-written.
// Compares the batched result against the plain unbatched reference (ground
// truth).
// ---------------------------------------------------------------------------
namespace {

// Copy the batch-annotation fields from a plain EvalExpr binarized tree
// onto the structurally-identical EvalExprTA tree (to_ta_node's tensor-branch
// reconstructor drops them). Both trees are built by the same binarize/
// transform_node post-order, so a lockstep walk aligns node-for-node.
void copy_batch_annotations(
    sequant::FullBinaryNode<sequant::EvalExpr> const& from,
    sequant::FullBinaryNode<sequant::EvalExprTA>& to) {
  const auto& mask = from->node_slice_mask();
  const_cast<sequant::EvalExprTA&>(*to).set_node_slice_mask(
      sequant::container::svector<
          std::pair<sequant::Index, sequant::BatchModeType>>(mask.begin(),
                                                             mask.end()));
  const auto& sliced = from->sliced_modes();
  const_cast<sequant::EvalExprTA&>(*to).set_sliced_modes(
      sequant::container::svector<sequant::Index>(sliced.begin(),
                                                  sliced.end()));
  const_cast<sequant::EvalExprTA&>(*to).set_batch_order_aware(
      from->batch_order_aware());
  const_cast<sequant::EvalExprTA&>(*to).set_batch_effective_count(
      from->batch_effective_count());
  if (!from.leaf()) {
    copy_batch_annotations(from.left(), to.left());
    copy_batch_annotations(from.right(), to.right());
  }
}

// Dump per-node annotations of a TA eval tree (result indices, node_slice_mask
// modes+kind, sliced_modes, order_aware) for eyeballing placement differences.
void dump_annotations(sequant::FullBinaryNode<sequant::EvalExprTA> const& n,
                      int depth = 0) {
  using sequant::BatchModeType;
  std::string pad(2 * depth, ' ');
  std::string res;
  for (auto const& ix : n->canon_indices()) {
    res += sequant::toUtf8(ix.full_label());
    res += ' ';
  }
  std::string bh;
  for (auto const& [ix, knd] : n->node_slice_mask()) {
    bh += sequant::toUtf8(ix.full_label());
    bh += (knd == BatchModeType::External ? ":ext " : ":con ");
  }
  std::string sm;
  for (auto const& ix : n->sliced_modes()) {
    sm += sequant::toUtf8(ix.full_label());
    sm += ' ';
  }
  std::string ident;
  try {
    ident = sequant::toUtf8(sequant::io::serialization::to_string(n->expr()));
  } catch (...) {
    ident = "?";
  }
  std::cerr << pad << (n.leaf() ? "LEAF " : "NODE ") << "res=[" << res
            << "] batched=[" << bh << "] sliced=[" << sm
            << "] oa=" << (n->batch_order_aware() ? 1 : 0)
            << " | expr=" << ident << "\n";
  if (!n.leaf()) {
    dump_annotations(n.left(), depth + 1);
    dump_annotations(n.right(), depth + 1);
  }
}

}  // namespace

TEST_CASE("dp_external_placement_correctness",
          "[eval][batched-external][node-placement]") {
  using namespace sequant;
  using TArrayD = TA::TSpArrayD;
  using node_t = FullBinaryNode<EvalExprTA>;
  using cache_t = CacheManager<node_t>;

  auto& world = TA::get_default_world();

  auto ibr = get_default_context().index_basis_registry();
  auto occ = ibr->retrieve_space(L"i");
  auto virt = ibr->retrieve_space(L"a");
  auto aux = ibr->retrieve_space(L"x");

  // Actual arrays: occ multi-tiled (40 in tiles of 4 -> 10 blocks) so the
  // external occ mode partitions into > 1 batch, and large vs virt/aux (4) so
  // an i_2-carrying intermediate (~i*i) is the peak and slicing i_2 lowers it,
  // forcing the DP to open i_2's external loop.
  rand_tensor_yield<double, TA::SparsePolicy> yield_{world, /*nocc=*/12,
                                                     /*nvirt=*/4, /*naux=*/4};
  yield_.set_max_tile(4);

  // Term with an external occ spectator (i_2, free on the result and carried on
  // internal nodes) plus a subtree invariant to it. Non-left-deep so the right
  // branch (w*p) carries NO i_2 (invariant), the left branch (g*h) carries i_2.
  // Flat 4-tensor product (no inner parens: those would parse as a nested
  // Product of sub-products, which optimize() treats as < 3 factors and leaves
  // unfactorized => no batch annotations). The optimizer factorizes it. i_2,i_6
  // are occ externals (free on the result, carried across internal nodes);
  // x_1,x_2 are contracted aux.
  // Particle-conserving ladder: virtuals alternate bra<->ket down the chain
  // (a_1: g ket / h bra, a_2: h ket / w bra, a_3: w ket / p bra) and the
  // contracted aux x_1,x_2 sit in the AUX slot (3-slot bra;ket;aux ->
  // Aux-origin, exempt from the strict-braket check). Topology, externals
  // (i_2,i_6) and per-space sizes are unchanged.
  std::string const expr_str =
      "g{i_2;a_1;x_1} * h{a_1;a_2;x_1} * w{a_2;a_3;x_2} * p{a_3;i_6;x_2}";
  std::string const target = "i_2,i_6";

  auto const is_occ = [occ](Index const& ix) { return ix.space() == occ; };
  auto const is_aux = [aux](Index const& ix) { return ix.space() == aux; };

  // idx_to_extent inflated so the modeled peak far exceeds peak_threshold and
  // the DP opens the external occ loop.
  std::function<std::size_t(Index const&)> idxsz =
      [occ, virt, aux](Index const& ix) -> std::size_t {
    if (ix.space() == occ) return 12;
    if (ix.space() == virt) return 4;
    if (ix.space() == aux) return 4;
    return 1;
  };
  std::function<std::size_t(Index const&)> bts =
      [occ](Index const& ix) -> std::size_t {
    return ix.space() == occ ? std::size_t{4} : std::size_t{60};
  };

  auto run = [&]() -> TArrayD {
    auto const expr =
        deserialize<ExprPtr>(std::wstring(expr_str.begin(), expr_str.end()));

    auto axes_map = std::make_shared<std::unordered_map<
        Expr const*, container::vector<NodeBatchAnnotation>>>();

    OptimizeOptions opts;
    opts.objective_function = ObjectiveFunction::DenseTimeSpaceBatched;
    opts.reorder = ReorderSum::NoReorder;
    opts.idx_to_extent = idxsz;
    opts.inner_pow = [](Index const&, std::size_t) -> double { return 1.0; };
    opts.batch_policy.is_batchable_contracted_index = is_aux;
    opts.batch_policy.is_batchable_external_index = is_occ;
    opts.batch_policy.batch_target_size = bts;
    opts.batch_policy.batch_spectator_indices = true;
    opts.batch_policy.peak_threshold = 1.0e3;  // tiny budget => force external
    opts.batch_policy.infeasible_is_error = false;  // may be over budget
    opts.term_batch_axes = axes_map;

    auto optimized = optimize(expr, opts);
    REQUIRE(optimized);
    auto it = axes_map->find(optimized.get());
    REQUIRE(it != axes_map->end());
    auto const& node_axes = it->second;
    std::cerr << "[diag] node_axes.size()=" << node_axes.size() << "\n";
    for (std::size_t j = 0; j < node_axes.size(); ++j) {
      std::cerr << "  node_axes[" << j << "] oa=" << node_axes[j].order_aware
                << " naxes=" << node_axes[j].axes.size() << " :";
      for (auto const& [ix, knd] : node_axes[j].axes)
        std::cerr << " " << toUtf8(ix.full_label())
                  << (knd == BatchModeType::External ? ":ext" : ":con");
      std::cerr << "\n";
    }

    BinarizationOptions bopts;
    bopts.node_batch_axes = node_axes;
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto enode = binarize(optimized, {}, bopts);
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END

    auto ta_node = to_ta_node(enode);
    copy_batch_annotations(enode, ta_node);

    // Populate sliced_modes from the External stamps (the real pipeline does
    // this inside CacheManager; the manual empty cache below does not).
    std::array<node_t, 1> forest{ta_node};
    sequant::eval::stamp_lifetime_masks(forest);
    ta_node = forest[0];

    std::cerr << "\n=== DP-opened external loops ===\n";
    dump_annotations(ta_node);

    auto accept = opts.batch_policy.is_batchable_index();
    auto cache = cache_t::empty();
    auto aops = yield_.array_ops();
    cache.set_array_ops(&aops);
    cache.set_custom_evaluator(make_batched_custom_evaluator(
        yield_, bts, accept, make_no_scope_guard{}, never_volatile{}));
    return evaluate(ta_node, target, yield_, cache)->get<TArrayD>();
  };

  // Ground truth: plain unbatched evaluation (no custom evaluator, batch
  // stamps ignored). Also fills yield_'s random leaf cache first.
  auto const ref_expr =
      deserialize<ExprPtr>(std::wstring(expr_str.begin(), expr_str.end()));
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
  auto ref_node = to_ta_node(binarize(ref_expr));
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
  auto const ref = evaluate(ref_node, target, yield_)->get<TArrayD>();
  REQUIRE(TA::norm2(ref) > 0.0);

  auto const res = run();

  TArrayD d;
  d("i,j") = ref("i,j") - res("i,j");
  double const err = TA::norm2(d);
  std::cerr << "\n[ext-placement] |ref| = " << TA::norm2(ref)
            << "  err = " << err << "\n";

  CHECK(err < 1e-10);
}

// A value whose canonicalization carries a PHASE (an antisymmetric tensor
// written with an odd permutation of its bra/ket: EvalExpr canonicalizes it
// in place and records canon_phase() == -1) must come out of the ordered
// executor with the same sign as it does out of forest descent.
//
// The orientation convention is the scope cache's, unchanged since before the
// cell table existed: the STORE holds the CANONICAL value and every READER
// applies the node's own phase once (eval.hpp's `cache.store_and_access(node,
// apply_phase(node, rb))` against the readers' `apply_phase`, and the root
// combine's `mult_by_phase`). Storage now lives in the cell registry, so its
// productions have to convert the same way -- evaluate_impl hands back the
// ORIENTED result and the phase is an involution, so a production multiplies
// it back in (ordered_executor.hpp's `canonical`, and eval.hpp's record_leaf).
// Store the oriented value instead and every read double-applies the phase:
// every value canonicalized through an odd permutation comes back with the
// wrong sign. Real TA data is required -- the DryRun backend carries no
// numbers and cannot witness a sign.
TEST_CASE("ordered-canon-phase-once", "[eval][ordered-executor]") {
  using sequant::evaluate;
  using sequant::eval::analyze_legality;
  using sequant::eval::build_ordered_schedule;
  using sequant::eval::compute_dag_boulevard;
  using sequant::eval::evaluate_ordered_schedule;
  using sequant::eval::OrderedSchedule;
  using sequant::eval::RichSchedule;
  using TArrayD = TA::TSpArrayD;
  using node_t = sequant::FullBinaryNode<sequant::EvalExprTA>;

  auto& world = TA::get_default_world();
  rand_tensor_yield<double, TA::SparsePolicy> yield_{world, 4, 6, 12};
  yield_.set_max_tile(4);

  // g is ANTISYMMETRIC and written with its two bra indices out of canonical
  // order, so its own leaf node canonicalizes to g{a_1,a_2;...} and carries
  // canon_phase() == -1; t is written canonically and carries +1.
  sequant::io::serialization::DeserializationOptions const antisymm{
      .def_perm_symm = sequant::Symmetry::Antisymm};
  auto const t1 = sequant::deserialize<sequant::ExprPtr>(
      L"g{a_2,a_1;i_3,i_4} * t{i_3,i_4;i_1,i_2}", antisymm);
  std::vector<node_t> forest{eval_node(t1)};
  std::string const target = "a_1,a_2,i_1,i_2";

  // Non-vacuous: at least one node of the forest really does carry a phase.
  std::size_t n_phased = 0;
  auto count_phases = [&](auto&& self, node_t const& n) -> void {
    if (n->canon_phase() != 1) ++n_phased;
    if (n.leaf()) return;
    self(self, n.left());
    self(self, n.right());
  };
  for (auto const& n : forest) count_phases(count_phases, n);
  // MEASURED on this fixture: the phase sits on the PRODUCT node (binarize
  // gives a contraction the phase of its canonicalized tensor operand), which
  // is this forest's root -- so it is the root production and the root
  // combine's own mult_by_phase that have to agree about the orientation the
  // registry holds. Undo the conversion at the root production and this case
  // fails with a relative difference of exactly 2, an outright sign flip.
  REQUIRE(n_phased > 0);

  auto const ref = evaluate(forest, target, yield_)->get<TArrayD>();

  sequant::BatchPolicy policy;
  policy.is_batchable_contracted_index = [](sequant::Index const&) {
    return false;
  };
  policy.is_batchable_external_index = [](sequant::Index const&) {
    return false;
  };

  sequant::eval::dryrun::SizeRegime const regime;
  sequant::eval::dryrun::CostModel const cm{regime};
  auto const block_of = [](sequant::Index const&) -> std::size_t { return 4; };
  RichSchedule const rich = compute_dag_boulevard(forest, cm, block_of);
  auto const legality = analyze_legality(rich, forest, policy);
  OrderedSchedule const ordered =
      build_ordered_schedule(rich, legality, policy, {});

  auto cache = sequant::CacheManager<node_t>::empty();
  std::function<std::size_t(sequant::Index const&)> const target_batch =
      [](sequant::Index const&) -> std::size_t { return 4; };
  auto const got = evaluate_ordered_schedule(forest, ordered, rich, target,
                                             yield_, cache, target_batch)
                       ->get<TArrayD>();

  // Same sign, same values: a double-applied phase shows up here as a
  // relative difference of 2 (the sign flip), not as FP noise.
  TArrayD diff;
  diff(target) = got(target) - ref(target);
  double const rel = TA::norm2(diff) / TA::norm2(ref);
  INFO("relative L2 diff (ordered vs forest descent, phased value) = " << rel);
  CHECK(rel < 1e-12);
}

// The root results the ordered executor hands back must never ALIAS the
// storage it still holds. The root combine (forest_combine.hpp) accumulates
// the forest's roots IN PLACE into the first one, and a root cell that is
// persistent is also published to the cache handle's PersistentValueStore --
// so handing out the registry's own buffer leaves root0 + root1 in the store,
// and the NEXT evaluation of the same schedule (which cache-halt serves
// straight from the store) starts from the corrupted value.
//
// Only a DEFAULT-constructed layout exposes it: with a real layout the combine
// permutes each root first, which allocates a fresh buffer. Real TA data is
// required -- the DryRun backend carries no numbers and cannot witness a
// double-counted addend.
TEST_CASE("ordered-root-no-alias", "[eval][ordered]") {
  using sequant::evaluate;
  using sequant::eval::analyze_legality;
  using sequant::eval::build_ordered_schedule;
  using sequant::eval::compute_dag_boulevard;
  using sequant::eval::evaluate_ordered_schedule;
  using sequant::eval::OrderedSchedule;
  using sequant::eval::RichSchedule;
  using TArrayD = TA::TSpArrayD;
  using node_t = sequant::FullBinaryNode<sequant::EvalExprTA>;

  auto& world = TA::get_default_world();
  rand_tensor_yield<double, TA::SparsePolicy> yield_{world, 4, 6, 12};
  yield_.set_max_tile(4);

  // Two INDEPENDENT roots over the same externals {i_1,i_2}. Nothing here is
  // volatile, so each root cell is persistent -- a root has no consumer at
  // all, which is exactly the "constant term" case cache-halt must serve from
  // the store instead of recomputing.
  auto const t1 =
      sequant::deserialize<sequant::ExprPtr>(L"g{a_1;i_1} * h{i_2;a_1}");
  auto const t2 =
      sequant::deserialize<sequant::ExprPtr>(L"p{a_1;i_1} * q{i_2;a_1}");
  std::vector<node_t> forest{eval_node(t1), eval_node(t2)};
  REQUIRE(forest.front()->annot() == forest.back()->annot());
  std::string const annot = forest.front()->annot();
  // The stored value is the CANONICAL orientation; keeping both roots at
  // phase +1 lets the store's content be compared to the evaluated root
  // directly.
  REQUIRE(forest.front()->canon_phase() == 1);
  // DEFAULT layout: combine_forest_roots does not permute, so root 0's own
  // buffer IS the accumulator.
  std::string const layout{};

  // Ground truth by forest descent: the whole sum, and root 0 alone (the leaf
  // yield caches by label, so every evaluation below sees the same numbers).
  auto const ref_sum =
      TA::clone(evaluate(forest, layout, yield_)->get<TArrayD>());
  std::vector<node_t> const root0_only{forest.front()};
  auto const ref_root0 =
      TA::clone(evaluate(root0_only, layout, yield_)->get<TArrayD>());
  REQUIRE(TA::norm2(ref_root0) > 0.0);

  sequant::BatchPolicy policy;
  policy.is_batchable_contracted_index = [](sequant::Index const&) {
    return false;
  };
  policy.is_batchable_external_index = [](sequant::Index const&) {
    return false;
  };

  sequant::eval::dryrun::SizeRegime const regime;
  sequant::eval::dryrun::CostModel const cm{regime};
  auto const block_of = [](sequant::Index const&) -> std::size_t { return 4; };
  RichSchedule const rich = compute_dag_boulevard(forest, cm, block_of);
  auto const legality = analyze_legality(rich, forest, policy);
  OrderedSchedule const ordered =
      build_ordered_schedule(rich, legality, policy, {});

  // ONE cache handle across both evaluations -- the persistent value store
  // lives on it, and is what carries a value from one evaluation to the next.
  auto cache = sequant::CacheManager<node_t>::empty();
  std::function<std::size_t(sequant::Index const&)> const target_batch =
      [](sequant::Index const&) -> std::size_t { return 4; };

  auto const first =
      TA::clone(evaluate_ordered_schedule(forest, ordered, rich, layout, yield_,
                                          cache, target_batch)
                    ->get<TArrayD>());
  {
    TArrayD diff;
    diff(annot) = first(annot) - ref_sum(annot);
    double const rel = TA::norm2(diff) / TA::norm2(ref_sum);
    INFO("relative L2 diff (ordered, first evaluation vs forest descent) = "
         << rel);
    CHECK(rel < 1e-12);
  }

  // The store still holds root 0 ALONE. With an aliased root result the
  // combine's add_inplace would have made this root0 + root1.
  auto const& store = cache.persistent_values();
  REQUIRE(store.holds(forest.front()->hash_value()));
  {
    auto const stored = store.get(forest.front()->hash_value());
    REQUIRE(stored);
    TArrayD diff;
    diff(annot) = stored->get<TArrayD>()(annot) - ref_root0(annot);
    double const rel = TA::norm2(diff) / TA::norm2(ref_root0);
    INFO(
        "relative L2 diff (persistent store's root 0 after the combine vs "
        "root 0 alone) = "
        << rel);
    CHECK(rel < 1e-12);
  }

  // The second evaluation is served from the store (both roots are persistent
  // and held, so cache-halt skips their production) and must reproduce the
  // first exactly -- with an aliased root 0 it would come out as
  // root0 + root1 + root1.
  auto const second =
      TA::clone(evaluate_ordered_schedule(forest, ordered, rich, layout, yield_,
                                          cache, target_batch)
                    ->get<TArrayD>());
  {
    TArrayD diff;
    diff(annot) = second(annot) - first(annot);
    double const rel = TA::norm2(diff) / TA::norm2(first);
    INFO("relative L2 diff (second evaluation vs first, one cache handle) = "
         << rel);
    CHECK(rel < 1e-12);
  }
}
