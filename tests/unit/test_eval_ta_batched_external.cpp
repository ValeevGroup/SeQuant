#include "test_eval_ta_fixtures.hpp"

using namespace sequant::tests::eval_ta;

TEST_CASE("batched_eval_external_axis_scatter", "[eval][batched-external]") {
  // Exactness gate: a batch mode stamped EXTERNAL on a node is a
  // external index that survives FREE onto the node's result (a Hadamard /
  // batched-matmul mode present on every operand and the result, contracted at
  // no node). The runtime must SCATTER the per-block partials into DISJOINT
  // slices of a pre-sized result (write_into_slice), NOT accumulate them
  // (add_inplace, correct only for a contracted mode). Because the blocks tile
  // the mode without overlap, the scattered result equals the unbatched one
  // EXACTLY -- a memory schedule, never an approximation.
  //
  // The external here is an AUXILIARY index x_1 rather than an occupied one
  // only because canonicalize forbids a NON-auxiliary index shared among > 2
  // tensor slots (no well-defined bra/ket slot type); a high-order aux
  // hyperindex carried into the result IS supported. The runtime scatter
  // mechanism under test is mode-agnostic -- occ vs aux matters only at the
  // DP/cost level (covered by the [dryrun-occ-*] gates). Same rationale as
  // "eval_forest_over_external_occ nests intra-term aux batching".
  //
  // The scatter destination starts with an all-zero shape (see
  // make_ta_array_ops), so this also checks that the scattered tiles land in
  // it.
  using sequant::evaluate;
  using sequant::make_batched_custom_evaluator;
  using TArrayD = TA::DistArray<TA::Tensor<double>, TA::SparsePolicy>;
  using node_t = sequant::FullBinaryNode<sequant::EvalExprTA>;
  using cache_t = sequant::CacheManager<node_t>;

  auto& world = TA::get_default_world();
  // virtual single-tiled (4); aux x multi-tiled (12 in tiles of 4 -> 3 tiles)
  // so the external mode genuinely partitions into > 1 block.
  rand_tensor_yield<double, TA::SparsePolicy> yield_{world, /*nocc=*/4,
                                                     /*nvirt=*/4,
                                                     /*naux=*/12};
  yield_.set_max_tile(4);

  // A flat, fully binary tree whose result carries aux x_1 as a pure external:
  //   inner = g*h  -> {a_1, a_3, x_1}   (contracts a_2)
  //   root  = inner*p -> {a_1, a_4, x_1} (contracts a_3)
  // x_1 is present on g, h, p AND the result, contracted at no node; a_1, a_4
  // are ordinary external modes. The ResultExpr LHS pins x_1 EXTERNAL (the
  // ExprPtr binarize would instead infer it contracted); stamping x_1 External
  // on the root and the inner node exercises the scatter branch at both levels.
  auto const res_expr = sequant::deserialize<sequant::ResultExpr>(
      L"R{a_1;a_4;x_1} = ((g{a_1;a_2;x_1} * h{a_2;a_3;x_1}) * p{a_3;a_4;x_1})");
  auto node = eval_node(res_expr);  // mutable: External mode stamped below
  std::string const target = node->annot();

  auto const aux_space =
      sequant::get_default_context().index_basis_registry()->retrieve_space(
          L"x");
  auto accept_aux = [aux_space](sequant::Index const& ix) {
    return ix.space() == aux_space;
  };

  // The external mode x_1, taken from the node's free/result outer modes so
  // the Index object matches the parsed one exactly (a plain outer mode, no
  // protos).
  sequant::Index mode;
  for (auto const& ix : node->canon_indices())
    if (accept_aux(ix) && !ix.has_proto_indices()) {
      mode = ix;
      break;
    }
  REQUIRE(mode.nonnull());
  REQUIRE(sequant::index_position(node, mode).has_value());  // free on result

  // Stamp External on the root AND the inner node (the mode appears free at
  // both levels).
  node->set_node_slice_mask({{mode, sequant::BatchModeType::External}});
  node_t* inner = nullptr;
  auto find_inner = [&](auto&& self, node_t& n) -> void {
    if (n.leaf()) return;
    if (&n != &node && sequant::index_position(n, mode).has_value()) inner = &n;
    self(self, n.left());
    self(self, n.right());
  };
  find_inner(find_inner, node);
  REQUIRE(inner != nullptr);
  REQUIRE(sequant::index_position(*inner, mode).has_value());
  (*inner)->set_node_slice_mask({{mode, sequant::BatchModeType::External}});

  // Reference: plain unbatched evaluation (node_slice_mask are ignored without
  // a custom evaluator). Computed first so yield_'s random leaves are generated
  // and cached, then reused by the batched run below.
  auto const ref = evaluate(node, target, yield_)->get<TArrayD>();
  REQUIRE(TA::norm2(ref) > 0.0);  // guard: reference is nontrivially nonzero

  // Spy scope-guard: records the block count each time the evaluator FIRES
  // (picks a mode partitioning into > 1 block). A silently-ignored External
  // stamp would never fire -- this is the RED/GREEN witness that the ON path
  // actually blocked.
  std::vector<std::size_t> guard_calls;
  auto spy = [&guard_calls](std::size_t n) {
    guard_calls.push_back(n);
    return sequant::no_scope_guard{};
  };

  // Two block sizes, both strictly below x_1's extent (12): 4 -> 3 blocks,
  // 8 -> 2 blocks. Exactness must hold for every partition.
  for (std::size_t const target_batch_size : {std::size_t{4}, std::size_t{8}}) {
    guard_calls.clear();
    auto cache = cache_t::empty();
    auto aops = yield_.array_ops();
    cache.set_array_ops(&aops);
    cache.set_custom_evaluator(make_batched_custom_evaluator(
        yield_,
        [target_batch_size](sequant::Index const&) -> std::size_t {
          return target_batch_size;
        },
        accept_aux, spy, sequant::never_volatile{}));
    auto const res = evaluate(node, target, yield_, cache)->get<TArrayD>();

    // Exactness: scattering disjoint blocks reconstructs the whole result.
    TArrayD diff;
    diff("i,j,k") = ref("i,j,k") - res("i,j,k");
    REQUIRE(TA::norm2(diff) < 1e-12);

    // The ON path blocked: the root fired at least once, partitioning into
    // > 1 block (block size strictly less than x_1's extent).
    REQUIRE_FALSE(guard_calls.empty());
    for (auto const n : guard_calls) CHECK(n > 1);
    std::size_t const expected_root_blocks = target_batch_size == 4 ? 3 : 2;
    CHECK(guard_calls.front() == expected_root_blocks);
  }
}

TEST_CASE("ta_array_ops_sparse_zeros_hold_no_tiles",
          "[eval][batched-external]") {
  // A sparse flat scatter destination, and its slice to an enclosing batch, are
  // zero by shape rather than by zero-filled tiles: nothing is allocated before
  // the scatter writes.
  using TArrayD = TA::DistArray<TA::Tensor<double>, TA::SparsePolicy>;
  auto& world = TA::get_default_world();
  rand_tensor_yield<double, TA::SparsePolicy> yield_{world, /*nocc=*/4,
                                                     /*nvirt=*/4, /*naux=*/12};
  yield_.set_max_tile(4);
  auto const aops = yield_.array_ops();

  auto const count_nonzero_tiles = [](TArrayD const& arr) {
    std::size_t n = 0;
    for (std::size_t ord = 0; ord != arr.trange().tiles_range().volume(); ++ord)
      if (!arr.is_zero(ord)) ++n;
    return n;
  };

  auto const zeros =
      aops.make_zeros({sequant::Index{L"i_1"}, sequant::Index{L"x_1"}});
  auto const& full = zeros->get<TArrayD>();
  CHECK(full.trange().dim(0).extent() == 4);
  CHECK(full.trange().dim(1).extent() == 12);
  CHECK(count_nonzero_tiles(full) == 0);
  CHECK(TA::norm2(full) == 0.0);

  auto const sliced = zeros->slice_mode(1, 4, 8);
  auto const& sl = sliced->get<TArrayD>();
  CHECK(sl.trange().dim(1).elements_range().first == 4);
  CHECK(sl.trange().dim(1).elements_range().second == 8);
  CHECK(count_nonzero_tiles(sl) == 0);

  // A scatter write adds exactly the written tile to the shape.
  TArrayD src(world, full.trange());
  src.fill_random();
  world.gop.fence();
  auto const src_res =
      sequant::eval_result<sequant::ResultTensorTA<TArrayD>>(src);
  auto const blk = src_res->slice_mode(1, 4, 8);
  zeros->write_into_slice(*blk, 1, 4, 8);
  CHECK(count_nonzero_tiles(full) == 1);
  CHECK(full.shape().nnz() == 1);
  CHECK(TA::norm2(full) == Catch::Approx(TA::norm2(blk->get<TArrayD>())));
}

TEST_CASE("ta_array_ops_sparse_tot_zeros_hold_no_tiles",
          "[eval][batched-external]") {
  // The nested counterpart of ta_array_ops_sparse_zeros_hold_no_tiles. The
  // slice of the zero destination goes through slice_array_over_mode's
  // empty-inner branch, which must not allocate either.
  using ToTArray = TA::DistArray<tot_tile_t<double>, TA::SparsePolicy>;
  using ResultToT = sequant::ResultTensorOfTensorTA<ToTArray>;
  auto& world = TA::get_default_world();
  rand_tensor_yield<double, TA::SparsePolicy> yield_{world, /*nocc=*/4,
                                                     /*nvirt=*/4, /*naux=*/12};
  yield_.set_max_tile(4);
  auto const aops = yield_.array_ops<ToTArray>();

  auto const count_nonzero_tiles = [](ToTArray const& arr) {
    std::size_t n = 0;
    for (std::size_t ord = 0; ord != arr.trange().tiles_range().volume(); ++ord)
      if (!arr.is_zero(ord)) ++n;
    return n;
  };

  // outer (i_1, x_1); a_1<i_1> is an inner mode
  auto const zeros =
      aops.make_zeros({sequant::Index{L"i_1"}, sequant::Index{L"x_1"},
                       sequant::Index(L"a_1", {L"i_1"})});
  REQUIRE(zeros->is<ResultToT>());
  auto const& full = zeros->get<ToTArray>();
  CHECK(full.trange().rank() == 2);
  CHECK(full.trange().dim(1).extent() == 12);
  CHECK(count_nonzero_tiles(full) == 0);

  auto const sliced = zeros->slice_mode(1, 4, 8);
  auto const& sl = sliced->get<ToTArray>();
  CHECK(sl.trange().dim(1).elements_range().first == 4);
  CHECK(sl.trange().dim(1).elements_range().second == 8);
  CHECK(count_nonzero_tiles(sl) == 0);

  // A scatter write adds exactly the written tile to the shape.
  ToTArray src(world, full.trange());
  for (auto it = src.begin(); it != src.end(); ++it)
    if (src.is_local(it.index()))
      *it = random_tensor_of_tensor<double>(it.make_range(), TA::Range{4});
  world.gop.fence();
  auto const src_res = sequant::eval_result<ResultToT>(src);
  auto const blk = src_res->slice_mode(1, 4, 8);
  zeros->write_into_slice(*blk, 1, 4, 8);
  CHECK(count_nonzero_tiles(full) == 1);
  CHECK(full.shape().nnz() == 1);
  CHECK(sequant::detail::tot_inner_rank(full) == 1);
  CHECK(TA::norm2(full) == Catch::Approx(TA::norm2(blk->get<ToTArray>())));
}

TEST_CASE("ta_array_ops_zeros_dense", "[eval][batched-external]") {
  // A dense array has no shape to express an empty destination with, so
  // make_zeros must allocate every tile, zero-filled (flat) or with empty
  // inner tensors (nested). (A sparse array never takes this path: TiledArray
  // keeps the SparseShape threshold positive, so its zero shape holds no
  // tiles.)
  using TArrayD = TA::TArrayD;
  using ToTArray = TA::DistArray<tot_tile_t<double>, TA::DensePolicy>;
  using ResultToT = sequant::ResultTensorOfTensorTA<ToTArray>;

  auto& world = TA::get_default_world();
  rand_tensor_yield<double, TA::DensePolicy> yield_{world, /*nocc=*/4,
                                                    /*nvirt=*/4, /*naux=*/12};
  yield_.set_max_tile(4);
  auto const aops = yield_.array_ops<ToTArray>();

  auto const all_tiles_nonzero = [](auto const& arr) {
    for (std::size_t ord = 0; ord != arr.trange().tiles_range().volume(); ++ord)
      if (arr.is_zero(ord)) return false;
    return true;
  };

  SECTION("flat") {
    auto const zeros =
        aops.make_zeros({sequant::Index{L"i_1"}, sequant::Index{L"x_1"}});
    REQUIRE(zeros->is<sequant::ResultTensorTA<TArrayD>>());
    auto const& full = zeros->get<TArrayD>();
    CHECK(all_tiles_nonzero(full));
    CHECK(TA::norm2(full) == 0.0);

    TArrayD src(world, full.trange());
    src.fill_random();
    world.gop.fence();
    auto const blk =
        sequant::eval_result<sequant::ResultTensorTA<TArrayD>>(src)->slice_mode(
            1, 4, 8);
    zeros->write_into_slice(*blk, 1, 4, 8);
    CHECK(TA::norm2(full) == Catch::Approx(TA::norm2(blk->get<TArrayD>())));
  }

  SECTION("nested") {
    auto const zeros =
        aops.make_zeros({sequant::Index{L"i_1"}, sequant::Index{L"x_1"},
                         sequant::Index(L"a_1", {L"i_1"})});
    REQUIRE(zeros->is<ResultToT>());
    auto const& full = zeros->get<ToTArray>();
    CHECK(all_tiles_nonzero(full));
    CHECK(sequant::detail::tot_inner_rank(full) == 0);

    // slice_array_over_mode's empty-inner branch takes the same path
    auto const sliced = zeros->slice_mode(1, 4, 8);
    auto const& sl = sliced->get<ToTArray>();
    CHECK(all_tiles_nonzero(sl));
    CHECK(sequant::detail::tot_inner_rank(sl) == 0);

    ToTArray src(world, full.trange());
    for (auto it = src.begin(); it != src.end(); ++it)
      if (src.is_local(it.index()))
        *it = random_tensor_of_tensor<double>(it.make_range(), TA::Range{4});
    world.gop.fence();
    auto const blk = sequant::eval_result<ResultToT>(src)->slice_mode(1, 4, 8);
    zeros->write_into_slice(*blk, 1, 4, 8);
    CHECK(sequant::detail::tot_inner_rank(full) == 1);
    CHECK(TA::norm2(full) == Catch::Approx(TA::norm2(blk->get<ToTArray>())));
  }
}

TEST_CASE("batched-external-hoist-slice-on-use", "[eval][batched-external]") {
  // This test used to pin the `has_demoted_external`
  // demotion veto: in an external SCATTER over a mode e, a DESCENDANT
  // intermediate M that itself CARRIES e free onto its result (an External
  // `node_slice_mask()` stamp absent from `sliced_modes()` -- a meet-demoted
  // carrier, the cross-pair giants in the C60 story) was NEVER hoisted, so it
  // was rebuilt SLICED once per e-block (h requested nblocks times). The veto
  // is now DELETED: placement is purely seed-based, so M is HOISTED to its
  // seed home. Its residency (`sliced_modes`) is empty, so the seed home is
  // the chain root: M is built exactly ONCE (h requested once), not nblocks
  // times.
  //
  // Crucially this changes PLACEMENT only, never the RESULT: M is cached FULL
  // at the root, and when R consumes it inside the e-scatter the batched
  // Enter-stage slice-on-use slices the ancestor-scope value to the current
  // e-block (the `[eval][slice-on-use]` witness below), so the scatter still
  // reassembles the whole result exactly (assertion 1, unchanged): seed
  // placement is what the runtime uses.
  //
  // The sibling INV (order-aware, NO External stamp) was already hoisted under
  // the veto and still is: both siblings now hoist, so h_evals and v_evals are
  // both 1. The test remains a real witness -- it proves the demoted carrier
  // now hoists (was nblocks, now 1) AND that the numerical result is untouched.
  //
  // Network:
  //   R{a_1;a_4;x_1} = ((g{a_1;a_2;x_1} * h{a_2;a_3;x_1})
  //                      * (u{a_3;a_5} * v{a_5;a_4}))
  //  - M = g*h contracts a_2 -> {a_1;a_3;x_1}: CARRIES the external x_1.
  //  - INV = u*v contracts a_5 -> {a_3;a_4}: a plain sibling of M, carrying
  //    NO aux mode at all (a drop-in structural replacement for the leaf `p`
  //    the earlier version of this test used, so R's own shape is unchanged).
  //  - R = M*INV contracts a_3 -> {a_1;a_4;x_1}: the x_1 scatter trigger.
  // M is stamped External on x_1 + order_aware(true), with sliced_modes left
  // empty (default). INV is stamped order_aware(true) only -- no
  // node_slice_mask entries at all. Both are collected + hoisted to the
  // term/root cache (built once, regardless of e-blocks).
  using sequant::evaluate;
  using sequant::make_batched_custom_evaluator;
  using TArrayD = TA::TSpArrayD;
  using node_t = sequant::FullBinaryNode<sequant::EvalExprTA>;
  using cache_t = sequant::CacheManager<node_t>;

  auto& world = TA::get_default_world();
  // virtual single-tiled (4); aux x multi-tiled (12 in tiles of 4 -> 3 tiles)
  // so the external mode genuinely partitions into > 1 block.
  rand_tensor_yield<double, TA::SparsePolicy> yield_{world, /*nocc=*/4,
                                                     /*nvirt=*/4, /*naux=*/12};
  yield_.set_max_tile(4);

  auto const res_expr = sequant::deserialize<sequant::ResultExpr>(
      L"R{a_1;a_4;x_1} = ((g{a_1;a_2;x_1} * h{a_2;a_3;x_1})"
      L" * (u{a_3;a_5} * v{a_5;a_4}))");
  auto node = eval_node(res_expr);
  std::string const target = node->annot();

  auto const aux_space =
      sequant::get_default_context().index_basis_registry()->retrieve_space(
          L"x");
  auto accept_aux = [aux_space](sequant::Index const& ix) {
    return ix.space() == aux_space;
  };

  // The external mode x_1 (a plain outer aux mode, no protos), taken from the
  // node's free/result modes so the Index matches the parsed one exactly.
  sequant::Index mode;
  for (auto const& ix : node->canon_indices())
    if (accept_aux(ix) && !ix.has_proto_indices()) {
      mode = ix;
      break;
    }
  REQUIRE(mode.nonnull());
  REQUIRE(sequant::index_position(node, mode).has_value());  // free on R

  // Stamp External on the root R (the scatter trigger over x_1).
  node->set_node_slice_mask({{mode, sequant::BatchModeType::External}});

  // Locate M (carries x_1 free) and INV (the sibling that does not) -- the
  // only two non-root, non-leaf nodes in this network.
  node_t* m = nullptr;
  node_t* inv = nullptr;
  auto find_mn = [&](auto&& self, node_t& n) -> void {
    if (n.leaf()) return;
    if (&n != &node) {
      if (sequant::index_position(n, mode).has_value())
        m = &n;
      else
        inv = &n;
    }
    self(self, n.left());
    self(self, n.right());
  };
  find_mn(find_mn, node);
  REQUIRE(m != nullptr);
  REQUIRE(inv != nullptr);
  REQUIRE(m != inv);
  REQUIRE(sequant::index_position(*m, mode).has_value());  // M carries x_1
  REQUIRE_FALSE(
      sequant::index_position(*inv, mode).has_value());  // INV does not

  // M: External stamp, order-aware, EMPTY sliced_modes (default) -- a
  // meet-demoted external carrier. With the veto deleted, M is collected +
  // hoisted to its seed home (the chain root, since its residency is empty),
  // built exactly once; slice-on-use slices the hoisted-full M per e-block.
  (*m)->set_node_slice_mask({{mode, sequant::BatchModeType::External}});
  (*m)->set_batch_order_aware(true);

  // INV: order-aware, no node_slice_mask entries at all (no External stamp) ->
  // collected + hoisted to the term/root cache, built exactly once (as before
  // -- INV was never demoted).
  (*inv)->set_batch_order_aware(true);

  // Reference: plain unbatched evaluation (fills yield_'s leaf cache too).
  auto const ref = evaluate(node, target, yield_)->get<TArrayD>();
  REQUIRE(TA::norm2(ref) > 0.0);  // guard: reference is nontrivially nonzero

  // Counting leaf evaluator. We count leaf h -- the RIGHT operand of M = g*h.
  // find_leaf_carrying(., x_1) returns the LEFTMOST carrier (g), so the
  // sliceability PROBE (mode_batches on the carrier leaf) only ever requests g;
  // h is requested purely when M is actually MATERIALIZED (the g*h
  // contraction). Counting h therefore isolates M's build count from probe
  // noise: hoisted-full => M built once at full e (h requested once, the RED
  // bug); rebuilt sliced per e-block => M built once per block (h requested
  // nblocks times, each at block-e extent, the fix). Leaf v (INV's right
  // operand) is never touched by any sliceability probe (INV carries no aux
  // mode at all), so it isolates INV's build count the same way.
  for (std::size_t const target_batch_size : {std::size_t{4}, std::size_t{8}}) {
    // nblocks (e-blocks at this batch size) = 3 for batch 4, 2 for batch 8
    int h_evals = 0;
    int v_evals = 0;
    auto counting_yield = [&yield_, &h_evals,
                           &v_evals](node_t const& leaf) -> sequant::ResultPtr {
      if (leaf->is_tensor()) {
        if (leaf->as_tensor().label() == L"h") ++h_evals;
        if (leaf->as_tensor().label() == L"v") ++v_evals;
      }
      return yield_(leaf);
    };

    auto cache = cache_t::empty();
    auto aops = yield_.array_ops();
    cache.set_array_ops(&aops);
    cache.set_custom_evaluator(make_batched_custom_evaluator(
        counting_yield,
        [target_batch_size](sequant::Index const&) -> std::size_t {
          return target_batch_size;
        },
        accept_aux, sequant::make_no_scope_guard{}, sequant::never_volatile{}));
    auto const res =
        evaluate(node, target, counting_yield, cache)->get<TArrayD>();

    // 1. Slice-exactness: scattering disjoint blocks reconstructs the whole
    //    result exactly (a memory schedule, never an approximation).
    TArrayD diff;
    diff("i,j,k") = ref("i,j,k") - res("i,j,k");
    REQUIRE(TA::norm2(diff) < 1e-12);

    // 2. M (the meet-demoted external carrier) is now HOISTED to its seed
    //    home and built exactly ONCE (h requested once), regardless of
    //    nblocks -- the veto that used to force nblocks sliced rebuilds is
    //    deleted. The result is still exact (assertion 1) because slice-on-use
    //    slices the hoisted-full M to each e-block at the consumer.
    CHECK(h_evals == 1);

    // 3. INV (no External stamp) is hoisted once to the term/root cache too,
    //    regardless of nblocks -- unchanged by the veto deletion.
    CHECK(v_evals == 1);
  }
}

TEST_CASE("batched-cached-slice-on-use", "[eval][slice-on-use]") {
  // Slice-on-use RED/GREEN witness for the LOAD-BEARING fix. A shared
  // intermediate M that carries an external mode e can be cached FULL at an
  // OUTER scope (correct where it is used unbatched / peak-acceptable there)
  // and then consumed inside a NESTED e-external loop -- the C60 summand-46
  // `s.C` pattern (`I(mu~; a<i1,i2>)` stored
  // full at scope -1, served WHOLE to the external use, contracting the
  // full 8.7 GB into a 2930 GB giant). The fix: the Enter-stage slice-on-use
  // slices a value fetched from an ANCESTOR scope (access_at hops > 0) to the
  // current block for the loops the fetch crossed, so M is realized at BLOCK-e
  // extent at the consumer, never full.
  //
  // Network: R{a_1;a_4;x_1} = (g{a_1;a_2;x_1} * h{a_2;a_3;x_1}) * p{a_3;a_4}.
  //  - M = g*h contracts a_2 -> {a_1;a_3;x_1}: carries the external x_1.
  //  - R = M*p contracts a_3 -> {a_1;a_4;x_1}: consumes M under an x_1 block.
  // M is pre-stored FULL in an OUTER cache; an INNER child cache carries a
  // batch_context of one x_1 block. Evaluating R on the inner cache fetches M
  // from the outer scope (hops == 1) and MUST slice it to the block.
  using sequant::evaluate;
  using TArrayD = TA::TSpArrayD;
  using node_t = sequant::FullBinaryNode<sequant::EvalExprTA>;
  using cache_t = sequant::CacheManager<node_t>;
  using hasher_t = sequant::TreeNodeHasher<node_t>;
  using comp_t = sequant::TreeNodeEqualityComparator<node_t>;

  auto& world = TA::get_default_world();
  // virtual single-tiled (4); aux x multi-tiled (12 in tiles of 4 -> 3 tiles).
  rand_tensor_yield<double, TA::SparsePolicy> yield_{world, /*nocc=*/4,
                                                     /*nvirt=*/4, /*naux=*/12};
  yield_.set_max_tile(4);

  auto const res_expr = sequant::deserialize<sequant::ResultExpr>(
      L"R{a_1;a_4;x_1} = ((g{a_1;a_2;x_1} * h{a_2;a_3;x_1}) * p{a_3;a_4})");
  auto node = eval_node(res_expr);  // R
  std::string const target = node->annot();

  auto const aux_space =
      sequant::get_default_context().index_basis_registry()->retrieve_space(
          L"x");
  // The external mode x_1 (a plain outer aux mode, no protos).
  sequant::Index mode;
  for (auto const& ix : node->canon_indices())
    if (ix.space() == aux_space && !ix.has_proto_indices()) {
      mode = ix;
      break;
    }
  REQUIRE(mode.nonnull());
  REQUIRE(sequant::index_position(node, mode).has_value());  // free on R

  // Locate M: the unique non-root internal node that carries x_1 free on its
  // result (g*h).
  node_t* m = nullptr;
  auto find_m = [&](auto&& self, node_t& n) -> void {
    if (n.leaf()) return;
    if (&n != &node && sequant::index_position(n, mode).has_value()) m = &n;
    self(self, n.left());
    self(self, n.right());
  };
  find_m(find_m, node);
  REQUIRE(m != nullptr);
  REQUIRE(sequant::index_position(*m, mode).has_value());  // M carries x_1
  // Keep the canonical phase trivial so the pre-stored value convention is
  // unambiguous (store V, Enter-stage hit applies phase == 1).
  REQUIRE((*m)->canon_phase() == 1);

  // The x_1 block under test: the first tile, [0, 4) of the extent-12 mode
  // (tile-aligned, as slice_mode requires).
  std::size_t const blk_lo = 0, blk_hi = 4;

  // Reference (full x_1): plain unbatched evaluation. Also fills yield_'s leaf
  // cache so the later runs reuse the same random leaves.
  auto const ref = evaluate(node, target, yield_)->get<TArrayD>();
  REQUIRE(TA::norm2(ref) > 0.0);
  auto const ref_vol = ref.trange().elements_range().volume();

  // Independent BLOCK reference: R computed with every x_1-carrying leaf sliced
  // to [0, 4) (the old le_g-style leaf slicing), i.e. the correct block-x_1
  // result -- built WITHOUT the slice-on-use path.
  auto block_ref_leaf = [&](node_t const& ln) -> sequant::ResultPtr {
    sequant::ResultPtr r = yield_(ln);
    if (auto const p = sequant::index_position(ln, mode))
      return r->slice_mode(*p, blk_lo, blk_hi);
    return r;
  };
  auto const block_ref = evaluate(node, target, block_ref_leaf)->get<TArrayD>();
  auto const block_vol = block_ref.trange().elements_range().volume();
  REQUIRE(block_vol > 0);
  REQUIRE(block_vol * 3 == ref_vol);  // one of three x_1 tiles

  // --- Slice-on-use run. ---
  // OUTER cache holds M FULL; a counting leaf evaluator proves M's children
  // (g, h) are NOT re-evaluated during the run (M is served from cache, sliced
  // on use -- not rebuilt).
  int gh_evals = 0;
  auto counting_yield = [&yield_,
                         &gh_evals](node_t const& leaf) -> sequant::ResultPtr {
    if (leaf->is_tensor() && (leaf->as_tensor().label() == L"g" ||
                              leaf->as_tensor().label() == L"h"))
      ++gh_evals;
    return yield_(leaf);
  };

  auto const m_full = evaluate(*m, (*m)->annot(), counting_yield);
  REQUIRE(m_full->get<TArrayD>().trange().elements_range().volume() ==
          block_ref.trange().elements_range().volume() * 3);  // full x_1

  std::unordered_map<node_t, std::size_t, hasher_t, comp_t> outer_counts;
  outer_counts.emplace(*m, 100);  // generous life; M stays alive across reads
  auto outer = cache_t(std::move(outer_counts));
  (void)outer.store_and_access(*m, m_full);
  REQUIRE(outer.alive(*m));

  auto inner = cache_t::empty();
  inner.set_parent(&outer);
  // DAG-scope level is a placeholder. This test exercises the exact-axis
  // resolution (forest-descent style: exact_axis == axis), which
  // slice_to_use resolves via index_position(nd, *exact_axis) -- byte-
  // identical to the `.axis` / `.range`-driven path this test was written
  // against.
  inner.set_batch_context(
      {{mode,
        sequant::DagScopeLevel{1, std::wstring(mode.space().base_key())},
        {blk_lo, blk_hi},
        mode}});

  gh_evals = 0;  // count only the slice-on-use run
  auto const res =
      evaluate(node, target, counting_yield, inner)->get<TArrayD>();

  // (a) SIZE PROBE: M is realized at BLOCK-x_1 extent at the consumer, so R is
  //     block-sized -- NOT the full-x_1 giant. (RED: served whole -> full-x_1
  //     product -> res_vol == ref_vol.)
  auto const res_vol = res.trange().elements_range().volume();
  REQUIRE(res_vol == block_vol);
  REQUIRE(res_vol * 3 == ref_vol);

  // M's children were never re-evaluated: M was fetched from the outer cache
  // and sliced on use, not rebuilt from g, h.
  REQUIRE(gh_evals == 0);

  // (b) EXACTNESS: the sliced-on-use result equals the independent block
  //     reference exactly (a memory schedule, never an approximation).
  TArrayD diff;
  diff("i,j,k") = block_ref("i,j,k") - res("i,j,k");
  REQUIRE(TA::norm2(diff) < 1e-12);
}

TEST_CASE("batched_scratch_no_seed_external", "[eval][batched-external]") {
  // RED/GREEN witness for the SEED DECISION in
  // make_batched_scratch. When an EXTERNAL mode is batched, a
  // persistent, alive intermediate that is INVARIANT under the batch's
  // *contracted* member mode but CARRIES that batched External mode must NOT be
  // seeded from the real cache: its full, unsliced-external value is wrong
  // under the external slice, so it must be recomputed sliced. Before this fix
  // the per-node signature and the `!e.sig && persistent && alive` seed rule
  // knew only the contracted member mode, so such a node (contracted signature
  // == absent) was wrongly seeded FULL -- the "seed full i" defect and the root
  // of the original w8-occbatch abort. After the fix the batch's External modes
  // enter the signature: a node carrying one is non-invariant and excluded.
  using node_t = sequant::FullBinaryNode<sequant::EvalExprTA>;
  using sequant::evaluate;
  using sequant::Index;

  auto& world = TA::get_default_world();
  rand_tensor_yield<double, TA::SparsePolicy> yield_{world, /*nocc=*/4,
                                                     /*nvirt=*/4, /*naux=*/12};
  yield_.set_max_tile(4);

  // P = g*h carries the aux external x_1 (the External mode); the volatile
  // head t makes P a PERSISTENT final. P does not carry the contracted batch
  // mode played by `mu` below.
  auto const expr = sequant::deserialize<sequant::ExprPtr>(
      L"(g{a_2;i_1;x_1} * h{i_3;a_2}) * t{i_1;i_3}");
  auto n = eval_node(expr);
  REQUIRE_FALSE(n.leaf());
  auto const& P = n.left();  // g*h, carries x_1
  REQUIRE_FALSE(P.leaf());

  // the External mode carried by P (a plain aux outer mode)
  auto const aux_space =
      sequant::get_default_context().index_basis_registry()->retrieve_space(
          L"x");
  Index x_ext;
  for (auto const& ix : P->canon_indices())
    if (ix.space() == aux_space && !ix.has_proto_indices()) {
      x_ext = ix;
      break;
    }
  REQUIRE(x_ext.nonnull());
  REQUIRE(sequant::index_position(P, x_ext).has_value());  // P carries it

  // real cache: P classifies persistent (volatile head t); store a value so it
  // is ALIVE (the state that makes a node a seed candidate). Built with n's
  // node_slice_mask() still empty (no External stamp anywhere yet), so the
  // gated cache_manager()'s own internal stamp_lifetime_masks() call (which
  // stamps every caller's forest, not just build_dryrun_cache's) sees no
  // External mode and leaves every mask all-full -- the mask veto stays inert
  // and P is persistent here exactly as before that change.
  auto is_volatile_t = [](node_t const& k) {
    return k.leaf() && k->is_tensor() && k->as_tensor().label() == L"t";
  };
  auto real = sequant::cache_manager(std::vector{n}, is_volatile_t);
  REQUIRE(real.persistent(P));
  (void)real.store_and_access(P, evaluate(P, P->annot(), yield_));
  REQUIRE(real.alive(P));

  // NOW stamp x_ext External on the member root, as the optimizer would; this
  // is the batch's External-mode set that make_batched_scratch partitions
  // out below. This is a DIFFERENT mechanism from the cache's lifetime mask
  // (make_batched_scratch reads node_slice_mask() directly, never the mask),
  // and is applied only after `real` is built so it does not also feed the
  // cache's own veto -- keeping this test isolated to the seed-decision guard
  // under test.
  n->set_node_slice_mask({{x_ext, sequant::BatchModeType::External}});

  // member mode = a contracted mode P does NOT carry (mu-analog): its per-node
  // signature over P is 'absent', so on the contracted mode alone P looks
  // seedable. The External mode is what must veto the seed.
  Index const mu(L"i_9");  // absent from the tree -> null contracted signature
  REQUIRE_FALSE(sequant::index_position(P, mu).has_value());
  std::vector<std::pair<node_t const*, Index>> const members{{&n, mu}};

  sequant::TreeNodeEqualityComparator<node_t> const eq;
  auto seeds_P = [&](auto const& bs) {
    return std::any_of(bs.seeds.begin(), bs.seeds.end(),
                       [&](node_t const* s) { return eq(*s, P); });
  };

  // POST-FIX: P carries the batched External mode, so it is NOT seeded (it will
  // be recomputed under the external slice). This is the seeding discipline
  // make_batched_scratch has -- the ordered (table-driven) executor does not
  // come through here at all.
  auto const bs = sequant::detail::make_batched_scratch(members, real);
  REQUIRE_FALSE(seeds_P(bs));

  // Control: with NO External mode stamped (a purely contracted batch), the
  // SAME invariant persistent P IS seeded -- proving the exclusion is driven by
  // the External stamp and that the contracted-only path is byte-identical.
  n->set_node_slice_mask({});
  auto const bs2 = sequant::detail::make_batched_scratch(members, real);
  REQUIRE(seeds_P(bs2));
}

TEST_CASE("batched_scratch_tot_presize_scatter", "[eval][batched-external]") {
  // The ToT scatter primitives (write_into_slice -> write_array_into_mode)
  // must reassemble a full-extent destination exactly. CSV/PNO-CCSD residuals
  // carry ToT tiles, so the external-mode scatter needs a nested destination
  // as well as a flat one. The runtime sequence is driven here directly:
  // build the zero destination over the carrier's whole outer tiling, then
  // write_into_slice every disjoint block. The reassembled ToT must equal the
  // original.
  using sequant::eval_result;
  using sequant::ResultTensorOfTensorTA;
  using sequant::slice_array_over_mode;
  using ToTArray = TA::DistArray<tot_tile_t<double>, TA::SparsePolicy>;
  using ResultToT = ResultTensorOfTensorTA<ToTArray>;
  auto& world = TA::get_default_world();

  // outer mode 0 (the external mode): two size-3 tiles (tile-spanning
  // split); outer mode 1: one size-2 tile. Inner tensors are 2x2 with
  // position-dependent values so a mis-scattered block changes the norm.
  TA::TiledRange const outer_tr{{0, 3, 6}, {0, 2}};
  TA::Range const inner_r(std::array<std::size_t, 2>{2, 2});

  auto build = [&world, inner_r](TA::TiledRange const& otr,
                                 bool zero) -> ToTArray {
    auto tile_fn = [inner_r, zero](TA::Range const& orng) {
      return make_tot_tile<double>(
          orng, [&](std::size_t) { return inner_r; },
          [&](std::size_t o, std::size_t k) {
            if (zero) return 0.0;
            double base = 0.0;
            for (auto c : orng.idx(o))
              base = base * 37.0 + static_cast<double>(c + 1);
            return base * 100.0 + static_cast<double>(k + 1);
          });
    };
    ToTArray arr{world, otr};
    for (auto it = arr.begin(); it != arr.end(); ++it)
      if (arr.is_local(it.index()))
        *it = world.taskq.add(tile_fn, it.make_range());
    world.gop.fence();
    return arr;
  };

  ToTArray R = build(outer_tr, /*zero=*/false);

  // gather the disjoint, tile-aligned blocks over outer mode 0 (the external
  // mode), exactly as the scatter branch's per-block evaluate would produce.
  auto const b0 = slice_array_over_mode(R, 0, 0, 1);  // outer elements [0,3)
  auto const b1 = slice_array_over_mode(R, 0, 1, 2);  // outer elements [3,6)

  // Build the full-extent zero ToT destination directly from the carrier R's
  // trange -- exactly the empty-inner zero ToT the scatter branch's make_zeros
  // produces (every outer tile a well-formed empty-inner tensor, overwritten by
  // the scattered blocks below).
  ToTArray dest_arr{world, R.trange()};
  for (auto it = dest_arr.begin(); it != dest_arr.end(); ++it)
    if (dest_arr.is_local(it.index()))
      *it = ToTArray::value_type{it.make_range()};
  world.gop.fence();
  sequant::ResultPtr const rdest = eval_result<ResultToT>(dest_arr);
  // the destination spans the FULL external mode (mode 0: {0,3,6})
  REQUIRE(rdest->get<ToTArray>().trange().dim(0).tile_extent() == 2);
  REQUIRE(rdest->get<ToTArray>().trange().dim(1).tile_extent() == 1);

  // SCATTER every block into its disjoint slice of the pre-sized destination.
  rdest->write_into_slice(*eval_result<ResultToT>(b0), 0, 0, 3);
  rdest->write_into_slice(*eval_result<ResultToT>(b1), 0, 3, 6);

  auto const& out = rdest->get<ToTArray>();
  ToTArray diff;
  diff("i,j;a,b") = out("i,j;a,b") - R("i,j;a,b");
  REQUIRE(Catch::Approx(diff("i,j;a,b").dot(diff("i,j;a,b"))) == 0.0);
  // guard against a trivial all-zero pass (e.g. a silently no-op'd scatter).
  REQUIRE(out("i,j;a,b").dot(out("i,j;a,b")) > 0.0);
}

TEST_CASE("batched_eval_external_proto_occ_scatter",
          "[eval][batched-external]") {
  // D2.2 regression: the external-mode scatter locates and slices an external
  // occupied index that a CSV/PNO-CCk giant carries ONLY as a protoindex of a
  // composite PNO leg `a<i,j>` -- never written as a plain top-level slot on
  // the giant or on the operands that build it. The giant here is the
  // particle-particle-ladder-shaped intermediate
  //   W{a1<i,j>,a2<i,j>} = (g{m1;m2} * C{a1<i,j>;m2}) * C{a2<i,j>;m1}
  // whose operand legs carry the occ i_1,i_2 only inside the composite `a`.
  //
  // WHY NO PROTO-AWARE LOCATOR IS NEEDED: occupied indices are GUARANTEED to be
  // OUTER modes of a nested (tensor-of-tensor) DistArray -- only the CSV/PNO
  // modes (the composite's base index `a`) are inner. Canonicalization promotes
  // such a proto occ to a plain, top-level, NON-proto entry of canon_indices()
  // (annot() = "i_1,i_2;a_1i_1i_2,a_2i_1i_2"), so index_position() already
  // finds it at its physical outer mode and the scatter slices it as an
  // ordinary outer mode -- there is no "proto-only" blind spot. D1 (the DP) now
  // emits BatchModeType::External for exactly these occ; THIS scatter is what
  // realizes the DP's modeled memory bound at runtime. The test pins that
  // end-to-end: the scatter fires over the occ, the assembled result
  // reconstructs the unbatched reference EXACTLY, and the per-block footprint
  // scales ~block/extent. It also guards the 2026-07-11 defect-1 over-throw (a
  // plain PAO leg alongside composite proto legs, or a leaf that does not carry
  // the occ at all): locating the occ returns nullopt there instead of
  // asserting.
  //
  // The scatter destination starts with an all-zero shape (see
  // make_ta_array_ops), so this also checks that the scattered tiles land in
  // it.
  using sequant::evaluate;
  using sequant::Index;
  using sequant::index_position;
  using sequant::make_batched_custom_evaluator;
  using node_t = sequant::FullBinaryNode<sequant::EvalExprTA>;
  using cache_t = sequant::CacheManager<node_t>;
  using ToTArray = TA::DistArray<tot_tile_t<double>, TA::SparsePolicy>;

  auto& world = TA::get_default_world();
  // occupied multi-tiled (12 in tiles of 4 -> 3 tiles) so the proto occ mode
  // genuinely partitions into > 1 block; virtual/pao single-tiled.
  rand_tensor_yield<double, TA::SparsePolicy> yield_{world, /*nocc=*/12,
                                                     /*nvirt=*/4,
                                                     /*naux=*/8};
  yield_.set_max_tile(4);

  // g is a flat (mu~) operand carrying NO occ; the two C legs are composite
  // (a<i_1,i_2>) carrying the occ only as protos.
  auto const expr = sequant::deserialize<sequant::ExprPtr>(
      L"(g{m_1;m_2} * C{m_2;a1<i_1,i_2>}) * C{a2<i_1,i_2>;m_1}");
  auto node = eval_node(expr);  // mutable: External mode stamped below
  std::string const target = node->annot();

  auto const occ =
      sequant::get_default_context().index_basis_registry()->retrieve_space(
          L"i");
  auto accept_occ = [occ](Index const& ix) {
    return ix.space() == occ && !ix.has_proto_indices();
  };

  // The external occ mode, taken from the giant's own (promoted) outer canon
  // indices so the Index matches exactly.
  Index mode;
  for (auto const& ix : node->canon_indices())
    if (accept_occ(ix)) {
      mode = ix;
      break;
    }
  REQUIRE(mode.nonnull());
  // The proto occ IS locatable: promoted to a plain outer canon index, so
  // index_position finds it at its physical outer mode -- no blind spot.
  REQUIRE(index_position(node, mode).has_value());

  // Mixed-leaf guard (2026-07-11 defect-1): the C leaf carries a plain PAO
  // outer leg (m_1/m_2) ALONGSIDE the composite a<i,j> whose proto is the occ;
  // the g leaf carries NO occ at all. Locating the occ must resolve on the C
  // leaf and return nullopt on the g leaf -- never over-throw on a leaf that
  // does not carry the target.
  node_t const* c_leaf = nullptr;
  node_t const* g_leaf = nullptr;
  auto find_leaves = [&](auto&& self, node_t const& n) -> void {
    if (n.leaf()) {
      if (n->is_tensor()) {
        auto const lbl = n->as_tensor().label();
        if (lbl == L"C" && !c_leaf) c_leaf = &n;
        if (lbl == L"g" && !g_leaf) g_leaf = &n;
      }
      return;
    }
    self(self, n.left());
    self(self, n.right());
  };
  find_leaves(find_leaves, node);
  REQUIRE(c_leaf != nullptr);
  REQUIRE(g_leaf != nullptr);
  // mixed leaf carries the occ as an outer mode (no over-throw)...
  CHECK(index_position(*c_leaf, mode).has_value());
  // ...and find_leaf_carrying descends past the g leaf to it, returning a
  // physical outer mode (the mode fed to slice_mode / mode_batches).
  auto const carrier = sequant::find_leaf_carrying(node, mode);
  REQUIRE(carrier.has_value());
  CHECK(index_position(carrier->first, mode) == carrier->second);
  // a leaf that does not carry the occ resolves to nullopt (no over-throw).
  CHECK_FALSE(index_position(*g_leaf, mode).has_value());

  // Stamp External on every node whose result carries the occ (the root and the
  // inner g*C product), as the optimizer would for a forest-level external
  // mode.
  node->set_node_slice_mask({{mode, sequant::BatchModeType::External}});
  auto stamp_carriers = [&](auto&& self, node_t& n) -> void {
    if (n.leaf()) return;
    if (&n != &node && index_position(n, mode).has_value())
      n->set_node_slice_mask({{mode, sequant::BatchModeType::External}});
    self(self, n.left());
    self(self, n.right());
  };
  stamp_carriers(stamp_carriers, node);

  // Reference: plain unbatched evaluation (node_slice_mask ignored without a
  // custom evaluator). Computed first so yield_'s random leaves are cached and
  // reused.
  auto const ref = evaluate(node, target, yield_)->get<ToTArray>();
  REQUIRE(TA::norm2(ref) > 0.0);  // guard: reference is nontrivially nonzero

  // Footprint scaling: the occ outer mode spans the full extent (12) in the
  // unbatched result, and a single tile-aligned block spans only the block
  // width (4) -- the ~block/extent shrink the scatter buys per block.
  auto const dest_mode = *index_position(node, mode);
  REQUIRE(ref.trange().dim(dest_mode).extent() == 12);

  // Spy scope-guard: records the block count each time the scatter fires.
  std::vector<std::size_t> guard_calls;
  auto spy = [&guard_calls](std::size_t n) {
    guard_calls.push_back(n);
    return sequant::no_scope_guard{};
  };

  auto cache = cache_t::empty();
  auto aops = yield_.template array_ops<ToTArray>();
  cache.set_array_ops(&aops);
  cache.set_custom_evaluator(make_batched_custom_evaluator(
      yield_, [](Index const&) -> std::size_t { return 4; }, accept_occ, spy,
      sequant::never_volatile{}));
  auto const res = evaluate(node, target, yield_, cache)->get<ToTArray>();

  // Exactness: scattering the disjoint occ blocks reconstructs the whole ToT.
  ToTArray diff;
  diff(target) = ref(target) - res(target);
  REQUIRE(diff(target).dot(diff(target)) == Catch::Approx(0.0).margin(1e-20));
  // guard against a trivial all-zero pass.
  REQUIRE(res(target).dot(res(target)) > 0.0);

  // The scatter genuinely blocked the proto occ: fired, > 1 block (occ extent
  // 12, block width 4 -> 3 blocks at the root).
  REQUIRE_FALSE(guard_calls.empty());
  for (auto const n : guard_calls) CHECK(n > 1);
  CHECK(guard_calls.front() == 3);
}

TEST_CASE("batched_eval_external_hadamard", "[eval][batched-external]") {
  // Genuine Hadamard external: guards the "slice EVERY carrying
  // operand to the same block" contract. i_3 is shared ELEMENTWISE (not
  // summed) between BOTH operands of ONE product node and survives to the
  // result -- a true Hadamard mode, not merely a free index that happens to
  // pass through a single operand untouched (that weaker shape is what the
  // free-only fixtures, e.g. batched_scratch_no_seed_external's g/h/t tree,
  // already exercise). A bug that sliced only ONE of the two operands
  // carrying i_3 (leaving the other at its full extent) would make the
  // per-block elementwise product over i_3 mismatched/wrong for THIS
  // fixture, whereas a free-only fixture (where only one operand ever
  // carries the mode) cannot distinguish "slice every carrying operand" from
  // "slice the first carrying operand" -- there is only ever one.
  //
  // R{i_3;a_1,a_2} = (g{i_3;a_1} * h{i_3;a_2}): i_3 appears on g, on h, AND
  // on the result, contracted at no node (a pure elementwise/outer product
  // over i_3 combined with a_1 from g and a_2 from h). Because i_3 is
  // repeated across exactly 2 operand tensors (the is_valid cap is > 2), it
  // is syntactically legal directly on the occupied space; plain binarize
  // would default to CONTRACTING a repeated-and-otherwise-unlisted index, so
  // the ResultExpr LHS pins i_3 into R's free/output indices (the same
  // technique batched_eval_external_axis_scatter uses for its aux
  // external).
  using sequant::evaluate;
  using sequant::make_batched_custom_evaluator;
  using TArrayD = TA::TSpArrayD;
  using node_t = sequant::FullBinaryNode<sequant::EvalExprTA>;
  using cache_t = sequant::CacheManager<node_t>;

  auto& world = TA::get_default_world();
  // occupied multi-tiled (12 in tiles of 4 -> 3 tiles) so the Hadamard
  // external i_3 genuinely partitions into > 1 block; virtual single-tiled
  // (4).
  rand_tensor_yield<double, TA::SparsePolicy> yield_{world, /*nocc=*/12,
                                                     /*nvirt=*/4};
  yield_.set_max_tile(4);

  auto const res_expr = sequant::deserialize<sequant::ResultExpr>(
      L"R{;a_1,a_2;i_3} = (g{;a_1;i_3} * h{;a_2;i_3})");
  auto node = eval_node(res_expr);  // mutable: External mode stamped below
  std::string const target = node->annot();

  auto const occ =
      sequant::get_default_context().index_basis_registry()->retrieve_space(
          L"i");
  auto accept_occ = [occ](sequant::Index const& ix) {
    return ix.space() == occ;
  };

  // The Hadamard external i_3, taken from the node's own free/result outer
  // modes.
  sequant::Index mode;
  for (auto const& ix : node->canon_indices())
    if (accept_occ(ix)) {
      mode = ix;
      break;
    }
  REQUIRE(mode.nonnull());
  REQUIRE(sequant::index_position(node, mode).has_value());  // free on result

  // Confirm i_3 is genuinely carried by BOTH operands (not just one) -- the
  // structural property that makes this fixture "genuinely Hadamard" rather
  // than free-only.
  REQUIRE(sequant::index_position(node.left(), mode).has_value());
  REQUIRE(sequant::index_position(node.right(), mode).has_value());

  node->set_node_slice_mask({{mode, sequant::BatchModeType::External}});

  // Reference: plain unbatched evaluation (the OFF path).
  auto const ref = evaluate(node, target, yield_)->get<TArrayD>();
  REQUIRE(TA::norm2(ref) > 0.0);  // guard: reference is nontrivially nonzero

  // Firing-witness spy, as in batched_eval_external_axis_scatter: records the
  // block count each time the evaluator fires. A silently-ignored External
  // stamp (or a stamp that degraded to a no-op because only one operand got
  // sliced and the runtime bailed) would never fire this.
  std::vector<std::size_t> guard_calls;
  auto spy = [&guard_calls](std::size_t n) {
    guard_calls.push_back(n);
    return sequant::no_scope_guard{};
  };

  // Two block sizes, both strictly below i_3's extent (12): 4 -> 3 blocks,
  // 8 -> 2 blocks (the ON configuration for both).
  for (std::size_t const target_batch_size : {std::size_t{4}, std::size_t{8}}) {
    guard_calls.clear();
    auto cache = cache_t::empty();
    auto aops = yield_.array_ops();
    cache.set_array_ops(&aops);
    cache.set_custom_evaluator(make_batched_custom_evaluator(
        yield_,
        [target_batch_size](sequant::Index const&) -> std::size_t {
          return target_batch_size;
        },
        accept_occ, spy, sequant::never_volatile{}));
    auto const res = evaluate(node, target, yield_, cache)->get<TArrayD>();

    // Exactness: if the runtime had sliced only ONE of the two Hadamard
    // operands, the per-block elementwise product over i_3 would combine
    // mismatched slices and this would fail (or the two mismatched ranges
    // would not even align into a well-formed contraction).
    TArrayD diff;
    diff("i,j,k") = ref("i,j,k") - res("i,j,k");
    REQUIRE(TA::norm2(diff) < 1e-12);

    // The ON path blocked: the evaluator fired at least once, partitioning
    // into > 1 block.
    REQUIRE_FALSE(guard_calls.empty());
    for (auto const n : guard_calls) CHECK(n > 1);
    std::size_t const expected_blocks = target_batch_size == 4 ? 3 : 2;
    CHECK(guard_calls.front() == expected_blocks);
  }
}

TEST_CASE("batched_eval_external_nested_contracted",
          "[eval][batched-external]") {
  // An External mode on an
  // OUTER node nested with a Contracted mode on an INNER node -- the mirror
  // image of "eval_batched_custom_evaluator nests inner mode" (which nests
  // Contracted-over-Contracted). Here the SCATTER branch at the
  // root must re-enter and let a plain accumulate (Contracted) fire WITHIN
  // each external block, i.e. "for external-block: for contracted-block:
  // accumulate, then scatter the block". This exercises the opposite
  // composition order from batched_eval_external_axis_scatter (which nests
  // External outside External at two DIFFERENT nodes, never mixing kinds).
  //
  // inner = u{i_1;i_2;x_1,x_2} * v{i_1;i_3;x_2}: contracts i_1 and x_2
  // (both present on both operands, absent from the term's declared free
  // set), leaving free x_1, i_2, i_3.
  // root = inner * p{i_3;i_6;x_1}: contracts i_3 (present on inner and p,
  // summed), leaving free x_1, i_2, i_6. x_1 is present on BOTH root
  // operands (inner and p) and is declared part of R's free indices, so it
  // is EXTERNAL (scattered) at the root -- a plain ExprPtr would instead
  // default to contracting a doubly-repeated, undeclared index, so (as in
  // batched_eval_external_axis_scatter) the ResultExpr LHS pins x_1 free.
  using sequant::evaluate;
  using sequant::make_batched_custom_evaluator;
  using TArrayD = TA::TSpArrayD;
  using node_t = sequant::FullBinaryNode<sequant::EvalExprTA>;
  using cache_t = sequant::CacheManager<node_t>;

  auto& world = TA::get_default_world();
  // occ/virt single-tiled (4, unused as batch modes here); aux multi-tiled
  // (12 in tiles of 4 -> 3 tiles) so both x_1 (outer, External) and x_2
  // (inner, Contracted) partition into > 1 block.
  rand_tensor_yield<double, TA::SparsePolicy> yield_{world, 4, 4, 12};
  yield_.set_max_tile(4);

  auto const res_expr = sequant::deserialize<sequant::ResultExpr>(
      L"R{i_2;i_6;x_1} = ((u{i_2;i_1;x_1,x_2} * v{i_1;i_3;x_2})"
      L" * p{i_3;i_6;x_1})");
  auto node = eval_node(res_expr);  // mutable: batch modes stamped below
  std::string const target = node->annot();

  auto const aux_space =
      sequant::get_default_context().index_basis_registry()->retrieve_space(
          L"x");
  auto accept_aux = [aux_space](sequant::Index const& ix) {
    return ix.space() == aux_space;
  };

  // The External mode x_1, taken from the ROOT's own free/result outer
  // modes (a plain outer mode, no protos).
  sequant::Index x1_axis;
  for (auto const& ix : node->canon_indices())
    if (accept_aux(ix) && !ix.has_proto_indices()) {
      x1_axis = ix;
      break;
    }
  REQUIRE(x1_axis.nonnull());
  REQUIRE(sequant::index_position(node, x1_axis).has_value());  // free on R

  // The inner node (u*v) carries a DIFFERENT aux mode (x_2), contracted
  // there -- locate it the same way "nests inner mode" does, by its
  // contracted aux index rather than tree position (robust to binarize's
  // operand ordering).
  node_t* inner = nullptr;
  std::optional<sequant::Index> x2_axis;
  auto find_inner = [&](auto&& self, node_t& n) -> void {
    if (n.leaf()) return;
    if (&n != &node) {
      if (auto ax = sequant::batch_axis(n, accept_aux)) {
        inner = &n;
        x2_axis = *ax;
      }
    }
    self(self, n.left());
    self(self, n.right());
  };
  find_inner(find_inner, node);
  REQUIRE(inner != nullptr);
  REQUIRE(x2_axis.has_value());
  REQUIRE(*x2_axis != x1_axis);  // the two nested modes differ

  // Stamp the root's external mode External, the inner's contracted mode
  // Contracted.
  node->set_node_slice_mask({{x1_axis, sequant::BatchModeType::External}});
  (*inner)->set_node_slice_mask(
      {{*x2_axis, sequant::BatchModeType::Contracted}});

  // Reference: plain unbatched evaluation (the OFF path). Computed first so
  // yield_'s random leaf arrays are generated and cached, then reused below.
  auto const ref = evaluate(node, target, yield_)->get<TArrayD>();
  REQUIRE(TA::norm2(ref) > 0.0);  // guard: reference is nontrivially nonzero

  // Firing-witness spy: records the block count each time EITHER level
  // fires (the root's scatter or the inner's accumulate). Distinct target
  // sizes for x_1 (outer, External) vs x_2 (inner, Contracted) make the two
  // levels' partitions unambiguous: x_1 -> target 4 over extent 12 = 3
  // blocks; x_2 -> target 8 = 2 blocks.
  auto target_batch_size = [x1_axis](sequant::Index const& ix) -> std::size_t {
    return ix == x1_axis ? std::size_t{4} : std::size_t{8};
  };
  std::vector<std::size_t> guard_calls;
  auto spy = [&guard_calls](std::size_t n) {
    guard_calls.push_back(n);
    return sequant::no_scope_guard{};
  };

  auto cache = cache_t::empty();
  auto aops = yield_.array_ops();
  cache.set_array_ops(&aops);
  cache.set_custom_evaluator(make_batched_custom_evaluator(
      yield_, target_batch_size, accept_aux, spy, sequant::never_volatile{}));
  auto const res = evaluate(node, target, yield_, cache)->get<TArrayD>();

  // Exactness: External-scatter-outside-Contracted-accumulate composes to
  // the same result as the unbatched contraction.
  TArrayD diff;
  diff("i,j,k") = ref("i,j,k") - res("i,j,k");
  REQUIRE(TA::norm2(diff) < 1e-10);

  // Depth-2 nesting engaged: the root fires once over x_1 (3 blocks) and,
  // per x_1 block, the inner fires once over x_2 (2 blocks) -- 1 + 3 = 4
  // firings, counts a permutation of {3, 2, 2, 2}. This is the
  // RED/GREEN witness that the scatter branch's per-block reinstall lets an
  // inner Contracted mode nest exactly as it does under an outer Contracted
  // mode.
  REQUIRE(guard_calls.size() == 4);
  CHECK(std::count(guard_calls.begin(), guard_calls.end(), 3u) == 1);
  CHECK(std::count(guard_calls.begin(), guard_calls.end(), 2u) == 3);
}

TEST_CASE("batched_eval_contracted_nested_external",
          "[eval][batched-external][forest_descent]") {
  // Issue #762: the mirror of batched_eval_external_nested_contracted -- a
  // Contracted loop on the ROOT with an External scatter on an INNER node
  // that also carries the root's loop mode. Every inner partial is sliced to
  // the root's current block, so the inner scatter destination must be too.
  //
  // inner = u{i_1;i_2;x_1,x_2} * v{i_2;i_3}: contracts i_2, leaving x_1 and
  // x_2 free on the inner node. root = inner * p{i_3;i_4;x_1,x_2}: contracts
  // i_3 and both aux modes. One aux mode is the root's Contracted loop; the
  // other is the inner node's External loop, opened inside it.
  using sequant::evaluate;
  using sequant::make_batched_custom_evaluator;
  using TArrayD = TA::TSpArrayD;
  using node_t = sequant::FullBinaryNode<sequant::EvalExprTA>;
  using cache_t = sequant::CacheManager<node_t>;

  auto& world = TA::get_default_world();
  // aux multi-tiled (12 in tiles of 4 -> 3 tiles), so both loops partition
  // on tile boundaries other than 0
  rand_tensor_yield<double, TA::SparsePolicy> yield_{world, 4, 4, 12};
  yield_.set_max_tile(4);

  auto const res_expr = sequant::deserialize<sequant::ResultExpr>(
      L"R{i_1;i_4} = ((u{i_1;i_2;x_1,x_2} * v{i_2;i_3})"
      L" * p{i_3;i_4;x_1,x_2})");
  auto node = eval_node(res_expr);  // mutable: batch modes stamped below
  std::string const target = node->annot();

  auto const aux_space =
      sequant::get_default_context().index_basis_registry()->retrieve_space(
          L"x");
  auto accept_aux = [aux_space](sequant::Index const& ix) {
    return ix.space() == aux_space;
  };

  // the inner node is the non-leaf child of the root
  REQUIRE_FALSE(node.leaf());
  node_t& inner = node.left().leaf() ? node.right() : node.left();
  REQUIRE_FALSE(inner.leaf());

  // the two aux modes, as the inner node carries them; the root carries
  // neither
  sequant::container::svector<sequant::Index> aux;
  for (auto const& ix : inner->canon_indices())
    if (accept_aux(ix)) aux.push_back(ix);
  REQUIRE(aux.size() == 2);
  sequant::Index const outer_axis = aux[0];  // root, Contracted
  sequant::Index const inner_axis = aux[1];  // inner node, External
  REQUIRE_FALSE(sequant::index_position(node, outer_axis).has_value());
  REQUIRE_FALSE(sequant::index_position(node, inner_axis).has_value());

  node->set_node_slice_mask({{outer_axis, sequant::BatchModeType::Contracted}});
  inner->set_node_slice_mask({{outer_axis, sequant::BatchModeType::Contracted},
                              {inner_axis, sequant::BatchModeType::External}});

  // Reference: plain unbatched evaluation.
  auto const ref = evaluate(node, target, yield_)->get<TArrayD>();
  REQUIRE(TA::norm2(ref) > 0.0);

  // outer -> target 4 over extent 12 = 3 blocks; inner -> target 8 = 2
  // blocks, [0,8) and [8,12)
  auto target_batch_size =
      [outer_axis](sequant::Index const& ix) -> std::size_t {
    return ix == outer_axis ? std::size_t{4} : std::size_t{8};
  };
  std::vector<std::size_t> guard_calls;
  auto spy = [&guard_calls](std::size_t n) {
    guard_calls.push_back(n);
    return sequant::no_scope_guard{};
  };

  auto cache = cache_t::empty();
  auto aops = yield_.array_ops();
  cache.set_array_ops(&aops);
  cache.set_custom_evaluator(make_batched_custom_evaluator(
      yield_, target_batch_size, accept_aux, spy, sequant::never_volatile{}));
  auto const res = evaluate(node, target, yield_, cache)->get<TArrayD>();

  TArrayD diff;
  diff("i,j") = ref("i,j") - res("i,j");
  REQUIRE(TA::norm2(diff) < 1e-10);

  // the root fires once over the outer mode (3 blocks) and, per outer block,
  // the inner node scatters over its mode (2 blocks): {3, 2, 2, 2}
  REQUIRE(guard_calls.size() == 4);
  CHECK(guard_calls.front() == 3u);
  CHECK(std::count(guard_calls.begin(), guard_calls.end(), 2u) == 3);
}

TEST_CASE("batched-nested-scope-guards", "[eval]") {
  // Prove that screening relaxation composes over nested batch
  // levels. The re-entrant inner evaluator is built by threading
  // make_scope_guard (along with accept/is_volatile/persistent_only/
  // target_batch_size) INTO the nested make_batched_custom_evaluator call
  // unchanged -- so the inner level constructs its own guard via
  // make_scope_guard(inner_batches), not a no-op. The outer scope_guard is a
  // local RAII variable held for the outer level's entire batch loop, which
  // includes the per-batch evaluate() calls that trigger the inner re-entry;
  // so when the inner guard is constructed, the outer guard is STILL ALIVE.
  // A backend guard that relaxes block-sparse screening scaled by its own
  // level's batch count therefore composes MULTIPLICATIVELY across nesting:
  // net relaxation = outer_batches * inner_batches, matching "a contribution
  // significant over the full product of batch modes must not be screened
  // away in an individual (outer-cell, inner-cell) cell."
  //
  // This test does not exercise REAL block-sparse screening -- every tile of
  // its yields is nonzero, so there is nothing to screen. Instead it proves the
  // STRUCTURAL composition: a custom ScopeGuardFactory whose RAII guard
  // records, on construction/destruction, the batch count it was built with
  // against a shared "currently alive" stack. If both guards are ever alive at
  // once with the stack holding {outer_n, inner_n}, that is exactly the
  // multiplicative-composition invariant a real backend guard would exploit.
  // Numeric validation against an actual screening threshold is left to an
  // end-to-end MPQC run, where a real TiledArray-SparseShape-backed guard
  // exists to relax.
  using sequant::evaluate;
  using sequant::make_batched_custom_evaluator;
  using TArrayD = TA::TSpArrayD;
  using node_t = sequant::FullBinaryNode<sequant::EvalExprTA>;
  using cache_t = sequant::CacheManager<node_t>;

  auto& world = TA::get_default_world();
  // Same depth-2 nesting shape as the "nests inner mode" test above: occ
  // single-tiled (4), aux multi-tiled (12 in tiles of 4 -> 3 tiles) so both
  // batch modes realize 3 batches each.
  rand_tensor_yield<double, TA::SparsePolicy> yield_{world, 4, 4, 12};
  yield_.set_max_tile(4);

  auto const expr = sequant::deserialize<sequant::ExprPtr>(
      L"(((u{i_1;i_2;x_1,x_2} * v{i_3;i_1;x_2}) * w{i_2;i_5;x_1})"
      L" * p{i_6;i_3;x_1})");
  std::string const target = "i_5,i_6";
  auto node = eval_node(expr);

  auto const aux_space =
      sequant::get_default_context().index_basis_registry()->retrieve_space(
          L"x");
  auto accept_aux = [aux_space](sequant::Index const& ix) {
    return ix.space() == aux_space;
  };

  auto const root_axis = sequant::batch_axis(node, accept_aux);
  REQUIRE(root_axis.has_value());
  node->set_node_slice_mask({{*root_axis, sequant::BatchModeType::Contracted}});

  node_t* inner = nullptr;
  std::optional<sequant::Index> inner_axis;
  auto find_inner = [&](auto&& self, node_t& n) -> void {
    if (n.leaf()) return;
    if (&n != &node) {
      if (auto ax = sequant::batch_axis(n, accept_aux)) {
        inner = &n;
        inner_axis = *ax;
      }
    }
    self(self, n.left());
    self(self, n.right());
  };
  find_inner(find_inner, node);
  REQUIRE(inner != nullptr);
  REQUIRE(inner_axis.has_value());
  REQUIRE(*inner_axis != *root_axis);
  (*inner)->set_node_slice_mask(
      {{*inner_axis, sequant::BatchModeType::Contracted}});

  auto const ref = evaluate(node, target, yield_)->get<TArrayD>();

  // Shared recording state: a stack of currently-alive guards' batch counts
  // (pushed on construction, popped on destruction). Whenever a guard is
  // constructed while another is already alive (stack depth becomes 2), the
  // product of the two counts is recorded -- that product is the net
  // relaxation factor a composing backend guard would apply at that instant.
  struct GuardState {
    std::vector<std::size_t> live;
    std::size_t max_depth = 0;
    std::vector<std::size_t> products_at_depth2;
  } state;

  struct TrackingGuard {
    GuardState* st;
    TrackingGuard(GuardState* s, std::size_t n) : st(s) {
      st->live.push_back(n);
      st->max_depth = std::max(st->max_depth, st->live.size());
      if (st->live.size() == 2)
        st->products_at_depth2.push_back(st->live[0] * st->live[1]);
    }
    TrackingGuard(TrackingGuard const&) = delete;
    TrackingGuard& operator=(TrackingGuard const&) = delete;
    ~TrackingGuard() { st->live.pop_back(); }
  };
  auto make_tracking_guard = [&state](std::size_t n) {
    return TrackingGuard(&state, n);
  };

  auto cache = cache_t::empty();
  auto aops = yield_.array_ops();
  cache.set_array_ops(&aops);
  cache.set_custom_evaluator(make_batched_custom_evaluator(
      yield_,
      [](sequant::Index const&) -> std::size_t { return std::size_t{4}; },
      accept_aux, make_tracking_guard, sequant::never_volatile{}));
  auto const res = evaluate(node, target, yield_, cache)->get<TArrayD>();

  // Exactness: nesting the guard-instrumented evaluator changes no numerics.
  TArrayD diff;
  diff("i,j") = ref("i,j") - res("i,j");
  auto const err = TA::norm2(diff);
  REQUIRE(err < 1e-10);

  // The real deliverable: both guards were alive simultaneously (depth 2
  // reached) at least once, and every such simultaneity records the product
  // of outer_batches (3, over x_1) and inner_batches (3, over x_2) -- 9, the
  // net relaxation a composing backend would apply. This happens once per
  // outer batch (3 outer batches), matching the "nests inner mode" test's
  // guard_calls count of 4 total firings (1 outer + 3 inner).
  REQUIRE(state.max_depth == 2);
  REQUIRE(state.products_at_depth2.size() == 3);
  for (auto const p : state.products_at_depth2) CHECK(p == 9);
  // All guards are popped by the end: no leaked / mismatched RAII scope.
  CHECK(state.live.empty());
}

TEST_CASE("make_evaluator BatchPolicy adapter", "[eval]") {
  // Strong equivalence: make_evaluator(policy, yielder) must produce the same
  // numerical result as a hand-built make_batched_custom_evaluator with the
  // matching Tensor->EvalNode volatile lift, on the same expression.
  using sequant::evaluate;
  using sequant::make_batched_custom_evaluator;
  using sequant::make_evaluator;
  using TArrayD = TA::TSpArrayD;
  using node_t = sequant::FullBinaryNode<sequant::EvalExprTA>;
  using cache_t = sequant::CacheManager<node_t>;

  auto& world = TA::get_default_world();
  // Multi-tile unoccupied space so batching actually engages (3 tiles of <=4).
  rand_tensor_yield<double, TA::SparsePolicy> yield_{world, 4, 12};
  yield_.set_max_tile(4);

  // Contracts a1,a2 (unoccupied) -> batch mode is an unoccupied index (3
  // tiles). The subtree contains a "t" leaf, which the policy marks volatile.
  auto const expr = sequant::deserialize<sequant::ExprPtr>(
      L"g_{i1,i2}^{a1,a2} * t_{a1,a2}^{i3,i4}",
      {.def_braket_symm = sequant::Hermiticity::NonHermitian});
  std::string const target = "i_1,i_2,i_3,i_4";
  auto const node = eval_node(expr);

  // Reference (no batching).
  auto const ref = evaluate(node, target, yield_)->get<TArrayD>();

  // Build a BatchPolicy: batch over any index at size 4; "t" is volatile;
  // persistent_only=true turns the volatility gate ON (default is across the
  // board).
  sequant::BatchPolicy policy{
      .is_batchable_contracted_index =
          [](sequant::Index const&) { return true; },
      .batch_target_size = [](sequant::Index const&) -> std::size_t {
        return 4;
      },
      .is_volatile_leaf =
          [](sequant::Tensor const& t) { return t.label() == L"t"; },
      .persistent_only = true};

  // (1) make_evaluator with volatile "t" gate (persistent_only=true): batching
  // should be DECLINED (same as the persistence gate test above for
  // make_batched_custom_evaluator).
  {
    bool batched = false;
    auto spy = [&batched](std::size_t) {
      batched = true;
      return sequant::no_scope_guard{};
    };
    auto cache = cache_t::empty();
    cache.set_custom_evaluator(make_evaluator(policy, yield_, spy));
    auto const res = evaluate(node, target, yield_, cache)->get<TArrayD>();
    REQUIRE_FALSE(batched);  // volatile subtree -> not batched
    REQUIRE(equal_tarrays<Loose>(res, ref));
  }

  // (2) Strong equivalence: make_evaluator (no volatile gate: empty predicate)
  // == hand-built make_batched_custom_evaluator, same result as reference.
  {
    sequant::BatchPolicy policy_nv{
        .is_batchable_contracted_index = policy.is_batchable_contracted_index,
        .batch_target_size = policy.batch_target_size,
        .is_volatile_leaf = {}  // no volatile gate
    };

    // make_evaluator path
    auto cache_me = cache_t::empty();
    cache_me.set_custom_evaluator(make_evaluator(policy_nv, yield_));
    auto const res_me =
        evaluate(node, target, yield_, cache_me)->get<TArrayD>();

    // hand-built path
    auto hand_nv = [](node_t const&) -> bool { return false; };
    auto cache_hb = cache_t::empty();
    cache_hb.set_custom_evaluator(make_batched_custom_evaluator(
        yield_, policy_nv.batch_target_size, policy_nv.is_batchable_index(),
        sequant::make_no_scope_guard{}, hand_nv));
    auto const res_hb =
        evaluate(node, target, yield_, cache_hb)->get<TArrayD>();

    // Both must match the reference (batched summation -> loose FP tolerance).
    REQUIRE(equal_tarrays<Loose>(res_me, ref));
    REQUIRE(equal_tarrays<Loose>(res_hb, ref));
    // And agree with each other to the same loose tolerance: both come from
    // the same batched-summation algorithm, whose accumulation order is
    // thread-non-deterministic, so the two independent evaluations agree only
    // up to FP noise (a Tight/exact compare here flakes by a few ULPs).
    REQUIRE(equal_tarrays<Loose>(res_me, res_hb));
  }
}
