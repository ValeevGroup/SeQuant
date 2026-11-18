#include "test_eval_ta_fixtures.hpp"

using namespace sequant::tests::eval_ta;

TEST_CASE("eval_custom_evaluator", "[eval]") {
  using sequant::evaluate;
  using sequant::ResultPtr;
  using TArrayD = TA::TSpArrayD;
  using node_t = sequant::FullBinaryNode<sequant::EvalExprTA>;
  using cache_t = sequant::CacheManager<node_t>;

  auto& world = TA::get_default_world();
  const size_t nocc = 2, nvirt = 20;
  auto yield_ = rand_tensor_yield<double, TA::SparsePolicy>{world, nocc, nvirt};

  // a multi-product expression: several non-leaf nodes in the eval tree.
  auto const expr = sequant::deserialize<sequant::ExprPtr>(
      L"-1/4 * g_{i3,i4}^{a3,a4} * t_{a2,a4}^{i1,i2} * t_{a1,a3}^{i3,i4}",
      {.def_perm_symm = sequant::Symmetry::Antisymm});
  std::string const target = "a_1,a_2,i_1,i_2";
  auto const node = eval_node(expr);

  // standard-scheme reference (no custom evaluator)
  auto const ref = evaluate(node, target, yield_)->get<TArrayD>();

  SECTION("declining evaluator defers to the standard scheme") {
    // A custom evaluator that always returns null is consulted at every
    // non-leaf node and the standard scheme produces the result.
    int consulted = 0;
    auto cache = cache_t::empty();
    cache.set_custom_evaluator(
        [&consulted](node_t const&, cache_t&) -> ResultPtr {
          ++consulted;
          return nullptr;
        });
    auto const res = evaluate(node, target, yield_, cache)->get<TArrayD>();
    REQUIRE(equal_tarrays(res, ref));
    REQUIRE(consulted > 1);  // multiple non-leaf nodes, all declined
  }

  SECTION("intercepting evaluator takes over a subtree") {
    // A custom evaluator that takes over the first (root) node it is consulted
    // on -- here by re-evaluating that subtree via the standard scheme on a
    // scratch cache. The result must still match, and the non-null return must
    // short-circuit the recursion, so the evaluator fires exactly once.
    int consulted = 0;
    auto cache = cache_t::empty();
    cache.set_custom_evaluator(
        [&consulted, &yield_](node_t const& n, cache_t&) -> ResultPtr {
          ++consulted;
          auto scratch = cache_t::empty();
          return evaluate(n, yield_, scratch);
        });
    auto const res = evaluate(node, target, yield_, cache)->get<TArrayD>();
    REQUIRE(equal_tarrays(res, ref));
    REQUIRE(consulted == 1);
  }
}

TEST_CASE("eval_batch_axis", "[eval]") {
  using sequant::batch_axis;
  using sequant::contracted_indices;

  auto node_of = [](std::wstring_view xpr) {
    return eval_node(sequant::deserialize<sequant::ExprPtr>(xpr));
  };

  SECTION("single contracted index") {
    // R_{a1}^{i1,i3} * f_{i3}^{i2} sums over i3.
    auto const node = node_of(L"R_{a1}^{i1,i3} * f_{i3}^{i2}");
    auto const c = contracted_indices(node);
    REQUIRE(c.size() == 1);
    auto const mode = batch_axis(node);
    REQUIRE(mode.has_value());
    REQUIRE(mode.value() == c.front());
  }

  SECTION("two contracted indices") {
    // g_{i1,i2}^{a1,a2} * t_{a1,a2}^{i3,i4} sums over a1,a2.
    auto const node = node_of(L"g_{i1,i2}^{a1,a2} * t_{a1,a2}^{i3,i4}");
    auto const c = contracted_indices(node);
    REQUIRE(c.size() == 2);
    auto const mode = batch_axis(node);
    REQUIRE(mode.has_value());
    REQUIRE(ranges::contains(c, mode.value()));
  }

  SECTION("leaf has no contracted index") {
    auto const node = node_of(L"f_{i1}^{a1}");
    REQUIRE(contracted_indices(node).empty());
    REQUIRE_FALSE(batch_axis(node).has_value());
  }

  SECTION("a sum is not a contraction") {
    auto const node = node_of(L"f_{i1}^{a1} + t_{a1}^{i1}");
    REQUIRE(contracted_indices(node).empty());
    REQUIRE_FALSE(batch_axis(node).has_value());
  }

  SECTION("predicate scopes the batch mode") {
    // contracts a1,a2 (unoccupied)
    auto const node = node_of(L"g_{i1,i2}^{a1,a2} * t_{a1,a2}^{i3,i4}");
    auto const c = contracted_indices(node);
    REQUIRE(c.size() == 2);

    // accept exactly one contracted index -> batch_axis returns it
    auto const only_first = [&c](sequant::Index const& ix) {
      return ix == c.front();
    };
    REQUIRE(batch_axis(node, only_first) == c.front());

    // accept none -> nullopt
    REQUIRE_FALSE(batch_axis(node, [](sequant::Index const&) {
                    return false;
                  }).has_value());

    // scope batching to a specific IndexSpace
    auto const unocc = c.front().space();
    auto const in_unocc = [&unocc](sequant::Index const& ix) {
      return ix.space() == unocc;
    };
    REQUIRE(batch_axis(node, in_unocc).has_value());

    // a node whose only contracted index is in a different (occupied) space
    auto const node_occ = node_of(L"R_{a1}^{i1,i3} * f_{i3}^{i2}");  // sums i3
    REQUIRE_FALSE(batch_axis(node_occ, in_unocc).has_value());
  }
}

TEST_CASE("eval_slice_array_over_mode", "[eval]") {
  using sequant::slice_array_over_mode;
  auto& world = TA::get_default_world();

  // arr: a(2 tiles), b(3 tiles), c(1 tile); bb shares b's TiledRange1.
  TA::TSpArrayD arr(world, TA::TiledRange{{0, 2, 4}, {0, 3, 6, 9}, {0, 5}});
  TA::TSpArrayD bb(world, TA::TiledRange{{0, 3, 6, 9}, {0, 7}});
  arr.fill_random();
  bb.fill_random();
  world.gop.fence();

  SECTION("trange of the sliced mode") {
    auto const s = slice_array_over_mode(arr, 1, 1, 3);  // b-tiles [1,3)
    REQUIRE(s.trange().dim(1).tile_extent() == 2);
    REQUIRE(s.trange().dim(0).tile_extent() ==
            arr.trange().dim(0).tile_extent());
    REQUIRE(s.trange().dim(2).tile_extent() ==
            arr.trange().dim(2).tile_extent());
  }

  SECTION("blocked-contraction-sliced-mode") {
    // full contraction over b
    TA::TSpArrayD full;
    full("a,c,d") = arr("a,b,c") * bb("b,d");

    // split b's 3 tiles into [0,1) and [1,3), contract each, sum
    auto const a0 = slice_array_over_mode(arr, 1, 0, 1);
    auto const a1 = slice_array_over_mode(arr, 1, 1, 3);
    auto const b0 = slice_array_over_mode(bb, 0, 0, 1);
    auto const b1 = slice_array_over_mode(bb, 0, 1, 3);
    TA::TSpArrayD p0, p1, summed;
    p0("a,c,d") = a0("a,b,c") * b0("b,d");
    p1("a,c,d") = a1("a,b,c") * b1("b,d");
    summed("a,c,d") = p0("a,c,d") + p1("a,c,d");

    REQUIRE(equal_tarrays(summed, full));
  }

  SECTION("Result::slice_mode takes tile-aligned element bounds") {
    // mode 1 (b) tiles {0,3,6,9}; element range [3,9) (tile-aligned, as
    // mode_batches produces) corresponds to tiles [1,3).
    sequant::ResultPtr const r =
        sequant::eval_result<sequant::ResultTensorTA<TA::TSpArrayD>>(arr);
    auto const via_result = r->slice_mode(1, 3, 9)->get<TA::TSpArrayD>();
    auto const direct = slice_array_over_mode(arr, 1, 1, 3);
    REQUIRE(equal_tarrays(via_result, direct));
  }

  SECTION("mode-batches-of-trange1") {
    using batches_t = sequant::container::svector<std::pair<size_t, size_t>>;
    // The TA realization of axis_batches: partition a TiledRange1 into tile-
    // grouped element ranges. mode 1 (b) has 3 tiles of 3 elements each (extent
    // 9). target_batch_size is an UPPER BOUND: each batch is the largest whole-
    // tile group whose total size does not exceed the target, with a floor of
    // one tile (so a target below the tile size still yields one tile per
    // batch). A batch must never exceed the target except via that one-tile
    // floor.
    auto const tr1 = arr.trange().dim(1);
    // (extra parens: compare as a single bool so Catch2 needn't stringify
    // pairs) target >= extent -> a single batch (whole mode).
    REQUIRE((sequant::mode_batches_of_trange1(tr1, 100) == batches_t{{0, 9}}));
    // target == 2 tiles -> two-tile batches (6 <= 6).
    REQUIRE((sequant::mode_batches_of_trange1(tr1, 6) ==
             batches_t{{0, 6}, {6, 9}}));
    // target just above the tile size (but below 2 tiles) -> ONE tile per
    // batch: a 2-tile batch (6) would exceed the target. Regression: the old
    // `acc >= target` rule rounded UP to a 2-tile batch, so any target a hair
    // above the tile size doubled the realized batch -- the aux_target_size
    // 236->243 crash (236-wide K tiles, 243 target -> 472-wide batch).
    REQUIRE((sequant::mode_batches_of_trange1(tr1, 5) ==
             batches_t{{0, 3}, {3, 6}, {6, 9}}));
    REQUIRE((sequant::mode_batches_of_trange1(tr1, 4) ==
             batches_t{{0, 3}, {3, 6}, {6, 9}}));
    // target == tile size -> one tile per batch.
    REQUIRE((sequant::mode_batches_of_trange1(tr1, 3) ==
             batches_t{{0, 3}, {3, 6}, {6, 9}}));
    // target below the tile size -> still one tile per batch (the floor).
    REQUIRE((sequant::mode_batches_of_trange1(tr1, 1) ==
             batches_t{{0, 3}, {3, 6}, {6, 9}}));
  }
}

// Numeric-invariance for Result::write_into_slice on the TA backend: the union
// of disjoint, tile-aligned blocks scattered into a pre-sized destination must
// reconstruct the whole array EXACTLY. write_into_slice() is the inverse of
// slice_array_over_mode() (which GATHERS a block out): here we gather disjoint
// blocks of a whole array R and scatter each back into a fresh, pre-sized
// destination, then require the reassembled destination == R elementwise.
TEST_CASE("eval_write_into_slice", "[eval]") {
  using sequant::eval_result;
  using sequant::ResultPtr;
  using sequant::ResultTensorOfTensorTA;
  using sequant::ResultTensorTA;
  using sequant::slice_array_over_mode;
  auto& world = TA::get_default_world();

  SECTION("flat: disjoint tile-aligned blocks reassemble exactly") {
    // mode 0 has 3 multi-element tiles (size 3 each; occ_tile_size>1 analog),
    // mode 1 a single tile. Split mode 0 into element ranges [0,3) and [3,9)
    // (tiles [0,1) and [1,3)): disjoint, tile-aligned, and tile-SPANNING.
    TA::TSpArrayD R(world, TA::TiledRange{{0, 3, 6, 9}, {0, 5}});
    R.fill_random();
    world.gop.fence();

    auto const b0 = slice_array_over_mode(R, 0, 0, 1);  // elements [0,3)
    auto const b1 = slice_array_over_mode(R, 0, 1, 3);  // elements [3,9)

    // destination pre-sized to R's full shape, then filled block by block.
    TA::TSpArrayD dest(world, R.trange());
    dest.fill_local(0.0);
    world.gop.fence();

    auto rdest = eval_result<ResultTensorTA<TA::TSpArrayD>>(dest);
    rdest->write_into_slice(*eval_result<ResultTensorTA<TA::TSpArrayD>>(b0), 0,
                            0, 3);
    rdest->write_into_slice(*eval_result<ResultTensorTA<TA::TSpArrayD>>(b1), 0,
                            3, 9);

    REQUIRE(equal_tarrays<Tight>(rdest->get<TA::TSpArrayD>(), R));
  }

  SECTION("flat-nonzero-lobound") {
    // mode 0 has element lobound 2 (a frozen-core-like offset) and two size-3
    // tiles; split into [2,5) and [5,8). The block bounds honor the lobound.
    TA::TSpArrayD R(world, TA::TiledRange{{2, 5, 8}, {0, 4}});
    R.fill_random();
    world.gop.fence();

    auto const b0 = slice_array_over_mode(R, 0, 0, 1);  // elements [2,5)
    auto const b1 = slice_array_over_mode(R, 0, 1, 2);  // elements [5,8)
    // the gathered blocks keep the source element lobound
    REQUIRE(b0.trange().dim(0).elements_range().first == 2);

    TA::TSpArrayD dest(world, R.trange());
    dest.fill_local(0.0);
    world.gop.fence();

    auto rdest = eval_result<ResultTensorTA<TA::TSpArrayD>>(dest);
    rdest->write_into_slice(*eval_result<ResultTensorTA<TA::TSpArrayD>>(b0), 0,
                            2, 5);
    rdest->write_into_slice(*eval_result<ResultTensorTA<TA::TSpArrayD>>(b1), 0,
                            5, 8);

    auto const& out = rdest->get<TA::TSpArrayD>();
    REQUIRE(out.trange().dim(0).elements_range().first == 2);  // lobound kept
    REQUIRE(equal_tarrays<Tight>(out, R));
  }

  SECTION("tot: disjoint tile-aligned blocks reassemble exactly") {
    using ToTArray = TA::DistArray<tot_tile_t<double>, TA::SparsePolicy>;
    using ResultToT = ResultTensorOfTensorTA<ToTArray>;

    // outer mode 0: two size-3 tiles (multi-element, tile-spanning split);
    // outer mode 1: one size-2 tile. Inner tensors are 2x2. Inner values are
    // position-dependent (derived from the outer coordinate) so a mis-scattered
    // block -- wrong offset -- would change the reassembled norm.
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

    auto const b0 = slice_array_over_mode(R, 0, 0, 1);  // outer elements [0,3)
    auto const b1 = slice_array_over_mode(R, 0, 1, 2);  // outer elements [3,6)

    // destination pre-sized to R's full shape with well-formed (zero) inners.
    ToTArray dest = build(outer_tr, /*zero=*/true);

    auto rdest = eval_result<ResultToT>(dest);
    rdest->write_into_slice(*eval_result<ResultToT>(b0), 0, 0, 3);
    rdest->write_into_slice(*eval_result<ResultToT>(b1), 0, 3, 6);

    // ToT elementwise invariance: norm of the difference (dot of the diff with
    // itself) must vanish. TA::norm2 does not support ToT tiles, so use dot.
    auto const& out = rdest->get<ToTArray>();
    ToTArray diff;
    diff("i,j;a,b") = out("i,j;a,b") - R("i,j;a,b");
    REQUIRE(Catch::Approx(diff("i,j;a,b").dot(diff("i,j;a,b"))) == 0.0);
    // and the reassembled result is nonzero (guards against a trivial all-zero
    // pass, e.g. if both blocks silently no-op'd).
    REQUIRE(out("i,j;a,b").dot(out("i,j;a,b")) > 0.0);
  }
}

// The on-disk form of a TiledArray result (Result::on_disk): the value lives
// in a TA::DiskArray, comes back whole (in_memory) or by tile block
// (slice_mode), and can be assembled by blocks (write_into_slice).
TEST_CASE("eval_on_disk_ta", "[eval][disk]") {
  using sequant::eval_result;
  using sequant::ResultPtr;
  using sequant::ResultTensorTA;
  using sequant::slice_array_over_mode;
  auto& world = TA::get_default_world();
  auto const disk_reads = [&world] {
    world.gop.fence();
    std::size_t n = TA::disk_io_stats().tiles_read;
    world.gop.sum(n);
    return n;
  };

  SECTION("flat") {
    using Res = ResultTensorTA<TA::TSpArrayD>;
    TA::TSpArrayD R(world, TA::TiledRange{{0, 3, 6, 9}, {0, 5}});
    R.fill_random();
    world.gop.fence();
    auto const r = eval_result<Res>(R);

    ResultPtr const d = r->on_disk();
    REQUIRE(d);
    CHECK(d->is_on_disk());
    CHECK(!r->is_on_disk());
    CHECK(d->size_in_bytes() >= 9 * 5 * sizeof(double));

    // whole
    ResultPtr const m = d->in_memory();
    CHECK(!m->is_on_disk());
    CHECK(equal_tarrays<Tight>(m->get<TA::TSpArrayD>(), R));

    // by block: equal to the slice of the in-memory value, and only the
    // block's tiles are read
    TA::reset_disk_io_stats();
    world.gop.fence();
    ResultPtr const sl = d->slice_mode(0, 3, 9);
    CHECK(!sl->is_on_disk());
    CHECK(equal_tarrays<Tight>(sl->get<TA::TSpArrayD>(),
                               slice_array_over_mode(R, 0, 1, 3)));
    CHECK(disk_reads() == 2);

    // a pending view goes to disk as the value it stands for
    ResultPtr const dv = r->slice_mode(0, 0, 3)->on_disk();
    REQUIRE(dv);
    CHECK(equal_tarrays<Tight>(dv->in_memory()->get<TA::TSpArrayD>(),
                               slice_array_over_mode(R, 0, 0, 1)));

    // the on-disk form is for storage, not arithmetic
    CHECK_THROWS(d->mult_by_phase(-1));
    CHECK_THROWS(d->add_inplace(*r));

    // assembled on disk block by block
    TA::TSpArrayD zero(world, R.trange());
    zero.fill_local(0.0);
    world.gop.fence();
    ResultPtr const dest = eval_result<Res>(zero)->on_disk();
    REQUIRE(dest);
    dest->write_into_slice(*eval_result<Res>(slice_array_over_mode(R, 0, 0, 1)),
                           0, 0, 3);
    dest->write_into_slice(*eval_result<Res>(slice_array_over_mode(R, 0, 1, 3)),
                           0, 3, 9);
    world.gop.fence();
    CHECK(dest->is_on_disk());
    CHECK(equal_tarrays<Tight>(dest->in_memory()->get<TA::TSpArrayD>(), R));
  }

  SECTION("nested") {
    using ToTArray = TA::DistArray<tot_tile_t<double>, TA::SparsePolicy>;
    using ResultToT = sequant::ResultTensorOfTensorTA<ToTArray>;
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
    auto const same = [](ToTArray const& a, ToTArray const& b) {
      ToTArray diff;
      diff("i,j;a,b") = a("i,j;a,b") - b("i,j;a,b");
      return diff("i,j;a,b").dot(diff("i,j;a,b")).get() == 0.0 &&
             a("i,j;a,b").dot(a("i,j;a,b")).get() > 0.0;
    };
    ToTArray const R = build(outer_tr, /*zero=*/false);
    auto const r = eval_result<ResultToT>(R);

    ResultPtr const d = r->on_disk();
    REQUIRE(d);
    CHECK(d->is_on_disk());
    CHECK(same(d->in_memory()->get<ToTArray>(), R));

    TA::reset_disk_io_stats();
    world.gop.fence();
    ResultPtr const sl = d->slice_mode(0, 3, 6);
    CHECK(!sl->is_on_disk());
    CHECK(same(sl->get<ToTArray>(), slice_array_over_mode(R, 0, 1, 2)));
    CHECK(disk_reads() == 1);

    ResultPtr const dest =
        eval_result<ResultToT>(build(outer_tr, /*zero=*/true))->on_disk();
    REQUIRE(dest);
    dest->write_into_slice(
        *eval_result<ResultToT>(slice_array_over_mode(R, 0, 0, 1)), 0, 0, 3);
    dest->write_into_slice(
        *eval_result<ResultToT>(slice_array_over_mode(R, 0, 1, 2)), 0, 3, 6);
    world.gop.fence();
    CHECK(same(dest->in_memory()->get<ToTArray>(), R));
  }
}

TEST_CASE("eval_batched_custom_evaluator", "[eval]") {
  using sequant::evaluate;
  using sequant::make_batched_custom_evaluator;
  using TArrayD = TA::TSpArrayD;
  using node_t = sequant::FullBinaryNode<sequant::EvalExprTA>;
  using cache_t = sequant::CacheManager<node_t>;

  auto& world = TA::get_default_world();
  // multi-tile arrays so batching over a contracted index actually engages:
  // unoccupied extent 12 in tiles of <=4 -> 3 tiles.
  rand_tensor_yield<double, TA::SparsePolicy> yield_{world, 4, 12};
  yield_.set_max_tile(4);

  // contracts a1,a2 (unoccupied) -> batch mode is an unoccupied index (3 tiles)
  auto const expr = sequant::deserialize<sequant::ExprPtr>(
      L"g_{i1,i2}^{a1,a2} * t_{a1,a2}^{i3,i4}");
  std::string const target = "i_1,i_2,i_3,i_4";
  auto node = eval_node(expr);

  // The runtime batches STRICTLY on the optimizer's annotations (there is no
  // heuristic fallback), so this test must state the mode it wants batched.
  // batch_axis() picks exactly what the removed depth-0 heuristic used to pick.
  auto const ax = sequant::batch_axis(node);
  REQUIRE(ax.has_value());
  node->set_node_slice_mask({{*ax, sequant::BatchModeType::Contracted}});

  // Reference first, so yield_'s (random) leaf arrays are generated and cached;
  // the batched evaluator below copies yield_ and thus reuses the same arrays.
  auto const ref = evaluate(node, target, yield_)->get<TArrayD>();

  // Batched evaluation must reproduce the reference for any target batch size,
  // since sum_K = sum_{batches} sum_{K in batch}. The batch mode is unoccupied
  // (extent 12, tiles of 4 -> 3 tiles). target_batch_size is in *elements*:
  // 100 -> 1 batch (no-op), 8 -> 2 batches ([0,8),[8,12)), 4 -> 3 batches, and
  // 1 -> 3 batches (each tile its own batch).
  for (std::size_t target_batch_size :
       {std::size_t{100}, std::size_t{8}, std::size_t{4}, std::size_t{1}}) {
    auto cache = cache_t::empty();
    auto aops = yield_.array_ops();
    cache.set_array_ops(&aops);
    cache.set_custom_evaluator(make_batched_custom_evaluator(
        yield_, [target_batch_size](sequant::Index const&) -> std::size_t {
          return target_batch_size;
        }));
    auto const res = evaluate(node, target, yield_, cache)->get<TArrayD>();
    // batched summation reorders the contraction, so allow a looser FP margin
    REQUIRE(equal_tarrays<Loose>(res, ref));
  }
}

TEST_CASE("eval_batched_custom_evaluator persistence gate", "[eval]") {
  using sequant::evaluate;
  using sequant::make_batched_custom_evaluator;
  using TArrayD = TA::TSpArrayD;
  using node_t = sequant::FullBinaryNode<sequant::EvalExprTA>;
  using cache_t = sequant::CacheManager<node_t>;

  auto& world = TA::get_default_world();
  rand_tensor_yield<double, TA::SparsePolicy> yield_{world, 4, 12};
  yield_.set_max_tile(4);

  // Contracts a1,a2 (unoccupied, 3 tiles) -> batchable over an unoccupied mode.
  // The subtree contains a "t" leaf, which we treat as volatile.
  auto const expr = sequant::deserialize<sequant::ExprPtr>(
      L"g_{i1,i2}^{a1,a2} * t_{a1,a2}^{i3,i4}");
  std::string const target = "i_1,i_2,i_3,i_4";
  auto node = eval_node(expr);
  // The runtime batches STRICTLY on the optimizer's annotations (there is no
  // heuristic fallback), so this test must state the mode it wants batched.
  // batch_axis() picks exactly what the removed depth-0 heuristic used to pick.
  auto const ax = sequant::batch_axis(node);
  REQUIRE(ax.has_value());
  node->set_node_slice_mask({{*ax, sequant::BatchModeType::Contracted}});
  auto const ref = evaluate(node, target, yield_)->get<TArrayD>();

  // Volatile-leaf predicate: the amplitude "t".
  auto is_volatile_t = [](node_t const& n) {
    return n.leaf() && n->is_tensor() && n->as_tensor().label() == L"t";
  };

  // (1) Gate ON (persistent_only=true): the subtree contains a volatile "t"
  // leaf, so batching is DECLINED -- the spy scope-guard (invoked only once a
  // node passes every gate and yields >1 batch) is never called, yet the
  // standard scheme still gives the correct result.
  {
    bool batched = false;
    auto spy = [&batched](std::size_t) {
      batched = true;
      return sequant::no_scope_guard{};
    };
    auto cache = cache_t::empty();
    auto aops = yield_.array_ops();
    cache.set_array_ops(&aops);
    cache.set_custom_evaluator(make_batched_custom_evaluator(
        yield_,
        [](sequant::Index const&) -> std::size_t { return std::size_t{4}; },
        sequant::accept_any_index{}, spy, is_volatile_t,
        /*persistent_only=*/true));
    auto const res = evaluate(node, target, yield_, cache)->get<TArrayD>();
    REQUIRE_FALSE(batched);  // gated volatile subtree -> not batched
    REQUIRE(equal_tarrays<Loose>(res, ref));
  }

  // (2) No gate (default never_volatile): the SAME node DOES batch -- confirms
  // (1)'s decline is due to the volatility gate, not the index/tiling setup.
  {
    bool batched = false;
    auto spy = [&batched](std::size_t) {
      batched = true;
      return sequant::no_scope_guard{};
    };
    auto cache = cache_t::empty();
    auto aops = yield_.array_ops();
    cache.set_array_ops(&aops);
    cache.set_custom_evaluator(make_batched_custom_evaluator(
        yield_,
        [](sequant::Index const&) -> std::size_t { return std::size_t{4}; },
        sequant::accept_any_index{}, spy, sequant::never_volatile{}));
    auto const res = evaluate(node, target, yield_, cache)->get<TArrayD>();
    REQUIRE(batched);  // no gate -> batched as before
    REQUIRE(equal_tarrays<Loose>(res, ref));
  }
}

TEST_CASE("eval_batched_custom_evaluator_tot", "[eval]") {
  using sequant::evaluate;
  using sequant::make_batched_custom_evaluator;
  using node_t = sequant::FullBinaryNode<sequant::EvalExprTA>;
  using cache_t = sequant::CacheManager<node_t>;

  auto& world = TA::get_default_world();
  // Multi-tile the occupied space so the contracted occupied index i3 spans
  // more than one tile and batching over it actually engages: extent 8 in
  // tiles of <=4 -> 2 tiles. The inner (virtual) space is left single-tiled.
  rand_tensor_yield<double> yield{world, /*nocc=*/8, /*nvirt=*/3};
  yield.set_max_tile(4);

  using ArrayToT = typename decltype(yield)::array_tot_type;

  // ToT * ToT -> ToT (the same expression as the "tot" section above). The
  // contracted indices are the occupied i3 (an *outer* mode of both leaves)
  // and the inner virtual a4. Scoping the batch mode to the occupied space
  // selects i3, so the batched partials slice ToT leaves over i3
  // (slice_array_over_mode) and sum ToT partials (add_inplace) -- the two
  // annotation-free ToT array operations that must emit an "outer;inner"
  // annotation rather than a flat one (else DistArray's is_tot_index() trips).
  // i_3 is explicitly contracted (ket of I, bra of s) together with the
  // composite a_4<i_2,i_3> it is a proto of, so it leaves the result: a
  // genuinely contracted occupied index of a ToT x ToT -> ToT product.
  auto const expr = sequant::deserialize<sequant::ExprPtr>(
      L"I{a4<i2,i3>,a1<i1,i2>;i1,i3} * s{i3;a4<i2,i3>}");
  auto node = eval_node(expr);
  std::string const target = node->annot();
  // Reference first (non-batched), so yield's random leaf arrays are generated
  // and cached; the batched evaluator reuses the same arrays.
  auto const ref = evaluate(node, target, yield)->get<ArrayToT>();

  auto const occ =
      sequant::get_default_context().index_basis_registry()->retrieve_space(
          L"i");
  auto accept_occ = [occ](sequant::Index const& ix) {
    return ix.space() == occ;
  };

  // The runtime batches STRICTLY on the optimizer's annotations (there is no
  // heuristic fallback), so this test must state the mode it wants batched.
  // batch_axis() picks exactly what the removed depth-0 heuristic used to pick.
  auto const ax = sequant::batch_axis(node, accept_occ);
  REQUIRE(ax.has_value());
  node->set_node_slice_mask({{*ax, sequant::BatchModeType::Contracted}});

  // TA::norm2 is unsupported for tensor-of-tensor tiles, so compare via the
  // self-dot of each array (a scalar norm^2); reordering the contraction over
  // i3 must not change it.
  auto self_dot = [](auto const& arr) {
    return arr("i,j;a").dot(arr("i,j;a"));
  };
  auto const ref_dot = self_dot(ref);

  for (std::size_t target_batch_size :
       {std::size_t{100}, std::size_t{4}, std::size_t{1}}) {
    auto cache = cache_t::empty();
    auto aops = yield.array_ops();
    cache.set_array_ops(&aops);
    cache.set_custom_evaluator(make_batched_custom_evaluator(
        yield,
        [target_batch_size](sequant::Index const&) -> std::size_t {
          return target_batch_size;
        },
        accept_occ));
    auto const res = evaluate(node, target, yield, cache)->get<ArrayToT>();
    REQUIRE(self_dot(res) == Catch::Approx(ref_dot));
  }
}

TEST_CASE("ta_tot_conj_complex", "[eval]") {
  // Exercises TiledArray's .conj() on a complex tensor-of-tensors through the
  // expression engine — the capability SeQuant's ToT adjoint() relies on after
  // teaching TA to recurse conj into nested tiles. With genuinely complex data
  // a missing/incorrect inner conj would flip the wrong imaginary sign.
  using namespace sequant;
  auto& world = TA::get_default_world();
  size_t const nocc = 2, nvirt = 3;
  rand_tensor_yield<std::complex<double>, TA::SparsePolicy> yield{world, nocc,
                                                                  nvirt};
  using ArrayToT = typename decltype(yield)::array_tot_type;

  std::string const annot{"i_2,i_3;a_3i_2i_3,a_4i_2i_3"};
  auto const t = deserialize<sequant::ExprPtr>(L"t{a3<i2,i3>,a4<i2,i3>;i2,i3}");
  auto const& src = yield(t->as<sequant::Tensor>())->get<ArrayToT>();

  ArrayToT conjd;
  conjd(annot) = src(annot).conj();
  ArrayToT::wait_for_lazy_cleanup(world);

  // elementwise: conjd inner == conj(src inner); identical tiling/layout (no
  // permutation in the assignment) lets us walk local tiles in lockstep.
  auto it_s = src.begin();
  auto it_c = conjd.begin();
  for (; it_s != src.end(); ++it_s, ++it_c) {
    auto const& souter = it_s->get();
    auto const& couter = it_c->get();
    REQUIRE(souter.size() == couter.size());
    for (std::size_t o = 0; o < souter.size(); ++o) {
      auto const& sinner = souter[o];
      auto const& cinner = couter[o];
      if (sinner.empty()) continue;
      for (std::size_t k = 0; k < sinner.size(); ++k) {
        CHECK(cinner[k].real() == Catch::Approx(sinner[k].real()));
        CHECK(cinner[k].imag() == Catch::Approx(-sinner[k].imag()));
      }
    }
  }
}

TEST_CASE("eval_batched_scratch", "[eval]") {
  using node_t = sequant::FullBinaryNode<sequant::EvalExprTA>;
  using cache_t = sequant::CacheManager<node_t>;
  using sequant::Index;

  // W-analog of the PNO-CCSD PPL intermediate: two canonically-equal internal
  // siblings, both carrying the auxiliary batch mode x_1 (free at the
  // children, contracted at the root; an aux-aux edge, like the DF index K).
  // Every orbital contraction pairs a bra with a ket.
  auto const expr = sequant::deserialize<sequant::ExprPtr>(
      L"(g{a_2;i_1;x_1} * h{i_3;a_2}) * (g{a_3;i_2;x_1} * h{i_4;a_3})");
  auto const node = eval_node(expr);
  REQUIRE_FALSE(node.leaf());
  REQUIRE_FALSE(node.left().leaf());
  REQUIRE_FALSE(node.right().leaf());
  // structural precondition: the two children ARE canonically equal
  {
    auto cm = sequant::cache_manager(ranges::views::single(node), 2);
    REQUIRE(cm.max_life(node.left()) == 2);
  }
  auto const x1 = [&] {  // the contracted (batch) mode
    auto modes = sequant::contracted_indices(node);
    REQUIRE(modes.size() == 1);
    return modes[0];
  }();

  auto real = cache_t::empty();

  SECTION("registers consistent repeats with in-pass counts") {
    std::vector<std::pair<node_t const*, Index>> const members{{&node, x1}};
    auto bs = sequant::detail::make_batched_scratch(members, real);
    REQUIRE(bs.cache.exists(node.left()));
    REQUIRE(bs.cache.max_life(node.left()) == 2);
    REQUIRE_FALSE(bs.cache.exists(node));  // member roots are not registered
    REQUIRE(bs.seeds.empty());
  }

  SECTION("signature-inconsistent subnodes are not registered") {
    // the same subtree appears under two members, but the second member's
    // mode (an index the shared subnode does not carry) gives it signature
    // 'absent' while the first gives a position -> inconsistent -> unshared
    auto const expr2 = sequant::deserialize<sequant::ExprPtr>(
        L"(g{a_2;i_1;x_1} * h{i_3;a_2}) * p{i_5;i_6;x_1}");
    auto const node2 = eval_node(expr2);
    auto const bogus_axis = Index(L"i_9");
    std::vector<std::pair<node_t const*, Index>> const members{
        {&node, x1}, {&node2, bogus_axis}};
    auto bs = sequant::detail::make_batched_scratch(members, real);
    REQUIRE_FALSE(bs.cache.exists(node.left()));
  }

  SECTION("descends through inconsistently-sliced re-encounters") {
    // M1 (mode x_1) = X * D2 with X = (D * u): D's two occurrences (inside X
    // and as the root's sibling D2) are visited with the same signature, so D
    // alone would be registered. M2 (bogus mode) re-encounters X with
    // signature 'absent' -- inconsistent, unshared. The walk must descend
    // through that re-encounter: under it D's signature is also 'absent', so
    // sharing D would serve M2's (per-occurrence) evaluation of X a wrongly
    // sliced value. D must end up unregistered.
    auto const expr_m1 = sequant::deserialize<sequant::ExprPtr>(
        L"((g{a_2;i_1;x_1} * h{i_3;a_2}) * u{i_5;i_3}) * "
        L"(g{a_3;i_2;x_1} * h{i_4;a_3})");
    auto const m1 = eval_node(expr_m1);
    auto const expr_m2 = sequant::deserialize<sequant::ExprPtr>(
        L"((g{a_2;i_1;x_1} * h{i_3;a_2}) * u{i_5;i_3}) * p{i_6;i_7;x_1}");
    auto const m2 = eval_node(expr_m2);
    // structural preconditions: m1 = X * D2, X = D * u, D2 == D == m2's X
    // child canonically
    sequant::TreeNodeEqualityComparator<node_t> const eq;
    REQUIRE_FALSE(m1.left().leaf());
    REQUIRE_FALSE(m1.left().left().leaf());
    auto const& D = m1.left().left();
    REQUIRE(eq(D, m1.right()));
    REQUIRE(eq(m1.left(), m2.left()));

    // positive control: M1 alone registers D (count 2, consistent signature)
    {
      std::vector<std::pair<node_t const*, Index>> const members{{&m1, x1}};
      auto bs = sequant::detail::make_batched_scratch(members, real);
      REQUIRE(bs.cache.exists(D));
      REQUIRE(bs.cache.max_life(D) == 2);
    }
    // adding M2 makes X inconsistent AND must expose D's 'absent' signature
    // beneath X's second occurrence
    {
      auto const bogus_axis = Index(L"i_9");
      std::vector<std::pair<node_t const*, Index>> const members{
          {&m1, x1}, {&m2, bogus_axis}};
      auto bs = sequant::detail::make_batched_scratch(members, real);
      REQUIRE_FALSE(bs.cache.exists(m1.left()));  // X: inconsistent
      REQUIRE_FALSE(bs.cache.exists(D));  // D: inconsistent via pruned branch
    }
  }
}

TEST_CASE("eval_batched_custom_evaluator nests inner mode", "[eval]") {
  // Exactness gate: two batchable modes annotated at DIFFERENT nodes
  // of one tree must batch by nesting -- the outer node slices mode A at the
  // top and, WITHIN each A batch, the evaluator re-enters on the per-batch
  // scratch to slice mode B at an inner node (`for A-batch: for B-batch:
  // replay`). The inner slice must compose on top of the outer one; the nested
  // sum equals the unbatched contraction because `sum_K = sum_b sum_{K in b}`
  // applies per mode.
  using sequant::evaluate;
  using sequant::make_batched_custom_evaluator;
  using TArrayD = TA::TSpArrayD;
  using node_t = sequant::FullBinaryNode<sequant::EvalExprTA>;
  using cache_t = sequant::CacheManager<node_t>;

  auto& world = TA::get_default_world();
  // occ single-tiled (4); aux multi-tiled (12 in tiles of 4 -> 3 tiles) so both
  // batch modes slice into 3 batches each and depth-2 nesting genuinely
  // engages.
  rand_tensor_yield<double, TA::SparsePolicy> yield_{world, 4, 4, 12};
  yield_.set_max_tile(4);

  // Root (outer) contracts the aux mode x_1; the inner product (u*v) contracts
  // a DIFFERENT aux mode x_2. Crucially, the inner leaf u carries BOTH x_1 and
  // x_2, so the inner re-derivation MUST slice x_1 to the outer batch (compose
  // le_g) -- reusing the original leaf evaluator would leave u at the full x_1
  // extent and the x_1 contraction above would then mismatch. Every product is
  // written fully binary so binarize preserves the nesting.
  auto const expr = sequant::deserialize<sequant::ExprPtr>(
      L"(((u{i_1;i_2;x_1,x_2} * v{i_3;i_1;x_2}) * w{i_2;i_5;x_1})"
      L" * p{i_6;i_3;x_1})");
  std::string const target = "i_5,i_6";
  auto node = eval_node(expr);  // mutable: batch modes are annotated below

  auto const aux_space =
      sequant::get_default_context().index_basis_registry()->retrieve_space(
          L"x");
  auto accept_aux = [aux_space](sequant::Index const& ix) {
    return ix.space() == aux_space;
  };

  // Annotate the root with its aux mode (x_1) and the unique descendant that
  // contracts an aux mode (the inner u*v, x_2). Locating the inner node by its
  // contracted aux mode keeps this robust to binarize's operand ordering.
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
  REQUIRE(*inner_axis != *root_axis);  // the two nested modes differ
  (*inner)->set_node_slice_mask(
      {{*inner_axis, sequant::BatchModeType::Contracted}});

  // Reference: plain (unbatched) evaluation. Computed first so yield_'s random
  // leaf arrays are generated and cached; the batched evaluator reuses them.
  auto const ref = evaluate(node, target, yield_)->get<TArrayD>();

  // Spy scope-guard: records the batch count each time an evaluator FIRES (i.e.
  // picks a mode partitioning into >1 batch). Without the re-entrant scratch
  // only the root fires (one record); with nesting the inner evaluator fires
  // once per outer batch -- this vector is what proves depth-2 nesting engaged.
  std::vector<std::size_t> guard_calls;
  auto spy = [&guard_calls](std::size_t n) {
    guard_calls.push_back(n);
    return sequant::no_scope_guard{};
  };

  auto cache = cache_t::empty();
  auto aops = yield_.array_ops();
  cache.set_array_ops(&aops);
  cache.set_custom_evaluator(make_batched_custom_evaluator(
      yield_,
      [](sequant::Index const&) -> std::size_t { return std::size_t{4}; },
      accept_aux, spy, sequant::never_volatile{}));
  auto const res = evaluate(node, target, yield_, cache)->get<TArrayD>();

  // Exactness gate: the nested batched result equals the unbatched one.
  TArrayD diff;
  diff("i,j") = ref("i,j") - res("i,j");
  auto const err = TA::norm2(diff);
  REQUIRE(err < 1e-10);

  // Depth-2 nesting engaged: the root fires once over x_1 (3 batches) and, per
  // x_1 batch, the inner evaluator fires over x_2 (3 batches) -- 1 + 3 = 4
  // firings, each partitioning into 3 batches. Before the re-entrant scratch,
  // only the root fired (guard_calls == {3}); size > 1 is the RED/GREEN gate.
  REQUIRE(guard_calls.size() > 1);
  CHECK(guard_calls.size() == 4);
  for (auto n : guard_calls) CHECK(n == 3);
}

TEST_CASE("eval_batched_custom_evaluator nests two modes on one node",
          "[eval]") {
  // Finding N2 regression gate: the single-term-opt DP can stamp MORE THAN ONE
  // batch mode on a SINGLE eval node (aprime is a bitmask;
  // reconstruct_batched_modes pushes one Index per set bit into that node's
  // node_slice_mask()). The runtime must slice ALL of them by nesting -- `for
  // K-batch: for mu1-batch: replay` at the SAME node -- otherwise the modeled
  // peak (which priced BOTH modes sliced) is a lie. Here ONE product node
  // carries two aux batch modes x_1 and x_2, both present on both leaves; the
  // evaluator must slice x_1 at depth 0 and, WITHIN each x_1 batch, re-enter on
  // the SAME node to slice x_2 at depth 1. Before the sliceability-aware pick
  // the runtime re-picked x_1 (the first annotated mode) on the re-entry, found
  // it already fully sliced (one batch), and DECLINED -- x_2 was never sliced
  // (max_depth stayed 1).
  using sequant::evaluate;
  using sequant::make_batched_custom_evaluator;
  using TArrayD = TA::TSpArrayD;
  using node_t = sequant::FullBinaryNode<sequant::EvalExprTA>;
  using cache_t = sequant::CacheManager<node_t>;

  auto& world = TA::get_default_world();
  // occ single-tiled (4); aux multi-tiled (12 in tiles of 4 -> 3 tiles).
  rand_tensor_yield<double, TA::SparsePolicy> yield_{world, 4, 4, 12};
  yield_.set_max_tile(4);

  // ONE product contracting i_1 (occ) plus BOTH aux modes x_1 and x_2; each
  // leaf carries x_1 and x_2, so slicing either aux mode slices both leaves.
  // Result free indices are i_5, i_6.
  auto const expr = sequant::deserialize<sequant::ExprPtr>(
      L"(g{i_1;i_5;x_1,x_2} * h{i_6;i_1;x_1,x_2})");
  std::string const target = "i_5,i_6";
  auto node = eval_node(expr);  // mutable: both batch modes annotated below

  auto const aux_space =
      sequant::get_default_context().index_basis_registry()->retrieve_space(
          L"x");
  auto accept_aux = [aux_space](sequant::Index const& ix) {
    return ix.space() == aux_space;
  };

  // Stamp BOTH aux contracted indices on the single node (the two-set-bit
  // aprime case). contracted_indices lists them; keep the aux ones.
  sequant::container::svector<sequant::Index> two_axes;
  for (sequant::Index const& ix : sequant::contracted_indices(node))
    if (accept_aux(ix)) two_axes.push_back(ix);
  REQUIRE(two_axes.size() == 2);
  REQUIRE(two_axes[0] != two_axes[1]);
  sequant::container::svector<std::pair<sequant::Index, sequant::BatchModeType>>
      two_axes_typed;
  for (sequant::Index const& ix : two_axes)
    two_axes_typed.push_back({ix, sequant::BatchModeType::Contracted});
  node->set_node_slice_mask(two_axes_typed);

  // Reference: plain (unbatched) evaluation; also generates yield_'s random
  // leaf arrays that the batched evaluator reuses.
  auto const ref = evaluate(node, target, yield_)->get<TArrayD>();

  // Distinct per-mode target sizes so the two nested levels realize DIFFERENT
  // batch counts. Batches snap to whole tiles of 4 and the target is an upper
  // bound: x_1 -> target 4 over extent 12 = 3 single-tile batches; x_2 ->
  // target 8 = 2 batches (two 4-tiles fit, then the last). A depth-2 product of
  // 3*2 = 6 (rather than 3*3 or 2*2) proves the outer and inner levels sliced
  // DISTINCT modes with distinct partitions.
  auto const& x1 = two_axes[0];
  auto target_batch_size = [x1](sequant::Index const& ix) -> std::size_t {
    return ix == x1 ? std::size_t{4} : std::size_t{8};
  };

  // Tracking guard (as in the compose test): records the live-guard stack so
  // depth-2 simultaneity and the per-instant product are observable.
  struct GuardState {
    std::vector<std::size_t> live;
    std::size_t max_depth = 0;
    std::vector<std::size_t> counts;
    std::vector<std::size_t> products_at_depth2;
  } state;
  struct TrackingGuard {
    GuardState* st;
    TrackingGuard(GuardState* s, std::size_t n) : st(s) {
      st->counts.push_back(n);
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
      yield_, target_batch_size, accept_aux, make_tracking_guard,
      sequant::never_volatile{}));
  auto const res = evaluate(node, target, yield_, cache)->get<TArrayD>();

  // Exactness: nesting two modes on one node changes no numerics.
  TArrayD diff;
  diff("i,j") = ref("i,j") - res("i,j");
  auto const err = TA::norm2(diff);
  REQUIRE(err < 1e-10);

  // Depth 2 reached on ONE node: x_1 sliced at depth 0, x_2 at depth 1. Before
  // the fix the node re-picked x_1 and declined -> max_depth would stay 1.
  REQUIRE(state.max_depth == 2);
  // The outer level fires once over x_1 (3 batches); per outer batch the inner
  // level fires once over x_2 (2 batches) -> 1 + 3 = 4 firings, counts a
  // permutation of {3, 2, 2, 2}.
  REQUIRE(state.counts.size() == 4);
  CHECK(std::count(state.counts.begin(), state.counts.end(), 3u) == 1);
  CHECK(std::count(state.counts.begin(), state.counts.end(), 2u) == 3);
  // Every depth-2 instant multiplies 3 (outer x_1) by 2 (inner x_2) = 6, i.e.
  // two DISTINCT modes with DISTINCT partitions were live together.
  REQUIRE(state.products_at_depth2.size() == 3);
  for (auto const p : state.products_at_depth2) CHECK(p == 6);
  CHECK(state.live.empty());
}

TEST_CASE("batched-hoist-loop-invariant", "[eval]") {
  // Hoisting gate: a network batched over an OUTER mode (x_1) whose
  // subtree contains an intermediate INVARIANT to that mode (I2 = a*b contracts
  // x_2 and carries no x_1) must build the invariant ONCE at its emitted scope
  // level and reuse it across the x_1 batches -- NOT rebuild it per batch. The
  // root R contracts x_1 (the outer batch loop, depth 0); I2 is a descendant
  // with scope_level == -1 (invariant to the whole nest), so it hoists to the
  // real cache and the scope chain serves it to every x_1 batch body.
  using sequant::evaluate;
  using sequant::make_batched_custom_evaluator;
  using TArrayD = TA::TSpArrayD;
  using node_t = sequant::FullBinaryNode<sequant::EvalExprTA>;
  using cache_t = sequant::CacheManager<node_t>;

  auto& world = TA::get_default_world();
  // occ single-tiled (4); aux multi-tiled (12 in tiles of 4 -> 3 tiles), so
  // x_1 slices into 3 batches.
  rand_tensor_yield<double, TA::SparsePolicy> yield_{world, 4, 4, 12};
  yield_.set_max_tile(4);

  // R = (((g{i_1;i_2;x_2} * h{i_1;i_3;x_2}) * w{i_2;i_5;x_1}) *
  // p{i_3;i_6;x_1}), a left-deep chain (same shape as the two-mode-nesting
  // test) whose INNERMOST pair is the x_1-invariant intermediate:
  // - I2 = g*h contracts i_1 and x_2 -> {i_2;i_3}: carries NO aux, so it is
  //   invariant to the outer x_1 loop.
  // - M = I2*w contracts i_2 -> {i_3;i_5;x_1}: carries x_1 (not hoistable).
  // - R = M*p contracts i_3 and x_1 -> {i_5;i_6}: the x_1 batch trigger.
  auto const expr = sequant::deserialize<sequant::ExprPtr>(
      L"(((g{i_1;i_2;x_2} * h{i_3;i_1;x_2}) * w{i_2;i_5;x_1})"
      L" * p{i_6;i_3;x_1})");
  std::string const target = "i_5,i_6";
  auto node = eval_node(expr);

  auto const aux_space =
      sequant::get_default_context().index_basis_registry()->retrieve_space(
          L"x");
  auto accept_aux = [aux_space](sequant::Index const& ix) {
    return ix.space() == aux_space;
  };

  // Root batches over its contracted aux mode x_1 (outer loop). Its own
  // scope_level is irrelevant (the trigger is never hoisted); leave it unset.
  auto const root_axis = sequant::batch_axis(node, accept_aux);
  REQUIRE(root_axis.has_value());
  node->set_node_slice_mask({{*root_axis, sequant::BatchModeType::Contracted}});

  // Locate I2 = the unique NON-root node that contracts an aux mode (x_2).
  node_t* i2 = nullptr;
  std::optional<sequant::Index> i2_axis;
  auto find_i2 = [&](auto&& self, node_t& n) -> void {
    if (n.leaf()) return;
    if (&n != &node) {
      if (auto ax = sequant::batch_axis(n, accept_aux)) {
        i2 = &n;
        i2_axis = *ax;
      }
    }
    self(self, n.left());
    self(self, n.right());
  };
  find_i2(find_i2, node);
  REQUIRE(i2 != nullptr);
  REQUIRE(i2_axis.has_value());
  REQUIRE(*i2_axis != *root_axis);  // I2 contracts a DIFFERENT aux mode
  // Annotate I2: batched over x_2, order-aware (emitted), and EMPTY residency
  // (empty sliced_modes) -- it is invariant to the whole enclosing nest, so
  // per-level placement hoists it above every batch loop to the real/term
  // (root) cache and builds it once. This drives the post-hoist_invariants
  // residency signals (batch_order_aware + sliced_modes), not a per-node scalar
  // placement level.
  (*i2)->set_node_slice_mask({{*i2_axis, sequant::BatchModeType::Contracted}});
  (*i2)->set_batch_order_aware(true);

  // Reference: plain (unbatched) evaluation; also fills yield_'s leaf cache.
  auto const ref = evaluate(node, target, yield_)->get<TArrayD>();

  // Counting leaf evaluator: counts how many times leaf `g` (which appears ONLY
  // inside I2) is requested. I2 built once => g requested once; I2 rebuilt per
  // x_1 (and per inner x_2) batch => requested many times. This is the
  // build-count proxy for I2.
  int g_evals = 0;
  auto counting_yield = [&yield_,
                         &g_evals](node_t const& leaf) -> sequant::ResultPtr {
    if (leaf->is_tensor() && leaf->as_tensor().label() == L"g") ++g_evals;
    return yield_(leaf);
  };

  auto cache = cache_t::empty();
  auto aops = yield_.array_ops();
  cache.set_array_ops(&aops);
  cache.set_custom_evaluator(make_batched_custom_evaluator(
      counting_yield,
      [](sequant::Index const&) -> std::size_t { return std::size_t{4}; },
      accept_aux, sequant::make_no_scope_guard{}, sequant::never_volatile{}));
  auto const res =
      evaluate(node, target, counting_yield, cache)->get<TArrayD>();

  // Exactness gate: hoisting changes no numerics.
  TArrayD diff;
  diff("i,j") = ref("i,j") - res("i,j");
  auto const err = TA::norm2(diff);
  REQUIRE(err < 1e-10);

  // Hoisting gate (RED/GREEN): the x_1-invariant intermediate I2 is built ONCE,
  // not once per x_1 batch. Before the hoist path, I2 was rebuilt (and itself
  // x_2-batched) inside every x_1 batch, so leaf `g` was evaluated many times.
  CHECK(g_evals == 1);
}

TEST_CASE("batched-hoist-intermediate-level", "[eval]") {
  // The runtime residency is sliced_modes(),
  // the all-batched-modes cross-occurrence meet -- which carries a node's
  // enclosing CONTRACTED (aux) modes directly, since a node variant to an outer
  // aux loop carries that aux free on a result slot. The pre-existing hoist
  // test only exercises the EMPTY-residency case (whole-nest invariant ->
  // root). This test exercises the NON-EMPTY case: a node P carries an OUTER
  // contracted mode (x_1) but not the INNER one (x_2) its containing trigger
  // batches over, so sliced_modes(P) == {x_1} and its home level is the OUTER
  // (x_1) scratch, not the term/root cache and not T's own inner (x_2) scratch.
  // P must therefore be built ONCE PER OUTER (x_1) BLOCK and reused across
  // every inner (x_2) batch within that block: rebuilt when x_1 advances,
  // untouched while only x_2 advances.
  //
  // Network (left-deep chain; verified by direct inspection of the binarized
  // tree -- canon_indices()/batch_axis() are LOCAL to a node's own top-level
  // pairing, not recomputed transitively, so the grouping below is exactly
  // what deserialize+binarize produce for this expression):
  //   D2 = q{i_1;a_9;x_1} * r{a_9;i_2}         -> {i_1;i_2;x_1}: carries the
  //        OUTER mode x_1, does NOT carry the INNER mode x_2 -- the witness
  //        node (non-empty sliced_modes = {x_1}).
  //   Q  = D2 * m{;;x_1,x_2}                    -> {i_1;i_2;x_2}: a drop-in
  //        structural replacement for a plain 3-index leaf (same free-index
  //        role as `g` in the "hoists loop-invariant descendant" test), so Q
  //        contracts x_1 internally and carries x_2 onward -- the rest of
  //        the network is unaffected by what is inside Q.
  //   T  = Q * h{i_1;i_3;x_2}                    -> {i_2;i_3;x_1}: contracts
  //        x_2 (T is the INNER trigger; D2 is T's grandchild).
  //   M  = T * w{i_2;i_5;x_1}                    -> {i_3;i_5;x_1}
  //   R  = M * p{i_3;i_6;x_1}                    -> {i_5;i_6}: contracts x_1
  //        (R is the OUTER trigger).
  using sequant::evaluate;
  using sequant::make_batched_custom_evaluator;
  using TArrayD = TA::TSpArrayD;
  using node_t = sequant::FullBinaryNode<sequant::EvalExprTA>;
  using cache_t = sequant::CacheManager<node_t>;

  auto& world = TA::get_default_world();
  // occ/virt single-tiled (4); aux multi-tiled (12 in tiles of 4 -> 3 blocks)
  // so both the outer (x_1) and inner (x_2) loops genuinely partition into
  // more than one block each.
  rand_tensor_yield<double, TA::SparsePolicy> yield_{world, 4, 4, 12};
  yield_.set_max_tile(4);

  auto const expr = sequant::deserialize<sequant::ExprPtr>(
      L"(((((q{i_1;a_9;x_1} * r{a_9;i_2}) * m{;;x_1,x_2}) * h{i_3;i_1;x_2}) * "
      L"w{i_2;i_5;x_1}) * p{i_6;i_3;x_1})");
  std::string const target = "i_5,i_6";
  auto node = eval_node(expr);

  auto const aux_space =
      sequant::get_default_context().index_basis_registry()->retrieve_space(
          L"x");
  auto accept_aux = [aux_space](sequant::Index const& ix) {
    return ix.space() == aux_space;
  };

  // Root's own aux mode (x_1, the OUTER trigger).
  auto const root_axis = sequant::batch_axis(node, accept_aux);
  REQUIRE(root_axis.has_value());
  node->set_node_slice_mask({{*root_axis, sequant::BatchModeType::Contracted}});

  // T = the unique non-root node contracting an aux mode different from
  // root's (x_2, the INNER trigger). Located via batch_axis (robust to
  // binarize's operand ordering), matching the existing hoist test's pattern.
  node_t* t = nullptr;
  std::optional<sequant::Index> t_axis;
  auto find_t = [&](auto&& self, node_t& n) -> void {
    if (n.leaf()) return;
    if (&n != &node) {
      if (auto ax = sequant::batch_axis(n, accept_aux)) {
        t = &n;
        t_axis = *ax;
      }
    }
    self(self, n.left());
    self(self, n.right());
  };
  find_t(find_t, node);
  REQUIRE(t != nullptr);
  REQUIRE(t_axis.has_value());
  REQUIRE(*t_axis != *root_axis);
  (*t)->set_node_slice_mask({{*t_axis, sequant::BatchModeType::Contracted}});

  // D2 = T's grandchild (T.left().left()), the q*r product. Reached by
  // direct structural navigation (confirmed empirically for this exact
  // expression, since batch_axis cannot locate D2 -- it contracts no aux
  // mode itself); guarded by a structural REQUIRE so a future
  // canonicalization change fails loudly here instead of silently testing
  // the wrong node.
  node_t& d2 = node.left().left().left().left();
  REQUIRE(sequant::index_position(d2, *root_axis).has_value());  // carries x_1
  REQUIRE_FALSE(sequant::index_position(d2, *t_axis).has_value());  // not x_2
  REQUIRE(d2.left().leaf());
  REQUIRE(d2.right().leaf());
  REQUIRE(d2.left()->as_tensor().label() == L"q");
  REQUIRE(d2.right()->as_tensor().label() == L"r");

  // Annotate D2: order-aware, with NON-EMPTY sliced_modes = {x_1} (the signal
  // under test -- the all-batched-modes meet carries the enclosing aux mode
  // x_1 that D2 holds free on its result) and no External node_slice_mask stamp
  // (so it is not demoted). Its residency is {x_1} only -- present at T's
  // firing (ectx = [x_1]) but NOT at root's own (ectx = []) -- so it is
  // correctly skipped at the outer pre-loop placement call and collected at
  // T's own (inner) placement call, homed at the OUTER (x_1) scratch level,
  // not the term/root cache and not T's own (x_2) scratch.
  d2->set_batch_order_aware(true);
  d2->set_sliced_modes({*root_axis});

  // Reference: plain (unbatched) evaluation; also fills yield_'s leaf cache.
  auto const ref = evaluate(node, target, yield_)->get<TArrayD>();

  // Counting leaf evaluator: `r` carries no aux index at all, so it is never
  // touched by the sliceability PROBE for either x_1 or x_2 (which only ever
  // requests a leaf that actually carries the probed mode) -- it is touched
  // ONLY when D2 is materialized. This isolates D2's build count, mirroring
  // the h_evals technique used for the external-carrier witness below.
  int r_evals = 0;
  auto counting_yield = [&yield_,
                         &r_evals](node_t const& leaf) -> sequant::ResultPtr {
    if (leaf->is_tensor() && leaf->as_tensor().label() == L"r") ++r_evals;
    return yield_(leaf);
  };

  auto cache = cache_t::empty();
  auto aops = yield_.array_ops();
  cache.set_array_ops(&aops);
  cache.set_custom_evaluator(make_batched_custom_evaluator(
      counting_yield,
      [](sequant::Index const&) -> std::size_t { return std::size_t{4}; },
      accept_aux, sequant::make_no_scope_guard{}, sequant::never_volatile{}));
  auto const res =
      evaluate(node, target, counting_yield, cache)->get<TArrayD>();

  // Exactness gate: placing D2 at the x_1 scratch (not term-global, not
  // rebuilt per x_2 sub-batch) changes no numerics -- it is genuinely
  // resliced to the correct x_1 block each time x_1 advances.
  TArrayD diff;
  diff("i,j") = ref("i,j") - res("i,j");
  REQUIRE(TA::norm2(diff) < 1e-10);

  // Placement gate: x_1 and x_2 both partition the 12-element aux extent into
  // 3 blocks of 4 with this target size, so a CORRECTLY-homed D2 is rebuilt
  // once per outer block: r_evals == 3. (Empirically checked, not shipped
  // here, since the fetch-time slice-on-use safety net keeps VALUES correct
  // either way -- this is a footprint/avoidable-work gate, not an exactness
  // one, matching the C60 story where a misplaced node is a peak/time bug,
  // not a wrong-answer bug: with d2->set_sliced_modes({*root_axis})
  // commented out, D2's residency is wrongly empty, it is misclassified as a
  // whole-nest invariant and hoisted to the term/root cache instead --
  // r_evals == 1, i.e. held at the FULL x_1 extent for the entire run
  // instead of resliced per block; with batch_order_aware(true) also removed
  // (no hoisting at all -- the pre-existing behavior before this residency
  // component existed), D2 is rebuilt inline once per (x_1, x_2) pair --
  // r_evals == 9.)
  CHECK(r_evals == 3);
}
