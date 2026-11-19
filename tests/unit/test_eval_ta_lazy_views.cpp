#include "test_eval_ta_fixtures.hpp"

using namespace sequant::tests::eval_ta;

// A flat TiledArray result may stand for a pending view of another result's
// array (detail::ArrayView): a batch slice, a mode permutation, a conjugation,
// a scaling by a constant. The view shares the parent's tiles until a consumer
// needs an array. These tests pin the contract: a view is an alias (no bytes
// of its own) until get<>() materializes it; the expression consumers (the
// products, dot, sums, accumulation, the scatter) take the view as is and
// leave it an alias; the array consumers (einsum, clone) materialize it
// exactly; writes into a view never reach the parent.

// The tiles a consuming TiledArray expression receives from a view: what
// with_expr() hands over, evaluated through its own engine in the view's
// stored mode order (so no permutation is applied), each local lazy tile
// converted as a product or a sum would convert it. Returns the evaluator's
// TiledRange and the (ordinal, tile) pairs of this rank's tiles.
template <typename ResultT>
auto view_consumer_tiles(ResultT const& r, std::string const& annot) {
  return r.with_expr(annot, [&](auto&& expr) {
    using Expr = std::decay_t<decltype(expr)>;
    typename TA::expressions::ExprTrait<Expr>::engine_type engine(expr);
    engine.init(TA::get_default_world(), {},
                TA::expressions::BipartiteIndexList(expr.annotation()));
    auto dist_eval = engine.make_dist_eval();
    dist_eval.eval();
    using Tile = typename decltype(dist_eval)::eval_type;
    std::vector<std::pair<std::size_t, Tile>> tiles;
    for (auto ord : *dist_eval.pmap())
      tiles.emplace_back(ord, static_cast<Tile>(dist_eval.get(ord).get()));
    dist_eval.wait();
    return std::make_pair(dist_eval.trange(), std::move(tiles));
  });
}

// true iff every tile of \p tiles (from view_consumer_tiles) that the parent
// array \p parent holds on this rank IS the parent's tile (same buffer); the
// view's tile (ordinal in \p view_trange) maps to the parent's at + \p lo
// (tile coordinates). A tile the parent holds elsewhere is a fetched copy by
// construction and is not inspected.
template <typename ArrayT, typename Tiles>
bool consumer_reads_parent_tiles(ArrayT const& parent,
                                 TA::TiledRange const& view_trange,
                                 Tiles const& tiles,
                                 std::vector<std::size_t> const& lo) {
  std::size_t checked = 0;
  for (auto const& [ord, tile] : tiles) {
    auto idx = view_trange.tiles_range().idx(ord);
    for (std::size_t d = 0; d < idx.size(); ++d) idx[d] += lo[d];
    if (!parent.is_local(idx)) continue;
    auto const parent_tile = parent.find_local(idx).get();
    if (static_cast<void const*>(tile.data()) !=
        static_cast<void const*>(parent_tile.data()))
      return false;
    ++checked;
  }
  // every rank holds some tiles of an array this small with no zero tiles
  return checked > 0;
}

TEST_CASE("ta_lazy_views", "[eval][ta][lazy]") {
  using sequant::DeNest;
  using sequant::eval_result;
  using sequant::ResultPtr;
  using sequant::ResultScalar;
  using sequant::ResultTensorTA;
  using TArrayD = TA::TSpArrayD;
  using TArrayZ = TA::TSpArrayZ;
  using ResultD = ResultTensorTA<TArrayD>;
  using ResultZ = ResultTensorTA<TArrayZ>;
  using Bounds = std::vector<std::size_t>;

  auto& world = TA::get_default_world();
  auto annot3 = [](std::string l, std::string r, std::string c) {
    return std::array<std::any, 3>{std::move(l), std::move(r), std::move(c)};
  };
  auto annot2 = [](std::string pre, std::string post) {
    return std::array<std::any, 2>{std::move(pre), std::move(post)};
  };
  // the eager slice: tiles [lo, hi) of every mode, element coordinates kept
  auto eager_block = [](auto const& arr, std::string const& ann, Bounds lo,
                        Bounds hi) {
    std::decay_t<decltype(arr)> out;
    out(ann) = arr(ann).block(lo, hi, TA::preserve_lobound);
    return out;
  };

  // A(a,b,c): tiles 3 x 4 x 2; B(c,d): tiles 2 x 2
  TA::TiledRange const trA{{0, 3, 6, 9}, {0, 2, 4, 6, 8}, {0, 5, 10}};
  TA::TiledRange const trB{{0, 5, 10}, {0, 4, 7}};
  TArrayD A(world, trA), B(world, trB);
  A.fill_random();
  B.fill_random();
  world.gop.fence();
  auto rA = eval_result<ResultD>(A);
  auto rB = eval_result<ResultD>(B);

  SECTION("the trace layout reports where the tiles live, on request") {
    EvalLogGuard const restore_log;
    auto& logger = sequant::Logger::instance();
    std::ostringstream os;
    os << A.trange();
    auto const layout = os.str();
    REQUIRE(rA->trange_annot() == layout);

    // the kind of the process map
    logger.eval.placement = 1;
    auto const kind = rA->trange_annot();
    REQUIRE(kind.starts_with(layout + " pmap="));
    REQUIRE(kind.ends_with("Pmap"));
    // a pending view reports its base's map, flagged, and stays pending
    auto s = rA->slice_mode(1, 2, 6);
    REQUIRE(s->trange_annot().ends_with("Pmap:view"));
    REQUIRE(s->is_buffer_alias());

    // also the bytes on each rank: max/mean, then one entry per rank
    logger.eval.placement = 2;
    auto const load = rA->trange_annot();
    REQUIRE(load.starts_with(kind + " load="));
    REQUIRE(load.ends_with("]"));
    auto const per_rank = load.substr(load.find('[', kind.size()));
    REQUIRE(std::count(per_rank.begin(), per_rank.end(), ',') ==
            world.size() - 1);
    if (world.size() == 1) REQUIRE(load.starts_with(kind + " load=1.00["));

    logger.eval.placement = 0;
    REQUIRE(rA->trange_annot() == layout);
  }

#ifdef SEQUANT_EVAL_TRACE
  SECTION("a batch movement record carries the backend's trace suffix") {
    std::ostringstream trace;  // outlives the guard, which unhooks it
    int calls = 0;             // likewise, for the hook capturing it
    EvalLogGuard const restore_log;
    auto& logger = sequant::Logger::instance();
    logger.eval.level = 1;
    logger.eval.stream = &trace;

    auto s = rA->slice_mode(1, 2, 6);
    (void)s->get<TArrayD>();  // materializes: one record, no suffix
    REQUIRE(trace.str().find("Batch | Materialize | ") == 0);
    REQUIRE(trace.str().find("suffix") == std::string::npos);

    trace.str({});
    logger.eval.heap_stats = [&calls] {
      ++calls;
      return std::string("suffix=1");
    };
    auto s2 = rA->slice_mode(1, 2, 6);
    (void)s2->get<TArrayD>();
    REQUIRE(calls == 1);
    REQUIRE(trace.str().find("B | suffix=1 | ") != std::string::npos);
  }
#endif

  SECTION("a batch slice is an alias until its value is asked for") {
    // mode 1, elements [2, 6) = tiles 1..2
    auto s = rA->slice_mode(1, 2, 6);
    REQUIRE(s->is_buffer_alias());
    REQUIRE(s->size_in_bytes() == 0);
    REQUIRE(s->measured_size_in_bytes() == 0);
    auto const exp = eager_block(A, "a,b,c", {0, 1, 0}, {3, 3, 2});
    // the layout it reports is the slice's, without materializing
    std::ostringstream os;
    os << exp.trange();
    REQUIRE(s->trange_annot() == os.str());
    REQUIRE(s->is_buffer_alias());
    // the value is the eager slice; asking for it materializes once
    REQUIRE(equal_tarrays(s->get<TArrayD>(), exp, "a,b,c"));
    REQUIRE_FALSE(s->is_buffer_alias());
    REQUIRE(s->measured_size_in_bytes() > 0);
    // A is untouched
    REQUIRE(A.trange() == trA);
  }

  SECTION("the consumer of a view reads the parent's tiles in place") {
    // what a product or a sum receives from a view is the parent's own tile,
    // not a copy: TiledArray hands the tiles of a const array's block that
    // keeps the array's element coordinates through, and with_expr() hands
    // the parent over as a const array
    {
      // a slice: mode 1, tiles 1..2
      auto s = rA->slice_mode(1, 2, 6);
      auto const [tr, tiles] = view_consumer_tiles(s->as<ResultD>(), "a,b,c");
      REQUIRE(consumer_reads_parent_tiles(A, tr, tiles, {0, 1, 0}));
      REQUIRE(s->is_buffer_alias());
    }
    {
      // a slice of a slice: mode 1 tiles 1..2, then mode 0 tile 2
      auto ss = rA->slice_mode(1, 2, 6)->slice_mode(0, 6, 9);
      auto const [tr, tiles] = view_consumer_tiles(ss->as<ResultD>(), "a,b,c");
      REQUIRE(consumer_reads_parent_tiles(A, tr, tiles, {2, 1, 0}));
    }
    {
      // a permuted slice: the permutation only relabels the parent's modes,
      // so the tiles are still the parent's
      auto ps = rA->slice_mode(1, 2, 6)->permute(annot2("a,b,c", "c,a,b"));
      auto const [tr, tiles] = view_consumer_tiles(ps->as<ResultD>(), "c,a,b");
      REQUIRE(consumer_reads_parent_tiles(A, tr, tiles, {0, 1, 0}));
    }
    {
      // a whole-array view (a permutation without a slice) too
      auto pv = rA->permute(annot2("a,b,c", "b,a,c"));
      auto const [tr, tiles] = view_consumer_tiles(pv->as<ResultD>(), "b,a,c");
      REQUIRE(consumer_reads_parent_tiles(A, tr, tiles, {0, 0, 0}));
    }
    {
      // a conjugated view of a complex array computes its tiles: a copy
      TArrayZ Z(world, trA);
      Z.fill_random();
      world.gop.fence();
      auto rZ = eval_result<ResultZ>(Z);
      auto cz = rZ->slice_mode(1, 2, 6)->adjoint(annot2("a,b,c", "a,b,c"));
      auto const [tr, tiles] = view_consumer_tiles(cz->as<ResultZ>(), "a,b,c");
      REQUIRE_FALSE(consumer_reads_parent_tiles(Z, tr, tiles, {0, 1, 0}));
    }
  }

  SECTION("a product takes a sliced operand as is") {
    // batch the contracted index c: slices of A (mode 2) and B (mode 0) over
    // the tile [5, 10)
    auto sA = rA->slice_mode(2, 5, 10);
    auto sB = rB->slice_mode(0, 5, 10);
    auto prod = sA->prod(*sB, annot3("a,b,c", "c,d", "a,b,d"), DeNest::False);
    // neither operand was materialized
    REQUIRE(sA->is_buffer_alias());
    REQUIRE(sB->is_buffer_alias());
    REQUIRE_FALSE(prod->is_buffer_alias());
    TArrayD exp;
    exp("a,b,d") = eager_block(A, "a,b,c", {0, 0, 1}, {3, 4, 2})("a,b,c") *
                   eager_block(B, "c,d", {1, 0}, {2, 2})("c,d");
    REQUIRE(equal_tarrays(prod->get<TArrayD>(), exp, "a,b,d"));

    // the batches add up to the full contraction
    auto s0A = rA->slice_mode(2, 0, 5);
    auto s0B = rB->slice_mode(0, 0, 5);
    auto p0 = s0A->prod(*s0B, annot3("a,b,c", "c,d", "a,b,d"), DeNest::False);
    p0->add_inplace(*prod);
    TArrayD full;
    full("a,b,d") = A("a,b,c") * B("c,d");
    REQUIRE(equal_tarrays(p0->get<TArrayD>(), full, "a,b,d"));
  }

  SECTION("a product with a Hadamard index takes the views as is") {
    // H(a,c): a appears on both operands and the result (a general product,
    // which the expression layer evaluates natively)
    TArrayD H(world, TA::TiledRange{{0, 3, 6, 9}, {0, 5, 10}});
    H.fill_random();
    world.gop.fence();
    auto rH = eval_result<ResultD>(H);
    auto sA = rA->slice_mode(2, 5, 10);
    auto sH = rH->slice_mode(1, 5, 10);
    auto prod = sA->prod(*sH, annot3("a,b,c", "a,c", "a,b"), DeNest::False);
    REQUIRE(sA->is_buffer_alias());
    REQUIRE(sH->is_buffer_alias());
    auto exp =
        TA::einsum(eager_block(A, "a,b,c", {0, 0, 1}, {3, 4, 2})("a,b,c"),
                   eager_block(H, "a,c", {0, 1}, {3, 2})("a,c"), "a,b");
    REQUIRE(equal_tarrays(prod->get<TArrayD>(), exp, "a,b"));
  }

  SECTION("a slice of a slice") {
    auto s = rA->slice_mode(1, 2, 6);
    auto ss = s->slice_mode(2, 5, 10);
    REQUIRE(ss->is_buffer_alias());
    REQUIRE(s->is_buffer_alias());
    auto const exp = eager_block(A, "a,b,c", {0, 1, 1}, {3, 3, 2});
    REQUIRE(equal_tarrays(ss->get<TArrayD>(), exp, "a,b,c"));
    // narrowing the already sliced mode
    auto sn = s->slice_mode(1, 4, 6);
    REQUIRE(equal_tarrays(sn->get<TArrayD>(),
                          eager_block(A, "a,b,c", {0, 2, 0}, {3, 3, 2}),
                          "a,b,c"));
    // slicing a view further never materializes it
    REQUIRE(s->is_buffer_alias());
    // a slice of a materialized result is an alias of that result
    REQUIRE(equal_tarrays(s->get<TArrayD>(),
                          eager_block(A, "a,b,c", {0, 1, 0}, {3, 3, 2}),
                          "a,b,c"));
    REQUIRE_FALSE(s->is_buffer_alias());
    auto sm = s->slice_mode(2, 0, 5);
    REQUIRE(sm->is_buffer_alias());
    REQUIRE(equal_tarrays(sm->get<TArrayD>(),
                          eager_block(A, "a,b,c", {0, 1, 0}, {3, 3, 1}),
                          "a,b,c"));
  }

  SECTION("a permutation is lazy and composes with a slice") {
    auto p = rA->permute(annot2("a,b,c", "c,a,b"));
    REQUIRE(p->is_buffer_alias());
    TArrayD expP;
    expP("c,a,b") = A("a,b,c");
    std::ostringstream os;
    os << expP.trange();
    REQUIRE(p->trange_annot() == os.str());
    // slicing the permuted view's mode 0 (= c of A)
    auto ps = p->slice_mode(0, 5, 10);
    REQUIRE(ps->is_buffer_alias());
    REQUIRE(p->is_buffer_alias());
    TArrayD expPS;
    expPS("c,a,b") = eager_block(A, "a,b,c", {0, 0, 1}, {3, 4, 2})("a,b,c");
    REQUIRE(equal_tarrays(ps->get<TArrayD>(), expPS, "c,a,b"));
    // the permuted slice as a product operand, relabelled once more
    auto ps2 = p->slice_mode(0, 5, 10);
    auto sB = rB->slice_mode(0, 5, 10);
    auto prod = ps2->prod(*sB, annot3("c,a,b", "c,d", "b,a,d"), DeNest::False);
    REQUIRE(ps2->is_buffer_alias());
    TArrayD exp;
    exp("b,a,d") = eager_block(A, "a,b,c", {0, 0, 1}, {3, 4, 2})("a,b,c") *
                   eager_block(B, "c,d", {1, 0}, {2, 2})("c,d");
    REQUIRE(equal_tarrays(prod->get<TArrayD>(), exp, "b,a,d"));
    // a permutation of a permutation
    auto pp = p->permute(annot2("c,a,b", "b,c,a"));
    REQUIRE(pp->is_buffer_alias());
    TArrayD expPP;
    expPP("b,c,a") = A("a,b,c");
    REQUIRE(equal_tarrays(pp->get<TArrayD>(), expPP, "b,c,a"));
    // the identity permutation is an alias too, and copies before a write:
    // accumulating into it leaves A alone
    auto const A0 = TA::clone(A);
    auto id = rA->permute(annot2("a,b,c", "a,b,c"));
    REQUIRE(id->is_buffer_alias());
    id->add_inplace(*eval_result<ResultD>(TA::clone(A)));
    REQUIRE_FALSE(id->is_buffer_alias());
    TArrayD twice;
    twice("a,b,c") = 2.0 * A0("a,b,c");
    REQUIRE(equal_tarrays(id->get<TArrayD>(), twice, "a,b,c"));
    REQUIRE(equal_tarrays(A, A0, "a,b,c"));
    REQUIRE(p->is_buffer_alias());
  }

  SECTION("the adjoint conjugates lazily") {
    TArrayZ Z(world, TA::TiledRange{{0, 3, 6}, {0, 5, 10}});
    TArrayZ W(world, TA::TiledRange{{0, 5, 10}, {0, 4}});
    Z.fill_random();
    W.fill_random();
    world.gop.fence();
    auto rZ = eval_result<ResultZ>(Z);
    auto rW = eval_result<ResultZ>(W);
    auto adj = rZ->adjoint(annot2("i,a", "a,i"));
    REQUIRE(adj->is_buffer_alias());
    TArrayZ expAdj;
    expAdj("a,i") = Z("i,a").conj();
    // consumed by a product: Z^dagger(a,i) W(a,j) -> (i,j), nothing
    // materialized
    auto prod = adj->prod(*rW, annot3("a,i", "a,j", "i,j"), DeNest::False);
    REQUIRE(adj->is_buffer_alias());
    TArrayZ exp;
    exp("i,j") = expAdj("a,i") * W("a,j");
    REQUIRE(equal_tarrays(prod->get<TArrayZ>(), exp, "i,j"));
    // a slice of the adjoint: mode 0 (= a) tile [5, 10)
    auto as = adj->slice_mode(0, 5, 10);
    REQUIRE(as->is_buffer_alias());
    TArrayZ expAS;
    expAS("a,i") = eager_block(Z, "i,a", {0, 1}, {2, 2})("i,a").conj();
    REQUIRE(equal_tarrays(as->get<TArrayZ>(), expAS, "a,i"));
    // the adjoint of the adjoint is a (trivial) view of the array
    auto adj2 = adj->adjoint(annot2("a,i", "i,a"));
    REQUIRE(adj2->is_buffer_alias());
    REQUIRE(equal_tarrays(adj2->get<TArrayZ>(), Z, "i,a"));
    REQUIRE_FALSE(adj2->is_buffer_alias());
    // and the value of the adjoint itself
    REQUIRE(equal_tarrays(adj->get<TArrayZ>(), expAdj, "a,i"));
  }

  SECTION("sum, scaling, dot and norm take views as is") {
    auto s1 = rA->slice_mode(2, 5, 10);
    auto s2 = rA->slice_mode(2, 5, 10);
    auto const blk = eager_block(A, "a,b,c", {0, 0, 1}, {3, 4, 2});
    auto sum = s1->sum(*s2, annot3("a,b,c", "a,b,c", "c,b,a"));
    REQUIRE(s1->is_buffer_alias());
    REQUIRE(s2->is_buffer_alias());
    TArrayD expSum;
    expSum("c,b,a") = blk("a,b,c") + blk("a,b,c");
    REQUIRE(equal_tarrays(sum->get<TArrayD>(), expSum, "c,b,a"));
    TArrayD unusedSum;
    unusedSum("a,b,c") = blk("a,b,c") + blk("a,b,c");
    expSum = unusedSum;

    auto scaled = s1->prod(*eval_result<ResultScalar<double>>(2.5),
                           annot3("a,b,c", "", "c,b,a"), DeNest::False);
    REQUIRE(s1->is_buffer_alias());
    TArrayD expScaled;
    expScaled("c,b,a") = 2.5 * blk("a,b,c");
    REQUIRE(equal_tarrays(scaled->get<TArrayD>(), expScaled, "c,b,a"));

    auto dot = s1->prod(*s2, annot3("a,b,c", "a,b,c", ""), DeNest::False);
    REQUIRE(s1->is_buffer_alias());
    double const expDot = TA::dot(blk("a,b,c"), blk("a,b,c"));
    REQUIRE(dot->get<double>() == Catch::Approx(expDot));

    REQUIRE(s1->norm2() == Catch::Approx(TA::norm2(blk)));
    REQUIRE(s1->is_buffer_alias());
  }

  SECTION("writes into a view never reach the parent") {
    auto const A0 = TA::clone(A);
    // accumulate into a slice: the slice becomes its own array
    auto s = rA->slice_mode(2, 5, 10);
    auto const blk = eager_block(A, "a,b,c", {0, 0, 1}, {3, 4, 2});
    auto other = eval_result<ResultD>(TA::clone(blk));
    s->add_inplace(*other);
    REQUIRE_FALSE(s->is_buffer_alias());
    TArrayD exp2;
    exp2("a,b,c") = 2.0 * blk("a,b,c");
    REQUIRE(equal_tarrays(s->get<TArrayD>(), exp2, "a,b,c"));
    REQUIRE(equal_tarrays(A, A0, "a,b,c"));
    // accumulating a view INTO a held value folds the view in
    auto held = eval_result<ResultD>(TA::clone(blk));
    held->add_inplace(*rA->slice_mode(2, 5, 10));
    REQUIRE(equal_tarrays(held->get<TArrayD>(), exp2, "a,b,c"));
    REQUIRE(equal_tarrays(A, A0, "a,b,c"));

    // the scatter: a view as the block, and a view as the destination
    TArrayD Zr(world, trA);
    Zr.fill(0.0);
    world.gop.fence();
    auto dest = eval_result<ResultD>(Zr);
    dest->write_into_slice(*rA->slice_mode(2, 5, 10), 2, 5, 10);
    TArrayD expZ;
    expZ("a,b,c") = Zr("a,b,c");
    expZ("a,b,c").block({0, 0, 1}, {3, 4, 2}, TA::preserve_lobound) =
        blk("a,b,c");
    REQUIRE(equal_tarrays(dest->get<TArrayD>(), expZ, "a,b,c"));
    auto destView = rA->slice_mode(1, 2, 6);  // a view as a destination
    destView->write_into_slice(*rA->slice_mode(1, 2, 6)->slice_mode(2, 5, 10),
                               2, 5, 10);
    REQUIRE_FALSE(destView->is_buffer_alias());
    REQUIRE(equal_tarrays(destView->get<TArrayD>(),
                          eager_block(A, "a,b,c", {0, 1, 0}, {3, 3, 2}),
                          "a,b,c"));
    REQUIRE(equal_tarrays(A, A0, "a,b,c"));
  }

  SECTION("clone of a view is a fresh array with the view applied") {
    auto s = rA->slice_mode(2, 5, 10);
    auto c = s->clone();
    REQUIRE(s->is_buffer_alias());
    REQUIRE_FALSE(c->is_buffer_alias());
    auto const blk = eager_block(A, "a,b,c", {0, 0, 1}, {3, 4, 2});
    REQUIRE(equal_tarrays(c->get<TArrayD>(), blk, "a,b,c"));
    // independent storage
    c->add_inplace(*eval_result<ResultD>(TA::clone(blk)));
    REQUIRE(equal_tarrays(s->get<TArrayD>(), blk, "a,b,c"));
  }

  SECTION("nested-array slices are views too") {
    using ToTArray = TA::DistArray<tot_tile_t<double>, TA::SparsePolicy>;
    using ResultToT = sequant::ResultTensorOfTensorTA<ToTArray>;
    // T(i,c;x,y) nested over tiles 2 x 2, inner 3 x 2
    TA::TiledRange const trT{{0, 2, 4}, {0, 5, 10}};
    TA::Range const inner({3, 2});
    ToTArray T(world, trT);
    for (auto it = T.begin(); it != T.end(); ++it)
      if (T.is_local(it.index()))
        *it = world.taskq.add(
            [inner](TA::Range const& orng) {
              return random_tensor_of_tensor<double>(orng, inner);
            },
            it.make_range());
    world.gop.fence();
    auto const T0 = TA::clone(T);
    auto rT = eval_result<ResultToT>(T);
    auto const eager = sequant::slice_array_over_mode(T, 1, 1, 2);

    // the slice is an alias, reports the slice's layout, materializes once
    auto sT = rT->slice_mode(1, 5, 10);
    REQUIRE(sT->is_buffer_alias());
    REQUIRE(sT->size_in_bytes() == 0);
    std::ostringstream os;
    os << eager.trange();
    REQUIRE(sT->trange_annot() == os.str());
    REQUIRE(equal_tarrays(sT->get<ToTArray>(), eager, "i,c;x,y"));
    REQUIRE_FALSE(sT->is_buffer_alias());

    // ToT * flat with both operands sliced on the contracted index: neither
    // is materialized; the batches add up to the full product
    auto sT1 = rT->slice_mode(1, 5, 10);
    auto sB1 = rB->slice_mode(0, 5, 10);
    auto p1 =
        sT1->prod(*sB1, annot3("i,c;x,y", "c,d", "i,d;x,y"), DeNest::False);
    REQUIRE(sT1->is_buffer_alias());
    REQUIRE(sB1->is_buffer_alias());
    ToTArray exp1;
    exp1("i,d;x,y") =
        eager("i,c;x,y") * eager_block(B, "c,d", {1, 0}, {2, 2})("c,d");
    REQUIRE(equal_tarrays(p1->get<ToTArray>(), exp1, "i,d;x,y"));
    auto p0 = rT->slice_mode(1, 0, 5)->prod(*rB->slice_mode(0, 0, 5),
                                            annot3("i,c;x,y", "c,d", "i,d;x,y"),
                                            DeNest::False);
    p0->add_inplace(*p1);
    ToTArray full;
    full("i,d;x,y") = T("i,c;x,y") * B("c,d");
    REQUIRE(equal_tarrays(p0->get<ToTArray>(), full, "i,d;x,y"));

    // a nested view as the operand of an accumulation, and as the block of a
    // scatter; neither reaches T
    auto acc = eval_result<ResultToT>(TA::clone(eager));
    acc->add_inplace(*rT->slice_mode(1, 5, 10));
    ToTArray twice;
    twice("i,c;x,y") = 2.0 * eager("i,c;x,y");
    REQUIRE(equal_tarrays(acc->get<ToTArray>(), twice, "i,c;x,y"));
    auto dest = eval_result<ResultToT>(TA::clone(T));
    dest->write_into_slice(*acc, 1, 5, 10);
    ToTArray expDest;
    expDest("i,c;x,y") = T("i,c;x,y");
    expDest("i,c;x,y").block({0, 1}, {2, 2}, TA::preserve_lobound) =
        twice("i,c;x,y");
    REQUIRE(equal_tarrays(dest->get<ToTArray>(), expDest, "i,c;x,y"));
    // scattering INTO a view materializes it; T is untouched throughout
    auto destView = rT->slice_mode(0, 0, 4);
    REQUIRE(destView->is_buffer_alias());
    destView->write_into_slice(*acc, 1, 5, 10);
    REQUIRE_FALSE(destView->is_buffer_alias());
    REQUIRE(equal_tarrays(destView->get<ToTArray>(), expDest, "i,c;x,y"));
    REQUIRE(equal_tarrays(T, T0, "i,c;x,y"));

    // a flat view (permuted and sliced) as the operand of a nested product
    auto pB = rB->permute(annot2("c,d", "d,c"))->slice_mode(1, 5, 10);
    auto p2 = rT->slice_mode(1, 5, 10)->prod(
        *pB, annot3("i,c;x,y", "d,c", "i,d;x,y"), DeNest::False);
    REQUIRE(pB->is_buffer_alias());
    REQUIRE(equal_tarrays(p2->get<ToTArray>(), exp1, "i,d;x,y"));
  }

  SECTION("a permutation or an adjoint of a slice is a view") {
    auto ps = rA->slice_mode(1, 2, 6)->permute(annot2("a,b,c", "c,a,b"));
    REQUIRE(ps->is_buffer_alias());
    TArrayD exp;
    exp("c,a,b") = eager_block(A, "a,b,c", {0, 1, 0}, {3, 3, 2})("a,b,c");
    REQUIRE(equal_tarrays(ps->get<TArrayD>(), exp, "c,a,b"));

    TArrayZ Z(world, TA::TiledRange{{0, 3, 6}, {0, 5, 10}});
    Z.fill_random();
    world.gop.fence();
    auto rZ = eval_result<ResultZ>(Z);
    auto as = rZ->slice_mode(1, 5, 10)->adjoint(annot2("i,a", "a,i"));
    REQUIRE(as->is_buffer_alias());
    TArrayZ expZ;
    expZ("a,i") = eager_block(Z, "i,a", {0, 1}, {2, 2})("i,a").conj();
    REQUIRE(equal_tarrays(as->get<TArrayZ>(), expZ, "a,i"));
  }

  SECTION("the array consumers of a view leave the parent alone") {
    auto const A0 = TA::clone(A);
    auto const blk = eager_block(A, "a,b,c", {0, 0, 1}, {3, 4, 2});
    // a phase: folded into the scaling, the view stays a view
    auto const vNeg = rA->slice_mode(2, 5, 10);
    auto neg = vNeg->mult_by_phase(-1);
    REQUIRE(vNeg->is_buffer_alias());
    TArrayD expNeg;
    expNeg("a,b,c") = -1.0 * blk("a,b,c");
    REQUIRE(equal_tarrays(neg->get<TArrayD>(), expNeg, "a,b,c"));
    REQUIRE(equal_tarrays(A, A0, "a,b,c"));
    // the symmetrizer (i <-> j with a <-> b): the same as on the
    // materialized slice; both members of a swapped pair are cut alike
    TArrayD M(world,
              TA::TiledRange{{0, 3, 6}, {0, 3, 6}, {0, 2, 4}, {0, 2, 4}});
    M.fill_random();
    world.gop.fence();
    auto rM = eval_result<ResultD>(M);
    auto symView = rM->slice_mode(0, 3, 6)->slice_mode(1, 3, 6)->symmetrize();
    auto symEager = eval_result<ResultD>(
                        eager_block(M, "i,j,a,b", {1, 1, 0, 0}, {2, 2, 2, 2}))
                        ->symmetrize();
    REQUIRE(equal_tarrays(symView->get<TArrayD>(), symEager->get<TArrayD>(),
                          "i,j,a,b"));
  }

  SECTION("a spectator mode with a nonzero element lobound keeps it") {
    // the modes of A2 and B2 start at element 4 (a frozen-core offset)
    TArrayD A2(world, TA::TiledRange{{4, 6, 8}, {0, 5, 10}});
    TArrayD B2(world, TA::TiledRange{{4, 6, 8}, {0, 4}});
    A2.fill_random();
    B2.fill_random();
    world.gop.fence();
    auto rA2 = eval_result<ResultD>(A2);
    auto rB2 = eval_result<ResultD>(B2);
    auto s = rA2->slice_mode(1, 5, 10);
    auto const blk = eager_block(A2, "a,b", {0, 1}, {2, 2});
    auto const& sv = s->as<ResultD>();
    REQUIRE(sv.logical_trange() == blk.trange());
    REQUIRE(sv.logical_trange().elements_range().lobound()[0] == 4);
    auto prod = s->prod(*rB2, annot3("a,b", "a,c", "b,c"), DeNest::False);
    REQUIRE(s->is_buffer_alias());
    TArrayD exp;
    exp("b,c") = blk("a,b") * B2("a,c");
    REQUIRE(equal_tarrays(prod->get<TArrayD>(), exp, "b,c"));
    // and through a permutation of the view
    auto pp = s->permute(annot2("a,b", "b,a"));
    REQUIRE(pp->as<ResultD>().logical_trange().elements_range().lobound()[1] ==
            4);
    TArrayD expP;
    expP("b,a") = blk("a,b");
    REQUIRE(equal_tarrays(pp->get<TArrayD>(), expP, "b,a"));
  }

  SECTION("scaling by a constant is a pending view") {
    auto const A0 = TA::clone(A);
    auto c = eval_result<ResultScalar<double>>(2.5);
    auto scaled = rA->prod(*c, annot3("a,b,c", "", "a,b,c"), DeNest::False);
    REQUIRE(scaled->is_buffer_alias());
    REQUIRE(scaled->measured_size_in_bytes() == 0);
    REQUIRE(scaled->as<ResultD>().logical_trange() == trA);
    // the scalar on the left (ResultScalar::prod forwards to the tensor)
    auto scaledL = c->prod(*rA, annot3("", "a,b,c", "a,b,c"), DeNest::False);
    REQUIRE(scaledL->is_buffer_alias());
    // the phase is a view as well, and the factors compose
    auto neg = scaled->mult_by_phase(-1);
    REQUIRE(neg->is_buffer_alias());
    REQUIRE(scaled->is_buffer_alias());
    TArrayD exp;
    exp("a,b,c") = 2.5 * A("a,b,c");
    REQUIRE(equal_tarrays(scaled->get<TArrayD>(), exp, "a,b,c"));
    REQUIRE_FALSE(scaled->is_buffer_alias());
    REQUIRE(equal_tarrays(scaledL->get<TArrayD>(), exp, "a,b,c"));
    TArrayD expNeg;
    expNeg("a,b,c") = -2.5 * A("a,b,c");
    REQUIRE(equal_tarrays(neg->get<TArrayD>(), expNeg, "a,b,c"));
    // the product annotation relabels the result: a permuted scaled view
    auto sp = rA->prod(*c, annot3("a,b,c", "", "c,a,b"), DeNest::False);
    REQUIRE(sp->is_buffer_alias());
    TArrayD expP;
    expP("c,a,b") = 2.5 * A("a,b,c");
    REQUIRE(equal_tarrays(sp->get<TArrayD>(), expP, "c,a,b"));
    // the factor survives a slice and a permutation of the scaled view
    auto const blk = eager_block(A, "a,b,c", {0, 0, 1}, {3, 4, 2});
    auto ss = rA->prod(*c, annot3("a,b,c", "", "a,b,c"), DeNest::False)
                  ->slice_mode(2, 5, 10)
                  ->permute(annot2("a,b,c", "c,a,b"));
    REQUIRE(ss->is_buffer_alias());
    TArrayD expSS;
    expSS("c,a,b") = 2.5 * blk("a,b,c");
    REQUIRE(equal_tarrays(ss->get<TArrayD>(), expSS, "c,a,b"));
    // and a scaling of a slice composes with its phase
    auto sc = rA->slice_mode(2, 5, 10)->mult_by_phase(-1)->prod(
        *c, annot3("a,b,c", "", "a,b,c"), DeNest::False);
    REQUIRE(sc->is_buffer_alias());
    TArrayD expSC;
    expSC("a,b,c") = -2.5 * blk("a,b,c");
    REQUIRE(equal_tarrays(sc->get<TArrayD>(), expSC, "a,b,c"));
    // a scaled view of a held value is a view of that value
    auto held = eval_result<ResultD>(TA::clone(blk));
    auto hs = held->prod(*c, annot3("a,b,c", "", "a,b,c"), DeNest::False);
    REQUIRE(hs->is_buffer_alias());
    REQUIRE_FALSE(held->is_buffer_alias());
    TArrayD expHS;
    expHS("a,b,c") = 2.5 * blk("a,b,c");
    REQUIRE(equal_tarrays(hs->get<TArrayD>(), expHS, "a,b,c"));
    REQUIRE(equal_tarrays(held->get<TArrayD>(), blk, "a,b,c"));
    // the parent is untouched throughout
    REQUIRE(equal_tarrays(A, A0, "a,b,c"));
    // a view that is the only holder of its parent's tiles owns their bytes:
    // a phase applied to a value nobody else kept
    auto const bytes =
        eval_result<ResultD>(TA::clone(A))->measured_size_in_bytes();
    REQUIRE(bytes > 0);
    auto sole = eval_result<ResultD>(TA::clone(A))->mult_by_phase(-1);
    REQUIRE(sole->is_buffer_alias());
    REQUIRE(sole->measured_size_in_bytes() == bytes);
    // not while the parent is held elsewhere
    REQUIRE(rA->mult_by_phase(-1)->measured_size_in_bytes() == 0);
  }

  SECTION("the consumers of a scaled view fold the factor in") {
    auto const A0 = TA::clone(A);
    auto c2 = eval_result<ResultScalar<double>>(2.0);
    auto c3 = eval_result<ResultScalar<double>>(-3.0);
    auto scale = [&](ResultPtr const& r, std::string const& ann,
                     ResultPtr const& c) {
      return r->prod(*c, annot3(ann, "", ann), DeNest::False);
    };
    auto const blkA = eager_block(A, "a,b,c", {0, 0, 1}, {3, 4, 2});
    auto const blkB = eager_block(B, "c,d", {1, 0}, {2, 2});
    auto sA = scale(rA->slice_mode(2, 5, 10), "a,b,c", c2);
    auto sB = scale(rB->slice_mode(0, 5, 10), "c,d", c3);
    // a contraction, both operands scaled
    auto prod = sA->prod(*sB, annot3("a,b,c", "c,d", "a,b,d"), DeNest::False);
    REQUIRE(sA->is_buffer_alias());
    REQUIRE(sB->is_buffer_alias());
    REQUIRE_FALSE(prod->is_buffer_alias());
    TArrayD exp;
    exp("a,b,d") = -6.0 * (blkA("a,b,c") * blkB("c,d"));
    REQUIRE(equal_tarrays(prod->get<TArrayD>(), exp, "a,b,d"));
    // a product with a Hadamard index
    TArrayD H(world, TA::TiledRange{{0, 3, 6, 9}, {0, 5, 10}});
    H.fill_random();
    world.gop.fence();
    auto sH = scale(eval_result<ResultD>(H)->slice_mode(1, 5, 10), "a,c", c3);
    auto hprod = sA->prod(*sH, annot3("a,b,c", "a,c", "a,b"), DeNest::False);
    REQUIRE(sH->is_buffer_alias());
    TArrayD expH;
    expH("a,b") =
        -6.0 * TA::einsum(blkA("a,b,c"),
                          eager_block(H, "a,c", {0, 1}, {3, 2})("a,c"),
                          "a,b")("a,b");
    REQUIRE(equal_tarrays(hprod->get<TArrayD>(), expH, "a,b"));
    // a dot product
    auto sA3 = scale(rA->slice_mode(2, 5, 10), "a,b,c", c3);
    auto dot = sA->prod(*sA3, annot3("a,b,c", "a,b,c", ""), DeNest::False);
    REQUIRE(sA->is_buffer_alias());
    REQUIRE(sA3->is_buffer_alias());
    REQUIRE(dot->get<double>() ==
            Catch::Approx(-6.0 * TA::dot(blkA("a,b,c"), blkA("a,b,c"))));
    // a sum, relabelled
    auto sum = sA->sum(*sA3, annot3("a,b,c", "a,b,c", "c,b,a"));
    REQUIRE(sA->is_buffer_alias());
    TArrayD expSum;
    expSum("c,b,a") = -1.0 * blkA("a,b,c");
    REQUIRE(equal_tarrays(sum->get<TArrayD>(), expSum, "c,b,a"));
    // an accumulation into a held value
    auto held = eval_result<ResultD>(TA::clone(blkA));
    held->add_inplace(*sA);
    REQUIRE(sA->is_buffer_alias());
    TArrayD exp3;
    exp3("a,b,c") = 3.0 * blkA("a,b,c");
    REQUIRE(equal_tarrays(held->get<TArrayD>(), exp3, "a,b,c"));
    // the norm
    REQUIRE(sA->norm2() == Catch::Approx(2.0 * TA::norm2(blkA)));
    REQUIRE(sA->is_buffer_alias());
    // the scatter of a scaled view
    TArrayD Zr(world, trA);
    Zr.fill(0.0);
    world.gop.fence();
    auto dest = eval_result<ResultD>(TA::clone(Zr));
    dest->write_into_slice(*sA, 2, 5, 10);
    REQUIRE(sA->is_buffer_alias());
    TArrayD expZ;
    expZ("a,b,c") = Zr("a,b,c");
    expZ("a,b,c").block({0, 0, 1}, {3, 4, 2}, TA::preserve_lobound) =
        2.0 * blkA("a,b,c");
    REQUIRE(equal_tarrays(dest->get<TArrayD>(), expZ, "a,b,c"));
    // the clone
    auto cl = sA->clone();
    REQUIRE(sA->is_buffer_alias());
    REQUIRE_FALSE(cl->is_buffer_alias());
    TArrayD expCl;
    expCl("a,b,c") = 2.0 * blkA("a,b,c");
    REQUIRE(equal_tarrays(cl->get<TArrayD>(), expCl, "a,b,c"));
    // writes into a scaled view never reach the parent
    auto sv = scale(rA, "a,b,c", c2);
    sv->add_inplace(*eval_result<ResultD>(TA::clone(A)));
    REQUIRE_FALSE(sv->is_buffer_alias());
    TArrayD expSV;
    expSV("a,b,c") = 3.0 * A("a,b,c");
    REQUIRE(equal_tarrays(sv->get<TArrayD>(), expSV, "a,b,c"));
    REQUIRE(equal_tarrays(A, A0, "a,b,c"));
  }

  SECTION("the adjoint of a scaled complex view conjugates the factor") {
    using Cplx = std::complex<double>;
    TArrayZ Z(world, TA::TiledRange{{0, 3, 6}, {0, 5, 10}});
    TArrayZ W(world, TA::TiledRange{{0, 5, 10}, {0, 4}});
    Z.fill_random();
    W.fill_random();
    world.gop.fence();
    auto rZ = eval_result<ResultZ>(Z);
    auto rW = eval_result<ResultZ>(W);
    Cplx const c{1.0, 2.0};
    auto rc = eval_result<ResultScalar<Cplx>>(c);
    auto cz = rZ->prod(*rc, annot3("i,a", "", "i,a"), DeNest::False);
    REQUIRE(cz->is_buffer_alias());
    REQUIRE(cz->norm2() == Catch::Approx(std::abs(c) * TA::norm2(Z)));
    // (c Z)^dagger = conj(c) Z^dagger
    auto adj = cz->adjoint(annot2("i,a", "a,i"));
    REQUIRE(adj->is_buffer_alias());
    REQUIRE(cz->is_buffer_alias());
    TArrayZ expAdj;
    expAdj("a,i") = std::conj(c) * Z("i,a").conj();
    // consumed by a product: nothing materialized
    auto prod = adj->prod(*rW, annot3("a,i", "a,j", "i,j"), DeNest::False);
    REQUIRE(adj->is_buffer_alias());
    TArrayZ expP;
    expP("i,j") = expAdj("a,i") * W("a,j");
    REQUIRE(equal_tarrays(prod->get<TArrayZ>(), expP, "i,j"));
    REQUIRE(equal_tarrays(adj->get<TArrayZ>(), expAdj, "a,i"));
    // scaling after the adjoint: c Z^dagger, the factor as given
    auto adjc = rZ->adjoint(annot2("i,a", "a,i"))
                    ->prod(*rc, annot3("a,i", "", "a,i"), DeNest::False);
    REQUIRE(adjc->is_buffer_alias());
    TArrayZ expAdjC;
    expAdjC("a,i") = c * Z("i,a").conj();
    REQUIRE(equal_tarrays(adjc->get<TArrayZ>(), expAdjC, "a,i"));
    // the adjoint of the adjoint of a scaled view is the scaled view
    auto adj2 =
        cz->adjoint(annot2("i,a", "a,i"))->adjoint(annot2("a,i", "i,a"));
    REQUIRE(adj2->is_buffer_alias());
    TArrayZ expZ;
    expZ("i,a") = c * Z("i,a");
    REQUIRE(equal_tarrays(adj2->get<TArrayZ>(), expZ, "i,a"));
    // a complex phase folds into the scaled conjugated view
    auto negAdj = cz->adjoint(annot2("i,a", "a,i"))->mult_by_phase(-1);
    REQUIRE(negAdj->is_buffer_alias());
    TArrayZ expNeg;
    expNeg("a,i") = -std::conj(c) * Z("i,a").conj();
    REQUIRE(equal_tarrays(negAdj->get<TArrayZ>(), expNeg, "a,i"));
    // the materialized scaled view itself
    REQUIRE(equal_tarrays(cz->get<TArrayZ>(), expZ, "i,a"));
    REQUIRE_FALSE(cz->is_buffer_alias());
  }

  SECTION("block-sparse arrays: the view keeps the shape") {
    using SpArray = TA::TSpArrayD;
    using ResultS = ResultTensorTA<SpArray>;
    SpArray SA = TA::to_sparse(A);
    SpArray SB = TA::to_sparse(B);
    auto rSA = eval_result<ResultS>(SA);
    auto rSB = eval_result<ResultS>(SB);
    auto sA = rSA->slice_mode(2, 5, 10);
    auto sB = rSB->slice_mode(0, 5, 10);
    REQUIRE(sA->is_buffer_alias());
    auto prod = sA->prod(*sB, annot3("a,b,c", "c,d", "a,b,d"), DeNest::False);
    REQUIRE(sA->is_buffer_alias());
    SpArray exp;
    exp("a,b,d") = eager_block(SA, "a,b,c", {0, 0, 1}, {3, 4, 2})("a,b,c") *
                   eager_block(SB, "c,d", {1, 0}, {2, 2})("c,d");
    REQUIRE(equal_tarrays(prod->get<SpArray>(), exp, "a,b,d"));
    // the view's value is the sparse block
    auto const blk = eager_block(SA, "a,b,c", {0, 0, 1}, {3, 4, 2});
    REQUIRE(equal_tarrays(sA->get<SpArray>(), blk, "a,b,c"));
    // scatter of a view into a sparse destination
    auto dest = eval_result<ResultS>(TA::to_sparse(A));
    SpArray twice;
    twice("a,b,c") = 2.0 * SA("a,b,c");
    dest->write_into_slice(*eval_result<ResultS>(twice)->slice_mode(2, 5, 10),
                           2, 5, 10);
    SpArray expD;
    expD("a,b,c") = SA("a,b,c");
    expD("a,b,c").block({0, 0, 1}, {3, 4, 2}, TA::preserve_lobound) =
        eager_block(twice, "a,b,c", {0, 0, 1}, {3, 4, 2})("a,b,c");
    REQUIRE(equal_tarrays(dest->get<SpArray>(), expD, "a,b,c"));
  }
}

// The nested-array result's views: a block over the outer modes, consumed by
// every nested product (nested x flat, nested x nested, the inner contraction
// that flattens, the dot), by sums, scaling, accumulation and the scatter.
TEST_CASE("ta_lazy_views_nested", "[eval][ta][lazy]") {
  using sequant::DeNest;
  using sequant::eval_result;
  using sequant::ResultScalar;
  using sequant::ResultTensorOfTensorTA;
  using ToTArray = TA::DistArray<tot_tile_t<double>, TA::SparsePolicy>;
  using ResultToT = ResultTensorOfTensorTA<ToTArray>;

  auto& world = TA::get_default_world();
  auto annot3 = [](std::string l, std::string r, std::string c) {
    return std::array<std::any, 3>{std::move(l), std::move(r), std::move(c)};
  };
  auto annot2 = [](std::string pre, std::string post) {
    return std::array<std::any, 2>{std::move(pre), std::move(post)};
  };
  // outer tiles 2 x 2, the same for every array
  TA::TiledRange const tr{{0, 2, 4}, {0, 5, 10}};
  auto make_tot = [&](TA::Range const& inner, TA::TiledRange const& outer) {
    ToTArray t(world, outer);
    for (auto it = t.begin(); it != t.end(); ++it)
      if (t.is_local(it.index()))
        *it = world.taskq.add(
            [inner](TA::Range const& orng) {
              return random_tensor_of_tensor<double>(orng, inner);
            },
            it.make_range());
    world.gop.fence();
    return t;
  };
  // L(i,j;a,b) inner 3 x 2, R(i,j;b,c) inner 2 x 4
  auto const L = make_tot(TA::Range({3, 2}), tr);
  auto const R = make_tot(TA::Range({2, 4}), tr);
  auto const L0 = TA::clone(L);
  auto rL = eval_result<ResultToT>(L);
  auto rR = eval_result<ResultToT>(R);
  // the eager slices over j = elements [5, 10) (tile 1)
  auto const Le = sequant::slice_array_over_mode(L, 1, 1, 2);
  auto const Re = sequant::slice_array_over_mode(R, 1, 1, 2);
  auto sL = [&] { return rL->slice_mode(1, 5, 10); };
  auto sR = [&] { return rR->slice_mode(1, 5, 10); };

  SECTION("the consumer of a nested view reads the parent's tiles in place") {
    {
      auto const [tr, tiles] =
          view_consumer_tiles(sL()->as<ResultToT>(), "i,j;a,b");
      REQUIRE(consumer_reads_parent_tiles(L, tr, tiles, {0, 1}));
    }
    {
      auto ss = sL()->slice_mode(0, 2, 4);
      auto const [tr, tiles] =
          view_consumer_tiles(ss->as<ResultToT>(), "i,j;a,b");
      REQUIRE(consumer_reads_parent_tiles(L, tr, tiles, {1, 1}));
    }
  }

  SECTION("a slice of a slice, and the clone") {
    auto ss = sL()->slice_mode(0, 2, 4);
    REQUIRE(ss->is_buffer_alias());
    auto const exp = sequant::slice_array_over_mode(Le, 0, 1, 2);
    REQUIRE(equal_tarrays(ss->get<ToTArray>(), exp, "i,j;a,b"));
    auto v = sL();
    auto c = v->clone();
    REQUIRE(v->is_buffer_alias());
    REQUIRE_FALSE(c->is_buffer_alias());
    REQUIRE(equal_tarrays(c->get<ToTArray>(), Le, "i,j;a,b"));
    c->add_inplace(*sL());
    REQUIRE(equal_tarrays(v->get<ToTArray>(), Le, "i,j;a,b"));
    REQUIRE(equal_tarrays(L, L0, "i,j;a,b"));
  }

  SECTION("nested x nested products take the views as is") {
    // outer Hadamard, inner contraction over b -> (i,j;a,c)
    auto l = sL();
    auto r = sR();
    auto p =
        l->prod(*r, annot3("i,j;a,b", "i,j;b,c", "i,j;a,c"), DeNest::False);
    REQUIRE(l->is_buffer_alias());
    REQUIRE(r->is_buffer_alias());
    auto const exp = TA::einsum(Le("i,j;a,b"), Re("i,j;b,c"), "i,j;a,c");
    REQUIRE(equal_tarrays(p->get<ToTArray>(), exp, "i,j;a,c"));
    // the result wants its inner modes in the other order
    auto p2 = sL()->prod(*sR(), annot3("i,j;a,b", "i,j;b,c", "i,j;c,a"),
                         DeNest::False);
    auto const exp2 = TA::einsum(Le("i,j;a,b"), Re("i,j;b,c"), "i,j;c,a");
    REQUIRE(equal_tarrays(p2->get<ToTArray>(), exp2, "i,j;c,a"));
    // the inner contraction that flattens: (i,j;a,b) . (i,j;a,b) -> (i,j)
    auto l3 = sL();
    auto d = l3->prod(*sL(), annot3("i,j;a,b", "i,j;a,b", "i,j"), DeNest::True);
    REQUIRE(l3->is_buffer_alias());
    auto const expD =
        TA::einsum<TA::DeNest::True>(Le("i,j;a,b"), Le("i,j;a,b"), "i,j");
    REQUIRE(equal_tarrays(d->get<TA::TSpArrayD>(), expD, "i,j"));
    // the dot
    auto l4 = sL();
    auto dot = l4->prod(*sL(), annot3("i,j;a,b", "i,j;a,b", ""), DeNest::False);
    REQUIRE(l4->is_buffer_alias());
    REQUIRE(dot->get<double>() ==
            Catch::Approx(TA::dot(Le("i,j;a,b"), Le("i,j;a,b"))));
  }

  SECTION("sum, scaling, permutation and adjoint of nested views") {
    auto l = sL();
    auto sum = l->sum(*sL(), annot3("i,j;a,b", "i,j;a,b", "i,j;a,b"));
    REQUIRE(l->is_buffer_alias());
    ToTArray twice;
    twice("i,j;a,b") = 2.0 * Le("i,j;a,b");
    REQUIRE(equal_tarrays(sum->get<ToTArray>(), twice, "i,j;a,b"));

    auto scaled = l->prod(*eval_result<ResultScalar<double>>(2.0),
                          annot3("i,j;a,b", "", "j,i;a,b"), DeNest::False);
    REQUIRE(l->is_buffer_alias());
    ToTArray expS;
    expS("j,i;a,b") = twice("i,j;a,b");
    REQUIRE(equal_tarrays(scaled->get<ToTArray>(), expS, "j,i;a,b"));

    auto perm = l->permute(annot2("i,j;a,b", "j,i;b,a"));
    REQUIRE(l->is_buffer_alias());
    ToTArray expP;
    expP("j,i;b,a") = Le("i,j;a,b");
    REQUIRE(equal_tarrays(perm->get<ToTArray>(), expP, "j,i;b,a"));

    auto adj = l->adjoint(annot2("i,j;a,b", "j,i;b,a"));
    REQUIRE(equal_tarrays(adj->get<ToTArray>(), expP, "j,i;b,a"));

    auto const vNeg = sL();
    auto neg = vNeg->mult_by_phase(-1);
    REQUIRE(vNeg->is_buffer_alias());
    ToTArray expN;
    expN("i,j;a,b") = -1.0 * Le("i,j;a,b");
    REQUIRE(equal_tarrays(neg->get<ToTArray>(), expN, "i,j;a,b"));
    REQUIRE(equal_tarrays(L, L0, "i,j;a,b"));

    // the symmetrizer (i <-> j with a <-> b): as on the eager slice of a
    // square array, both outer modes cut alike
    auto const S =
        make_tot(TA::Range({2, 2}), TA::TiledRange{{0, 2, 4}, {0, 2, 4}});
    auto rS = eval_result<ResultToT>(S);
    auto symV = rS->slice_mode(0, 2, 4)->slice_mode(1, 2, 4)->symmetrize();
    auto const Se = sequant::slice_array_over_mode(
        sequant::slice_array_over_mode(S, 0, 1, 2), 1, 1, 2);
    auto symE = eval_result<ResultToT>(TA::clone(Se))->symmetrize();
    REQUIRE(
        equal_tarrays(symV->get<ToTArray>(), symE->get<ToTArray>(), "i,j;a,b"));
  }

  SECTION("scaling a nested view is a pending view") {
    auto c = eval_result<ResultScalar<double>>(2.0);
    auto scale = [&](sequant::ResultPtr const& r, std::string const& ann,
                     sequant::ResultPtr const& c) {
      return r->prod(*c, annot3(ann, "", ann), DeNest::False);
    };
    auto l = sL();
    auto scaled = scale(l, "i,j;a,b", c);
    REQUIRE(scaled->is_buffer_alias());
    REQUIRE(scaled->measured_size_in_bytes() == 0);
    REQUIRE(l->is_buffer_alias());
    auto neg = scaled->mult_by_phase(-1);
    REQUIRE(neg->is_buffer_alias());
    ToTArray twice;
    twice("i,j;a,b") = 2.0 * Le("i,j;a,b");
    REQUIRE(scaled->norm2() ==
            Catch::Approx(std::sqrt(TA::squared_norm(twice).get())));
    REQUIRE(scaled->is_buffer_alias());
    REQUIRE(equal_tarrays(scaled->get<ToTArray>(), twice, "i,j;a,b"));
    REQUIRE_FALSE(scaled->is_buffer_alias());
    ToTArray expN;
    expN("i,j;a,b") = -2.0 * Le("i,j;a,b");
    REQUIRE(equal_tarrays(neg->get<ToTArray>(), expN, "i,j;a,b"));
    // a slice of a scaled whole-array view keeps the factor
    auto ws = scale(rL, "i,j;a,b", c)->slice_mode(1, 5, 10);
    REQUIRE(ws->is_buffer_alias());
    REQUIRE(equal_tarrays(ws->get<ToTArray>(), twice, "i,j;a,b"));
    // a scaling that relabels is applied, as the permutation is
    auto sp = scale(sL(), "i,j;a,b", c);
    auto spp = sL()->prod(*c, annot3("i,j;a,b", "", "j,i;a,b"), DeNest::False);
    ToTArray expPP;
    expPP("j,i;a,b") = twice("i,j;a,b");
    REQUIRE(equal_tarrays(spp->get<ToTArray>(), expPP, "j,i;a,b"));
    // the consumers fold the factor in: nested x nested (outer Hadamard,
    // inner contraction), the flattening inner contraction, the dot, the
    // sum, the accumulation, the scatter
    auto r = scale(sR(), "i,j;b,c", eval_result<ResultScalar<double>>(-3.0));
    auto p =
        sp->prod(*r, annot3("i,j;a,b", "i,j;b,c", "i,j;a,c"), DeNest::False);
    REQUIRE(sp->is_buffer_alias());
    REQUIRE(r->is_buffer_alias());
    ToTArray expP;
    expP("i,j;a,c") =
        -6.0 * TA::einsum(Le("i,j;a,b"), Re("i,j;b,c"), "i,j;a,c")("i,j;a,c");
    REQUIRE(equal_tarrays(p->get<ToTArray>(), expP, "i,j;a,c"));
    auto d = sp->prod(*scale(sL(), "i,j;a,b", c),
                      annot3("i,j;a,b", "i,j;a,b", "i,j"), DeNest::True);
    REQUIRE(sp->is_buffer_alias());
    TA::TSpArrayD expD;
    expD("i,j") = 4.0 * TA::einsum<TA::DeNest::True>(
                            Le("i,j;a,b"), Le("i,j;a,b"), "i,j")("i,j");
    REQUIRE(equal_tarrays(d->get<TA::TSpArrayD>(), expD, "i,j"));
    auto dot = sp->prod(*scale(sL(), "i,j;a,b", c),
                        annot3("i,j;a,b", "i,j;a,b", ""), DeNest::False);
    REQUIRE(dot->get<double>() ==
            Catch::Approx(4.0 * TA::dot(Le("i,j;a,b"), Le("i,j;a,b"))));
    auto sum = sp->sum(*scale(sL(), "i,j;a,b", c),
                       annot3("i,j;a,b", "i,j;a,b", "i,j;a,b"));
    ToTArray four;
    four("i,j;a,b") = 4.0 * Le("i,j;a,b");
    REQUIRE(equal_tarrays(sum->get<ToTArray>(), four, "i,j;a,b"));
    auto acc = eval_result<ResultToT>(TA::clone(Le));
    acc->add_inplace(*sp);
    REQUIRE(sp->is_buffer_alias());
    ToTArray three;
    three("i,j;a,b") = 3.0 * Le("i,j;a,b");
    REQUIRE(equal_tarrays(acc->get<ToTArray>(), three, "i,j;a,b"));
    auto dest = eval_result<ResultToT>(TA::clone(L));
    dest->write_into_slice(*sp, 1, 5, 10);
    REQUIRE(sp->is_buffer_alias());
    ToTArray expDest;
    expDest("i,j;a,b") = L("i,j;a,b");
    expDest("i,j;a,b").block({0, 1}, {2, 2}, TA::preserve_lobound) =
        twice("i,j;a,b");
    REQUIRE(equal_tarrays(dest->get<ToTArray>(), expDest, "i,j;a,b"));
    // the permutation and the adjoint apply the factor
    auto perm = sp->permute(annot2("i,j;a,b", "j,i;b,a"));
    auto adj = sp->adjoint(annot2("i,j;a,b", "j,i;b,a"));
    REQUIRE(sp->is_buffer_alias());
    ToTArray expPerm;
    expPerm("j,i;b,a") = twice("i,j;a,b");
    REQUIRE(equal_tarrays(perm->get<ToTArray>(), expPerm, "j,i;b,a"));
    REQUIRE(equal_tarrays(adj->get<ToTArray>(), expPerm, "j,i;b,a"));
    // writes into a scaled view never reach the parent
    auto sv = scale(rL, "i,j;a,b", c);
    sv->add_inplace(*eval_result<ResultToT>(TA::clone(L)));
    REQUIRE_FALSE(sv->is_buffer_alias());
    ToTArray expSV;
    expSV("i,j;a,b") = 3.0 * L("i,j;a,b");
    REQUIRE(equal_tarrays(sv->get<ToTArray>(), expSV, "i,j;a,b"));
    REQUIRE(equal_tarrays(L, L0, "i,j;a,b"));
    // a view that is the only holder of its parent's tiles owns their bytes,
    // and not while the parent is held elsewhere
    auto const bytes =
        eval_result<ResultToT>(TA::clone(L))->measured_size_in_bytes();
    REQUIRE(bytes > 0);
    auto sole = eval_result<ResultToT>(TA::clone(L))->mult_by_phase(-1);
    REQUIRE(sole->is_buffer_alias());
    REQUIRE(sole->measured_size_in_bytes() == bytes);
    REQUIRE(rL->mult_by_phase(-1)->measured_size_in_bytes() == 0);
  }

  SECTION("the norm and the tile diff of a nested view read it in place") {
    auto const v = sL();
    REQUIRE(v->norm2() == Catch::Approx(std::sqrt(TA::squared_norm(Le).get())));
    // a view of a slice of the same parent: the same tiles
    auto const v2 = sL();
    REQUIRE(v->tile_diff(*v2).find("only_this=0 only_other=0") !=
            std::string::npos);
    REQUIRE(v->is_buffer_alias());
    REQUIRE(v2->is_buffer_alias());
  }

  SECTION("a nested view as the block of a scatter") {
    ToTArray twiceFull;
    twiceFull("i,j;a,b") = 2.0 * L("i,j;a,b");
    auto dest = eval_result<ResultToT>(TA::clone(L));
    auto blockView = eval_result<ResultToT>(twiceFull)->slice_mode(1, 5, 10);
    dest->write_into_slice(*blockView, 1, 5, 10);
    REQUIRE(blockView->is_buffer_alias());
    ToTArray exp;
    exp("i,j;a,b") = L("i,j;a,b");
    exp("i,j;a,b").block({0, 1}, {2, 2}, TA::preserve_lobound) =
        sequant::slice_array_over_mode(twiceFull, 1, 1, 2)("i,j;a,b");
    REQUIRE(equal_tarrays(dest->get<ToTArray>(), exp, "i,j;a,b"));
    REQUIRE(equal_tarrays(L, L0, "i,j;a,b"));
  }

  SECTION("a view of an all-empty nested array is zero over the slice") {
    ToTArray E(world, tr);
    for (auto it = E.begin(); it != E.end(); ++it)
      if (E.is_local(it.index())) *it = tot_tile_t<double>(it.make_range());
    world.gop.fence();
    auto rE = eval_result<ResultToT>(E);
    auto v = rE->slice_mode(1, 5, 10);
    REQUIRE(v->is_buffer_alias());
    REQUIRE(v->trange_annot() == sL()->trange_annot());
    auto const& val = v->get<ToTArray>();
    REQUIRE(val.trange() == Le.trange());
    REQUIRE(sequant::detail::tot_inner_rank(val) == 0);
  }
}

// The result-shape-constrained product (apply_shaped_product, the product the
// evaluator emits when a result-shape provider is set) takes views as is, for
// every operand kind it supports.
TEST_CASE("ta_lazy_views_shaped", "[eval][ta][lazy]") {
  using sequant::eval_result;
  using sequant::ResultScalar;
  using sequant::ResultTensorOfTensorTA;
  using sequant::ResultTensorTA;
  using FlatArray = TA::DistArray<TA::Tensor<double>, TA::SparsePolicy>;
  using ToTArray = TA::DistArray<tot_tile_t<double>, TA::SparsePolicy>;
  using ResultF = ResultTensorTA<FlatArray>;
  using ResultToT = ResultTensorOfTensorTA<ToTArray>;

  auto& world = TA::get_default_world();
  auto annot3 = [](std::string l, std::string r, std::string c) {
    return std::array<std::any, 3>{std::move(l), std::move(r), std::move(c)};
  };
  auto full_shape = [](TA::TiledRange const& tr) {
    return TA::SparseShape<float>(TA::Tensor<float>(tr.tiles_range(), 1.0f),
                                  tr);
  };
  auto make_tot = [&](TA::Range const& inner, TA::TiledRange const& outer) {
    ToTArray t(world, outer, full_shape(outer));
    for (auto it = t.begin(); it != t.end(); ++it)
      if (t.is_local(it.index()))
        *it = world.taskq.add(
            [inner](TA::Range const& orng) {
              return random_tensor_of_tensor<double>(orng, inner);
            },
            it.make_range());
    world.gop.fence();
    return t;
  };
  auto shaped = [&](sequant::Result const& l, sequant::Result const& r,
                    std::array<std::any, 3> const& annot, bool de_nest) {
    auto const tr =
        sequant::result_outer_trange_from_results<double, TA::SparsePolicy>(
            l, r, annot);
    // a shape that drops the result's first outer tile
    TA::Tensor<float> norms{tr.tiles_range(), 1.0f};
    norms[0] = 0.0f;
    return sequant::apply_shaped_product<double, TA::SparsePolicy>(
        l, r, annot, TA::SparseShape<float>(norms, tr), de_nest);
  };

  // the flat operands
  TA::TSpArrayD Ad(world,
                   TA::TiledRange{{0, 3, 6, 9}, {0, 2, 4, 6, 8}, {0, 5, 10}});
  TA::TSpArrayD Bd(world, TA::TiledRange{{0, 5, 10}, {0, 4, 7}});
  Ad.fill_random();
  Bd.fill_random();
  world.gop.fence();
  FlatArray const A = TA::to_sparse(Ad);
  FlatArray const B = TA::to_sparse(Bd);
  auto rA = eval_result<ResultF>(A);
  auto rB = eval_result<ResultF>(B);
  auto const Ae = sequant::slice_array_over_mode(A, 2, 1, 2);
  auto const Be = sequant::slice_array_over_mode(B, 0, 1, 2);
  // the nested operands: T(i,c;x,y) over c, L(i,j;a,b) and R(i,j;b,c) over j
  TA::TiledRange const trT{{0, 2, 4}, {0, 5, 10}};
  auto const T = make_tot(TA::Range({3, 2}), trT);
  auto const L = make_tot(TA::Range({3, 2}), trT);
  auto const R = make_tot(TA::Range({2, 4}), trT);
  auto rT = eval_result<ResultToT>(T);
  auto rL = eval_result<ResultToT>(L);
  auto rR = eval_result<ResultToT>(R);
  auto const Te = sequant::slice_array_over_mode(T, 1, 1, 2);
  auto const Le = sequant::slice_array_over_mode(L, 1, 1, 2);
  auto const Re = sequant::slice_array_over_mode(R, 1, 1, 2);

  SECTION("flat x flat") {
    auto sA = rA->slice_mode(2, 5, 10);
    auto sB = rB->slice_mode(0, 5, 10);
    auto const an = annot3("a,b,c", "c,d", "a,b,d");
    auto p = shaped(*sA, *sB, an, false);
    REQUIRE(p);
    REQUIRE(sA->is_buffer_alias());
    REQUIRE(sB->is_buffer_alias());
    auto e =
        shaped(*eval_result<ResultF>(Ae), *eval_result<ResultF>(Be), an, false);
    REQUIRE(p->get<FlatArray>().is_zero(0));
    REQUIRE(equal_tarrays(p->get<FlatArray>(), e->get<FlatArray>(), "a,b,d"));
  }

  SECTION("nested x flat, both orders") {
    auto sT = rT->slice_mode(1, 5, 10);
    auto sB = rB->slice_mode(0, 5, 10);
    auto const an = annot3("i,c;x,y", "c,d", "i,d;x,y");
    auto p = shaped(*sT, *sB, an, false);
    REQUIRE(p);
    REQUIRE(sT->is_buffer_alias());
    REQUIRE(sB->is_buffer_alias());
    auto e = shaped(*eval_result<ResultToT>(Te), *eval_result<ResultF>(Be), an,
                    false);
    REQUIRE(p->get<ToTArray>().is_zero(0));
    REQUIRE(equal_tarrays(p->get<ToTArray>(), e->get<ToTArray>(), "i,d;x,y"));
    auto const an2 = annot3("c,d", "i,c;x,y", "i,d;x,y");
    auto p2 = shaped(*sB, *sT, an2, false);
    REQUIRE(p2);
    REQUIRE(sT->is_buffer_alias());
    REQUIRE(equal_tarrays(p2->get<ToTArray>(), e->get<ToTArray>(), "i,d;x,y"));
    // scaled operands: the nested operand's factor is put on the flat one
    auto c2 = eval_result<ResultScalar<double>>(2.0);
    auto c3 = eval_result<ResultScalar<double>>(-3.0);
    auto scT =
        sT->prod(*c2, annot3("i,c;x,y", "", "i,c;x,y"), sequant::DeNest::False);
    auto scB = sB->prod(*c3, annot3("c,d", "", "c,d"), sequant::DeNest::False);
    REQUIRE(scT->is_buffer_alias());
    REQUIRE(scB->is_buffer_alias());
    auto ps = shaped(*scT, *scB, an, false);
    REQUIRE(ps);
    ToTArray e6;
    e6("i,d;x,y") = -6.0 * e->get<ToTArray>()("i,d;x,y");
    REQUIRE(equal_tarrays(ps->get<ToTArray>(), e6, "i,d;x,y"));
    auto ps2 = shaped(*scB, *scT, an2, false);
    REQUIRE(ps2);
    REQUIRE(equal_tarrays(ps2->get<ToTArray>(), e6, "i,d;x,y"));
    // the unshaped product takes the same operands
    auto pp = scT->prod(*scB, an, sequant::DeNest::False);
    auto eu = eval_result<ResultToT>(Te)->prod(*eval_result<ResultF>(Be), an,
                                               sequant::DeNest::False);
    ToTArray eu6;
    eu6("i,d;x,y") = -6.0 * eu->get<ToTArray>()("i,d;x,y");
    REQUIRE(equal_tarrays(pp->get<ToTArray>(), eu6, "i,d;x,y"));
    // a zero factor
    auto zT =
        sT->prod(*eval_result<ResultScalar<double>>(0.0),
                 annot3("i,c;x,y", "", "i,c;x,y"), sequant::DeNest::False);
    REQUIRE(zT->is_buffer_alias());
    REQUIRE(zT->norm2() == 0.0);
    auto pz = shaped(*zT, *sB, an, false);
    REQUIRE(pz);
    REQUIRE(pz->norm2() == 0.0);
    auto zB = sB->prod(*eval_result<ResultScalar<double>>(0.0),
                       annot3("c,d", "", "c,d"), sequant::DeNest::False);
    REQUIRE(zB->norm2() == 0.0);
    REQUIRE(zB->get<FlatArray>().trange() == Be.trange());
  }

  SECTION("nested x nested, nested and flattened results") {
    auto sL = rL->slice_mode(1, 5, 10);
    auto sR = rR->slice_mode(1, 5, 10);
    auto const an = annot3("i,j;a,b", "i,j;b,c", "i,j;a,c");
    auto p = shaped(*sL, *sR, an, false);
    REQUIRE(p);
    REQUIRE(sL->is_buffer_alias());
    REQUIRE(sR->is_buffer_alias());
    auto e = shaped(*eval_result<ResultToT>(Le), *eval_result<ResultToT>(Re),
                    an, false);
    REQUIRE(p->get<ToTArray>().is_zero(0));
    REQUIRE(equal_tarrays(p->get<ToTArray>(), e->get<ToTArray>(), "i,j;a,c"));

    auto sL2 = rL->slice_mode(1, 5, 10);
    auto const anD = annot3("i,j;a,b", "i,j;a,b", "i,j");
    auto d = shaped(*sL2, *rL->slice_mode(1, 5, 10), anD, true);
    REQUIRE(d);
    REQUIRE(sL2->is_buffer_alias());
    auto eD = shaped(*eval_result<ResultToT>(Le), *eval_result<ResultToT>(Le),
                     anD, true);
    REQUIRE(d->get<FlatArray>().is_zero(0));
    REQUIRE(equal_tarrays(d->get<FlatArray>(), eD->get<FlatArray>(), "i,j"));
  }
}

TEST_CASE("ta_lazy_view_cell_bytes_balance", "[eval][ta][lazy]") {
  // A pending view is charged 0 bytes when a cell stores it; a reader that
  // asks for its value materializes it in place. Releasing the cell must give
  // back what was charged, or live_bytes underflows.
  using ResultD = sequant::ResultTensorTA<TA::TSpArrayD>;
  auto& world = TA::get_default_world();
  TA::TSpArrayD arr(world, TA::TiledRange{{0, 2, 4}, {0, 3, 6}});
  arr.fill(1.0);
  world.gop.fence();
  auto const view = sequant::eval_result<ResultD>(arr)->slice_mode(0, 0, 2);
  REQUIRE(view->is_buffer_alias());

  sequant::eval::CellTable t;
  sequant::eval::TableCell cell;
  cell.value_id = 0;
  cell.production.kind = sequant::eval::ProductionKind::Build;
  t.cells.push_back(cell);
  sequant::eval::CellRegistry reg(t, {});

  auto& logger = sequant::Logger::instance();
  auto const prev_level = logger.eval.level;
  logger.eval.level = 1;  // size_in_bytes() is 0 untraced
  reg.set(0, view);
  CHECK(reg.live_bytes() == 0);
  auto const r = reg.read(0);
  CHECK(r->get<TA::TSpArrayD>().trange().dim(0).extent() == 2);
  CHECK_FALSE(r->is_buffer_alias());
  CHECK(r->size_in_bytes() > 0);
  CHECK(reg.take(0) == r);
  CHECK(reg.live_bytes() == 0);
  logger.eval.level = prev_level;
}
