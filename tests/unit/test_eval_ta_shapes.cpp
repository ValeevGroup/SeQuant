#include "test_eval_ta_fixtures.hpp"

using namespace sequant::tests::eval_ta;

TEST_CASE("shape_spike_tot_general_product", "[shape-spike]") {
  // De-risk spike: can (A(la) * B(ra)).set_shape(s) for a ToT general product
  // (a) evaluate without throwing, (b) honor the imposed SparseShape (zeroed
  // tiles stay zero), and (c) produce the same inner data as TA::einsum for
  // the surviving tiles?
  //
  // Uses the same ToT*ToT->ToT annotation as the ToT_times_ToT_to_ToT section
  // of eval_with_tiledarray, but with SparsePolicy and a multi-tile outer
  // TiledRange so there are enough outer tiles to zero one.
  //
  // Annotation (same structure as the existing ToT*ToT test):
  //   lhs: "i_1,i_2,i_3;a_4i_2i_3,a_1i_1i_2"
  //   rhs: "i_1,i_2,i_3;a_2i_1i_2,a_4i_2i_3"
  //   result: "i_2,i_1;a_1i_1i_2,a_2i_1i_2"
  // Contraction: outer i_3, inner a_4.  Fused outer (Hadamard): i_1, i_2.
  // Result outer: (i_2, i_1), each with 2 tiles -> 4 outer tiles total.

  using ToTArray = TA::DistArray<tot_tile_t<double>, TA::SparsePolicy>;

  auto& world = TA::get_default_world();

  // 2 tiles per occ mode; 1 tile per virt mode.
  // occ extent = 4, tile size = 2; virt extent = 3.
  TA::TiledRange1 const tr1_occ{0, 2, 4};
  std::size_t const virt_ext = 3;

  // Outer TiledRange for lhs and rhs: 3D (i1, i2, i3) all occ.
  TA::TiledRange const outer_tr{tr1_occ, tr1_occ, tr1_occ};

  // Inner tensor Range:
  //   lhs inner = (a4, a1), each virt_ext = (3, 3)
  //   rhs inner = (a2, a4), each virt_ext = (3, 3)
  TA::Range const inner_r{virt_ext, virt_ext};

  // Build a sparse ToT array with all outer tiles populated (norm > threshold).
  auto make_tot = [&](TA::TiledRange const& otr,
                      TA::Range const& ir) -> ToTArray {
    // "Dense" SparseShape: every tile has scaled norm 1.0.
    TA::SparseShape<float> all_nonzero{1.0f, otr};
    ToTArray arr{world, otr, all_nonzero};
    for (auto it = arr.begin(); it != arr.end(); ++it) {
      if (arr.is_local(it.index())) {
        // Capture ir by value: the lambda must not capture it by reference
        // since the taskq may run after this scope.
        TA::Range ir_copy = ir;
        *it = world.taskq.add(
            [ir_copy](TA::Range const& orng) {
              return random_tensor_of_tensor<double>(orng, ir_copy);
            },
            it.make_range());
      }
    }
    world.gop.fence();
    return arr;
  };

  auto lhs = make_tot(outer_tr, inner_r);
  auto rhs = make_tot(outer_tr, inner_r);

  // Matching the existing ToT*ToT test annotations.
  std::string const la{"i_1,i_2,i_3;a_4i_2i_3,a_1i_1i_2"};
  std::string const ra{"i_1,i_2,i_3;a_2i_1i_2,a_4i_2i_3"};
  std::string const ca{"i_2,i_1;a_1i_1i_2,a_2i_1i_2"};

  // Baseline: einsum.
  ToTArray C0 = TA::einsum(lhs(la), rhs(ra), ca);
  world.gop.fence();
  REQUIRE_FALSE(C0.is_initialized() == false);

  // Result outer TiledRange: 2D (i_2, i_1), both with {0,2,4} -> 2x2=4 tiles.
  TA::TiledRange const res_tr{tr1_occ, tr1_occ};
  REQUIRE(res_tr == C0.trange());

  // Construct SparseShape that zeros outer tile (0,0) of the result.
  // tile_norms uses do_not_scale=true: the values are treated as already-scaled
  // (so a value of 0.0 is exactly zero and >0 is nonzero).
  TA::Tensor<float> tile_norms{res_tr.tiles_range(), 1.0f};
  tile_norms[{0, 0}] = 0.0f;
  TA::SparseShape<float> s{tile_norms, res_tr, /*do_not_scale=*/true};
  // s must outlive the assignment below (set_shape stores &s internally).

  // Standard-layer evaluation with imposed shape.
  ToTArray C1;
  REQUIRE_NOTHROW(C1(ca) = (lhs(la) * rhs(ra)).set_shape(s));
  world.gop.fence();

  // Outcome (1): shape constraint honored -- tile (0,0) must be zero/absent.
  REQUIRE(C1.is_zero({0, 0}));

  // Outcome (2): surviving outer tiles must match the einsum baseline.
  // Walk C1's local tiles; for each nonzero outer tile, compare inner data
  // element-wise with C0.
  for (auto it = C1.begin(); it != C1.end(); ++it) {
    if (!C1.is_local(it.index()) || C1.is_zero(it.index())) continue;
    auto const& idx = it.index();
    REQUIRE_FALSE(C0.is_zero(idx));

    auto const& outer0 = C0.find(idx).get();  // tot_tile_t<double>
    auto const& outer1 = it->get();
    REQUIRE(outer0.range() == outer1.range());

    for (std::size_t o = 0; o < outer0.size(); ++o) {
      auto const& inner0 = outer0[o];
      auto const& inner1 = outer1[o];
      REQUIRE(inner0.range() == inner1.range());
      for (std::size_t k = 0; k < inner0.size(); ++k) {
        CHECK(inner0[k] == Catch::Approx(inner1[k]));
      }
    }
  }
}

TEST_CASE("shape_spike_T_times_ToT_general_product", "[shape-spike]") {
  // Case A: T x ToT -> ToT (mixed operand). Models a dressed DF integral g
  // (flat T) contracted with a PNO coefficient C (ToT) to yield a ToT result.
  // Confirms (T_op(la) * ToT_op(ra)).set_shape(s) for this mixed general
  // product (a) evaluates, (b) honors the imposed SparseShape, (c) matches the
  // TA::einsum baseline on surviving tiles.
  //
  // Annotation mirrors the existing T_times_ToT_to_ToT section:
  //   T_op (flat):  "i_3,i_1"
  //   ToT_op (ToT): "i_2,i_3;a_3i_2i_3,a_4i_2i_3"
  //   result (ToT): "i_1,i_2;a_3i_2i_3,a_4i_2i_3"
  // Shared/contracted outer mode: i_3.  Free outer: i_1 (from T), i_2 (from
  // ToT).  Result outer (i_1, i_2), each with 2 tiles -> 4 outer tiles.

  using FlatArray = TA::DistArray<TA::Tensor<double>, TA::SparsePolicy>;
  using ToTArray = TA::DistArray<tot_tile_t<double>, TA::SparsePolicy>;

  auto& world = TA::get_default_world();

  TA::TiledRange1 const tr1_occ{0, 2, 4};  // 2 tiles, occ extent 4
  std::size_t const virt_ext = 3;

  // T_op outer: (i_3, i_1), both occ.
  TA::TiledRange const t_outer_tr{tr1_occ, tr1_occ};
  // ToT_op outer: (i_2, i_3), both occ; inner (a_3, a_4) each virt.
  TA::TiledRange const tot_outer_tr{tr1_occ, tr1_occ};
  TA::Range const inner_r{virt_ext, virt_ext};

  // Build a flat sparse T array with all tiles populated.
  auto make_flat = [&](TA::TiledRange const& otr) -> FlatArray {
    TA::SparseShape<float> all_nonzero{1.0f, otr};
    FlatArray arr{world, otr, all_nonzero};
    for (auto it = arr.begin(); it != arr.end(); ++it) {
      if (arr.is_local(it.index())) {
        *it = world.taskq.add(
            [](TA::Range const& rng) { return random_tensor<double>(rng); },
            it.make_range());
      }
    }
    world.gop.fence();
    return arr;
  };

  // Build a sparse ToT array with all outer tiles populated.
  auto make_tot = [&](TA::TiledRange const& otr,
                      TA::Range const& ir) -> ToTArray {
    TA::SparseShape<float> all_nonzero{1.0f, otr};
    ToTArray arr{world, otr, all_nonzero};
    for (auto it = arr.begin(); it != arr.end(); ++it) {
      if (arr.is_local(it.index())) {
        TA::Range ir_copy = ir;
        *it = world.taskq.add(
            [ir_copy](TA::Range const& orng) {
              return random_tensor_of_tensor<double>(orng, ir_copy);
            },
            it.make_range());
      }
    }
    world.gop.fence();
    return arr;
  };

  auto T_op = make_flat(t_outer_tr);
  auto ToT_op = make_tot(tot_outer_tr, inner_r);

  std::string const la{"i_3,i_1"};
  std::string const ra{"i_2,i_3;a_3i_2i_3,a_4i_2i_3"};
  std::string const ca{"i_1,i_2;a_3i_2i_3,a_4i_2i_3"};

  // Baseline: einsum.
  ToTArray C0 = TA::einsum(T_op(la), ToT_op(ra), ca);
  world.gop.fence();

  // Result outer TiledRange: 2D (i_1, i_2) -> 2x2 = 4 tiles.
  TA::TiledRange const res_tr{tr1_occ, tr1_occ};
  REQUIRE(res_tr == C0.trange());

  // SparseShape zeroing outer tile (0,0).
  TA::Tensor<float> tile_norms{res_tr.tiles_range(), 1.0f};
  tile_norms[{0, 0}] = 0.0f;
  TA::SparseShape<float> s{tile_norms, res_tr, /*do_not_scale=*/true};

  // Standard-layer evaluation with imposed shape.
  ToTArray C1;
  REQUIRE_NOTHROW(C1(ca) = (T_op(la) * ToT_op(ra)).set_shape(s));
  world.gop.fence();

  // (1) shape constraint honored.
  REQUIRE(C1.is_zero({0, 0}));

  // (2) surviving tiles match einsum baseline.
  for (auto it = C1.begin(); it != C1.end(); ++it) {
    if (!C1.is_local(it.index()) || C1.is_zero(it.index())) continue;
    auto const& idx = it.index();
    REQUIRE_FALSE(C0.is_zero(idx));

    auto const& outer0 = C0.find(idx).get();
    auto const& outer1 = it->get();
    REQUIRE(outer0.range() == outer1.range());

    for (std::size_t o = 0; o < outer0.size(); ++o) {
      auto const& inner0 = outer0[o];
      auto const& inner1 = outer1[o];
      REQUIRE(inner0.range() == inner1.range());
      for (std::size_t k = 0; k < inner0.size(); ++k) {
        CHECK(inner0[k] == Catch::Approx(inner1[k]));
      }
    }
  }
}

TEST_CASE("shape_spike_ToT_inner_contraction_to_flat_T", "[shape-spike]") {
  // Case B: ToT x ToT inner-contraction -> flat T (the dot_inner / DeNest::True
  // path; see result.hpp:581 TA::einsum<TA::DeNest::True>). Two ToT operands
  // share their outer (occ) modes (Hadamard) and fully contract their inner
  // (composite) modes, denesting to a plain tensor-of-scalars result.
  //
  // The standard-expression-layer equivalent of einsum<DeNest::True> is the
  // .dot_inner() expression (TA einsum/tiledarray.h:603):
  //   C(c) = A(a + inner.a).dot_inner(B(b + inner.b));
  // DotInnerExpr derives from Expr, so it exposes set_shape(); the spike here
  // is whether that override is honored for the denested flat-T result.
  //
  // Annotation:
  //   lhs (ToT): "i_1,i_2;a_1i_1i_2,a_2i_1i_2"
  //   rhs (ToT): "i_1,i_2;a_1i_1i_2,a_2i_1i_2"
  //   result (flat T): "i_1,i_2"   (inner a_1,a_2 dotted away)
  // Outer i_1,i_2 are Hadamard (survive); result outer (i_1,i_2) -> 2x2 tiles.

  using FlatArray = TA::DistArray<TA::Tensor<double>, TA::SparsePolicy>;
  using ToTArray = TA::DistArray<tot_tile_t<double>, TA::SparsePolicy>;

  auto& world = TA::get_default_world();

  TA::TiledRange1 const tr1_occ{0, 2, 4};  // 2 tiles, occ extent 4
  std::size_t const virt_ext = 3;

  // Both operands: outer (i_1, i_2) occ; inner (a_1, a_2) virt.
  TA::TiledRange const outer_tr{tr1_occ, tr1_occ};
  TA::Range const inner_r{virt_ext, virt_ext};

  auto make_tot = [&](TA::TiledRange const& otr,
                      TA::Range const& ir) -> ToTArray {
    TA::SparseShape<float> all_nonzero{1.0f, otr};
    ToTArray arr{world, otr, all_nonzero};
    for (auto it = arr.begin(); it != arr.end(); ++it) {
      if (arr.is_local(it.index())) {
        TA::Range ir_copy = ir;
        *it = world.taskq.add(
            [ir_copy](TA::Range const& orng) {
              return random_tensor_of_tensor<double>(orng, ir_copy);
            },
            it.make_range());
      }
    }
    world.gop.fence();
    return arr;
  };

  auto lhs = make_tot(outer_tr, inner_r);
  auto rhs = make_tot(outer_tr, inner_r);

  // Outer annotations (a/b) and inner annotations (inner.a/inner.b) for
  // .dot_inner(); the result is the bare outer annotation (no inner part).
  std::string const a_outer{"i_1,i_2"};
  std::string const a_inner{";a_1i_1i_2,a_2i_1i_2"};
  std::string const b_outer{"i_1,i_2"};
  std::string const b_inner{";a_1i_1i_2,a_2i_1i_2"};
  std::string const ca{"i_1,i_2"};

  // Baseline: einsum<DeNest::True>. The full ToT annotation carries the inner
  // part; the result annotation is outer-only (denested).
  FlatArray C0 = TA::einsum<TA::DeNest::True>(lhs(a_outer + a_inner),
                                              rhs(b_outer + b_inner), ca);
  world.gop.fence();

  // Result outer TiledRange: 2D (i_1, i_2) -> 2x2 = 4 tiles.
  TA::TiledRange const res_tr{tr1_occ, tr1_occ};
  REQUIRE(res_tr == C0.trange());

  // SparseShape zeroing outer tile (0,0).
  TA::Tensor<float> tile_norms{res_tr.tiles_range(), 1.0f};
  tile_norms[{0, 0}] = 0.0f;
  TA::SparseShape<float> s{tile_norms, res_tr, /*do_not_scale=*/true};

  // Standard-layer denesting product via .dot_inner(), with imposed shape.
  FlatArray C1;
  REQUIRE_NOTHROW(C1(ca) = (lhs(a_outer + a_inner))
                               .dot_inner(rhs(b_outer + b_inner))
                               .set_shape(s));
  world.gop.fence();

  // (1) shape constraint honored.
  REQUIRE(C1.is_zero({0, 0}));

  // (2) surviving flat tiles match the einsum<DeNest::True> baseline.
  for (auto it = C1.begin(); it != C1.end(); ++it) {
    if (!C1.is_local(it.index()) || C1.is_zero(it.index())) continue;
    auto const& idx = it.index();
    REQUIRE_FALSE(C0.is_zero(idx));

    // Copy the tiles out BY VALUE: C0.find(idx) and *it are temporary Futures
    // whose payload get() only references, so binding `auto const&` would
    // dangle once those temporaries die at the semicolon (ASan
    // stack-use-after-scope on the tile's Range). A TA::Tensor copy is a cheap
    // shallow handle that keeps its own Range alive for the loop body.
    TA::Tensor<double> const t0 = C0.find(idx).get();
    TA::Tensor<double> const t1 = it->get();
    REQUIRE(t0.range() == t1.range());
    for (std::size_t k = 0; k < t0.size(); ++k) {
      CHECK(t0[k] == Catch::Approx(t1[k]));
    }
  }
}

TEST_CASE("shape_provider_general_product", "[shape-provider]") {
  // A provider returning a REAL SparseShape for a ToT*ToT->ToT general
  // product (the `(lhs(la) * rhs(ra)).set_shape(s)` emission form) must
  // constrain the eval result: the zeroed outer tiles are is_zero and the
  // surviving tiles equal the unshaped (einsum) baseline. Also covers the
  // no-op (full-ones shape) and nullopt (decline) cases.
  //
  // Graph: ToT * ToT -> ToT, single internal Product node. SparsePolicy so a
  // SparseShape is meaningful; multi-tiled occ so >=1 outer tile can be zeroed.
  using sequant::CacheManager;
  using sequant::evaluate;
  using sequant::TAEvalContext;
  using node_t = sequant::FullBinaryNode<sequant::EvalExprTA>;
  using cache_t = CacheManager<node_t>;

  auto& world = TA::get_default_world();
  // occ extent 4 in tiles of <=2 -> 2 outer tiles per occ mode.
  rand_tensor_yield<double, TA::SparsePolicy> yield{world, /*nocc=*/4,
                                                    /*nvirt=*/3};
  yield.set_max_tile(2);

  using ToTArray = typename decltype(yield)::array_tot_type;

  constexpr std::wstring_view expr_str =
      L"I{a4<i2,i3>,a1<i1,i2>;i1,i2}"
      L" * "
      L"s{a2<i1,i2>;a4<i2,i3>}";
  auto const node = eval_node(sequant::deserialize<sequant::ExprPtr>(expr_str));
  // i_3 (pair basis of the contracted composite) is an outer mode of the
  // result; the node's canonical layout is the target.
  std::string const target = node->annot();

  // Unshaped reference (no hook).
  auto const ref = evaluate(node, target, yield)->get<ToTArray>();
  TA::TiledRange const res_tr = ref.trange();
  // The result's OUTER modes are i_1, i_2 and i_3 -- the contracted composite
  // a4<i2,i3> is summed over, but its pair-basis index i_3 survives as a free
  // outer mode (see the comment above `target`). So an outer tile index is
  // 3-dimensional and the result's TA annotation is `target`, not a hand-made
  // two-outer-mode string. Pin the rank so a layout change fails here loudly
  // rather than as an opaque TA_ASSERT inside a tile lookup.
  REQUIRE(res_tr.tiles_range().rank() == 3);

  auto make_ctx = [](auto shape_fn) {
    TAEvalContext ctx;
    ctx.result_shape_provider =
        [shape_fn](
            sequant::FullBinaryNode<sequant::EvalExprTA> const&,
            TA::TiledRange const& tr) -> std::optional<TA::SparseShape<float>> {
      return shape_fn(tr);
    };
    return ctx;
  };

  SECTION("real shape: zeroed tile is_zero, survivors match baseline") {
    // Zero outer tile (0,0,0); keep the rest.
    auto ctx = make_ctx([](TA::TiledRange const& tr) {
      TA::Tensor<float> norms{tr.tiles_range(), 1.0f};
      norms[{0, 0, 0}] = 0.0f;
      return TA::SparseShape<float>{norms, tr, /*do_not_scale=*/true};
    });
    auto cache = cache_t::empty();
    cache.set_shaped_product_hook(
        TAEvalContext::make_hook<double, TA::SparsePolicy>(ctx));
    auto const res = evaluate(node, target, yield, cache)->get<ToTArray>();

    // (1) zeroed tile gone.
    REQUIRE(res.is_zero({0, 0, 0}));
    // (2) survivors equal the unshaped baseline.
    for (auto it = res.begin(); it != res.end(); ++it) {
      if (!res.is_local(it.index()) || res.is_zero(it.index())) continue;
      auto const& idx = it.index();
      REQUIRE_FALSE(ref.is_zero(idx));
      auto const& o0 = ref.find(idx).get();
      auto const& o1 = it->get();
      REQUIRE(o0.range() == o1.range());
      for (std::size_t o = 0; o < o0.size(); ++o) {
        auto const& in0 = o0[o];
        auto const& in1 = o1[o];
        REQUIRE(in0.range() == in1.range());
        for (std::size_t k = 0; k < in0.size(); ++k)
          CHECK(in0[k] == Catch::Approx(in1[k]));
      }
    }
  }

  SECTION("no-op full-ones shape: result equals unshaped") {
    auto ctx = make_ctx([](TA::TiledRange const& tr) {
      return TA::SparseShape<float>{1.0f, tr};
    });
    auto cache = cache_t::empty();
    cache.set_shaped_product_hook(
        TAEvalContext::make_hook<double, TA::SparsePolicy>(ctx));
    auto const res = evaluate(node, target, yield, cache)->get<ToTArray>();
    std::string const& annot = target;
    REQUIRE(Catch::Approx(ref(annot).dot(ref(annot))) ==
            res(annot).dot(res(annot)));
  }

  SECTION("nullopt: hook declines, falls through to unshaped prod()") {
    auto ctx = make_ctx(
        [](TA::TiledRange const&) -> std::optional<TA::SparseShape<float>> {
          return std::nullopt;
        });
    auto cache = cache_t::empty();
    cache.set_shaped_product_hook(
        TAEvalContext::make_hook<double, TA::SparsePolicy>(ctx));
    auto const res = evaluate(node, target, yield, cache)->get<ToTArray>();
    std::string const& annot = target;
    REQUIRE(Catch::Approx(ref(annot).dot(ref(annot))) ==
            res(annot).dot(res(annot)));
  }

  SECTION("no hook installed => unshaped behavior unchanged") {
    auto const res = evaluate(node, target, yield)->get<ToTArray>();
    std::string const& annot = target;
    REQUIRE(Catch::Approx(ref(annot).dot(ref(annot))) ==
            res(annot).dot(res(annot)));
  }
}

TEST_CASE("shape_provider_denest_to_flat", "[shape-provider]") {
  // The DeNest::True path (ToT * ToT -> flat T), emitted as
  // `lhs(la).dot_inner(rhs(ra)).set_shape(s)`. A ResultExpr pins the flat head
  // (bare occ indices, no protos), so node->tot() is false while both operands
  // are ToT => de_nest == True at the binary site. A provider returning a real
  // shape must zero the designated outer tile and leave survivors matching the
  // unshaped (einsum<DeNest::True>) baseline.
  using sequant::CacheManager;
  using sequant::evaluate;
  using sequant::TAEvalContext;
  using node_t = sequant::FullBinaryNode<sequant::EvalExprTA>;
  using cache_t = CacheManager<node_t>;

  auto& world = TA::get_default_world();
  rand_tensor_yield<double, TA::SparsePolicy> yield{world, /*nocc=*/4,
                                                    /*nvirt=*/3};
  yield.set_max_tile(2);

  using FlatArray = typename decltype(yield)::array_type;

  // Result head D{i1,i2} (flat): both operands carry composite virtuals
  // a1<i1,i2>,a2<i1,i2> that fully contract (inner), leaving bare i1,i2.
  auto const res_expr = sequant::deserialize<sequant::ResultExpr>(
      L"D{i1,i2} = "
      L"R{a1<i1,i2>,a2<i1,i2>;i1,i2} * g{i1,i2;a1<i1,i2>,a2<i1,i2>}",
      {.def_braket_symm = sequant::Hermiticity::NonHermitian});
  auto const node = eval_node(res_expr);
  // Confirm this really is the denest (DeNest::True) path.
  REQUIRE(node.left()->tot());
  REQUIRE(node.right()->tot());
  REQUIRE_FALSE(node->tot());

  std::string const target{"i_1,i_2"};
  auto const ref = evaluate(node, target, yield)->get<FlatArray>();
  TA::TiledRange const res_tr = ref.trange();

  SECTION("denest-real-shape") {
    TAEvalContext ctx;
    ctx.result_shape_provider =
        [](sequant::FullBinaryNode<sequant::EvalExprTA> const&,
           TA::TiledRange const& tr) -> std::optional<TA::SparseShape<float>> {
      TA::Tensor<float> norms{tr.tiles_range(), 1.0f};
      norms[{0, 0}] = 0.0f;
      return TA::SparseShape<float>{norms, tr, /*do_not_scale=*/true};
    };
    auto cache = cache_t::empty();
    cache.set_shaped_product_hook(
        TAEvalContext::make_hook<double, TA::SparsePolicy>(ctx));
    auto const res = evaluate(node, target, yield, cache)->get<FlatArray>();

    REQUIRE(res.is_zero({0, 0}));
    for (auto it = res.begin(); it != res.end(); ++it) {
      if (!res.is_local(it.index()) || res.is_zero(it.index())) continue;
      auto const& idx = it.index();
      REQUIRE_FALSE(ref.is_zero(idx));
      auto const& t0 = ref.find(idx).get();
      auto const& t1 = it->get();
      REQUIRE(t0.range() == t1.range());
      for (std::size_t k = 0; k < t0.size(); ++k)
        CHECK(t0[k] == Catch::Approx(t1[k]));
    }
  }

  SECTION("nullopt on the denest result: unchanged") {
    TAEvalContext ctx;
    ctx.result_shape_provider =
        [](sequant::FullBinaryNode<sequant::EvalExprTA> const&,
           TA::TiledRange const&) -> std::optional<TA::SparseShape<float>> {
      return std::nullopt;
    };
    auto cache = cache_t::empty();
    cache.set_shaped_product_hook(
        TAEvalContext::make_hook<double, TA::SparsePolicy>(ctx));
    auto const res = evaluate(node, target, yield, cache)->get<FlatArray>();
    REQUIRE(equal_tarrays(res, ref, "i,j"));
  }
}
