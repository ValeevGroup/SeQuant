#include "test_eval_ta_fixtures.hpp"

using namespace sequant::tests::eval_ta;

// Regression tests for tensor-of-tensor (ToT) ops on arrays whose inner tiles
// are all empty (e.g. a fully screened CSV residual). tot_inner_rank() reads
// the rank from a populated inner tile, so it returns 0 here; downstream ToT
// ops must not turn that 0 into a degenerate (non-ToT) annotation. Reproduces
// the mpqc4 PAO-CSV abort fixed in add_inplace and column_symmetrize_ta.
TEST_CASE("tot_all_empty_inner", "[eval][tot]") {
  using sequant::eval_result;
  using sequant::ResultPtr;
  using sequant::ResultTensorOfTensorTA;
  using sequant::detail::tot_inner_rank;
  using ToTArray = TA::DistArray<tot_tile_t<double>, TA::SparsePolicy>;
  using ResultToT = ResultTensorOfTensorTA<ToTArray>;

  auto& world = TA::get_default_world();

  // Build a ToT array; `empty_inner` leaves every inner tensor empty.
  auto build = [&world](unsigned outer_rank, bool empty_inner,
                        TA::TiledRange1 outer_tr1 =
                            TA::TiledRange1{0, 3}) -> ToTArray {
    std::vector<TA::TiledRange1> tr1s(outer_rank, outer_tr1);
    TA::TiledRange outer_tr(tr1s.begin(), tr1s.end());
    std::vector<std::size_t> inner_ext(outer_rank, 4);
    TA::Range inner_r(inner_ext);

    auto tile_fn = [inner_r, empty_inner](TA::Range const& orng) {
      return make_tot_tile<double>(
          orng,
          [&](std::size_t) { return empty_inner ? TA::Range{} : inner_r; },
          [](std::size_t, std::size_t) { return 1.0; });
    };

    ToTArray arr{world, outer_tr};
    for (auto it = arr.begin(); it != arr.end(); ++it)
      if (arr.is_local(it.index()))
        *it = world.taskq.add(tile_fn, it.make_range());
    world.gop.fence();
    return arr;
  };

  SECTION("add-inplace-empty-populated") {
    auto empty = build(1, /*empty_inner=*/true);
    auto full = build(1, /*empty_inner=*/false);
    REQUIRE(tot_inner_rank(empty) == 0);
    REQUIRE(tot_inner_rank(full) == 1);

    auto r_empty = eval_result<ResultToT>(empty);
    auto r_full = eval_result<ResultToT>(full);
    REQUIRE_NOTHROW(r_empty->add_inplace(*r_full));
    auto const& res = r_empty->get<ToTArray>();
    REQUIRE(tot_inner_rank(res) == 1);
    // the populated operand's data (all 1.0) must actually have landed.
    for (auto it = res.begin(); it != res.end(); ++it)
      if (res.is_local(it.index()))
        for (auto const& inner : it->get()) {
          REQUIRE_FALSE(inner.empty());
          for (std::size_t k = 0; k < inner.size(); ++k)
            REQUIRE(inner[k] == 1.0);
        }
  }

  SECTION("add_inplace: populated += empty is a no-op") {
    auto full = build(1, false);
    auto empty = build(1, true);
    auto r_full = eval_result<ResultToT>(full);
    auto r_empty = eval_result<ResultToT>(empty);
    REQUIRE_NOTHROW(r_full->add_inplace(*r_empty));
    REQUIRE(tot_inner_rank(r_full->get<ToTArray>()) == 1);
  }

  SECTION("add_inplace: empty += empty is an identity no-op") {
    auto e1 = build(1, true);
    auto e2 = build(1, true);
    auto r1 = eval_result<ResultToT>(e1);
    auto r2 = eval_result<ResultToT>(e2);
    REQUIRE_NOTHROW(r1->add_inplace(*r2));
    REQUIRE(tot_inner_rank(r1->get<ToTArray>()) == 0);
  }

  SECTION("symmetrize: all-empty ToT is the identity") {
    auto empty = build(1, true);
    REQUIRE(tot_inner_rank(empty) == 0);
    auto r_empty = eval_result<ResultToT>(empty);
    ResultPtr sym;
    REQUIRE_NOTHROW(sym = r_empty->symmetrize());
    REQUIRE(tot_inner_rank(sym->get<ToTArray>()) == 0);
  }

  SECTION("slice_array_over_mode: all-empty ToT slices to a zero ToT") {
    auto empty = build(1, true);
    REQUIRE(tot_inner_rank(empty) == 0);
    ToTArray sliced;
    REQUIRE_NOTHROW(sliced = sequant::slice_array_over_mode(empty, 0, 0, 1));
    REQUIRE(tot_inner_rank(sliced) == 0);
    REQUIRE(sliced.trange().dim(0).tile_extent() == 1);
  }

  SECTION("slice-rebase-multitile") {
    // Slice across two tiles at a non-zero base; hand-built trange must match
    // what the populated path gets from TA's block().
    TA::TiledRange1 const tr1{0, 3, 6, 9};
    auto empty = build(1, /*empty_inner=*/true, tr1);
    auto full = build(1, /*empty_inner=*/false, tr1);
    auto const s_empty = sequant::slice_array_over_mode(empty, 0, 1, 3);
    auto const s_full = sequant::slice_array_over_mode(full, 0, 1, 3);
    REQUIRE(tot_inner_rank(s_empty) == 0);
    REQUIRE(s_empty.trange() == s_full.trange());
  }
}

// Result::pow on a ToT: an absent (empty) inner tile is exactly zero, which a
// positive exponent keeps; a non-positive exponent would make it inf or 1,
// which an empty tile cannot hold, so it is refused.
TEST_CASE("tot_elementwise_pow", "[eval][tot]") {
  using sequant::rational;
  using sequant::detail::ta_elementwise_pow;
  using ToTArray = TA::DistArray<tot_tile_t<double>, TA::SparsePolicy>;

  auto& world = TA::get_default_world();
  TA::Range const inner_r(std::vector<std::size_t>{2, 2});

  // a single outer tile of 3 inner tiles; inner element j of inner tile i is
  // 1 + i + j/4 (never zero). With `empty_middle` inner tile 1 is empty.
  auto build = [&world, inner_r](bool empty_middle) -> ToTArray {
    auto tile_fn = [inner_r, empty_middle](TA::Range const& orng) {
      return make_tot_tile<double>(
          orng,
          [&](std::size_t i) {
            return empty_middle && i == 1 ? TA::Range{} : inner_r;
          },
          [](std::size_t i, std::size_t j) {
            return 1.0 + static_cast<double>(i) + j / 4.0;
          });
    };
    ToTArray arr{world, TA::TiledRange{TA::TiledRange1{0, 3}}};
    for (auto it = arr.begin(); it != arr.end(); ++it)
      if (arr.is_local(it.index()))
        *it = world.taskq.add(tile_fn, it.make_range());
    world.gop.fence();
    return arr;
  };

  // checks that every inner element of `got` is f(element of `in`), and that
  // the empty inner tiles of `in` are exactly those of `got`
  auto check = [](ToTArray const& in, ToTArray const& got, auto const& f) {
    for (auto it = in.begin(); it != in.end(); ++it) {
      auto const t_in = it->get();
      auto const t_got = got.find_local(it.index()).get();
      REQUIRE(t_got.size() == t_in.size());
      for (std::size_t i = 0; i < t_in.size(); ++i) {
        REQUIRE(t_got.data()[i].empty() == t_in.data()[i].empty());
        for (std::size_t j = 0; j < t_in.data()[i].size(); ++j)
          CHECK(t_got.data()[i].data()[j] ==
                Catch::Approx(f(t_in.data()[i].data()[j])).margin(1e-12));
      }
    }
  };

  SECTION("positive exponent keeps an absent inner tile absent") {
    auto const arr = build(/*empty_middle=*/true);
    check(arr, ta_elementwise_pow(arr, rational{2}),
          [](double x) { return x * x; });
  }

  SECTION("negative exponent without absent inner tiles") {
    auto const arr = build(/*empty_middle=*/false);
    check(arr, ta_elementwise_pow(arr, rational{-1}),
          [](double x) { return 1 / x; });
  }

  SECTION("non-positive exponent of an absent inner tile is refused") {
    if (sequant::assert_behavior() != sequant::AssertBehavior::Abort) {
      auto const arr = build(/*empty_middle=*/true);
      REQUIRE_THROWS_AS(ta_elementwise_pow(arr, rational{-1}),
                        sequant::Exception);
      REQUIRE_THROWS_AS(ta_elementwise_pow(arr, rational{0}),
                        sequant::Exception);
    }
  }
}

TEST_CASE("measured size counts a shared nested cell once", "[eval][tot]") {
  // A pair-mirrored leaf stores its K <= L tiles and gives each L > K tile
  // cells that are shallow views of its partner's cells (MPQC's
  // mirror_pair_tiles), so the array holds the upper tiles' bytes only.
  // TA::size_of counts a shared cell once per tile it appears in; the
  // measurement must count the storage once.
  using ToTArray = TA::DistArray<tot_tile_t<double>, TA::SparsePolicy>;
  auto& world = TA::get_default_world();
  if (world.size() != 1) return;  // a mirror tile must be local to its partner
  TA::TiledRange const tr{{0, 2, 4}, {0, 2, 4}};
  auto const& tiles = tr.tiles_range();
  auto const inner = [](std::size_t) { return TA::Range{3}; };
  auto const fill = [](std::size_t o, std::size_t k) {
    return 1.0 + static_cast<double>(o) + static_cast<double>(k);
  };
  auto const owning = [&](std::size_t ord) {
    return make_tot_tile<double>(tr.make_tile_range(ord), inner, fill);
  };
  auto build = [&](bool mirrored) {
    ToTArray arr(world, tr);
    auto const upper = owning(tiles.ordinal({0, 1}));
    arr.set(tiles.ordinal({0, 0}), owning(tiles.ordinal({0, 0})));
    arr.set(tiles.ordinal({1, 1}), owning(tiles.ordinal({1, 1})));
    arr.set(tiles.ordinal({0, 1}), upper);
    auto const lower_ord = tiles.ordinal({1, 0});
    if (!mirrored) {
      arr.set(lower_ord, owning(lower_ord));
    } else {
      // null cells, each rebound to the transposed cell of the partner
      tot_tile_t<double> lower(tr.make_tile_range(lower_ord));
      for (auto const& e : lower.range()) {
        std::array<std::size_t, 2> const m{static_cast<std::size_t>(e[1]),
                                           static_cast<std::size_t>(e[0])};
        lower(e) = upper(m);
      }
      arr.set(lower_ord, std::move(lower));
    }
    world.gop.fence();
    return arr;
  };
  auto const ta_bytes = [](ToTArray const& arr) {
    return TA::size_of<TA::MemorySpace::Host>(arr);
  };
  auto const measured = [](ToTArray const& arr) {
    return sequant::eval_result<sequant::ResultTensorOfTensorTA<ToTArray>>(arr)
        ->measured_size_in_bytes();
  };
  // no sharing: the measurement is TA's count
  auto const plain = build(false);
  REQUIRE(measured(plain) == ta_bytes(plain));
  // the four cells of the lower tile view the upper tile's cells: TA counts
  // their storage twice, the measurement once
  auto const mirrored = build(true);
  REQUIRE(ta_bytes(mirrored) - measured(mirrored) ==
          4 * tot_inner_t<double>::cell_size(3));
  // and the mirrored array holds less than the plain one by the same amount
  REQUIRE(measured(plain) - measured(mirrored) ==
          4 * tot_inner_t<double>::cell_size(3));
}

TEST_CASE("measured-size-trace-gate", "[eval]") {
  // size_in_bytes() is collective for a DistArray and so switched off
  // untraced; the schedule document's run record measures through
  // measured_size_in_bytes() instead
  auto& world = TA::get_default_world();
  TA::TSpArrayD arr(world, TA::TiledRange{{0, 2, 4}, {0, 3, 6}});
  arr.fill(1.0);
  world.gop.fence();
  auto const res =
      sequant::eval_result<sequant::ResultTensorTA<TA::TSpArrayD>>(arr);
  auto& logger = sequant::Logger::instance();
  auto const prev_level = logger.eval.level;
  logger.eval.level = 0;
  REQUIRE(res->size_in_bytes() == 0);
  REQUIRE(res->measured_size_in_bytes() >= 4 * 6 * sizeof(double));
  logger.eval.level = 1;
  REQUIRE(res->size_in_bytes() == res->measured_size_in_bytes());
  logger.eval.level = prev_level;
}

TEST_CASE("canon-phase-of-a-cached-view", "[eval]") {
  // A cache entry memoizes its bytes when first measured. An entry that is
  // itself a phase view, the only holder of its parent's tiles, owns those
  // bytes; the read that applies the node's phase again must not make it
  // look shared before the cache measures it.
  using node_t = sequant::FullBinaryNode<sequant::EvalExprTA>;
  using ResultD = sequant::ResultTensorTA<TA::TArrayD>;
  auto& world = TA::get_default_world();
  // the phase of the antisymmetric g, written out of canonical order, sits on
  // the product (see ordered-canon-phase-once)
  auto const node = eval_node(sequant::deserialize<sequant::ExprPtr>(
      L"g{a_2,a_1;i_3,i_4} * t{i_3,i_4;i_1,i_2}",
      {.def_perm_symm = sequant::Symmetry::Antisymm}));
  REQUIRE(node->canon_phase() == -1);
  TA::TArrayD arr(world,
                  TA::TiledRange{{0, 2, 4}, {0, 2, 4}, {0, 2, 4}, {0, 2, 4}});
  arr.fill(1.0);
  world.gop.fence();

  auto& level = sequant::Logger::instance().eval.level;
  struct Restore {
    std::size_t& level;
    std::size_t prev;
    ~Restore() { level = prev; }
  } const restore{level, level};
  level = 1;  // size_in_bytes() measures only traced

  auto const bytes = sequant::eval_result<ResultD>(arr)->size_in_bytes();
  REQUIRE(bytes > 0);
  sequant::CacheManager<node_t> cache{
      std::vector<std::pair<node_t, std::size_t>>{{node, 2}}};
  auto const stored = cache.store_and_access(
      node, sequant::eval_result<ResultD>(TA::clone(arr))->mult_by_phase(-1));
  REQUIRE(stored->is_buffer_alias());
  auto const read =
      sequant::apply_canon_phase<sequant::Trace::On>(node, stored, cache);
  REQUIRE(read->is_buffer_alias());
  REQUIRE(cache.size_in_bytes() == bytes);
}

TEST_CASE("eval_with_tiledarray", "[eval]") {
  // Reproducer for the mpqc4 cck real-field NaN regression. The eval-graph
  // head's bra/ket split is purely positional: each external index ends up in
  // whichever slot (bra or ket) it occupied in its source tensor, summed
  // across the term. Two orientations of a tensor that are equivalent under
  // bra<->ket-swap (Symm braket_symmetry) produce different head bra_rank,
  // which breaks downstream code (e.g. mpqc's jacobi_update) that assumes a
  // conventional 2:2 (vir,vir;occ,occ) layout for a CCSD T2 residual head.
  // Bug is independent of scalar Field — fires under both Conjugate and Symm
  // whenever the canonical orientation puts an external on the "wrong" side.
  SECTION("head-braket-positional") {
    using sequant::deserialize;
    using sequant::EvalExprTA;
    using sequant::ExprPtr;

    // The test expressions below put an internal index (e.g. a_3) in the bra
    // of two factors (g.bra and t.bra), which violates the default-context
    // strict bra↔ket-symmetry policy (each internal must appear at most once
    // per side under Conjugate). Disable the policy for this scope; we are
    // probing the eval-graph head's bra/ket layout, not the canonicalizer's
    // covariance assumptions.
    auto ctx_resetter = sequant::set_scoped_default_context(
        sequant::Context{sequant::get_default_context()}.set(
            sequant::AssertStrictBraKetSymmetry::No));

    auto report = [](sequant::ExprPtr const& e, std::string const& label) {
      auto node = eval_node(e);
      auto const& head = node->as_tensor();
      std::wstring head_str = sequant::to_latex(head);
      std::string head_str8{head_str.begin(), head_str.end()};
      INFO(label + " head: " + head_str8 +
           "  bra_rank=" + std::to_string(head.bra_rank()) +
           "  ket_rank=" + std::to_string(head.ket_rank()));
      return std::make_pair(head.bra_rank(), head.ket_rank());
    };

    // Representative SF R2 residual term (h2o-cck-2-631g-pvdz) with externals
    // {a_1, a_2, i_1, i_2}. For a CCSD T2 the head ought to be I{a_1,a_2; i_1,
    // i_2} (vir,vir;occ,occ) so downstream Jacobi-style updates index orbital
    // energies correctly.

    // (A) g written with the occ external i_2 in its bra slot: regardless of
    // braket_symmetry, the head ends up I{i_2,a_1,a_2; i_1} — head bra/ket is
    // assigned by external SLOT, not by space.
    auto expr_bra_external_symm = deserialize(
        L"2 g{i_2,a_3;i_3,i_4}:N-S-S * t{a_3;i_3}:N-N-S "
        L"* t{a_1,a_2;i_1,i_4}:N-N-S");
    REQUIRE(expr_bra_external_symm);
    auto [br_symm, kr_symm] =
        report(expr_bra_external_symm, "g{i_2,a_3;...}:N-S-S");

    auto expr_bra_external_conj = deserialize(
        L"2 g{i_2,a_3;i_3,i_4}:N-C-S * t{a_3;i_3}:N-N-S "
        L"* t{a_1,a_2;i_1,i_4}:N-N-S");
    REQUIRE(expr_bra_external_conj);
    auto [br_conj, kr_conj] =
        report(expr_bra_external_conj, "g{i_2,a_3;...}:N-C-S");

    // (B) Same expression with g's bra/ket pre-swapped (mathematically
    // equivalent under Symm braket_symmetry). Braket-symmetric tensors now
    // canonicalize by bra/ket *color* (index spaces), so g canonicalizes to the
    // same orientation as in (A) regardless of how it is written: the head
    // comes out the same 3:1 as (A). Pre-swapping no longer steers the head
    // layout — use the ResultExpr API (C) to pin it.
    auto expr_ket_external_swap = deserialize(
        L"2 g{i_3,i_4;i_2,a_3}:N-S-S * t{a_3;i_3}:N-N-S "
        L"* t{a_1,a_2;i_1,i_4}:N-N-S");
    REQUIRE(expr_ket_external_swap);
    auto [br_swap, kr_swap] =
        report(expr_ket_external_swap, "g{i_3,i_4;i_2,a_3}:N-S-S (pre-swap)");

    // (C) Caller-supplied head layout via ResultExpr: the right shape of the
    // public API for this. The caller writes the LHS with the bra/ket layout
    // it wants the head to have; binarize(ResultExpr) at eval_expr.hpp:435
    // overwrites the eval-tree root's tensor with res.result_as_tensor(),
    // making the IR's positional choice irrelevant.
    auto res_explicit_layout = sequant::deserialize<sequant::ResultExpr>(
        L"R2{a_1,a_2;i_1,i_2}:N-N-S = "
        L"2 g{i_2,a_3;i_3,i_4}:N-S-S * t{a_3;i_3}:N-N-S "
        L"* t{a_1,a_2;i_1,i_4}:N-N-S");
    auto node_explicit = eval_node(res_explicit_layout);
    auto const& head_explicit = node_explicit->as_tensor();
    std::wstring head_explicit_str = sequant::to_latex(head_explicit);
    std::string head_explicit_str8{head_explicit_str.begin(),
                                   head_explicit_str.end()};
    INFO("ResultExpr R2{a_1,a_2;i_1,i_2} head: " + head_explicit_str8 +
         "  bra_rank=" + std::to_string(head_explicit.bra_rank()) +
         "  ket_rank=" + std::to_string(head_explicit.ket_rank()));

    // Document the current behavior:
    // (A) same 3:1 split regardless of Symm vs Conjugate — head bra/ket is
    // assigned positionally (by external slot), independent of scalar Field.
    CHECK(br_symm == 3);
    CHECK(kr_symm == 1);
    CHECK(br_conj == 3);
    CHECK(kr_conj == 1);
    // (B) the pre-swapped form canonicalizes (by bra/ket color) to the same
    // orientation as (A), so the head is the same 3:1 — pre-swapping is no
    // longer a lever on head layout.
    CHECK(br_swap == 3);
    CHECK(kr_swap == 1);
    // (C) ResultExpr API gives the caller exact control; head matches the LHS
    // verbatim no matter how the RHS factors are oriented internally.
    CHECK(head_explicit.bra_rank() == 2);
    CHECK(head_explicit.ket_rank() == 2);
    CHECK(head_explicit.bra().at(0).label() == L"a_1");
    CHECK(head_explicit.bra().at(1).label() == L"a_2");
    CHECK(head_explicit.ket().at(0).label() == L"i_1");
    CHECK(head_explicit.ket().at(1).label() == L"i_2");
  }

  SECTION("real") {
    using ranges::views::transform;
    using sequant::EvalExprTA;
    using sequant::evaluate;
    using sequant::evaluate_antisymm;
    using sequant::evaluate_symm;

    using TArrayD = TA::TSpArrayD;

    auto parse_antisymm = [](auto const& xpr) {
      return sequant::deserialize<sequant::ExprPtr>(
          xpr, {.def_perm_symm = sequant::Symmetry::Antisymm});
    };

    auto& world = TA::get_default_world();

    const size_t nocc = 2, nvirt = 20;
    auto yield_ =
        rand_tensor_yield<double, TA::SparsePolicy>{world, nocc, nvirt};
    auto yield = [&yield_](std::wstring_view lbl) -> TA::TSpArrayD const& {
      return yield_(lbl)->get<TA::TSpArrayD>();
    };

    auto yield_d = [&yield_](std::wstring_view lbl) ->
        typename TA::TSpArrayD::numeric_type {
          return yield_(lbl)->get<typename TA::TSpArrayD::numeric_type>();
        };

    auto eval = [&yield_](sequant::ExprPtr const& expr,
                          std::string const& target_labels) {
      return evaluate(eval_node(expr), target_labels, yield_)
          ->get<TA::TSpArrayD>();
    };

    auto eval_symm = [&yield_](sequant::ExprPtr const& expr,
                               std::string const& target_labels) {
      return evaluate_symm(eval_node(expr), target_labels, yield_)
          ->get<TA::TSpArrayD>();
    };

    auto eval_antisymm = [&yield_](sequant::ExprPtr const& expr,
                                   std::string const& target_labels) {
      return evaluate_antisymm(eval_node(expr), target_labels, yield_)
          ->get<TA::TSpArrayD>();
    };

    auto eval_biorthogonal_nns_project_with =
        [](decltype(yield_)& yield, sequant::ExprPtr const& expr,
           std::string const& target_labels) {
          auto result = evaluate(eval_node(expr), target_labels, yield);
          return sequant::mbpt::biorthogonal_nns_project(
              result->get<TA::TSpArrayD>(),
              eval_node(expr)->as_tensor().bra_rank());
        };

    auto eval_biorthogonal_nns_project = [&](sequant::ExprPtr const& expr,
                                             std::string const& target_labels) {
      return eval_biorthogonal_nns_project_with(yield_, expr, target_labels);
    };

    SECTION("summation") {
      auto expr1 = parse_antisymm(L"t_{a1}^{i1} + f_{i1}^{a1}");

      auto sum1_eval = eval(expr1, "i_1,a_1");

      auto sum1_man = TArrayD{};
      sum1_man("i1,a1") =
          yield(L"t{a1;i1}")("a1,i1") + yield(L"f{i1;a1}")("i1,a1");

      REQUIRE(equal_tarrays(sum1_eval, sum1_man));

      auto expr2 = parse_antisymm(L"2 * t_{a1}^{i1} + 3/2 * f_{i1}^{a1}");

      auto sum2_eval = eval(expr2, "i_1,a_1");

      auto sum2_man = TArrayD{};
      sum2_man("i1,a1") =
          2 * yield(L"t{a1;i1}")("a1,i1") + 1.5 * yield(L"f{i1;a1}")("i1,a1");
      REQUIRE(equal_tarrays(sum2_eval, sum2_man));
    }

    SECTION("product") {
      auto expr1 =
          parse_antisymm(L"1/2 * g_{i2,i4}^{a2,a4} * t_{a1,a2}^{i1,i2}");
      auto prod1_eval = eval(expr1, "i_4,a_1,a_4,i_1");

      TArrayD prod1_man{};
      prod1_man("i4,a1,a4,i1") = 1 / 2.0 *
                                 yield(L"g{i2,i4;a2,a4}")("i2,i4,a2,a4") *
                                 yield(L"t{a1,a2;i1,i2}")("a1,a2,i1,i2");

      REQUIRE(equal_tarrays(prod1_eval, prod1_man));

      // the layout of a head: its indices as written, sorted or not
      const auto head =
          sequant::deserialize<sequant::ExprPtr>(L"X{i4,a1;a4,i1}")
              ->as<sequant::Tensor>();
      const auto in_head = evaluate(sequant::binarize<sequant::EvalExprTA>(
                                        sequant::ResultExpr(head, expr1)),
                                    sequant::EvalExprTA::layout(head), yield_)
                               ->get<TA::TSpArrayD>();
      REQUIRE(equal_tarrays(in_head, prod1_man));

      // the indices with protoindices are the modes of the inner tensors
      REQUIRE(
          sequant::EvalExprTA::layout(sequant::deserialize<sequant::ExprPtr>(
                                          L"X{a4<i2,i3>,i2;a1<i1,i2>,i1}")
                                          ->as<sequant::Tensor>()) ==
          "i_2,i_1;a_4i_2i_3,a_1i_1i_2");

      auto expr2 = parse_antisymm(
          L"-1/4 * g_{i3,i4}^{a3,a4} * t_{a2,a4}^{i1,i2} * t_{a1,a3}^{ i3, "
          L"i4}");
      auto prod2_eval = eval(expr2, "a_1,a_2,i_1,i_2");

      auto prod2_man = TArrayD{};
      prod2_man("a1,a2,i1,i2") = -1 / 4.0 *
                                 yield(L"g{i3,i4;a3,a4}")("i3,i4,a3,a4") *
                                 yield(L"t{a2,a4;i1,i2}")("a2,a4,i1,i2") *
                                 yield(L"t{a1,a3;i3,i4}")("a1,a3,i3,i4");

      // Unlike the two-tensor cases above, this three-tensor chain does NOT
      // reproduce the hand-written TA expression bit for bit: eval lays each
      // intermediate out in the node's canonical mode order, TA lays this
      // expression's intermediate out in the order the annotations impose, so
      // the two gemms accumulate in different orders. With nvirt = 20 the
      // final reduction runs over 80 terms of O(1) random values, so the
      // resulting last-bit disagreement is O(1e-13) in ABSOLUTE norm -- above
      // equal_tarrays' fixed margin of 100*eps (2.2e-14), which only ever
      // passed here because the two orders happened to coincide on some
      // platforms (it does not on Linux/g++ + reference BLAS). Compare
      // RELATIVELY, which is what the check is actually about: a wrong
      // contraction or a wrong permutation shows up at O(0.1), not O(1e-13).
      TArrayD prod2_diff;
      prod2_diff("a1,a2,i1,i2") =
          prod2_eval("a1,a2,i1,i2") - prod2_man("a1,a2,i1,i2");
      double const prod2_rel = TA::norm2(prod2_diff) / TA::norm2(prod2_man);
      INFO("prod2 relative L2 diff = " << prod2_rel);
      REQUIRE(prod2_rel < 1e-12);

      auto expr3 = sequant::deserialize<sequant::ExprPtr>(
          L"R_{a1}^{i1,i3} * f_{i3}^{i2}");
      auto prod3_eval = eval(expr3, "a_1,i_1,i_2");
      auto prod3_man = TArrayD{};
      prod3_man("a1,i1,i2") =
          yield(L"R{a1;i1,i3}")("a1,i1,i3") * yield(L"f{i3;i2}")("i3,i2");
      REQUIRE(equal_tarrays(prod3_eval, prod3_man));

      auto expr4 = sequant::deserialize<sequant::ExprPtr>(
          L"1/4 * R_{a1,a2,a3}^{i2,i3} * g_{i2,i3}^{i1,a3}");
      auto prod4_eval = eval(expr4, "i_1,a_1,a_2");
      auto prod4_man = TArrayD{};
      prod4_man("i1,a1,a2") = 1 / 4.0 *
                              yield(L"R{a1,a2,a3;i2,i3}")("a1,a2,a3,i2,i3") *
                              yield(L"g{i2,i3;i1,a3}")("i2,i3,i1,a3");
      REQUIRE(equal_tarrays(prod4_eval, prod4_man));
    }

    SECTION("sum and product") {
      auto expr1 = parse_antisymm(
          L"-1/4 * g_{i3,i4}^{a3,a4} * t_{a2,a4}^{i1,i2} * t_{a1,a3}^{i3,i4}"
          " + "
          " 1/16 * g_{i3,i4}^{a3,a4} * t_{a1,a2}^{i3,i4} * t_{a3,a4}^{i1,i2}");
      auto eval1 = eval(expr1, "a_1,a_2,i_1,i_2");

      auto man1 = TArrayD{};
      man1("a1,a2,i1,i2") = -1.0 / 4 * yield(L"g{i3,i4;a3,a4}")("i3,i4,a3,a4") *
                                yield(L"t{a2,a4;i1,i2}")("a2,a4,i1,i2") *
                                yield(L"t{a1,a3;i3,i4}")("a1,a3,i3,i4") +
                            1.0 / 16 * yield(L"g{i3,i4;a3,a4}")("i3,i4,a3,a4") *
                                yield(L"t{a1,a2;i3,i4}")("a1,a2,i3,i4") *
                                yield(L"t{a3,a4;i1,i2}")("a3,a4,i1,i2");
      // the canonical intermediate layouts reorder the summation vs the
      // hand-written reference: agreement to rounding, not bitwise
      REQUIRE(equal_tarrays<Loose>(eval1, man1));

      auto expr2 = sequant::deserialize<sequant::ExprPtr>(
          L"1/4 * R_{a1,a2,a3}^{i2,i3} * g_{i2,i3}^{i1,a3} + R_{a1,a3}^{i1} * "
          L"f_{i2}^{a3} * t_{a2}^{i2}");
      auto eval2 = eval(expr2, "i_1,a_1,a_2");

      auto man2 = TArrayD{};
      man2("i1,a1,a2") =
          1 / 4.0 * yield(L"R{a1,a2,a3;i2,i3}")("a1,a2,a3,i2,i3") *
              yield(L"g{i2,i3;i1,a3}")("i2,i3,i1,a3") +
          yield(L"R{a1,a3;i1}")("a1,a3,i1") * yield(L"f{i2;a3}")("i2,a3") *
              yield(L"t{a2;i2}")("a2,i2");
      REQUIRE(equal_tarrays(eval2, man2));
    }

    SECTION("variable at leaves") {
      auto expr2 =
          parse_antisymm(L"(α * 2 * t_{a1}^{i1} * β) + (3/2 * f_{i1}^{a1})");

      auto sum2_eval = eval(expr2, "i_1,a_1");

      auto sum2_man = TArrayD{};
      sum2_man("i1,a1") =
          yield_d(L"α") * 2 * yield(L"t{a1;i1}")("a1,i1") * yield_d(L"β") +
          1.5 * yield(L"f{i1;a1}")("i1,a1");

      REQUIRE(equal_tarrays(sum2_eval, sum2_man));
    }

    SECTION("power at leaves") {
      using sequant::Constant;
      using sequant::ex;
      using sequant::Power;
      using sequant::rational;
      using sequant::Variable;

      {
        auto pw = ex<Power>(ex<Constant>(2), rational{1, 2});
        auto t = parse_antisymm(L"t_{a1}^{i1}");
        auto expr = pw * t;

        auto eval1 = eval(expr, "i_1,a_1");

        auto man1 = TArrayD{};
        man1("i1,a1") = std::sqrt(2.0) * yield(L"t{a1;i1}")("a1,i1");
        REQUIRE(equal_tarrays(eval1, man1));
      }

      {
        auto pw = ex<Power>(ex<Variable>(L"α"), rational{2, 1});
        auto f = parse_antisymm(L"f_{i1}^{a1}");
        auto expr = pw * f;
        auto eval1 = eval(expr, "i_1,a_1");

        auto alpha_val = yield_d(L"α");
        auto man1 = TArrayD{};
        man1("i1,a1") = (alpha_val * alpha_val) * yield(L"f{i1;a1}")("i1,a1");
        REQUIRE(equal_tarrays(eval1, man1));
      }

      {
        auto pw = ex<Power>(ex<Constant>(2), rational{1, 2});
        auto t = parse_antisymm(L"t_{a1}^{i1}");
        auto alpha = ex<Variable>(L"α");
        auto f = parse_antisymm(L"f_{i1}^{a1}");

        auto expr = pw * t + alpha * f;

        auto eval1 = eval(expr, "i_1,a_1");

        auto man1 = TArrayD{};
        man1("i1,a1") = std::sqrt(2.0) * yield(L"t{a1;i1}")("a1,i1") +
                        yield_d(L"α") * yield(L"f{i1;a1}")("i1,a1");
        REQUIRE(equal_tarrays(eval1, man1));
      }

      {
        auto pw = ex<Power>(ex<Constant>(rational{1, 3}), rational{2, 1});
        auto g = parse_antisymm(L"g_{i1, i2}^{a1, a2}");
        auto expr = pw * g;

        auto eval1 = eval(expr, "i_1,i_2,a_1,a_2");

        auto man1 = TArrayD{};
        man1("i1,i2,a1,a2") =
            std::pow(1.0 / 3.0, 2.0) * yield(L"g{i1,i2;a1,a2}")("i1,i2,a1,a2");
        REQUIRE(equal_tarrays(eval1, man1));
      }
    }

    SECTION("Antisymmetrization") {
      auto expr1 = parse_antisymm(L"g_{i1, i2}^{a1, a2}");
      auto eval1 = eval_antisymm(expr1, "i_1,i_2,a_1,a_2");
      auto const& arr1 = yield(L"g{i1,i2;a1,a2}");

      auto man1 = TArrayD{};
      man1("0,1,2,3") =
          arr1("0,1,2,3") - arr1("1,0,2,3") + arr1("1,0,3,2") - arr1("0,1,3,2");

      man1("0,1,2,3") = 0.25 * man1("0,1,2,3");

      REQUIRE(equal_tarrays(eval1, man1));

      // odd-ranked tensor
      auto expr2 = parse_antisymm(L"g_{i1, i2, i3}^{a1, a2}");
      auto eval2 = eval_antisymm(expr2, "i_1,i_2,i_3,a_1,a_2");
      auto const& arr2 = yield(L"g{i1,i2,i3;a1,a2}");

      auto man2 = TArrayD{};
      man2("0,1,2,3,4") =
          arr2("0,1,2,3,4") - arr2("1,0,2,3,4") + arr2("1,2,0,3,4") -
          arr2("2,1,0,3,4") + arr2("2,0,1,3,4") - arr2("0,2,1,3,4") -
          arr2("0,1,2,4,3") + arr2("1,0,2,4,3") - arr2("1,2,0,4,3") +
          arr2("2,1,0,4,3") - arr2("2,0,1,4,3") + arr2("0,2,1,4,3");

      REQUIRE(equal_tarrays(eval2, man2));

      auto expr3 = parse_antisymm(L"R_{a1,a2}^{}");
      auto eval3 = eval_antisymm(expr3, "a_1,a_2");
      auto const& arr3 = yield(L"R{a1,a2;}");
      auto man3 = TArrayD{};
      man3("0,1") = arr3("0,1") - arr3("1,0");
      man3("0,1") = 0.5 * man3("0,1");

      REQUIRE(equal_tarrays(eval3, man3));
    }

    SECTION("Symmetrization") {
      auto expr1 = parse_antisymm(L"g_{i1, i2}^{a1, a2}");
      auto eval1 = eval_symm(expr1, "i_1,i_2,a_1,a_2");
      auto const& arr1 = yield(L"g{i1,i2;a1,a2}");

      auto man1 = TArrayD{};
      man1("0,1,2,3") = arr1("0,1,2,3") + arr1("1,0,3,2");
      man1("0,1,2,3") = 0.5 * man1("0,1,2,3");

      REQUIRE(equal_tarrays(eval1, man1));

      auto expr2 = parse_antisymm(L"g_{i1,i2,i3}^{a1,a2,a3}");

      auto eval2 = eval_symm(expr2, "i_1,i_2,i_3,a_1,a_2,a_3");
      auto const& arr2 = yield(L"g{i1,i2,i3;a1,a2,a3}");
      TArrayD man2;
      man2("0,1,2,3,4,5") = arr2("0,1,2,3,4,5") + arr2("0,2,1,3,5,4") +
                            arr2("2,0,1,5,3,4") + arr2("2,1,0,5,4,3") +
                            arr2("1,2,0,4,5,3") + arr2("1,0,2,4,3,5");
      man2("0,1,2,3,4,5") = (1.0 / 6.0) * man2("0,1,2,3,4,5");
      REQUIRE(equal_tarrays(eval2, man2));
    }

    SECTION("Biorthogonal Cleanup") {
      // low-rank residuals: skip nns
      auto expr1 = parse_antisymm(L"R_{a1, a2}^{i1, i2}");
      auto eval1 = eval_biorthogonal_nns_project(expr1, "a_1,a_2,i_1,i_2");
      auto const& arr1 = yield(L"R{a1,a2;i1,i2}");

      auto man1 = TArrayD{};
      man1("0,1,2,3") = arr1("0,1,2,3");

      REQUIRE(equal_tarrays(eval1, man1));

      // for rank 3 residual, nns applies:
      // result = NNS_P * sum_of_ket_permutations
      auto expr2 = parse_antisymm(L"R_{a1, a2, a3}^{i1, i2, i3}");
      auto eval2 =
          eval_biorthogonal_nns_project(expr2, "a_1,a_2,a_3,i_1,i_2,i_3");
      auto const& arr2 = yield(L"R{a1,a2,a3;i1,i2,i3}");

      auto man2 = TArrayD{};
      man2("0,1,2,3,4,5") =
          (5.0 / 6.0) * arr2("0,1,2,3,4,5") -
          (1.0 / 6.0) *
              (arr2("0,1,2,3,5,4") + arr2("0,1,2,4,3,5") + arr2("0,1,2,4,5,3") +
               arr2("0,1,2,5,3,4") + arr2("0,1,2,5,4,3"));

      REQUIRE(equal_tarrays(eval2, man2));

      // for rank 4 residual, nns applies:
      // result = NNS_P * sum_of_ket_permutations
      // (own extents: nocc >= rank, else every ket tuple repeats an index and
      // the antisymmetric part of the projector goes unchecked; small nvirt,
      // since the projection acts on the kets only)
      auto yield4_ = rand_tensor_yield<double, TA::SparsePolicy>{world, 4, 2};
      auto expr3 = parse_antisymm(L"R_{a1, a2, a3, a4}^{i1, i2, i3, i4}");
      auto eval3 = eval_biorthogonal_nns_project_with(
          yield4_, expr3, "a_1,a_2,a_3,a_4,i_1,i_2,i_3,i_4");
      auto const& arr3 =
          yield4_(L"R{a1,a2,a3,a4;i1,i2,i3,i4}")->get<TA::TSpArrayD>();

      auto man3 = TArrayD{};
      man3("0,1,2,3,4,5,6,7") = 14.0 / 24.0 * arr3("0,1,2,3,4,5,6,7") +
                                -4.0 / 24.0 * arr3("0,1,2,3,4,5,7,6") +
                                -4.0 / 24.0 * arr3("0,1,2,3,4,6,5,7") +
                                -1.0 / 24.0 * arr3("0,1,2,3,4,6,7,5") +
                                -1.0 / 24.0 * arr3("0,1,2,3,4,7,5,6") +
                                -4.0 / 24.0 * arr3("0,1,2,3,4,7,6,5") +
                                -4.0 / 24.0 * arr3("0,1,2,3,5,4,6,7") +
                                2.0 / 24.0 * arr3("0,1,2,3,5,4,7,6") +
                                -1.0 / 24.0 * arr3("0,1,2,3,5,6,4,7") +
                                2.0 / 24.0 * arr3("0,1,2,3,5,6,7,4") +
                                2.0 / 24.0 * arr3("0,1,2,3,5,7,4,6") +
                                -1.0 / 24.0 * arr3("0,1,2,3,5,7,6,4") +
                                -1.0 / 24.0 * arr3("0,1,2,3,6,4,5,7") +
                                2.0 / 24.0 * arr3("0,1,2,3,6,4,7,5") +
                                -4.0 / 24.0 * arr3("0,1,2,3,6,5,4,7") +
                                -1.0 / 24.0 * arr3("0,1,2,3,6,5,7,4") +
                                2.0 / 24.0 * arr3("0,1,2,3,6,7,4,5") +
                                2.0 / 24.0 * arr3("0,1,2,3,6,7,5,4") +
                                2.0 / 24.0 * arr3("0,1,2,3,7,4,5,6") +
                                -1.0 / 24.0 * arr3("0,1,2,3,7,4,6,5") +
                                -1.0 / 24.0 * arr3("0,1,2,3,7,5,4,6") +
                                -4.0 / 24.0 * arr3("0,1,2,3,7,5,6,4") +
                                2.0 / 24.0 * arr3("0,1,2,3,7,6,4,5") +
                                2.0 / 24.0 * arr3("0,1,2,3,7,6,5,4");

      REQUIRE(equal_tarrays<Loose>(eval3, man3));
    }

    SECTION("Others") {
      using namespace std::string_literals;
      auto expr1 = parse_antisymm(
          L"-1/4 * g_{i3,i4}^{a3,a4} * t_{a2,a4}^{i1,i2} * t_{a1,a3}^{i3,i4}"
          " + "
          " 1/16 * g_{i3,i4}^{a3,a4} * t_{a1,a2}^{i3,i4} * t_{a3,a4}^{i1,i2}");

      auto eval1 = evaluate(eval_node(expr1), "i_1,i_2,a_1,a_2"s, yield_)
                       ->get<TA::TSpArrayD>();

      auto nodes1 = *expr1 | ranges::views::transform([](auto&& x) {
        return eval_node(x);
      }) | ranges::to_vector;

      auto eval2 =
          evaluate(nodes1, "i_1,i_2,a_1,a_2"s, yield_)->get<TA::TSpArrayD>();

      REQUIRE(equal_tarrays(eval1, eval2));
    }

    SECTION("non-covariant indices") {
      using sequant::deserialize;
      using sequant::EvalExprTA;
      using sequant::evaluate;

      auto& world = TA::get_default_world();
      const size_t nocc = 2, nvirt = 4, naux = 12;

      auto yield_ =
          rand_tensor_yield<double, TA::SparsePolicy>{world, nocc, nvirt, naux};
      auto yield = [&yield_](std::wstring_view lbl) -> TA::TSpArrayD const& {
        return yield_(lbl)->get<TA::TSpArrayD>();
      };

      auto eval = [&yield_](sequant::ExprPtr const& expr,
                            std::string const& target_labels) {
        return evaluate(eval_node(expr), target_labels, yield_)
            ->get<TA::TSpArrayD>();
      };

      auto expr1 = deserialize(
          L"((X{a1;;x1} X{;a2;x1}) Y{;;x1,x2})(X{a3;;x2} X{;a4;x2})");
      auto eval1 = eval(expr1, "a_1,a_2,a_3,a_4");
      auto man1 = [&]() {
        auto X1 = yield(L"X{a1;;x1}");
        REQUIRE(X1.trange().elements_range().extent(0) == nvirt);
        REQUIRE(X1.trange().elements_range().extent(1) == naux);
        auto X2 = yield(L"X{;a2;x1}");
        REQUIRE(X2.trange().elements_range().extent(0) == nvirt);
        REQUIRE(X2.trange().elements_range().extent(1) == naux);
        auto X3 = yield(L"X{a3;;x2}");
        REQUIRE(X3.trange().elements_range().extent(0) == nvirt);
        REQUIRE(X3.trange().elements_range().extent(1) == naux);
        auto X4 = yield(L"X{;a4;x2}");
        REQUIRE(X4.trange().elements_range().extent(0) == nvirt);
        REQUIRE(X4.trange().elements_range().extent(1) == naux);
        auto Y = yield(L"Y{;;x1,x2}");
        REQUIRE(Y.trange().elements_range().extent(0) == naux);
        REQUIRE(Y.trange().elements_range().extent(1) == naux);
        auto X12 = TA::einsum("ax,bx->abx", X1, X2);
        REQUIRE(X12.trange().elements_range().extent(0) == nvirt);
        REQUIRE(X12.trange().elements_range().extent(1) == nvirt);
        REQUIRE(X12.trange().elements_range().extent(2) == naux);
        auto X12Y = TA::einsum("abx,xy->aby", X12, Y);
        auto X34 = TA::einsum("cy,dy->cdy", X3, X4);
        return TA::einsum("aby,cdy->abcd", X12Y, X34);
      }();
      REQUIRE(equal_tarrays(eval1, man1, "a1,a2,a3,a4"));

      // cluster-specific RDM: γ{a2;a1;i1,i2} = t{i1,i2;a1,a3} T2{a2,a3;i1,i2}
      // i1,i2 form standard bra-ket pairs across factors but are external
      // (in aux of result); binarize(ResultExpr) must keep them uncontracted
      {
        auto res = deserialize<sequant::ResultExpr>(
            L"GAM{a2;a1;i1,i2} = t{i1,i2;a1,a3} T2{a2,a3;i1,i2}");
        auto node = eval_node(res);
        auto eval_rdm = evaluate(node, std::string("a_2,a_1,i_1,i_2"), yield_)
                            ->get<TA::TSpArrayD>();
        auto man_rdm = [&]() {
          auto t = yield(L"t{i1,i2;a1,a3}");
          auto T2 = yield(L"T2{a2,a3;i1,i2}");
          return TA::einsum("ijab,cbij->caij", t, T2);
        }();
        REQUIRE(equal_tarrays(eval_rdm, man_rdm, "a2,a1,i1,i2"));
      }

      {  // multiple bra or ket indices require AssertStrictBraKetSymmetry::No
        auto ctx_resetter = sequant::set_scoped_default_context(
            sequant::Context{sequant::get_default_context()}.set(
                sequant::AssertStrictBraKetSymmetry::No));

        // hyperindex i1 in ket slots of 3 tensors
        auto expr2 = deserialize(L"T{a1;i1} T{a2;i1} T{a3;i1}");
        auto eval2 = eval(expr2, "a_1,a_2,a_3");
        auto man2 = [&]() {
          auto T1 = yield(L"T{a1;i1}");
          auto T2 = yield(L"T{a2;i1}");
          auto T3 = yield(L"T{a3;i1}");
          auto T12 = TA::einsum("ai,bi->abi", T1, T2);
          return TA::einsum("abi,ci->abc", T12, T3);
        }();
        REQUIRE(equal_tarrays(eval2, man2, "a1,a2,a3"));

        // hyperindex a1 in bra slots of 3 tensors
        auto expr3 = deserialize(L"T{a1;i1} T{a1;i2} T{a1;i3}");
        auto eval3 = eval(expr3, "i_1,i_2,i_3");
        auto man3 = [&]() {
          auto T1 = yield(L"T{a1;i1}");
          auto T2 = yield(L"T{a1;i2}");
          auto T3 = yield(L"T{a1;i3}");
          auto T12 = TA::einsum("ai,aj->aij", T1, T2);
          return TA::einsum("aij,ak->ijk", T12, T3);
        }();
        REQUIRE(equal_tarrays(eval3, man3, "i1,i2,i3"));
      }

      {  // high-order aux hyperindex carried into the result: aux index x1 is
         // shared by 3 tensors AND external (in the result aux slot). This is
         // the >2-tensor named-hyperedge case the TNv3 canonicalize_slots fix
         // unlocks; TA evaluates the resulting batched-over-x1 contraction
         // natively (einsum keeps an index common to both operands and result).
         //   R{;;x1} = A{;a1;x1} B{a1;a2;x1} C{a2;;x1}
         //   R[x] = sum_{a1,a2} A[a1,x] B[a1,a2,x] C[a2,x]
        auto res = deserialize<sequant::ResultExpr>(
            L"R{;;x1} = A{;a1;x1} B{a1;a2;x1} C{a2;;x1}");
        auto evalR = evaluate(eval_node(res), std::string("x_1"), yield_)
                         ->get<TA::TSpArrayD>();
        auto manR = [&]() {
          auto A = yield(L"A{;a1;x1}");    // [a1, x]
          auto B = yield(L"B{a1;a2;x1}");  // [a1, a2, x]
          auto C = yield(L"C{a2;;x1}");    // [a2, x]
          auto AB = TA::einsum("ax,abx->bx", A, B);
          return TA::einsum("bx,bx->x", AB, C);
        }();
        REQUIRE(equal_tarrays(evalR, manR, "x1"));
      }
    }  // multiple bra or ket indices
  }

  SECTION("complex") {
    using TArrayC =
        TA::DistArray<TA::Tensor<std::complex<double>>, TA::SparsePolicy>;

    const size_t nocc = 2, nvirt = 20;
    auto& world = TA::get_default_world();

    auto yield_ = rand_tensor_yield<std::complex<double>, TA::SparsePolicy>{
        world, nocc, nvirt};

    auto yield = [&yield_](std::wstring_view lbl) -> TArrayC const& {
      return yield_(lbl)->get<TArrayC>();
    };

    auto eval = [&yield_](sequant::ExprPtr const& expr,
                          std::string const& target_labels) {
      return evaluate(eval_node(expr), target_labels, yield_)->get<TArrayC>();
    };

    auto eval_symm = [&yield_](sequant::ExprPtr const& expr,
                               std::string const& target_labels) {
      return evaluate_symm(eval_node(expr), target_labels, yield_)
          ->get<TArrayC>();
    };

    auto eval_antisymm = [&yield_](sequant::ExprPtr const& expr,
                                   std::string const& target_labels) {
      return evaluate_antisymm(eval_node(expr), target_labels, yield_)
          ->get<TArrayC>();
    };

    using namespace sequant;
    using namespace std::string_literals;

    SECTION("summation") {
      auto expr1 = deserialize<sequant::ExprPtr>(L"t_{a1}^{i1} + f_{i1}^{a1}");

      auto sum1_eval = eval(expr1, "i_1,a_1");

      auto sum1_man = TArrayC{};
      sum1_man("i1,a1") =
          yield(L"t{a1;i1}")("a1,i1") + yield(L"f{i1;a1}")("i1,a1");

      REQUIRE(equal_tarrays(sum1_eval, sum1_man));

      auto expr2 =
          deserialize<sequant::ExprPtr>(L"2 * t_{a1}^{i1} + 3/2 * f_{i1}^{a1}");

      auto sum2_eval = eval(expr2, "i_1,a_1");

      auto sum2_man = TArrayC{};
      sum2_man("i1,a1") =
          std::complex<double>{2} * yield(L"t{a1;i1}")("a1,i1") +
          std::complex<double>{1.5} * yield(L"f{i1;a1}")("i1,a1");

      REQUIRE(equal_tarrays(sum2_eval, sum2_man));
    }

    SECTION("product") {
      auto expr1 = deserialize<sequant::ExprPtr>(
          L"1/2 * g_{i2,i4}^{a2,a4} * t_{a1,a2}^{i1,i2}");
      auto prod1_eval = eval(expr1, "i_4,a_1,a_4,i_1");

      TArrayC prod1_man{};
      prod1_man("i4,a1,a4,i1") = std::complex<double>{1 / 2.0} *
                                 yield(L"g{i2,i4;a2,a4}")("i2,i4,a2,a4") *
                                 yield(L"t{a1,a2;i1,i2}")("a1,a2,i1,i2");

      REQUIRE(equal_tarrays(prod1_eval, prod1_man));

      auto expr2 = deserialize<sequant::ExprPtr>(
          L"-1/4 * g_{i3,i4}^{a3,a4} * t_{a2,a4}^{i1,i2} * t_{a1,a3}^{ i3, "
          L"i4}");
      auto prod2_eval = eval(expr2, "a_1,a_2,i_1,i_2");

      auto prod2_man = TArrayC{};
      prod2_man("a1,a2,i1,i2") = std::complex<double>{-1 / 4.0} *
                                 yield(L"g{i3,i4;a3,a4}")("i3,i4,a3,a4") *
                                 yield(L"t{a2,a4;i1,i2}")("a2,a4,i1,i2") *
                                 yield(L"t{a1,a3;i3,i4}")("a1,a3,i3,i4");

      REQUIRE(equal_tarrays(prod2_eval, prod2_man));

      auto expr3 = sequant::deserialize<sequant::ExprPtr>(
          L"R_{a1}^{i1,i3} * f_{i3}^{i2}");
      auto prod3_eval = eval(expr3, "a_1,i_1,i_2");
      auto prod3_man = TArrayC{};
      prod3_man("a1,i1,i2") =
          yield(L"R{a1;i1,i3}")("a1,i1,i3") * yield(L"f{i3;i2}")("i3,i2");

      REQUIRE(equal_tarrays(prod3_eval, prod3_man));

      auto expr4 = sequant::deserialize<sequant::ExprPtr>(
          L"1/4 * R_{a1,a2,a3}^{i2,i3} * g_{i2,i3}^{i1,a3}");
      auto prod4_eval = eval(expr4, "i_1,a_1,a_2");
      auto prod4_man = TArrayC{};
      prod4_man("i1,a1,a2") = 1 / 4.0 *
                              yield(L"R{a1,a2,a3;i2,i3}")("a1,a2,a3,i2,i3") *
                              yield(L"g{i2,i3;i1,a3}")("i2,i3,i1,a3");
      REQUIRE(equal_tarrays(prod4_eval, prod4_man));
    }

    SECTION("sum and product") {
      auto expr1 = deserialize<sequant::ExprPtr>(
          L"-1/4 * g_{i3,i4}^{a3,a4} * t_{a2,a4}^{i1,i2} * t_{a1,a3}^{i3,i4}"
          " + "
          " 1/16 * g_{i3,i4}^{a3,a4} * t_{a1,a2}^{i3,i4} * t_{a3,a4}^{i1,i2}");
      auto eval1 = eval(expr1, "a_1,a_2,i_1,i_2");

      auto man1 = TArrayC{};
      man1("a1,a2,i1,i2") = std::complex<double>{-1.0 / 4} *
                                yield(L"g{i3,i4;a3,a4}")("i3,i4,a3,a4") *
                                yield(L"t{a2,a4;i1,i2}")("a2,a4,i1,i2") *
                                yield(L"t{a1,a3;i3,i4}")("a1,a3,i3,i4") +
                            std::complex<double>{1.0 / 16} *
                                yield(L"g{i3,i4;a3,a4}")("i3,i4,a3,a4") *
                                yield(L"t{a1,a2;i3,i4}")("a1,a2,i3,i4") *
                                yield(L"t{a3,a4;i1,i2}")("a3,a4,i1,i2");

      REQUIRE(equal_tarrays(eval1, man1));

      auto expr2 = sequant::deserialize<sequant::ExprPtr>(
          L"1/4 * R_{a1,a2,a3}^{i2,i3} * g_{i2,i3}^{i1,a3} + R_{a1,a3}^{i1} * "
          L"f_{i2}^{a3} * t_{a2}^{i2}");
      auto eval2 = eval(expr2, "i_1,a_1,a_2");

      auto man2 = TArrayC{};
      man2("i1,a1,a2") =
          1 / 4.0 * yield(L"R{a1,a2,a3;i2,i3}")("a1,a2,a3,i2,i3") *
              yield(L"g{i2,i3;i1,a3}")("i2,i3,i1,a3") +
          yield(L"R{a1,a3;i1}")("a1,a3,i1") * yield(L"f{i2;a3}")("i2,a3") *
              yield(L"t{a2;i2}")("a2,i2");
      REQUIRE(equal_tarrays(eval2, man2));
    }

    SECTION("Antisymmetrization") {
      auto expr1 = deserialize<sequant::ExprPtr>(L"g_{i1, i2}^{a1, a2}");
      auto eval1 = eval_antisymm(expr1, "i_1,i_2,a_1,a_2");
      auto const& arr1 = yield(L"g{i1,i2;a1,a2}");

      auto man1 = TArrayC{};
      man1("0,1,2,3") =
          arr1("0,1,2,3") - arr1("1,0,2,3") + arr1("1,0,3,2") - arr1("0,1,3,2");

      man1("0,1,2,3") = std::complex<double>{0.25} * man1("0,1,2,3");

      REQUIRE(equal_tarrays(eval1, man1));

      // odd-ranked tensor
      auto expr2 = deserialize<sequant::ExprPtr>(L"g_{i1, i2, i3}^{a1, a2}");
      auto eval2 = eval_antisymm(expr2, "i_1,i_2,i_3,a_1,a_2");
      auto const& arr2 = yield(L"g{i1,i2,i3;a1,a2}");

      auto man2 = TArrayC{};
      man2("0,1,2,3,4") =
          arr2("0,1,2,3,4") - arr2("1,0,2,3,4") + arr2("1,2,0,3,4") -
          arr2("2,1,0,3,4") + arr2("2,0,1,3,4") - arr2("0,2,1,3,4") -
          arr2("0,1,2,4,3") + arr2("1,0,2,4,3") - arr2("1,2,0,4,3") +
          arr2("2,1,0,4,3") - arr2("2,0,1,4,3") + arr2("0,2,1,4,3");

      REQUIRE(equal_tarrays(eval2, man2));

      auto expr3 = deserialize<sequant::ExprPtr>(L"R_{a1,a2}^{}");
      auto eval3 = eval_antisymm(expr3, "a_1,a_2");
      auto const& arr3 = yield(L"R{a1,a2;}");
      auto man3 = TArrayC{};
      man3("0,1") = arr3("0,1") - arr3("1,0");
      man3("0,1") = std::complex<double>{0.5} * man3("0,1");

      REQUIRE(equal_tarrays(eval3, man3));
    }

    SECTION("Symmetrization") {
      auto expr1 = deserialize<sequant::ExprPtr>(L"g_{i1, i2}^{a1, a2}");
      auto eval1 = eval_symm(expr1, "i_1,i_2,a_1,a_2");
      auto const& arr1 = yield(L"g{i1,i2;a1,a2}");

      auto man1 = TArrayC{};
      man1("0,1,2,3") = arr1("0,1,2,3") + arr1("1,0,3,2");
      man1("0,1,2,3") = 0.5 * man1("0,1,2,3");

      REQUIRE(equal_tarrays(eval1, man1));

      auto expr2 = deserialize<sequant::ExprPtr>(L"g_{i1,i2,i3}^{a1,a2,a3}");

      auto eval2 = eval_symm(expr2, "i_1,i_2,i_3,a_1,a_2,a_3");
      auto const& arr2 = yield(L"g{i1,i2,i3;a1,a2,a3}");
      TArrayC man2;
      man2("0,1,2,3,4,5") = arr2("0,1,2,3,4,5") + arr2("0,2,1,3,5,4") +
                            arr2("2,0,1,5,3,4") + arr2("2,1,0,5,4,3") +
                            arr2("1,2,0,4,5,3") + arr2("1,0,2,4,3,5");
      man2("0,1,2,3,4,5") = (1.0 / 6.0) * man2("0,1,2,3,4,5");

      REQUIRE(equal_tarrays(eval2, man2));
    }

    SECTION("Others") {
      using namespace std::string_literals;
      auto expr1 = deserialize<sequant::ExprPtr>(
          L"-1/4 * g_{i3,i4}^{a3,a4} * t_{a2,a4}^{i1,i2} * t_{a1,a3}^{i3,i4}"
          " + "
          " 1/16 * g_{i3,i4}^{a3,a4} * t_{a1,a2}^{i3,i4} * t_{a3,a4}^{i1,i2}");

      auto eval1 = evaluate(eval_node(expr1), "i_1,i_2,a_1,a_2"s, yield_)
                       ->get<TArrayC>();

      auto nodes1 = *expr1 | ranges::views::transform([](auto&& x) {
        return eval_node(x);
      }) | ranges::to_vector;

      auto eval2 = evaluate(nodes1, "i_1,i_2,a_1,a_2"s, yield_)->get<TArrayC>();

      REQUIRE(equal_tarrays(eval1, eval2));
    }
  }

  SECTION("tot") {
    using namespace sequant;

    //
    // eg: approx_equal("i,j;a,b", arr1, arr2)
    // - arr1 and arr2 are DistArrays with equal TiledRange and matching Range
    // for
    //   the inner tensors at the corresponding tile positions.
    // - 'i', 'j', 'a', and 'b' are dummy indices that annotate the modes of
    // outer
    //   and inner tensors. Why? Because TA::norm2 function is not supported for
    //   tensor-of-tensor tiles
    //
    auto approx_equal = [](std::string const& annot, auto const& lhs,
                           auto const& rhs) -> bool {
      return Catch::Approx(lhs(annot).dot(lhs(annot))) ==
             rhs(annot).dot(rhs(annot));
    };

    auto& world = TA::get_default_world();

    size_t const nocc = 2;
    size_t const nvirt = 3;

    rand_tensor_yield<double> yield{world, nocc, nvirt};

    using ArrayT = typename decltype(yield)::array_type;
    using ArrayToT = typename decltype(yield)::array_tot_type;
    using NumericT = typename decltype(yield)::numeric_type;

    SECTION("T_times_ToT_to_ToT") {
      constexpr std::wstring_view expr_str =
          L"3"
          L" * "
          L"f{i3;i1}"
          L" * "
          L"t{a3<i2,i3>,a4<i2,i3>;i2,i3}";
      auto const node = eval_node(deserialize<sequant::ExprPtr>(expr_str));
      std::string const target_layout{"i_1,i_2,i_3;a_3i_2i_3,a_4i_2i_3"};
      auto result = evaluate(node, target_layout, yield)->get<ArrayToT>();
      ArrayToT ref;
      {
        auto const& lhs = yield(L"f{i3;i1}")->get<ArrayT>();
        auto const& rhs =
            yield(L"t{a3<i2,i3>,a4<i2,i3>;i2,i3}")->get<ArrayToT>();
        ref = TA::einsum(lhs("i_3,i_1"), rhs("i_2,i_3;a_3i_2i_3,a_4i_2i_3"),
                         target_layout);
        ref(target_layout) = 3 * ref(target_layout);
      }
      REQUIRE(approx_equal("i,j,k;a,b", result, ref));
    }

    SECTION("ToT_times_ToT_to_ToT") {
      constexpr std::wstring_view expr_str =
          L"I{a4<i2,i3>,a1<i1,i2>;i1,i2}"
          L" * "
          L"s{a2<i1,i2>;a4<i2,i3>}";

      auto const node = eval_node(deserialize<sequant::ExprPtr>(expr_str));
      // The contracted composite a_4<i_2,i_3> leaves the result a function
      // of i_3 (the pair basis), so i_3 is an outer mode of the result: use
      // the node's own canonical layout, and keep i_3 in the reference.
      std::string const target_layout = node->annot();

      auto result = evaluate(node, target_layout, yield)->get<ArrayToT>();

      ArrayToT ref;
      {
        auto const& lhs =
            yield(L"I{a4<i2,i3>,a1<i1,i2>;i1,i2}")->get<ArrayToT>();
        auto const& rhs = yield(L"s{a2<i1,i2>;a4<i2,i3>}")->get<ArrayToT>();
        ref = TA::einsum(lhs("i_1,i_2,i_3;a_4i_2i_3,a_1i_1i_2"),
                         rhs("i_1,i_2,i_3;a_2i_1i_2,a_4i_2i_3"), target_layout);
      }
      REQUIRE(approx_equal("i,j,k;a,b", result, ref));
    }

    SECTION("ToT_times_ToT_to_Scalar") {
      constexpr std::wstring_view expr_str =
          L"I{a1<i1,i2>,a2<i1,i2>;i1,i2}"
          L" * "
          L"g{i1,i2;a2<i1,i2>,a1<i1,i2>}";
      auto const node = eval_node(deserialize<sequant::ExprPtr>(expr_str));

      auto result = evaluate(node, yield)->get<NumericT>();

      NumericT ref;
      {
        auto const& lhs =
            yield(L"I{a1<i1,i2>,a2<i1,i2>;i1,i2}")->get<ArrayToT>();
        auto const& rhs =
            yield(L"g{i1,i2;a2<i1,i2>,a1<i1,i2>}")->get<ArrayToT>();
        ref = TA::dot(lhs("i_1,i_2;a_1i_1i_2,a_2i_1i_2"),
                      rhs("i_1,i_2;a_2i_1i_2,a_1i_1i_2"));
      }
      REQUIRE(result == Catch::Approx(ref));
    }

    SECTION("symmetrize") {
      // n=2 (logical rank 4): outer = occupied (i1,i2), inner = virtual
      // (a1,a2). symmetrize() must produce S(i,j;a,b) = 1/2! * (R(i,j;a,b) +
      // R(j,i;b,a))
      // -- outer and inner modes permute in lockstep.
      {
        auto const Rnode = eval_node(
            deserialize<sequant::ExprPtr>(L"R{a1<i1,i2>,a2<i1,i2>;i1,i2}"));
        auto const Rres = yield(Rnode);
        auto const& R = Rres->get<ArrayToT>();

        auto const symm = Rres->symmetrize()->get<ArrayToT>();

        std::string const annot{"i_1,i_2;a_1i_1i_2,a_2i_1i_2"};
        ArrayToT ref;
        ref(annot) = R(annot) + R("i_2,i_1;a_2i_1i_2,a_1i_1i_2");
        ref(annot) = 0.5 * ref(annot);

        ArrayToT diff;
        diff(annot) = symm(annot) - ref(annot);
        // Norm-squared of the difference; an absolute margin is needed because
        // Approx(0.0) is a pure relative test. The reference sums the same
        // permutations in a different order than the production helper, so a
        // tiny rounding residual (~1e-30) is expected; a real mismatch is O(1).
        REQUIRE(diff(annot).dot(diff(annot)) ==
                Catch::Approx(0.0).margin(1e-12));
      }

      // n=3 (logical rank 6): exercises the 6-permutation / 1/3! path. The 6
      // particle permutations permute outer (i1,i2,i3) and inner (a1,a2,a3) in
      // lockstep; each inner index keeps its proto-suffix "i_1i_2i_3".
      {
        auto const Rnode = eval_node(deserialize<sequant::ExprPtr>(
            L"R{a1<i1,i2,i3>,a2<i1,i2,i3>,a3<i1,i2,i3>;i1,i2,i3}"));
        auto const Rres = yield(Rnode);
        auto const& R = Rres->get<ArrayToT>();

        auto const symm = Rres->symmetrize()->get<ArrayToT>();

        std::string const annot{
            "i_1,i_2,i_3;a_1i_1i_2i_3,a_2i_1i_2i_3,a_3i_1i_2i_3"};
        ArrayToT ref;
        ref(annot) = R("i_1,i_2,i_3;a_1i_1i_2i_3,a_2i_1i_2i_3,a_3i_1i_2i_3") +
                     R("i_1,i_3,i_2;a_1i_1i_2i_3,a_3i_1i_2i_3,a_2i_1i_2i_3") +
                     R("i_3,i_1,i_2;a_3i_1i_2i_3,a_1i_1i_2i_3,a_2i_1i_2i_3") +
                     R("i_3,i_2,i_1;a_3i_1i_2i_3,a_2i_1i_2i_3,a_1i_1i_2i_3") +
                     R("i_2,i_3,i_1;a_2i_1i_2i_3,a_3i_1i_2i_3,a_1i_1i_2i_3") +
                     R("i_2,i_1,i_3;a_2i_1i_2i_3,a_1i_1i_2i_3,a_3i_1i_2i_3");
        ref(annot) = (1.0 / 6.0) * ref(annot);

        ArrayToT diff;
        diff(annot) = symm(annot) - ref(annot);
        // Norm-squared of the difference; an absolute margin is needed because
        // Approx(0.0) is a pure relative test. The reference sums the same
        // permutations in a different order than the production helper, so a
        // tiny rounding residual (~1e-30) is expected; a real mismatch is O(1).
        REQUIRE(diff(annot).dot(diff(annot)) ==
                Catch::Approx(0.0).margin(1e-12));
      }
    }
  }
}

TEST_CASE("eval_tensor_power_ta", "[eval]") {
  using namespace sequant;
  auto& world = TA::get_default_world();
  const size_t nocc = 2, nvirt = 3;

  const TensorSymmetries dsym{.perm = Symmetry::Symm,
                              .braket = BraKetSymmetry::Symm,
                              .column = ColumnSymmetry::Symm};
  auto D = ex<Tensor>(L"D", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"}, dsym);
  auto M = ex<Tensor>(L"M", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"},
                      TensorSymmetries{.perm = Symmetry::Antisymm});
  std::string const target = "i_1,i_2,a_1,a_2";

  // the elements of an array in local-tile order; arrays evaluated into the
  // same layout from the same yield share a tiling, so these line up
  auto elements = [](auto const& arr) {
    using T = typename std::remove_cvref_t<decltype(arr)>::numeric_type;
    std::vector<T> out;
    for (auto it = arr.begin(); it != arr.end(); ++it) {
      auto const tile = it->get();
      out.insert(out.end(), tile.begin(), tile.end());
    }
    return out;
  };
  auto check_elementwise = [](auto const& got, auto const& expected) {
    REQUIRE(got.size() == expected.size());
    REQUIRE(!got.empty());
    for (std::size_t k = 0; k < got.size(); ++k) {
      CHECK(std::real(got[k]) ==
            Catch::Approx(std::real(expected[k])).margin(1e-12));
      CHECK(std::imag(got[k]) ==
            Catch::Approx(std::imag(expected[k])).margin(1e-12));
    }
  };

  SECTION("real") {
    auto yield_ =
        rand_tensor_yield<double, TA::SparsePolicy>{world, nocc, nvirt};
    auto eval_at = [&](ExprPtr const& expr) {
      return elements(
          evaluate(eval_node(expr), target, yield_)->get<TA::TSpArrayD>());
    };
    auto map = [](std::vector<double> v, auto const& f) {
      for (auto& x : v) x = f(x);
      return v;
    };

    auto node = eval_node(pow(D, rational{-1}));
    REQUIRE(node->op_type() == EvalOp::Power);
    auto const src =
        TA::clone(yield_(node.left())->get<TA::TSpArrayD>());  // pristine copy
    auto const d = eval_at(D);

    check_elementwise(eval_at(pow(D, rational{-1})),
                      map(d, [](double x) { return 1 / x; }));
    check_elementwise(eval_at(pow(D, rational{0})),
                      map(d, [](double) { return 1.0; }));
    check_elementwise(eval_at(pow(D, rational{3})),
                      map(d, [](double x) { return x * x * x; }));
    // the harness's data are signed: a fractional power of a negative
    // element is NaN
    {
      auto const got = eval_at(pow(D, rational{1, 2}));
      REQUIRE(got.size() == d.size());
      for (std::size_t k = 0; k < d.size(); ++k) {
        if (d[k] < 0)
          CHECK(std::isnan(got[k]));
        else
          CHECK(got[k] == Catch::Approx(std::sqrt(d[k])).margin(1e-12));
      }
    }

    // ℓ = M ⊙ D^-1, a Hadamard product
    auto const m = eval_at(M);
    ResultExpr ell(Tensor(L"ℓ", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"},
                          TensorSymmetries{.perm = Symmetry::Antisymm}),
                   M * pow(D, rational{-1}));
    {
      auto expected = m;
      for (std::size_t k = 0; k < m.size(); ++k) expected[k] /= d[k];
      check_elementwise(
          elements(
              evaluate(eval_node(ell), target, yield_)->get<TA::TSpArrayD>()),
          expected);
    }

    // E = L · (M ⊙ D^-1): every index appears on three tensors
    {
      auto ctx = get_default_context();
      ctx.set(AssertStrictBraKetSymmetry::No);
      auto _ = set_scoped_default_context(ctx);
      auto L = ex<Tensor>(L"L", bra{L"a_1", L"a_2"}, ket{L"i_1", L"i_2"},
                          TensorSymmetries{.perm = Symmetry::Antisymm});
      ResultExpr E(Variable(L"E"), L * M * pow(D, rational{-1}));
      auto const got =
          evaluate(eval_node(E), std::string{}, yield_)->get<double>();
      auto const l = eval_at(L);
      double expected = 0;
      for (std::size_t k = 0; k < l.size(); ++k) expected += l[k] * m[k] / d[k];
      CHECK(got == Catch::Approx(expected).margin(1e-10));
    }

    // M ⊙ D^-1 + D^-1, with D^-1 as either summand, evaluates correctly and
    // leaves the cached D untouched
    for (bool const power_first : {false, true}) {
      auto const hadamard = M * pow(D, rational{-1});
      auto const power = pow(D, rational{-1});
      ResultExpr sum(Tensor(L"ℓ", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"},
                            TensorSymmetries{.perm = Symmetry::Nonsymm}),
                     power_first ? power + hadamard : hadamard + power);
      auto expected = m;
      for (std::size_t k = 0; k < m.size(); ++k)
        expected[k] = (m[k] + 1) / d[k];
      check_elementwise(
          elements(
              evaluate(eval_node(sum), target, yield_)->get<TA::TSpArrayD>()),
          expected);
      check_elementwise(elements(yield_(node.left())->get<TA::TSpArrayD>()),
                        elements(src));
    }
  }

  SECTION("sparse") {
    using TSpArrayD = TA::DistArray<TA::Tensor<double>, TA::SparsePolicy>;
    auto yield_ =
        rand_tensor_yield<double, TA::SparsePolicy>{world, nocc, nvirt};
    auto const d =
        elements(evaluate(eval_node(D), target, yield_)->get<TSpArrayD>());
    auto expected = d;
    for (auto& x : expected) x = x * x;
    check_elementwise(
        elements(evaluate(eval_node(pow(D, rational{2})), target, yield_)
                     ->get<TSpArrayD>()),
        expected);

    // a zero tile of a sparse array would become inf under a negative
    // exponent, which a zero tile cannot hold
    if (assert_behavior() != AssertBehavior::Abort) {
      TA::TiledRange const tr{TA::TiledRange1{0, 2, 4}};
      TA::Tensor<float> norms(tr.tiles_range(), 1.0f);
      norms.at_ordinal(1) = 0.0f;
      TSpArrayD arr(world, tr, TA::SparseShape<float>(norms, tr));
      arr.fill(2.0);
      REQUIRE(arr.shape().sparsity() > 0);
      REQUIRE_THROWS_AS(detail::ta_elementwise_pow(arr, rational{-1}),
                        Exception);
    }
  }

  SECTION("complex") {
    using TArrayC =
        TA::DistArray<TA::Tensor<std::complex<double>>, TA::SparsePolicy>;
    auto yield_ = rand_tensor_yield<std::complex<double>, TA::SparsePolicy>{
        world, nocc, nvirt};
    auto const d =
        elements(evaluate(eval_node(D), target, yield_)->get<TArrayC>());
    auto square = d;
    auto inverse = d;
    for (std::size_t k = 0; k < d.size(); ++k) {
      square[k] = d[k] * d[k];
      inverse[k] = 1.0 / d[k];
    }
    check_elementwise(
        elements(evaluate(eval_node(pow(D, rational{2})), target, yield_)
                     ->get<TArrayC>()),
        square);
    check_elementwise(
        elements(evaluate(eval_node(pow(D, rational{-1})), target, yield_)
                     ->get<TArrayC>()),
        inverse);
  }
}

TEST_CASE("eval_broadcast_sum_ta", "[eval][broadcast]") {
  using namespace sequant;
  auto& world = TA::get_default_world();
  const size_t nocc = 2, nvirt = 3;
  auto yield_ = rand_tensor_yield<double, TA::SparsePolicy>{world, nocc, nvirt};
  // a Broadcast node shapes the modes it adds through the array ops
  auto eval_with = [](auto const& node, std::string const& layout,
                      auto const& yield) {
    auto const aops = yield.array_ops();
    auto cache = CacheManager<std::remove_cvref_t<decltype(node)>>::empty();
    cache.set_array_ops(&aops);
    return evaluate(node, layout, yield, cache);
  };
  auto eval_aops = [&](auto const& node, std::string const& layout) {
    return eval_with(node, layout, yield_);
  };

  // the elements of an array in row-major order of its (single-tile) layout
  auto elements = [](TA::TSpArrayD const& arr) {
    std::vector<double> out;
    for (auto it = arr.begin(); it != arr.end(); ++it) {
      auto const tile = it->get();
      out.insert(out.end(), tile.begin(), tile.end());
    }
    return out;
  };
  auto check_elementwise = [](std::vector<double> const& got,
                              std::vector<double> const& expected) {
    REQUIRE(got.size() == expected.size());
    REQUIRE(!got.empty());
    for (std::size_t k = 0; k < got.size(); ++k)
      CHECK(got[k] == Catch::Approx(expected[k]).margin(1e-12));
  };

  auto eps_a = [](wchar_t const* a) { return ex<Tensor>(L"ε", bra{}, ket{a}); };
  auto eps_i = [](wchar_t const* i) { return ex<Tensor>(L"ε", bra{i}, ket{}); };
  // the yielded ε vectors: every ket-ε occurrence shares one array by key,
  // likewise every bra-ε occurrence
  auto const ea =
      elements(yield_(eps_a(L"a_1")->as<Tensor>())->get<TA::TSpArrayD>());
  auto const ei =
      elements(yield_(eps_i(L"i_1")->as<Tensor>())->get<TA::TSpArrayD>());
  REQUIRE(ea.size() == nvirt);
  REQUIRE(ei.size() == nocc);

  SECTION("four-term denominator") {
    const TensorSymmetries dsym{.perm = Symmetry::Symm,
                                .braket = BraKetSymmetry::Symm,
                                .column = ColumnSymmetry::Symm};
    ResultExpr res(
        Tensor(L"D", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"}, dsym),
        eps_a(L"a_1") + eps_a(L"a_2") - eps_i(L"i_1") - eps_i(L"i_2"));
    std::vector<double> expected;
    for (std::size_t i1 = 0; i1 < nocc; ++i1)
      for (std::size_t i2 = 0; i2 < nocc; ++i2)
        for (std::size_t a1 = 0; a1 < nvirt; ++a1)
          for (std::size_t a2 = 0; a2 < nvirt; ++a2)
            expected.push_back(ea[a1] + ea[a2] - ei[i1] - ei[i2]);
    check_elementwise(
        elements(
            eval_aops(eval_node(res), "i_1,i_2,a_1,a_2")->get<TA::TSpArrayD>()),
        expected);
  }

  SECTION("asymmetric head") {
    // distinct leaves and a nonsymmetric head: the expected tensor changes
    // under any mode swap, so this pins the result layout
    // the yield's same-shaped random vectors can coincide, so the four leaves
    // get chosen, pairwise distinct data
    auto const e1 = ex<Tensor>(L"ε1", bra{}, ket{L"a_1"});
    auto const e2 = ex<Tensor>(L"ε2", bra{}, ket{L"a_2"});
    auto const e3 = ex<Tensor>(L"ε3", bra{L"i_1"}, ket{});
    auto const e4 = ex<Tensor>(L"ε4", bra{L"i_2"}, ket{});
    std::vector<double> const v1{0.5, 1.7, 2.3}, v2{-3.1, 4.9, 0.2},
        v3{0.37, -1.3}, v4{2.9, 0.61};
    REQUIRE(v1.size() == nvirt);
    REQUIRE(v3.size() == nocc);
    for (auto const& [e, v] : {std::pair{e1, v1}, std::pair{e2, v2},
                               std::pair{e3, v3}, std::pair{e4, v4}}) {
      TA::TSpArrayD arr(world, TA::TiledRange{yield_.tr1_for(L"", v.size())});
      arr.init_tiles([v](TA::Range const& r) {
        TA::Tensor<double> t(r);
        for (std::size_t k = 0; k < t.size(); ++k)
          t.at_ordinal(k) = v[r.lobound()[0] + k];
        return t;
      });
      yield_.set_leaf(e->as<Tensor>(),
                      eval_result<ResultTensorTA<TA::TSpArrayD>>(arr));
    }
    const TensorSymmetries nsym{.perm = Symmetry::Nonsymm,
                                .braket = BraKetSymmetry::Nonsymm,
                                .column = ColumnSymmetry::Nonsymm};
    ResultExpr res(
        Tensor(L"D'", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"}, nsym),
        e1 + e2 - e3 - e4);
    std::vector<double> expected;
    for (std::size_t i1 = 0; i1 < nocc; ++i1)
      for (std::size_t i2 = 0; i2 < nocc; ++i2)
        for (std::size_t a1 = 0; a1 < nvirt; ++a1)
          for (std::size_t a2 = 0; a2 < nvirt; ++a2)
            expected.push_back(v1[a1] + v2[a2] - v3[i1] - v4[i2]);
    check_elementwise(
        elements(
            eval_aops(eval_node(res), "i_1,i_2,a_1,a_2")->get<TA::TSpArrayD>()),
        expected);
  }

  SECTION("scalar summand") {
    ResultExpr res(Tensor(L"D", bra{L"i_1"}, ket{L"a_1"}),
                   eps_a(L"a_1") - eps_i(L"i_1") - ex<Variable>(L"E"));
    double const e = yield_(Variable(L"E"))->get<double>();
    std::vector<double> expected;
    for (std::size_t i1 = 0; i1 < nocc; ++i1)
      for (std::size_t a1 = 0; a1 < nvirt; ++a1)
        expected.push_back(ea[a1] - ei[i1] - e);
    check_elementwise(
        elements(eval_aops(eval_node(res), "i_1,a_1")->get<TA::TSpArrayD>()),
        expected);
  }

  SECTION("without array ops") {
    ResultExpr res(Tensor(L"D", bra{L"i_1"}, ket{L"a_1"}),
                   eps_a(L"a_1") - eps_i(L"i_1"));
    REQUIRE_THROWS_WITH(
        evaluate(eval_node(res), std::string{"i_1,a_1"}, yield_),
        Catch::Matchers::ContainsSubstring("needs BackendArrayOps"));
  }

  SECTION("sparse policy") {
    using TSpArrayD = TA::DistArray<TA::Tensor<double>, TA::SparsePolicy>;
    auto syield =
        rand_tensor_yield<double, TA::SparsePolicy>{world, nocc, nvirt};
    syield.set_max_tile(1);
    // single-element tiles; the second of the ket vector is zero, so the
    // sparse leaf screens it out
    auto const e1 = ex<Tensor>(L"ε1", bra{}, ket{L"a_1"});
    auto const e3 = ex<Tensor>(L"ε3", bra{L"i_1"}, ket{});
    std::vector<double> const v1{0.5, 0.0, 2.3}, v3{0.37, -1.3};
    auto sparse_of = [&](std::vector<double> const& v) {
      return TA::make_array<TSpArrayD>(
          world, TA::TiledRange{syield.tr1_for(L"", v.size())},
          [v](TA::Tensor<double>& t, TA::Range const& r) -> float {
            t = TA::Tensor<double>(r);
            for (std::size_t k = 0; k < t.size(); ++k)
              t.at_ordinal(k) = v[r.lobound()[0] + k];
            return t.norm();
          });
    };
    for (auto const& [e, v] : {std::pair{e1, v1}, std::pair{e3, v3}})
      syield.set_leaf(e->as<Tensor>(),
                      eval_result<ResultTensorTA<TSpArrayD>>(sparse_of(v)));

    // a result tile is zero iff the operand tile it repeats is
    {
      auto const aops = syield.array_ops();
      ResultPtr const operand = syield(e1->as<Tensor>());
      REQUIRE(operand->get<TSpArrayD>().is_zero({1}));
      auto const b =
          operand->broadcast({std::string{"a_1"}, std::string{"i_1,a_1"}},
                             {aops.make_zeros({Index{L"i_1"}})});
      auto const& arr = b->get<TSpArrayD>();
      REQUIRE(arr.trange().dim(0).tile_extent() == nocc);
      REQUIRE(arr.trange().dim(1).tile_extent() == nvirt);
      for (std::size_t i = 0; i < nocc; ++i)
        for (std::size_t a = 0; a < nvirt; ++a) {
          CHECK(arr.is_zero({i, a}) == (a == 1));
          if (a != 1)
            CHECK(arr.find({i, a}).get()(i, a) ==
                  Catch::Approx(v1[a]).margin(1e-12));
        }
    }

    // D{i;a} = ε1{;a} - ε3{i;}; the screened ε1 tile still gives a nonzero
    // result tile through ε3
    ResultExpr res(Tensor(L"D", bra{L"i_1"}, ket{L"a_1"}), e1 - e3);
    auto const d =
        eval_with(eval_node(res), "i_1,a_1", syield)->get<TSpArrayD>();
    REQUIRE(d.trange().dim(0).tile_extent() == nocc);
    REQUIRE(d.trange().dim(1).tile_extent() == nvirt);
    for (std::size_t i = 0; i < nocc; ++i)
      for (std::size_t a = 0; a < nvirt; ++a) {
        REQUIRE_FALSE(d.is_zero({i, a}));
        CHECK(d.find({i, a}).get()(i, a) ==
              Catch::Approx(v1[a] - v3[i]).margin(1e-12));
      }
  }
}

TEST_CASE("eval_inlined_system_ta", "[eval][inline]") {
  using namespace sequant;
  namespace tp = sequant::tests::tensor_power;
  namespace data = tp::data;
  using System = container::svector<ResultExpr>;
  using Flags = container::svector<bool>;
  auto& world = TA::get_default_world();
  // three occupied orbitals: with two, an antisymmetric pair leaves E{i_2;i_3}
  // diagonal
  const size_t nocc = 3, nvirt = 3;
  auto yield_ = rand_tensor_yield<double, TA::SparsePolicy>{world, nocc, nvirt};
  // a Broadcast node shapes the modes it adds through the array ops
  auto eval_with = [](auto const& node, std::string const& layout,
                      auto const& yield) {
    auto const aops = yield.array_ops();
    auto cache = CacheManager<std::remove_cvref_t<decltype(node)>>::empty();
    cache.set_array_ops(&aops);
    return evaluate(node, layout, yield, cache);
  };
  auto eval_aops = [&](auto const& node, std::string const& layout) {
    return eval_with(node, layout, yield_);
  };
  // the elements of an array in row-major order of its (single-tile) layout
  auto elements = [](TA::TSpArrayD const& arr) {
    std::vector<double> out;
    for (auto it = arr.begin(); it != arr.end(); ++it) {
      auto const tile = it->get();
      out.insert(out.end(), tile.begin(), tile.end());
    }
    return out;
  };
  auto check_elementwise = [](std::vector<double> const& got,
                              std::vector<double> const& expected) {
    REQUIRE(got.size() == expected.size());
    REQUIRE(!got.empty());
    for (std::size_t k = 0; k < got.size(); ++k)
      CHECK(got[k] == Catch::Approx(expected[k]).margin(1e-12));
  };
  auto pin = [&](ExprPtr const& tnsr, auto f) {
    yield_.pin(tnsr->as<Tensor>(), std::move(f));
  };
  auto const M = tp::antisymm(L"M", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"});
  auto const D = tp::power_base();
  auto const L = tp::contractor();
  auto pin_M = [&] {
    pin(M, [](std::vector<std::size_t> const& x) {
      return data::M(x[0], x[1], x[2], x[3]);
    });
  };
  auto pin_D = [&] {
    pin(D, [](std::vector<std::size_t> const& x) {
      return data::D(x[0], x[1], x[2], x[3]);
    });
  };
  auto pin_L = [&] {
    pin(L, [](std::vector<std::size_t> const& x) {
      return data::L(x[0], x[1], x[2], x[3]);
    });
  };
  // the value of the one root of a system
  auto eval_system = [&](System const& sys, std::string const& layout) {
    auto const res = binarize<EvalExprTA>(sys);
    REQUIRE(res.roots().size() == 1);
    return eval_aops(res.roots()[0].get(), layout);
  };
  // ℓ read at the given slots
  auto ell_use = [](auto b, auto k) {
    return ex<Tensor>(tp::relabelled(tp::ell_definition().result_as_tensor(),
                                     std::move(b), std::move(k)));
  };
  std::string const e_layout = "i_2,i_3";
  // E{i_2;i_3} of consumer_E(), row-major
  auto const expected_E = [&] {
    std::vector<double> result;
    for (std::size_t i2 = 0; i2 < nocc; ++i2)
      for (std::size_t i3 = 0; i3 < nocc; ++i3)
        result.push_back(data::contract_with_ell(data::L, i3, i2, nocc, nvirt));
    return result;
  };

  SECTION("inlined E equals the driver-chained E and the direct form") {
    pin_M();
    pin_D();
    pin_L();
    System const sys{tp::ell_definition(), tp::consumer_E()};
    REQUIRE(binarize(sys).consumed == Flags{true, false});
    auto const inlined =
        elements(eval_system(sys, e_layout)->get<TA::TSpArrayD>());
    check_elementwise(inlined, expected_E());

    // ℓ evaluated by itself, then served as the leaf E reads, in either
    // orientation (see pin)
    auto serve_ell = [&](ExprPtr const& use, std::string const& layout) {
      yield_.set_leaf(
          use->as<Tensor>(),
          evaluate(eval_node(tp::ell_definition()), layout, yield_));
    };
    serve_ell(ell_use(bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"}),
              "i_1,i_2,a_1,a_2");
    serve_ell(ell_use(bra{L"a_1", L"a_2"}, ket{L"i_1", L"i_2"}),
              "a_1,a_2,i_1,i_2");
    check_elementwise(inlined, elements(evaluate(eval_node(tp::consumer_E()),
                                                 e_layout, yield_)
                                            ->get<TA::TSpArrayD>()));

    // the direct form contracts L with M and D in one network, which the
    // strict bra-ket rule rejects
    auto _ = set_scoped_default_context(tp::lax_braket_context());
    check_elementwise(inlined,
                      elements(evaluate(eval_node(tp::hadamard_power_term()),
                                        e_layout, yield_)
                                   ->get<TA::TSpArrayD>()));
  }

  SECTION("three levels: ε into D into ℓ into E") {
    pin_M();
    pin_L();
    pin(ex<Tensor>(L"ε", bra{}, ket{L"a_1"}),
        [](std::vector<std::size_t> const& x) { return data::eps_virt(x[0]); });
    pin(ex<Tensor>(L"ε", bra{L"i_1"}, ket{}),
        [](std::vector<std::size_t> const& x) { return data::eps_occ(x[0]); });
    System const sys{tp::D_broadcast(), tp::ell_definition(), tp::consumer_E()};
    REQUIRE(binarize(sys).consumed == Flags{true, true, false});
    check_elementwise(
        elements(eval_system(sys, e_layout)->get<TA::TSpArrayD>()),
        expected_E());
  }

  SECTION("three levels batched over i_1, a mode the broadcasts add") {
    // one occupied orbital per tile, so the loop over i_1 runs one batch per
    // orbital; the Broadcast nodes of ε{;a_1}, ε{;a_2} and ε{i_2;} add i_1,
    // so each builds its batch's slice of i_1 from the array ops
    yield_.set_max_tile(1);
    pin_M();
    pin_L();
    pin(ex<Tensor>(L"ε", bra{}, ket{L"a_1"}),
        [](std::vector<std::size_t> const& x) { return data::eps_virt(x[0]); });
    pin(ex<Tensor>(L"ε", bra{L"i_1"}, ket{}),
        [](std::vector<std::size_t> const& x) { return data::eps_occ(x[0]); });
    System const sys{tp::D_broadcast(), tp::ell_definition(), tp::consumer_E()};
    auto const res = binarize<EvalExprTA>(sys);
    REQUIRE(res.roots().size() == 1);
    auto batched = res.roots()[0].get();
    tp::stamp_occ_loop(batched);
    auto policy = tp::occ_contracted_policy();
    auto const unbatched = elements(
        eval_aops(res.roots()[0].get(), e_layout)->get<TA::TSpArrayD>());
    check_elementwise(unbatched, expected_E());
    auto const aops = yield_.array_ops();
    for (auto const scheduler :
         {BatchScheduler::ordered, BatchScheduler::forest_descent}) {
      CAPTURE(static_cast<int>(scheduler));
      policy.scheduler = scheduler;
      auto cache =
          CacheManager<std::remove_cvref_t<decltype(batched)>>::empty();
      cache.set_array_ops(&aops);
      if (scheduler == BatchScheduler::forest_descent)
        cache.set_custom_evaluator(
            make_evaluator(policy, yield_, make_no_scope_guard{}));
      auto const slices_before = detail::slice_counter().load();
      auto const got =
          evaluate(std::vector{batched}, policy, e_layout, yield_, cache);
      check_elementwise(elements(got->get<TA::TSpArrayD>()), unbatched);
      // the ordered executor runs the loop: the operands and the added modes
      // of the Broadcast nodes are sliced per batch
      if (scheduler == BatchScheduler::ordered)
        CHECK(detail::slice_counter().load() > slices_before);
    }
  }

  SECTION("two uses of one definition") {
    pin_M();
    pin_D();
    pin_L();
    // G{i_2;i_3} = L{a_1,a_2;i_1,i_3} ℓ{i_1,i_2;a_1,a_2}
    //            + K{a_1,a_2;i_1,i_2} ℓ{i_1,i_3;a_1,a_2}
    auto const K = tp::antisymm(L"K", bra{L"a_1", L"a_2"}, ket{L"i_1", L"i_2"});
    pin(K, [](std::vector<std::size_t> const& x) {
      return data::K(x[0], x[1], x[2], x[3]);
    });
    ResultExpr const g(
        tp::head(L"G"),
        L * ell_use(bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"}) +
            K * ell_use(bra{L"i_1", L"i_3"}, ket{L"a_1", L"a_2"}));
    System const sys{tp::ell_definition(), g};
    REQUIRE(binarize(sys).consumed == Flags{true, false});
    std::vector<double> expected;
    for (std::size_t i2 = 0; i2 < nocc; ++i2)
      for (std::size_t i3 = 0; i3 < nocc; ++i3)
        expected.push_back(
            data::contract_with_ell(data::L, i3, i2, nocc, nvirt) +
            data::contract_with_ell(data::K, i2, i3, nocc, nvirt));
    check_elementwise(
        elements(eval_system(sys, e_layout)->get<TA::TSpArrayD>()), expected);
  }

  SECTION("a bra-ket-swapped use of a bra-ket-symmetric head") {
    pin(ex<Tensor>(L"ε", bra{}, ket{L"a_1"}),
        [](std::vector<std::size_t> const& x) { return data::eps_virt(x[0]); });
    pin(ex<Tensor>(L"ε", bra{L"i_1"}, ket{}),
        [](std::vector<std::size_t> const& x) { return data::eps_occ(x[0]); });
    auto const t = ex<Tensor>(L"t", bra{L"i_2"}, ket{L"a_2"});
    pin(t,
        [](std::vector<std::size_t> const& x) { return data::t(x[0], x[1]); });
    // Y{i_1;a_1} = t{i_2;a_2} D{a_1,a_2;i_1,i_2}, D read with bra and ket
    // swapped
    auto y_of = [&t](auto b, auto k) {
      return ResultExpr(
          Tensor(L"Y", bra{L"i_1"}, ket{L"a_1"}),
          t * ex<Tensor>(L"D", std::move(b), std::move(k), tp::symmetric()));
    };
    System const sys{tp::D_broadcast(),
                     y_of(bra{L"a_1", L"a_2"}, ket{L"i_1", L"i_2"})};
    REQUIRE(binarize(sys).consumed == Flags{true, false});
    auto const swapped =
        elements(eval_system(sys, "i_1,a_1")->get<TA::TSpArrayD>());
    std::vector<double> expected;
    for (std::size_t i1 = 0; i1 < nocc; ++i1)
      for (std::size_t a1 = 0; a1 < nvirt; ++a1) {
        double y = 0;
        for (std::size_t i2 = 0; i2 < nocc; ++i2)
          for (std::size_t a2 = 0; a2 < nvirt; ++a2)
            y += data::t(i2, a2) * data::D(i1, i2, a1, a2);
        expected.push_back(y);
      }
    check_elementwise(swapped, expected);

    // the use spelled like the head
    System const plain{tp::D_broadcast(),
                       y_of(bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"})};
    check_elementwise(
        swapped, elements(eval_system(plain, "i_1,a_1")->get<TA::TSpArrayD>()));
  }

  SECTION("a scalar head") {
    pin_M();
    pin_D();
    pin_L();
    // E = L{a_1,a_2;i_1,i_2} ℓ{i_1,i_2;a_1,a_2}
    ResultExpr const e(
        Variable(L"E"),
        tp::antisymm(L"L", bra{L"a_1", L"a_2"}, ket{L"i_1", L"i_2"}) *
            ell_use(bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"}));
    // F{i_1;a_1} = E t{i_1;a_1}
    auto const t = ex<Tensor>(L"t", bra{L"i_1"}, ket{L"a_1"});
    pin(t,
        [](std::vector<std::size_t> const& x) { return data::t(x[0], x[1]); });
    ResultExpr const f(Tensor(L"F", bra{L"i_1"}, ket{L"a_1"}),
                       ex<Variable>(L"E") * t);
    System const sys{tp::ell_definition(), e, f};
    REQUIRE(binarize(sys).consumed == Flags{true, true, false});

    // E through its own tree
    double const e_val =
        eval_system(System{tp::ell_definition(), e}, std::string{})
            ->get<double>();
    double e_ref = 0;
    for (std::size_t i2 = 0; i2 < nocc; ++i2)
      e_ref += data::contract_with_ell(data::L, i2, i2, nocc, nvirt);
    CHECK(e_val == Catch::Approx(e_ref).margin(1e-12));

    std::vector<double> expected;
    for (std::size_t i1 = 0; i1 < nocc; ++i1)
      for (std::size_t a1 = 0; a1 < nvirt; ++a1)
        expected.push_back(e_val * data::t(i1, a1));
    check_elementwise(
        elements(
            eval_system(sys, std::string{"i_1,a_1"})->get<TA::TSpArrayD>()),
        expected);
  }
}

// The spin-orbital MP2 energy of water (STO-3G, the data of the eval
// integration test), authored as a system of three equations:
//
//   D{i_1,i_2;a_1,a_2} = ε{;a_1} + ε{;a_2} - ε{i_1;} - ε{i_2;}
//   t{i_1,i_2;a_1,a_2} = g{i_1,i_2;a_1,a_2} ⊙ pow(D{i_1,i_2;a_1,a_2}, -1)
//   E = -1/4 g{a_1,a_2;i_1,i_2} t{i_1,i_2;a_1,a_2}
//
// t is the first-order doubles amplitude, here with the sign convention
// t = g / (ε_a + ε_b - ε_i - ε_j), hence the minus in E. D is a broadcast
// sum: each rank-1 ε is repeated over the modes of D it lacks (EvalOp::
// Broadcast, which shapes the added modes through the array ops of the
// cache); pow(D, -1) is an element-wise power (TensorPower). t reads D, and E
// reads t, through the head of the defining equation, so binarize(system)
// inlines D into t and t into E: only the tree of E is left to evaluate.
TEST_CASE("mp2-energy", "[eval][mp2]") {
  using namespace sequant;
  namespace tp = sequant::tests::tensor_power;
  using Flags = container::svector<bool>;

  // Fock matrix and antisymmetrized integrals <pq||rs> over all spin orbitals,
  // occupied ones first, row-major as the files store them
  std::string const dir = SEQUANT_EVAL_DATA_DIR;
  REQUIRE(std::ifstream(dir + "/fock_so.dat").good());
  REQUIRE(std::ifstream(dir + "/eri_so.dat").good());
  eval::DataInfo const info(dir + "/fock_so.dat", dir + "/eri_so.dat");
  std::size_t const nocc = info.nocc(), nvirt = info.nvirt(),
                    nobs = nocc + nvirt;
  REQUIRE(nocc == 10);
  REQUIRE(nvirt == 4);
  std::vector<double> fock(nobs * nobs), eri(nobs * nobs * nobs * nobs);
  eval::read_tensor(info.fock_file(), fock);
  eval::read_tensor(info.eri_file(), eri);
  // orbital energies: the diagonal of the (canonical) Fock matrix
  auto eps = [&](std::size_t p) { return fock[p * nobs + p]; };
  auto g = [&](std::size_t p, std::size_t q, std::size_t r, std::size_t s) {
    return eri[((p * nobs + q) * nobs + r) * nobs + s];
  };

  // leaves: ε as rank-1 tensors, g antisymmetric in each slot pair and
  // bra-ket symmetric
  auto& world = TA::get_default_world();
  auto yield_ = rand_tensor_yield<double, TA::SparsePolicy>{world, nocc, nvirt};
  auto pin_leaves = [&] {
    yield_.pin(
        Tensor(L"ε", bra{}, ket{L"a_1"}),
        [&](std::vector<std::size_t> const& x) { return eps(nocc + x[0]); });
    yield_.pin(Tensor(L"ε", bra{L"i_1"}, ket{}),
               [&](std::vector<std::size_t> const& x) { return eps(x[0]); });
    // pin serves g{a,b;i,j} with the values of g{i,j;a,b}: the data is
    // bra-ket symmetric (checked below)
    yield_.pin(Tensor(L"g", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"},
                      tp::antisymmetric()),
               [&](std::vector<std::size_t> const& x) {
                 return g(x[0], x[1], nocc + x[2], nocc + x[3]);
               });
  };

  // the reference: E = -1/4 Σ_{ijab} g_{ijab} g_{abij} / (ε_a+ε_b-ε_i-ε_j)
  double e_ref = 0;
  for (std::size_t i = 0; i < nocc; ++i)
    for (std::size_t j = 0; j < nocc; ++j)
      for (std::size_t a = nocc; a < nobs; ++a)
        for (std::size_t b = nocc; b < nobs; ++b) {
          REQUIRE(g(a, b, i, j) == Catch::Approx(g(i, j, a, b)).margin(1e-12));
          e_ref -= g(i, j, a, b) * g(a, b, i, j) /
                   (4 * (eps(a) + eps(b) - eps(i) - eps(j)));
        }
  // the eval integration program (tests/integration/eval/ta) prints the same
  // value as the CCSD energy of its second iteration, the energy of the
  // amplitudes after one update from zero
  CHECK(e_ref == Catch::Approx(-0.04914988716234519).margin(1e-12));

  // the system
  Tensor const D(L"D", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"},
                 TensorSymmetries{.perm = Symmetry::Symm,
                                  .braket = BraKetSymmetry::Symm,
                                  .column = ColumnSymmetry::Symm});
  Tensor const t(L"t", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"},
                 TensorSymmetries{.perm = Symmetry::Antisymm,
                                  .braket = BraKetSymmetry::Symm,
                                  .column = ColumnSymmetry::Symm});
  auto eps_a = [](wchar_t const* a) { return ex<Tensor>(L"ε", bra{}, ket{a}); };
  auto eps_i = [](wchar_t const* i) { return ex<Tensor>(L"ε", bra{i}, ket{}); };
  auto g_of = [](auto b, auto k) {
    return ex<Tensor>(L"g", std::move(b), std::move(k), tp::antisymmetric());
  };
  ResultExpr const D_eq(
      D, eps_a(L"a_1") + eps_a(L"a_2") - eps_i(L"i_1") - eps_i(L"i_2"));
  ResultExpr const t_eq(t, g_of(bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"}) *
                               pow(ex<Tensor>(D), rational{-1}));
  ResultExpr const E_eq(Variable(L"E"),
                        ex<Constant>(rational{-1, 4}) *
                            g_of(bra{L"a_1", L"a_2"}, ket{L"i_1", L"i_2"}) *
                            ex<Tensor>(t));
  container::svector<ResultExpr> const system{D_eq, t_eq, E_eq};

  // D and t are inlined into their uses; E alone remains
  auto const inlined = binarize<EvalExprTA>(system);
  REQUIRE(inlined.consumed == Flags{true, true, false});
  REQUIRE(inlined.roots().size() == 1);

  SECTION("unbatched") {
    pin_leaves();
    auto const aops = yield_.array_ops();
    auto const roots = inlined.roots();
    auto const& root = roots[0].get();
    auto cache = CacheManager<std::remove_cvref_t<decltype(root)>>::empty();
    cache.set_array_ops(&aops);
    double const e =
        evaluate(root, std::string{}, yield_, cache)->get<double>();
    CHECK(e == Catch::Approx(e_ref).margin(1e-12));
  }

  SECTION("batched over i_1") {
    // one occupied orbital per tile, so the loop over i_1 runs one batch per
    // orbital
    yield_.set_max_tile(1);
    pin_leaves();
    auto batched = inlined.roots()[0].get();
    // the root scales by -1/4 the contraction g·t, its left operand, which
    // sums over i_1
    tp::stamp_occ_loop(batched.left());
    auto policy = tp::occ_contracted_policy();
    auto const aops = yield_.array_ops();
    for (auto const scheduler :
         {BatchScheduler::ordered, BatchScheduler::forest_descent}) {
      CAPTURE(static_cast<int>(scheduler));
      policy.scheduler = scheduler;
      auto cache =
          CacheManager<std::remove_cvref_t<decltype(batched)>>::empty();
      cache.set_array_ops(&aops);
      if (scheduler == BatchScheduler::forest_descent)
        cache.set_custom_evaluator(
            make_evaluator(policy, yield_, make_no_scope_guard{}));
      auto const slices_before = detail::slice_counter().load();
      double const e =
          evaluate(std::vector{batched}, policy, std::string{}, yield_, cache)
              ->get<double>();
      CHECK(e == Catch::Approx(e_ref).margin(1e-12));
      // the ordered executor runs the loop, slicing the operands per batch
      if (scheduler == BatchScheduler::ordered)
        CHECK(detail::slice_counter().load() > slices_before);
    }
  }
}

// Minimal TA-level probe for the nested-tile permutation crash seen in the
// ToT culprit replica: permute a nested array on its outer modes only, on
// inner only, and on both, through TA expressions.
TEST_CASE("TA nested-tile permutation probe", "[tot][tot-perm-probe]") {
  using ToTArray = TA::DistArray<tot_tile_t<double>, TA::SparsePolicy>;
  auto& world = TA::get_default_world();
  TA::TiledRange tr{{0, 2, 4}, {0, 3, 6}};
  ToTArray in{world, tr};
  for (auto it = in.begin(); it != in.end(); ++it) {
    if (!in.is_local(it.index())) continue;
    // inner 2x3 tensors, a few left EMPTY (screened cells)
    *it = make_tot_tile<double>(
        it.make_range(),
        [](std::size_t o) {
          return o % 5 == 3 ? TA::Range{} : TA::Range{2, 3};
        },
        [](std::size_t o, std::size_t) {
          return 2.0 + static_cast<double>(o);
        });
  }
  world.gop.fence();
  ToTArray outer_only, inner_only, both;
  outer_only("j,i;a,b") = in("i,j;a,b");
  inner_only("i,j;b,a") = in("i,j;a,b");
  both("j,i;b,a") = in("i,j;a,b");
  world.gop.fence();
  CHECK(TA::squared_norm(outer_only) > 0.0);
  CHECK(TA::squared_norm(inner_only) > 0.0);
  CHECK(TA::squared_norm(both) > 0.0);
}

// Tile-level probe: outer transposition of a nested tile holding EMPTY inner
// cells, through the permuting constructor (the einsum lazy-tile cast path).
TEST_CASE("TA nested tile with empty inner cells: permuting constructor",
          "[tot][tot-tile-perm-probe]") {
  using Inner = TA::Tensor<double>;
  using Outer = TA::Tensor<Inner>;
  Outer tile(TA::Range{2, 3});
  std::size_t k = 0;
  for (auto& el : tile) {
    if ((k++ % 2) == 1) continue;  // leave every other inner cell EMPTY
    el = Inner(TA::Range{2, 2}, 1.0 + k);
  }
  TA::Permutation const outer_t{1, 0};
  TA::Permutation const inner_id{0, 1};
  std::cerr << "[probe] clone of tile with empty inners" << std::endl;
  Outer const c0 = tile.clone();
  CHECK(c0.range().extent(0) == 2);
  std::cerr << "[probe] outer_t * tile" << std::endl;
  Outer const t1 = outer_t * tile;  // outer transpose, inner unchanged
  CHECK(t1.range().extent(0) == 3);
  std::cerr << "[probe] Outer(tile, bp outer_t x id)" << std::endl;
  // the concatenated form: second partition offset by the first's size
  TA::BipartitePermutation const bp(TA::Permutation{1, 0, 2, 3}, 2);
  Outer const t2(tile, bp);
  CHECK(t2.range().extent(0) == 3);
  std::cerr << "[probe] Outer(tile, bp outer_t x inner_t)" << std::endl;
  TA::BipartitePermutation const bp2(TA::Permutation{1, 0, 3, 2}, 2);
  Outer const t3(tile, bp2);
  CHECK(t3.range().extent(0) == 3);
  std::cerr << "[probe] done" << std::endl;
}

// Roots that are results in their own right -- residuals of different
// ranks, an energy -- come back one each, in the forest's order, each what
// evaluating it alone gives; the summed entry refuses them.
TEST_CASE("evaluate_each_ta", "[eval]") {
  using namespace sequant;
  using Catch::Matchers::ContainsSubstring;
  auto& world = TA::get_default_world();
  const size_t nocc = 2, nvirt = 3;
  auto yield_ = rand_tensor_yield<double, TA::SparsePolicy>{world, nocc, nvirt};
  auto const aops = yield_.array_ops();

  auto root_of = [](std::wstring_view equation) {
    return binarize<EvalExprTA>(deserialize<ResultExpr>(equation));
  };
  using Node = EvalNode<EvalExprTA>;
  std::vector<Node> const forest{
      root_of(L"S{a_1;i_1}:N-N-N = f{a_1;a_2}:N-N-N * t{a_2;i_1}:N-N-N"),
      root_of(L"D{a_1,a_2;i_1,i_2}:N-N-N = g{a_1,a_2;a_3,a_4}:N-N-N * "
              L"t{a_3,a_4;i_1,i_2}:N-N-N"),
      root_of(L"E = f{i_1;a_1}:N-N-N * t{a_1;i_1}:N-N-N")};
  std::vector<std::string> const annots{"0,1", "0,1,2,3", ""};

  for (auto const scheduler :
       {BatchScheduler::ordered, BatchScheduler::forest_descent}) {
    CAPTURE(static_cast<int>(scheduler));
    BatchPolicy policy;
    policy.scheduler = scheduler;

    auto cache = CacheManager<Node>::empty();
    cache.set_array_ops(&aops);
    auto const results = evaluate_each(forest, policy, yield_, cache);
    REQUIRE(results.size() == forest.size());

    for (std::size_t i = 0; i != forest.size(); ++i) {
      CAPTURE(i);
      REQUIRE(results[i]);
      auto const alone = evaluate(forest[i], yield_);
      if (forest[i]->is_scalar()) {
        CHECK(results[i]->get<double>() == Catch::Approx(alone->get<double>()));
        continue;
      }
      auto const& got = results[i]->get<TA::TSpArrayD>();
      auto const& ref = alone->get<TA::TSpArrayD>();
      REQUIRE(got.trange() == ref.trange());
      TA::TSpArrayD diff;
      diff(annots[i]) = got(annots[i]) - ref(annots[i]);
      CHECK(TA::norm2(diff) == Catch::Approx(0).margin(1e-12));
      CHECK(TA::norm2(ref) > 0);
    }

    auto summed_cache = CacheManager<Node>::empty();
    summed_cache.set_array_ops(&aops);
    REQUIRE_THROWS_WITH(
        evaluate(forest, policy, std::string{}, yield_, summed_cache),
        ContainsSubstring("evaluate_each"));
  }
}

TEST_CASE("eval_sum_in_place_shared_operand_ta", "[eval]") {
  using namespace sequant;
  auto& world = TA::get_default_world();
  const size_t nocc = 2, nvirt = 3;
  auto yield_ = rand_tensor_yield<double, TA::SparsePolicy>{world, nocc, nvirt};
  auto const aops = yield_.array_ops();
  using Node = EvalNode<EvalExprTA>;
  auto node_of = [](std::wstring_view expr, Symmetry symm = Symmetry::Nonsymm) {
    return eval_node(
        deserialize<ExprPtr>(expr, {.def_perm_symm = symm,
                                    .def_braket_symm = BraKetSymmetry::Nonsymm,
                                    .def_col_symm = ColumnSymmetry::Nonsymm}));
  };
  auto cache_of = [&aops](std::vector<Node> const& nodes) {
    auto cache = cache_manager(nodes);
    cache.set_array_ops(&aops);
    return cache;
  };
  auto norm_of_difference = [](TA::TSpArrayD const& a, TA::TSpArrayD const& b,
                               std::string const& annot) {
    TA::TSpArrayD diff;
    diff(annot) = a(annot) - b(annot);
    return TA::norm2(diff);
  };

  SECTION("a pending operand of a product") {
    std::vector const nodes{
        node_of(L"(p{a_2;i_5} + q{a_2;i_5}) * (p{a_1;i_3} + q{a_1;i_3} + "
                L"g{a_1;i_3}) * W{i_3,i_5;a_1,a_2}")};
    auto const alone = evaluate(nodes[0], yield_)->get<double>();
    auto cache = cache_of(nodes);
    CHECK(evaluate(nodes[0], yield_, cache)->get<double>() ==
          Catch::Approx(alone).margin(1e-10));
  }

  SECTION("a result the caller holds") {
    std::vector const nodes{node_of(L"f{a_1;a_2} * t{a_2;i_1}"),
                            node_of(L"f{a_1;a_2} * t{a_2;i_1} + g{a_1;i_1}")};
    REQUIRE(nodes[1].left()->hash_value() == nodes[0]->hash_value());
    auto cache = cache_of(nodes);
    for (std::size_t k : {0, 1}) {
      CAPTURE(k);
      auto const alone = evaluate(nodes[k], yield_);
      auto const cached = evaluate(nodes[k], yield_, cache);
      CHECK(norm_of_difference(cached->get<TA::TSpArrayD>(),
                               alone->get<TA::TSpArrayD>(),
                               "0,1") == Catch::Approx(0).margin(1e-10));
    }
    // the first result, read again once the second root has been evaluated
    auto const first = evaluate(nodes[0], yield_);
    auto first_cache = cache_of(nodes);
    auto const held = evaluate(nodes[0], yield_, first_cache);
    (void)evaluate(nodes[1], yield_, first_cache);
    CHECK(norm_of_difference(held->get<TA::TSpArrayD>(),
                             first->get<TA::TSpArrayD>(),
                             "0,1") == Catch::Approx(0).margin(1e-10));
  }

  SECTION("a result the caller holds, read with the opposite phase") {
    // the two products are one value of opposite phases: what the Sum
    // accumulates into is a view of the tiles of the result held
    auto const plus = L"g{a_1,a_2;i_3,i_4} * t{i_3,i_4;i_1,i_2}";
    auto const minus = L"g{a_2,a_1;i_3,i_4} * t{i_3,i_4;i_1,i_2}";
    for (bool plus_first : {true, false}) {
      CAPTURE(plus_first);
      std::vector const nodes{
          node_of(plus_first ? plus : minus, Symmetry::Antisymm),
          node_of(std::wstring{plus_first ? minus : plus} +
                      L" + h{a_1,a_2;i_1,i_2}",
                  Symmetry::Antisymm)};
      REQUIRE(nodes[1].left()->hash_value() == nodes[0]->hash_value());
      REQUIRE(nodes[1].left()->canon_phase() == -nodes[0]->canon_phase());
      auto cache = cache_of(nodes);
      auto const held = evaluate(nodes[0], yield_, cache);
      TA::TSpArrayD const before = TA::clone(held->get<TA::TSpArrayD>());
      REQUIRE(TA::norm2(before) > 0);
      (void)evaluate(nodes[1], yield_, cache);
      CHECK(norm_of_difference(held->get<TA::TSpArrayD>(), before, "0,1,2,3") ==
            Catch::Approx(0).margin(1e-10));
    }
  }

  SECTION("the first root of a summed forest") {
    std::vector const nodes{node_of(L"f{a_1;a_2} * t{a_2;i_1} * Z{i_1;a_1}"),
                            node_of(L"g{a_1;i_1} * Z{i_1;a_1}"),
                            node_of(L"f{a_1;a_2} * t{a_2;i_1} * Z{i_1;a_1}")};
    auto const sum = 2 * evaluate(nodes[0], yield_)->get<double>() +
                     evaluate(nodes[1], yield_)->get<double>();
    auto cache = cache_of(nodes);
    CHECK(evaluate(nodes, std::string{}, yield_, cache)->get<double>() ==
          Catch::Approx(sum).margin(1e-10));
  }
}
