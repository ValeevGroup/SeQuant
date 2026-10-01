#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/wick.hpp>
#include <SeQuant/core/wick_extended.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>

#include <catch2/catch_test_macros.hpp>
#include "catch2_sequant.hpp"

TEST_CASE("wick_extended", "[algorithms][wick][valgrind_skip]") {
  using namespace sequant;
  using namespace sequant::tests;

  auto ctx = get_default_context();
  ctx.set(mbpt::make_mr_spaces());
  ctx.set(Vacuum::MultiProduct);
  auto ctx_resetter = set_scoped_default_context(ctx);

  // helper: standard Wick with all partial contractions + provenance from the
  // (uncanonicalized) input; adequate here because every index is external
  auto wick_partial = [](const FNOperatorSeq& nopseq, OpProvenance& prov) {
    prov.clear();
    std::size_t ord = 0;
    for (const auto& nop : nopseq) {
      for (const auto& op : nop) prov.emplace(op.index(), ord);
      ++ord;
    }
    FWickTheorem wick{std::make_shared<FNOperatorSeq>(nopseq)};
    return wick.full_contractions(false).compute();
  };

  SECTION("cumulant_expand: ⟨{a†a†a}{a}⟩ = κ2") {
    FNOperatorSeq in{FNOperator(cre({L"u_1", L"u_2"}), ann({L"u_3"})),
                     FNOperator(cre({}), ann({L"u_4"}))};
    OpProvenance prov;
    auto wick_out = wick_partial(in, prov);
    auto result = cumulant_expand<Statistics::FermiDirac>(wick_out, prov, {});
    REQUIRE_THAT(result, EquivalentTo(L"κ{u_4,u_3;u_1,u_2}:A-H-S"));
  }

  SECTION("cumulant_expand: ⟨{a†a}{a†a}⟩ = γη + κ2") {
    FNOperatorSeq in{FNOperator(cre({L"u_1"}), ann({L"u_2"})),
                     FNOperator(cre({L"u_3"}), ann({L"u_4"}))};
    OpProvenance prov;
    auto wick_out = wick_partial(in, prov);
    auto result = cumulant_expand<Statistics::FermiDirac>(wick_out, prov, {});
    REQUIRE_THAT(result, EquivalentTo(L"γ{u_4;u_1}:N-H-S * η{u_2;u_3}:N-H-S "
                                      L"+ κ{u_2,u_4;u_1,u_3}:A-H-S"));
  }

  SECTION("cumulant_expand: max_cumulant_rank = 1 means pairs only") {
    FNOperatorSeq in{FNOperator(cre({L"u_1"}), ann({L"u_2"})),
                     FNOperator(cre({L"u_3"}), ann({L"u_4"}))};
    OpProvenance prov;
    auto wick_out = wick_partial(in, prov);
    auto result = cumulant_expand<Statistics::FermiDirac>(
        wick_out, prov, {.max_cumulant_rank = 1});
    REQUIRE_THAT(result, EquivalentTo(L"γ{u_4;u_1}:N-H-S * η{u_2;u_3}:N-H-S"));
    auto result0 = cumulant_expand<Statistics::FermiDirac>(
        wick_out, prov, {.max_cumulant_rank = 0});
    REQUIRE(simplify(result - result0) == ex<Constant>(0));
  }

  SECTION("cumulant_expand: unbalanced survivors vanish") {
    FNOperatorSeq in{FNOperator(cre({L"u_1", L"u_2"}), ann({})),
                     FNOperator(cre({}), ann({L"u_3"}))};
    OpProvenance prov;
    auto wick_out = wick_partial(in, prov);
    auto result = cumulant_expand<Statistics::FermiDirac>(wick_out, prov, {});
    REQUIRE(result == ex<Constant>(0));
  }

  SECTION("cumulant_expand: a block within one operator vanishes") {
    // {a†_u1 a†_u2 a_u3 a_u4}{a†_u5 a_u6}: the 4 legs of nop 0 may not form
    // a block by themselves; every κ must involve u_5 or u_6
    FNOperatorSeq in{FNOperator(cre({L"u_1", L"u_2"}), ann({L"u_3", L"u_4"})),
                     FNOperator(cre({L"u_5"}), ann({L"u_6"}))};
    OpProvenance prov;
    auto wick_out = wick_partial(in, prov);
    auto result = cumulant_expand<Statistics::FermiDirac>(wick_out, prov, {});
    REQUIRE(result->is<Sum>());
    std::size_t nkappa2 = 0, nkappa3 = 0;
    for (const auto& term : *result) {
      const ExprPtrList single{term};
      for (const auto& f :
           term->is<Product>() ? term->as<Product>().factors() : single) {
        if (f->is<Tensor>() && f->as<Tensor>().label() == L"κ") {
          const auto& t = f->as<Tensor>();
          (t.bra_rank() == 2 ? nkappa2 : nkappa3) += 1;
          bool touches_nop1 = false;
          for (const auto& idx : t.const_braket())
            if (idx == Index(L"u_5") || idx == Index(L"u_6"))
              touches_nop1 = true;
          REQUIRE(touches_nop1);
        }
      }
    }
    REQUIRE(nkappa2 == 4);
    REQUIRE(nkappa3 == 1);
  }

  SECTION("cumulant_expand: 3-body truncation") {
    // {a†a†a}{a†a a} has a κ3 term (all six legs) that max_cumulant_rank=2
    // must drop, leaving everything else unchanged
    FNOperatorSeq in{FNOperator(cre({L"u_1", L"u_2"}), ann({L"u_3"})),
                     FNOperator(cre({L"u_4"}), ann({L"u_5", L"u_6"}))};
    OpProvenance prov;
    auto wick_out = wick_partial(in, prov);
    auto full = cumulant_expand<Statistics::FermiDirac>(wick_out, prov, {});
    auto trunc = cumulant_expand<Statistics::FermiDirac>(
        wick_out, prov, {.max_cumulant_rank = 2});
    auto diff = simplify(full - trunc);
    // exactly the κ3 term
    REQUIRE_THAT(diff, EquivalentTo(L"-κ{u_5,u_6,u_3;u_1,u_2,u_4}:A-H-S"));
  }

  SECTION("cumulant_expand: partial contractions leave a GNO remainder") {
    FNOperatorSeq in{FNOperator(cre({L"u_1"}), ann({L"u_2"})),
                     FNOperator(cre({L"u_3"}), ann({L"u_4"}))};
    OpProvenance prov;
    auto wick_out = wick_partial(in, prov);
    auto result = cumulant_expand<Statistics::FermiDirac>(
        wick_out, prov, {.full_contractions = false});
    // Eq. (ext. Wick, 1-body×1-body): the Wick output's 4 terms plus κ2;
    // nothing else, because a block needs ≥2 ops from ≥2 nops and the only
    // such balanced set is all four legs
    REQUIRE_THAT(result, EquivalentTo(L"ã{u_2,u_4;u_1,u_3} "
                                      L"- γ{u_4;u_1}:N-H-S * ã{u_2;u_3} "
                                      L"+ η{u_2;u_3}:N-H-S * ã{u_4;u_1} "
                                      L"+ γ{u_4;u_1}:N-H-S * η{u_2;u_3}:N-H-S "
                                      L"+ κ{u_2,u_4;u_1,u_3}:A-H-S"));
    for (const auto& term : *result)
      for (const auto& f : *term)
        if (f->is<FNOperator>())
          REQUIRE(f->as<FNOperator>().vacuum() == Vacuum::MultiProduct);
  }

  SECTION("cumulant_expand: an inactive survivor is not a cumulant leg") {
    // {a†_u1 a_i1}{a†_u2 a_u3}: the only balanced set of legs from both nops
    // includes the core annihilator i_1, so no κ forms
    FNOperatorSeq in{FNOperator(cre({L"u_1"}), ann({L"i_1"})),
                     FNOperator(cre({L"u_2"}), ann({L"u_3"}))};
    OpProvenance prov;
    auto wick_out = wick_partial(in, prov);
    REQUIRE(cumulant_expand<Statistics::FermiDirac>(wick_out, prov, {}) ==
            ex<Constant>(0));
    auto partial = cumulant_expand<Statistics::FermiDirac>(
        wick_out, prov, {.full_contractions = false});
    REQUIRE(simplify(partial - wick_out) == ex<Constant>(0));
  }

  SECTION("extended_wick: vacuum must be MultiProduct") {
    auto sr_ctx = get_default_context();
    sr_ctx.set(Vacuum::SingleProduct);
    auto sr_resetter = set_scoped_default_context(sr_ctx);
    auto in = ex<FNOperator>(cre({L"u_1"}), ann({L"u_2"})) *
              ex<FNOperator>(cre({L"u_3"}), ann({L"u_4"}));
    REQUIRE_THROWS_AS(extended_wick<Statistics::FermiDirac>(in), Exception);
  }

  SECTION("extended_wick: the input is not modified") {
    auto in = ex<FNOperator>(cre({L"u_1"}), ann({L"u_2"})) *
              ex<Tensor>(L"h", bra{L"u_2"}, ket{L"u_1"}, Symmetry::Nonsymm,
                         BraKetSymmetry::Nonsymm, ColumnSymmetry::Symm);
    extended_wick<Statistics::FermiDirac>(in);
    REQUIRE(in->as<Product>().factor(0)->is<FNOperator>());
  }

  SECTION("extended_wick: pure-active identities via the wrapper") {
    auto in = ex<FNOperator>(cre({L"u_1"}), ann({L"u_2"})) *
              ex<FNOperator>(cre({L"u_3"}), ann({L"u_4"}));
    auto result = extended_wick<Statistics::FermiDirac>(in);
    REQUIRE_THAT(result, EquivalentTo(L"γ{u_4;u_1}:N-H-S * η{u_2;u_3}:N-H-S "
                                      L"+ κ{u_2,u_4;u_1,u_3}:A-H-S"));
  }

  SECTION("extended_wick: general indices split into core δ + active γ") {
    // ⟨{a†_p1 a_p2}{a†_p3 a_p4}⟩ with p = M ∪ E = {o,i,u,a,g}:
    // cre·ann pair over R = M: δ on core O + γ on active u;
    // ann·cre pair over U = E: δ on virtual {a,g} + η on active u, where
    // {a,g} is not a registered space, so its δ splits into δ on a + δ on g;
    // plus κ2 on the all-active projection: 2 × 3 + 1 = 7 terms
    auto in = ex<FNOperator>(cre({L"p_1"}), ann({L"p_2"})) *
              ex<FNOperator>(cre({L"p_3"}), ann({L"p_4"}));
    auto result = extended_wick<Statistics::FermiDirac>(in);
    // every γ/η/κ index is active, also with partial contractions
    auto require_active_densities = [](const ExprPtr& expr) {
      REQUIRE(expr->is<Sum>());
      for (const auto& term : *expr)
        for (const auto& f : *term)
          if (f->is<Tensor>()) {
            const auto& t = f->as<Tensor>();
            if (t.label() == L"γ" || t.label() == L"η" || t.label() == L"κ")
              for (const auto& idx : t.const_braket())
                REQUIRE(idx.space() == Index(L"u_1").space());
          }
    };
    require_active_densities(result);
    require_active_densities(extended_wick<Statistics::FermiDirac>(
        in, {.full_contractions = false}));
    REQUIRE(result->size() == 7);
    // projecting is substituting a_p = Σ_x δ(p,x) a_x, so every term keeps
    // the + of γ{p_4;p_1}·η{p_2;p_3} + κ{p_2,p_4;p_1,p_3}, with γ over core
    // = δ and η over virtual = δ
    REQUIRE_THAT(
        result,
        EquivalentTo(
            L"κ{u_3,u_4;u_1,u_2}:A-C-S * δ{u_1;p_1}:N-C-S * "
            L"δ{u_2;p_3}:N-C-S * δ{p_2;u_3}:N-C-S * δ{p_4;u_4}:N-C-S "
            L"+ γ{u_2;u_1}:N-C-S * δ{u_1;p_1}:N-C-S * δ{a_1;p_3}:N-C-S * "
            L"δ{p_2;a_1}:N-C-S * δ{p_4;u_2}:N-C-S "
            L"+ γ{u_2;u_1}:N-C-S * δ{u_1;p_1}:N-C-S * δ{g_1;p_3}:N-C-S * "
            L"δ{p_2;g_1}:N-C-S * δ{p_4;u_2}:N-C-S "
            L"+ η{u_2;u_1}:N-C-S * δ{u_1;p_3}:N-C-S * δ{O_1;p_1}:N-C-S * "
            L"δ{p_2;u_2}:N-C-S * δ{p_4;O_1}:N-C-S "
            L"+ δ{a_1;p_3}:N-C-S * δ{O_1;p_1}:N-C-S * δ{p_2;a_1}:N-C-S * "
            L"δ{p_4;O_1}:N-C-S "
            L"+ δ{g_1;p_3}:N-C-S * δ{O_1;p_1}:N-C-S * δ{p_2;g_1}:N-C-S * "
            L"δ{p_4;O_1}:N-C-S "
            L"+ η{u_3;u_1}:N-C-S * γ{u_4;u_2}:N-C-S * δ{u_1;p_3}:N-C-S * "
            L"δ{u_2;p_1}:N-C-S * δ{p_2;u_3}:N-C-S * δ{p_4;u_4}:N-C-S"));
  }

  SECTION("extended_wick: dummy indices keep their provenance") {
    // a one-body h summed against its operator's indices, times an active
    // one-body operator: every op index of the first factor is a dummy
    const Index p1(L"p_1"), p2(L"p_2");
    auto in = ex<Tensor>(L"h", bra{p1}, ket{p2}, Symmetry::Nonsymm,
                         BraKetSymmetry::Conjugate, ColumnSymmetry::Symm) *
              ex<FNOperator>(cre({p1}), ann({p2})) *
              ex<FNOperator>(cre({L"u_3"}), ann({L"u_4"}));
    ExprPtr result;
    REQUIRE_NOTHROW(result = extended_wick<Statistics::FermiDirac>(in));
    // every γ/η/κ index is active
    for (const auto& term : *result)
      for (const auto& f : *term)
        if (f->is<Tensor>()) {
          const auto& t = f->as<Tensor>();
          if (t.label() == L"γ" || t.label() == L"η" || t.label() == L"κ")
            for (const auto& idx : t.const_braket())
              REQUIRE(idx.space() == Index(L"u_1").space());
        }
    REQUIRE_NOTHROW(extended_wick<Statistics::FermiDirac>(
        in, {.full_contractions = false}));
  }

  SECTION("extended_wick: δs over summed indices are applied") {
    // ⟨h^p_q a†_p a_q⟩: the δ binding each summed op index to its projection
    // is applied, the external ones (see "general indices split into core δ +
    // active γ") are kept
    auto h = [](const Index& b, const Index& k) {
      return ex<Tensor>(L"h", bra{b}, ket{k}, Symmetry::Nonsymm,
                        BraKetSymmetry::Conjugate, ColumnSymmetry::Symm);
    };
    const Index p1(L"p_1"), p2(L"p_2"), O1(L"O_1"), u1(L"u_1"), u2(L"u_2"),
        u3(L"u_3");
    auto in = h(p1, p2) * ex<FNOperator>(cre({p1}), ann({})) *
              ex<FNOperator>(cre({}), ann({p2}));
    auto result = extended_wick<Statistics::FermiDirac>(in);
    REQUIRE(simplify(result - h(O1, O1) -
                     h(u2, u1) * density::make_rdm(u1, u2)) == ex<Constant>(0));
    // ... also for a δ that η = δ - γ introduces
    auto in2 = h(u3, u2) * ex<FNOperator>(cre({L"u_1"}), ann({L"u_2"})) *
               ex<FNOperator>(cre({L"u_3"}), ann({L"u_4"}));
    auto result2 = extended_wick<Statistics::FermiDirac>(
        in2, {.eta_as_delta_minus_gamma = true});
    result2->visit(
        [](const ExprPtr& e) {
          if (e->is<Tensor>()) REQUIRE(e->as<Tensor>().label() != L"δ");
        },
        /*atoms_only=*/true);
  }

  SECTION("extended_wick: tensors commute with the theorem") {
    // operators carrying dummies of two different tensors: the result must
    // equal that of the same operators with external indices, times the
    // tensors
    auto tensors = ex<Tensor>(L"h", bra{L"u_1", L"u_2"}, ket{L"u_3", L"u_4"},
                              Symmetry::Nonsymm, BraKetSymmetry::Nonsymm,
                              ColumnSymmetry::Symm) *
                   ex<Tensor>(L"g", bra{L"u_5", L"u_6"}, ket{L"u_7", L"u_8"},
                              Symmetry::Nonsymm, BraKetSymmetry::Nonsymm,
                              ColumnSymmetry::Symm);
    auto ops = ex<FNOperator>(cre({L"u_1", L"u_2"}), ann({L"u_3", L"u_4"})) *
               ex<FNOperator>(cre({L"u_5", L"u_6"}), ann({L"u_7", L"u_8"})) *
               ex<FNOperator>(cre({L"u_9"}), ann({L"u_10"}));
    for (const bool full : {true, false}) {
      const ExtendedWickOptions opts{.full_contractions = full};
      auto lhs = extended_wick<Statistics::FermiDirac>(tensors * ops, opts);
      auto rhs =
          simplify(tensors * extended_wick<Statistics::FermiDirac>(ops, opts));
      REQUIRE(simplify(lhs - rhs) == ex<Constant>(0));
    }
  }

  SECTION("extended_wick: projected survivors with partial contractions") {
    // {a†_p1 a_p2}{a†_u3 a_u4}: projecting an op is substituting
    // a_p = Σ_x δ(p,x) a_x over the active, core and virtual parts of p, so
    // each term of the pure-active result (see "partial contractions leave a
    // GNO remainder") reappears with its sign, its indices projected
    auto in = ex<FNOperator>(cre({L"p_1"}), ann({L"p_2"})) *
              ex<FNOperator>(cre({L"u_3"}), ann({L"u_4"}));
    auto result =
        extended_wick<Statistics::FermiDirac>(in, {.full_contractions = false});
    // result contains `term` exactly once, with its sign
    auto contains = [&](std::wstring_view term) {
      return simplify(result - deserialize(term))->size() == result->size() - 1;
    };
    // a_p2·a†_u3 contracted (adjacent, so +), p_1 of the surviving
    // {a†_p1 a_u4} projected onto active and onto virtual a
    REQUIRE(
        contains(L"η{u_2;u_3}:N-C-S * δ{u_1;p_1}:N-C-S * δ{p_2;u_2}:N-C-S * "
                 L"ã{u_4;u_1}"));
    REQUIRE(
        contains(L"η{u_1;u_3}:N-C-S * δ{a_1;p_1}:N-C-S * δ{p_2;u_1}:N-C-S * "
                 L"ã{u_4;a_1}"));
    // nothing contracted, p_1 and p_2 projected onto active: the four legs
    // form +κ{u_2,u_4;u_1,u_3} = -κ{u_4,u_2;u_1,u_3}
    REQUIRE(
        contains(L"-1 κ{u_4,u_2;u_1,u_3}:A-C-S * δ{u_1;p_1}:N-C-S * "
                 L"δ{p_2;u_2}:N-C-S"));
  }

  SECTION("extended_wick: connectivity") {
    // partial contractions: with full ones, connecting 0 to 1 but not to 2
    // leaves 2 isolated, and no term survives
    auto in = ex<FNOperator>(cre({L"u_1"}), ann({L"u_2"})) *
              ex<FNOperator>(cre({L"u_3"}), ann({L"u_4"})) *
              ex<FNOperator>(cre({L"u_5"}), ann({L"u_6"}));
    auto all =
        extended_wick<Statistics::FermiDirac>(in, {.full_contractions = false});
    auto filtered = extended_wick<Statistics::FermiDirac>(
        in, {.full_contractions = false,
             .nop_connections = {{0, 1}},
             .nop_avoided_connections = {{0, 2}}});
    REQUIRE(filtered->is<Sum>());
    REQUIRE(filtered->size() > 0);
    REQUIRE(filtered->size() < all->size());
    auto ord = [](const Index& i) -> int {
      const auto l = i.label();
      if (l == L"u_1" || l == L"u_2") return 0;
      if (l == L"u_3" || l == L"u_4") return 1;
      return 2;
    };
    for (const auto& term : *filtered) {
      bool e01 = false, e02 = false;
      const ExprPtrList single{term};
      for (const auto& f :
           term->is<Product>() ? term->as<Product>().factors() : single) {
        if (!f->is<Tensor>()) continue;
        const auto& t = f->as<Tensor>();
        container::set<int> ords;
        for (const auto& idx : t.const_braket()) ords.insert(ord(idx));
        if (ords.contains(0) && ords.contains(1)) e01 = true;
        if (ords.contains(0) && ords.contains(2)) e02 = true;
      }
      REQUIRE(e01);
      REQUIRE(!e02);
    }
  }

  SECTION("extended_wick: connectivity through a cumulant") {
    // partial contractions: a 0-2 connection can be realized by a cumulant
    // block alone, which no pair of operators expresses
    auto in = ex<FNOperator>(cre({L"u_1"}), ann({L"u_2"})) *
              ex<FNOperator>(cre({L"u_3"}), ann({L"u_4"})) *
              ex<FNOperator>(cre({L"u_5"}), ann({L"u_6"}));
    auto all =
        extended_wick<Statistics::FermiDirac>(in, {.full_contractions = false});
    auto filtered = extended_wick<Statistics::FermiDirac>(
        in, {.full_contractions = false, .nop_connections = {{0, 2}}});
    REQUIRE(filtered->is<Sum>());
    REQUIRE(filtered->size() > 0);
    REQUIRE(filtered->size() < all->size());
    auto ord = [](const Index& i) -> int {
      const auto l = i.label();
      if (l == L"u_1" || l == L"u_2") return 0;
      if (l == L"u_3" || l == L"u_4") return 1;
      return 2;
    };
    bool kappa_only = false;
    for (const auto& term : *filtered) {
      bool by_pair = false, by_kappa = false;
      const ExprPtrList single{term};
      for (const auto& f :
           term->is<Product>() ? term->as<Product>().factors() : single) {
        if (!f->is<Tensor>()) continue;
        const auto& t = f->as<Tensor>();
        container::set<int> ords;
        for (const auto& idx : t.const_braket()) ords.insert(ord(idx));
        if (ords.contains(0) && ords.contains(2))
          (t.label() == L"κ" ? by_kappa : by_pair) = true;
      }
      REQUIRE((by_pair || by_kappa));
      if (by_kappa && !by_pair) kappa_only = true;
    }
    REQUIRE(kappa_only);
  }

  SECTION("extended_wick: a coefficient tensor is not a connection") {
    // h{u_1;u_6} spans operators 0 and 2, but only the factors the theorem
    // produces connect operators
    auto h = ex<Tensor>(L"h", bra{L"u_1"}, ket{L"u_6"}, Symmetry::Nonsymm,
                        BraKetSymmetry::Nonsymm, ColumnSymmetry::Symm);
    auto ops = ex<FNOperator>(cre({L"u_1"}), ann({L"u_2"})) *
               ex<FNOperator>(cre({L"u_3"}), ann({L"u_4"})) *
               ex<FNOperator>(cre({L"u_5"}), ann({L"u_6"}));
    for (const auto& opts : {ExtendedWickOptions{.full_contractions = false,
                            .nop_avoided_connections = {{0, 2}}},
                             ExtendedWickOptions{.full_contractions = false,
                             .nop_connections = {{0, 2}}}}) {
      auto without_h = extended_wick<Statistics::FermiDirac>(ops, opts);
      REQUIRE(without_h->size() > 0);
      auto lhs = extended_wick<Statistics::FermiDirac>(h * ops, opts);
      REQUIRE(simplify(lhs - simplify(h * without_h)) == ex<Constant>(0));
    }
  }

  SECTION("extended_wick: connection ordinals must name input operators") {
    auto in = ex<FNOperator>(cre({L"u_1"}), ann({L"u_2"})) *
              ex<FNOperator>(cre({L"u_3"}), ann({L"u_4"}));
    REQUIRE_THROWS_AS(extended_wick<Statistics::FermiDirac>(
                          in, {.nop_connections = {{0, 2}}}),
                      Exception);
    REQUIRE_THROWS_AS(extended_wick<Statistics::FermiDirac>(
                          in, {.nop_avoided_connections = {{2, 1}}}),
                      Exception);
  }

  SECTION("extended_wick: Sum input") {
    auto in = ex<FNOperator>(cre({L"u_1"}), ann({L"u_2"})) *
                  ex<FNOperator>(cre({L"u_3"}), ann({L"u_4"})) +
              ex<Constant>(2) * ex<FNOperator>(cre({L"u_1"}), ann({L"u_2"})) *
                  ex<FNOperator>(cre({L"u_3"}), ann({L"u_4"}));
    auto result = extended_wick<Statistics::FermiDirac>(in);
    REQUIRE_THAT(result, EquivalentTo(L"3 γ{u_4;u_1}:N-H-S * η{u_2;u_3}:N-H-S "
                                      L"+ 3 κ{u_2,u_4;u_1,u_3}:A-H-S"));
  }

  SECTION("extended_wick: η = δ - γ") {
    auto in = ex<FNOperator>(cre({L"u_1"}), ann({L"u_2"})) *
              ex<FNOperator>(cre({L"u_3"}), ann({L"u_4"}));
    auto result = extended_wick<Statistics::FermiDirac>(
        in, {.eta_as_delta_minus_gamma = true});
    REQUIRE_THAT(result, EquivalentTo(L"γ{u_4;u_1}:N-H-S * δ{u_2;u_3} "
                                      L"- γ{u_4;u_1}:N-H-S * γ{u_2;u_3}:N-H-S "
                                      L"+ κ{u_2,u_4;u_1,u_3}:A-H-S"));
    bool has_eta = false;
    result->visit(
        [&](const ExprPtr& e) {
          if (e->is<Tensor>() && e->as<Tensor>().label() == L"η")
            has_eta = true;
        },
        /*atoms_only=*/true);
    REQUIRE(!has_eta);
  }
}
