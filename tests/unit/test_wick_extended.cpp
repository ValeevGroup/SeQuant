#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/utility/indices.hpp>
#include <SeQuant/core/wick.hpp>
#include <SeQuant/core/wick_extended.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "catch2_sequant.hpp"

#include <array>

namespace sequant {

struct WickExtendedAccessor {};

/// the standard theorem's contractions under any vacuum, bypassing the
/// MultiProduct dispatch of compute()
template <>
template <>
struct WickTheorem<Statistics::FermiDirac>::access_by<WickExtendedAccessor> {
  ExprPtr compute_contractions(WickTheorem<Statistics::FermiDirac>& wick) {
    return wick.compute_contractions(/*count_only=*/false,
                                     /*skip_input_canonicalization=*/false);
  }
};

}  // namespace sequant

TEST_CASE("wick_extended", "[algorithms][wick][valgrind_skip]") {
  using namespace sequant;
  using namespace sequant::tests;

  auto ctx = get_default_context();
  ctx.set(mbpt::make_mr_spaces());
  ctx.set(Vacuum::MultiProduct);
  auto ctx_resetter = set_scoped_default_context(ctx);

  // helper: standard Wick with all partial contractions + provenance from the
  // (uncanonicalized) input; adequate here because every index is external
  auto wick_partial = [](const FNOperatorSeq& nopseq,
                         detail::OpProvenance& prov) {
    prov.clear();
    std::size_t ord = 0;
    for (const auto& nop : nopseq) {
      for (const auto& op : nop) prov.emplace(op.index(), ord);
      ++ord;
    }
    FWickTheorem wick{std::make_shared<FNOperatorSeq>(nopseq)};
    wick.full_contractions(false);
    return FWickTheorem::access_by<WickExtendedAccessor>{}.compute_contractions(
        wick);
  };

  // WickTheorem under the MultiProduct vacuum, configured by @p opts
  auto wick_mp = [](const ExprPtr& in,
                    const detail::ExtendedWickOptions& opts = {}) {
    FWickTheorem wick{in};
    wick.full_contractions(opts.full_contractions)
        .max_cumulant_rank(opts.max_cumulant_rank)
        .eta_as_delta_minus_gamma(opts.eta_as_delta_minus_gamma)
        .set_nop_connections(opts.nop_connections)
        .set_nop_avoided_connections(opts.nop_avoided_connections);
    return wick.compute();
  };

  SECTION("cumulant_expand: ⟨{a†a†a}{a}⟩ = κ2") {
    FNOperatorSeq in{FNOperator(cre({L"u_1", L"u_2"}), ann({L"u_3"})),
                     FNOperator(cre({}), ann({L"u_4"}))};
    detail::OpProvenance prov;
    auto wick_out = wick_partial(in, prov);
    auto result =
        detail::cumulant_expand<Statistics::FermiDirac>(wick_out, prov, {});
    REQUIRE_THAT(result, EquivalentTo(L"κ{u_4,u_3;u_1,u_2}"));
  }

  SECTION("cumulant_expand: ⟨{a†a}{a†a}⟩ = γη + κ2") {
    FNOperatorSeq in{FNOperator(cre({L"u_1"}), ann({L"u_2"})),
                     FNOperator(cre({L"u_3"}), ann({L"u_4"}))};
    detail::OpProvenance prov;
    auto wick_out = wick_partial(in, prov);
    auto result =
        detail::cumulant_expand<Statistics::FermiDirac>(wick_out, prov, {});
    REQUIRE_THAT(result, EquivalentTo(L"γ{u_4;u_1} * η{u_2;u_3} "
                                      L"+ κ{u_2,u_4;u_1,u_3}"));
  }

  SECTION("cumulant_expand: max_cumulant_rank = 1 means pairs only") {
    FNOperatorSeq in{FNOperator(cre({L"u_1"}), ann({L"u_2"})),
                     FNOperator(cre({L"u_3"}), ann({L"u_4"}))};
    detail::OpProvenance prov;
    auto wick_out = wick_partial(in, prov);
    auto result = detail::cumulant_expand<Statistics::FermiDirac>(
        wick_out, prov, {.max_cumulant_rank = 1});
    REQUIRE_THAT(result, EquivalentTo(L"γ{u_4;u_1} * η{u_2;u_3}"));
    auto result0 = detail::cumulant_expand<Statistics::FermiDirac>(
        wick_out, prov, {.max_cumulant_rank = 0});
    REQUIRE(simplify(result - result0) == ex<Constant>(0));
  }

  SECTION("cumulant_expand: input shapes") {
    // a c-number passes through unchanged
    REQUIRE(detail::cumulant_expand<Statistics::FermiDirac>(
                ex<Constant>(3), {}, {}) == ex<Constant>(3));
    REQUIRE(detail::cumulant_expand<Statistics::FermiDirac>(
                ex<Constant>(0), {}, {}) == ex<Constant>(0));
    // a bare NormalOperator is expanded like a Product holding it
    const auto nop =
        ex<FNOperator>(cre({L"u_1", L"u_2"}), ann({L"u_3", L"u_4"}));
    const detail::OpProvenance prov{{Index(L"u_1"), 0},
                                    {Index(L"u_3"), 0},
                                    {Index(L"u_2"), 1},
                                    {Index(L"u_4"), 1}};
    const auto bare =
        detail::cumulant_expand<Statistics::FermiDirac>(nop, prov, {});
    REQUIRE(bare != ex<Constant>(0));
    REQUIRE(bare == detail::cumulant_expand<Statistics::FermiDirac>(
                        ex<Product>(ExprPtrList{nop}), prov, {}));
    REQUIRE_THAT(bare, EquivalentTo(L"κ{u_3,u_4;u_1,u_2}"));
  }

  SECTION("cumulant_expand: unbalanced survivors vanish") {
    FNOperatorSeq in{FNOperator(cre({L"u_1", L"u_2"}), ann({})),
                     FNOperator(cre({}), ann({L"u_3"}))};
    detail::OpProvenance prov;
    auto wick_out = wick_partial(in, prov);
    auto result =
        detail::cumulant_expand<Statistics::FermiDirac>(wick_out, prov, {});
    REQUIRE(result == ex<Constant>(0));
  }

  SECTION("cumulant_expand: a block within one operator vanishes") {
    // {a†_u1 a†_u2 a_u3 a_u4}{a†_u5 a_u6}: the 4 legs of nop 0 may not form
    // a block by themselves; every κ must involve u_5 or u_6
    FNOperatorSeq in{FNOperator(cre({L"u_1", L"u_2"}), ann({L"u_3", L"u_4"})),
                     FNOperator(cre({L"u_5"}), ann({L"u_6"}))};
    detail::OpProvenance prov;
    auto wick_out = wick_partial(in, prov);
    auto result =
        detail::cumulant_expand<Statistics::FermiDirac>(wick_out, prov, {});
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
    detail::OpProvenance prov;
    auto wick_out = wick_partial(in, prov);
    auto full =
        detail::cumulant_expand<Statistics::FermiDirac>(wick_out, prov, {});
    auto trunc = detail::cumulant_expand<Statistics::FermiDirac>(
        wick_out, prov, {.max_cumulant_rank = 2});
    auto diff = simplify(full - trunc);
    // exactly the κ3 term
    REQUIRE_THAT(diff, EquivalentTo(L"-κ{u_5,u_6,u_3;u_1,u_2,u_4}"));
  }

  SECTION("cumulant_expand: partial contractions leave a GNO remainder") {
    FNOperatorSeq in{FNOperator(cre({L"u_1"}), ann({L"u_2"})),
                     FNOperator(cre({L"u_3"}), ann({L"u_4"}))};
    detail::OpProvenance prov;
    auto wick_out = wick_partial(in, prov);
    auto result = detail::cumulant_expand<Statistics::FermiDirac>(
        wick_out, prov, {.full_contractions = false});
    // Eq. (ext. Wick, 1-body×1-body): the Wick output's 4 terms plus κ2;
    // nothing else, because a block needs ≥2 ops from ≥2 nops and the only
    // such balanced set is all four legs
    REQUIRE_THAT(result, EquivalentTo(L"ã{u_2,u_4;u_1,u_3} "
                                      L"- γ{u_4;u_1} * ã{u_2;u_3} "
                                      L"+ η{u_2;u_3} * ã{u_4;u_1} "
                                      L"+ γ{u_4;u_1} * η{u_2;u_3} "
                                      L"+ κ{u_2,u_4;u_1,u_3}"));
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
    detail::OpProvenance prov;
    auto wick_out = wick_partial(in, prov);
    REQUIRE(detail::cumulant_expand<Statistics::FermiDirac>(
                wick_out, prov, {}) == ex<Constant>(0));
    auto partial = detail::cumulant_expand<Statistics::FermiDirac>(
        wick_out, prov, {.full_contractions = false});
    REQUIRE(simplify(partial - wick_out) == ex<Constant>(0));
  }

  SECTION("WickTheorem: the operators must use the MultiProduct vacuum") {
    // an operator normal-ordered relative to another vacuum is not a GNO
    // string; it is rejected rather than reinterpreted
    for (const auto vacuum : {Vacuum::Physical, Vacuum::SingleProduct}) {
      auto in = ex<FNOperator>(cre({L"u_1"}), ann({L"u_2"}), vacuum) *
                ex<FNOperator>(cre({L"u_3"}), ann({L"u_4"}), vacuum);
      REQUIRE_THROWS_AS(FWickTheorem{in}.compute(), Exception);
      REQUIRE_THROWS_AS(FWickTheorem{in}.full_contractions(false).compute(),
                        Exception);
    }
    // the elementary-operator spelling of the parser is Physical
    REQUIRE_THROWS_AS(
        FWickTheorem{deserialize(L"a{u_1;u_2} a{u_3;u_4}")}.compute(),
        Exception);
  }

  SECTION("WickTheorem: bosons are not supported") {
    auto in = ex<BNOperator>(cre({L"u_1"}), ann({L"u_2"}), Vacuum::Physical) *
              ex<BNOperator>(cre({L"u_3"}), ann({L"u_4"}), Vacuum::Physical);
    REQUIRE_THROWS_MATCHES(BWickTheorem{in}.compute(), Exception,
                           Catch::Matchers::MessageMatches(
                               Catch::Matchers::ContainsSubstring("bosons")));
  }

  SECTION("WickTheorem: count_only is not supported") {
    auto in = ex<FNOperator>(cre({L"u_1"}), ann({L"u_2"})) *
              ex<FNOperator>(cre({L"u_3"}), ann({L"u_4"}));
    REQUIRE_THROWS_MATCHES(
        FWickTheorem{in}.compute(/*count_only=*/true), Exception,
        Catch::Matchers::MessageMatches(
            Catch::Matchers::ContainsSubstring("count_only")));
  }

  SECTION("WickTheorem: spin-free operators are not supported") {
    auto sf_ctx = get_default_context();
    sf_ctx.set(SPBasis::Spinfree);
    auto sf_resetter = set_scoped_default_context(sf_ctx);
    auto in = ex<FNOperator>(cre({L"u_1"}), ann({L"u_2"})) *
              ex<FNOperator>(cre({L"u_3"}), ann({L"u_4"}));
    REQUIRE_THROWS_MATCHES(
        FWickTheorem{in}.compute(), Exception,
        Catch::Matchers::MessageMatches(
            Catch::Matchers::ContainsSubstring("spin-free")));
  }

  SECTION("WickTheorem: the input is not modified") {
    // the coefficient's bra is paired with the creator, as mbpt builds them
    auto in = ex<FNOperator>(cre({L"u_1"}), ann({L"u_2"})) *
              ex<Tensor>(L"h", bra{L"u_1"}, ket{L"u_2"}, Symmetry::Nonsymm,
                         BraKetSymmetry::Nonsymm, ColumnSymmetry::Symm);
    wick_mp(in);
    REQUIRE(in->as<Product>().factor(0)->is<FNOperator>());
  }

  SECTION("WickTheorem: pure-active identities via the wrapper") {
    auto in = ex<FNOperator>(cre({L"u_1"}), ann({L"u_2"})) *
              ex<FNOperator>(cre({L"u_3"}), ann({L"u_4"}));
    auto result = wick_mp(in);
    REQUIRE_THAT(result, EquivalentTo(L"γ{u_4;u_1} * η{u_2;u_3} "
                                      L"+ κ{u_2,u_4;u_1,u_3}"));
  }

  SECTION("WickTheorem: general indices split into core δ + active γ") {
    // ⟨{a†_p1 a_p2}{a†_p3 a_p4}⟩ with p = M ∪ E = {o,i,u,a,g}:
    // cre·ann pair over R = M: δ on core O + γ on active u;
    // ann·cre pair over U = E: δ on virtual {a,g} + η on active u, where
    // {a,g} is not a registered space, so its δ splits into δ on a + δ on g;
    // plus κ2 on the all-active projection: 2 × 3 + 1 = 7 terms
    auto in = ex<FNOperator>(cre({L"p_1"}), ann({L"p_2"})) *
              ex<FNOperator>(cre({L"p_3"}), ann({L"p_4"}));
    auto result = wick_mp(in);
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
    require_active_densities(wick_mp(in, {.full_contractions = false}));
    REQUIRE(result->size() == 7);
    // projecting is substituting a_p = Σ_x δ(p,x) a_x, so every term keeps
    // the + of γ{p_4;p_1}·η{p_2;p_3} + κ{p_2,p_4;p_1,p_3}, with γ over core
    // = δ and η over virtual = δ
    REQUIRE_THAT(
        result,
        EquivalentTo(
            L"κ{u_3,u_4;u_1,u_2} * δ{u_1;p_1}:N-C-S * "
            L"δ{u_2;p_3}:N-C-S * δ{p_2;u_3}:N-C-S * δ{p_4;u_4}:N-C-S "
            L"+ γ{u_2;u_1} * δ{u_1;p_1}:N-C-S * δ{a_1;p_3}:N-C-S * "
            L"δ{p_2;a_1}:N-C-S * δ{p_4;u_2}:N-C-S "
            L"+ γ{u_2;u_1} * δ{u_1;p_1}:N-C-S * δ{g_1;p_3}:N-C-S * "
            L"δ{p_2;g_1}:N-C-S * δ{p_4;u_2}:N-C-S "
            L"+ η{u_2;u_1} * δ{u_1;p_3}:N-C-S * δ{O_1;p_1}:N-C-S * "
            L"δ{p_2;u_2}:N-C-S * δ{p_4;O_1}:N-C-S "
            L"+ δ{a_1;p_3}:N-C-S * δ{O_1;p_1}:N-C-S * δ{p_2;a_1}:N-C-S * "
            L"δ{p_4;O_1}:N-C-S "
            L"+ δ{g_1;p_3}:N-C-S * δ{O_1;p_1}:N-C-S * δ{p_2;g_1}:N-C-S * "
            L"δ{p_4;O_1}:N-C-S "
            L"+ η{u_3;u_1} * γ{u_4;u_2} * δ{u_1;p_3}:N-C-S * "
            L"δ{u_2;p_1}:N-C-S * δ{p_2;u_3}:N-C-S * δ{p_4;u_4}:N-C-S"));
  }

  SECTION("WickTheorem: an input density over a wider space is split") {
    // a one-body γ (η) of the input is the reference (hole) density over
    // whatever space its indices range over: like the ones the theorem
    // produces it is δ on the core (virtual) space, γ (η) on the active one
    // and zero on the virtual (core) one
    auto h = [](const Index& b, const Index& k) {
      return ex<Tensor>(L"h", bra{b}, ket{k}, Symmetry::Nonsymm,
                        BraKetSymmetry::Conjugate, ColumnSymmetry::Symm);
    };
    const Index p1(L"p_1"), p2(L"p_2"), O1(L"O_1"), u5(L"u_5"), u6(L"u_6"),
        a1(L"a_1"), g1(L"g_1");
    auto ops = ex<FNOperator>(cre({L"u_1"}), ann({L"u_2"})) *
               ex<FNOperator>(cre({L"u_3"}), ann({L"u_4"}));
    const auto ops_av = wick_mp(ops);
    ExprPtr gamma_result;
    REQUIRE_NOTHROW(gamma_result =
                        wick_mp(density::make_rdm(p1, p2) * h(p2, p1) * ops));
    const auto gamma_coeff = h(O1, O1) + h(u6, u5) * density::make_rdm(u5, u6);
    REQUIRE(simplify(gamma_result - gamma_coeff * ops_av) == ex<Constant>(0));
    ExprPtr eta_result;
    REQUIRE_NOTHROW(
        eta_result = wick_mp(density::make_hole_rdm(p1, p2) * h(p2, p1) * ops));
    // the virtual space {a,g} is not registered, so its δ splits over a and g
    const auto eta_coeff =
        h(a1, a1) + h(g1, g1) + h(u6, u5) * density::make_hole_rdm(u5, u6);
    REQUIRE(simplify(eta_result - eta_coeff * ops_av) == ex<Constant>(0));
  }

  SECTION("WickTheorem: a split density keeps the bases of its indices") {
    // η{a_1<i_1>;a_2<i_2>}: its virtual block binds a_1 to an index in its
    // basis, so the result keeps the overlap between the two bases and every
    // δ is between indices of one basis
    const Index i1(L"i_1"), i2(L"i_2"), a1(L"a_1", {i1}), a2(L"a_2", {i2});
    auto in = density::make_hole_rdm(a1, a2) *
              ex<FNOperator>(cre({L"u_3"}), ann({L"u_4"})) *
              ex<FNOperator>(cre({L"u_5"}), ann({L"u_6"}));
    ExprPtr result;
    REQUIRE_NOTHROW(result = wick_mp(in));
    std::size_t noverlaps = 0;
    result->visit(
        [&](const ExprPtr& e) {
          if (!e->is<Tensor>()) return;
          const auto& t = e->as<Tensor>();
          INFO(toUtf8(to_latex(e)));
          if (t.label() == reserved::kronecker_label())
            REQUIRE(t.bra()[0].proto_indices() == t.ket()[0].proto_indices());
          if (t.label() == reserved::overlap_label()) {
            ++noverlaps;
            REQUIRE(container::set<Index>{t.bra()[0], t.ket()[0]} ==
                    container::set<Index>{a1, a2});
          }
        },
        /* atoms_only = */ true);
    REQUIRE(noverlaps > 0);

    // as for an operator, a protoindexed index of a density must not reach
    // the active space
    if (assert_behavior() == AssertBehavior::Throw) {
      const Index p1(L"p_1", {i1}), p2(L"p_2", {i2});
      REQUIRE_THROWS_WITH(wick_mp(density::make_rdm(p1, p2) *
                                  ex<FNOperator>(cre({L"u_3"}), ann({L"u_4"})) *
                                  ex<FNOperator>(cre({L"u_5"}), ann({L"u_6"}))),
                          Catch::Matchers::ContainsSubstring(
                              "must not reach the active space"));
    }
  }

  SECTION("WickTheorem: dummy indices keep their provenance") {
    // a one-body h summed against its operator's indices, times an active
    // one-body operator: every op index of the first factor is a dummy
    const Index p1(L"p_1"), p2(L"p_2");
    auto in = ex<Tensor>(L"h", bra{p1}, ket{p2}, Symmetry::Nonsymm,
                         BraKetSymmetry::Conjugate, ColumnSymmetry::Symm) *
              ex<FNOperator>(cre({p1}), ann({p2})) *
              ex<FNOperator>(cre({L"u_3"}), ann({L"u_4"}));
    ExprPtr result;
    REQUIRE_NOTHROW(result = wick_mp(in));
    // every γ/η/κ index is active
    for (const auto& term : *result)
      for (const auto& f : *term)
        if (f->is<Tensor>()) {
          const auto& t = f->as<Tensor>();
          if (t.label() == L"γ" || t.label() == L"η" || t.label() == L"κ")
            for (const auto& idx : t.const_braket())
              REQUIRE(idx.space() == Index(L"u_1").space());
        }
    REQUIRE_NOTHROW(wick_mp(in, {.full_contractions = false}));
  }

  SECTION("WickTheorem: δs over summed indices are applied") {
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
    auto result = wick_mp(in);
    REQUIRE(simplify(result - h(O1, O1) -
                     h(u2, u1) * density::make_rdm(u1, u2)) == ex<Constant>(0));
    // ... also for a δ that η = δ - γ introduces
    auto in2 = h(u3, u2) * ex<FNOperator>(cre({L"u_1"}), ann({L"u_2"})) *
               ex<FNOperator>(cre({L"u_3"}), ann({L"u_4"}));
    auto result2 = wick_mp(in2, {.eta_as_delta_minus_gamma = true});
    result2->visit(
        [](const ExprPtr& e) {
          if (e->is<Tensor>()) REQUIRE(e->as<Tensor>().label() != L"δ");
        },
        /*atoms_only=*/true);
  }

  SECTION("WickTheorem: an index shared by two operators") {
    // ⟨{a†_u1 a†_u2}{a_u2 a_u1}⟩ is ⟨{a†_u1 a†_u2}{a_u3 a_u4}⟩ with u_3 = u_2
    // and u_4 = u_1
    const Index u1(L"u_1"), u2(L"u_2"), u3(L"u_3"), u4(L"u_4");
    const container::map<Index, Index> identify{{u3, u2}, {u4, u1}};
    auto rename = [&](ExprPtr expr) {
      expr->visit(
          [&](const ExprPtr& e) {
            if (e->is<Tensor>()) {
              e->as<Tensor>().transform_indices(identify);
              e->as<Tensor>().reset_tags();
            } else if (e->is<FNOperator>()) {
              e->as<FNOperator>().transform_indices(identify);
              for (const auto& op : e->as<FNOperator>()) op.index().reset_tag();
            }
          },
          /*atoms_only=*/true);
      return simplify(expr);
    };
    // the shared indices are annihilators, then creators, of the later
    // operator: ⟨{a_u1 a_u2}{a†_u2 a†_u1}⟩ is ⟨{a_u1 a_u2}{a†_u3 a†_u4}⟩ with
    // the same identification
    const std::pair<ExprPtr, ExprPtr> shared_distinct[] = {
        {ex<FNOperator>(cre({u1, u2}), ann({})) *
             ex<FNOperator>(cre({}), ann({u2, u1})),
         ex<FNOperator>(cre({u1, u2}), ann({})) *
             ex<FNOperator>(cre({}), ann({u3, u4}))},
        {ex<FNOperator>(cre({}), ann({u1, u2})) *
             ex<FNOperator>(cre({u2, u1}), ann({})),
         ex<FNOperator>(cre({}), ann({u1, u2})) *
             ex<FNOperator>(cre({u3, u4}), ann({}))}};
    for (const auto& [shared, distinct] : shared_distinct) {
      for (const bool full : {true, false}) {
        const detail::ExtendedWickOptions opts{.full_contractions = full};
        auto result = wick_mp(shared, opts);
        auto expected = rename(wick_mp(distinct, opts));
        INFO("input: " << toUtf8(to_latex(shared)) << "\nfull=" << full
                       << "\nresult: " << toUtf8(to_latex(result))
                       << "\nexpected: " << toUtf8(to_latex(expected)));
        REQUIRE(simplify(result - expected) == ex<Constant>(0));
        bool has_kappa = false;
        result->visit(
            [&](const ExprPtr& e) {
              if (e->is<Tensor>() && e->as<Tensor>().label() == L"κ")
                has_kappa = true;
            },
            /*atoms_only=*/true);
        REQUIRE(has_kappa);
      }
    }
    // the shared index alone does not connect the operators
    auto two =
        ex<FNOperator>(cre({u1}), ann({})) * ex<FNOperator>(cre({}), ann({u1}));
    auto disconnected = wick_mp(
        two, {.full_contractions = false, .nop_avoided_connections = {{0, 1}}});
    REQUIRE(disconnected->is<FNOperator>());
    REQUIRE(disconnected->as<FNOperator>().size() == 2);
  }

  SECTION("WickTheorem: tensors commute with the theorem") {
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
      const detail::ExtendedWickOptions opts{.full_contractions = full};
      auto lhs = wick_mp(tensors * ops, opts);
      auto rhs = simplify(tensors * wick_mp(ops, opts));
      REQUIRE(simplify(lhs - rhs) == ex<Constant>(0));
    }
  }

  SECTION("WickTheorem: projected survivors with partial contractions") {
    // {a†_p1 a_p2}{a†_u3 a_u4}: projecting an op is substituting
    // a_p = Σ_x δ(p,x) a_x over the active, core and virtual parts of p, so
    // each term of the pure-active result (see "partial contractions leave a
    // GNO remainder") reappears with its sign, its indices projected
    auto in = ex<FNOperator>(cre({L"p_1"}), ann({L"p_2"})) *
              ex<FNOperator>(cre({L"u_3"}), ann({L"u_4"}));
    auto result = wick_mp(in, {.full_contractions = false});
    // result contains `term` exactly once, with its sign
    auto contains = [&](std::wstring_view term) {
      return simplify(result - deserialize(term))->size() == result->size() - 1;
    };
    // a_p2·a†_u3 contracted (adjacent, so +), p_1 of the surviving
    // {a†_p1 a_u4} projected onto active and onto virtual a
    REQUIRE(
        contains(L"η{u_2;u_3} * δ{u_1;p_1}:N-C-S * δ{p_2;u_2}:N-C-S * "
                 L"ã{u_4;u_1}"));
    REQUIRE(
        contains(L"η{u_1;u_3} * δ{a_1;p_1}:N-C-S * δ{p_2;u_1}:N-C-S * "
                 L"ã{u_4;a_1}"));
    // nothing contracted, p_1 and p_2 projected onto active: the four legs
    // form +κ{u_2,u_4;u_1,u_3} = -κ{u_4,u_2;u_1,u_3}
    REQUIRE(
        contains(L"-1 κ{u_4,u_2;u_1,u_3} * δ{u_1;p_1}:N-C-S * "
                 L"δ{p_2;u_2}:N-C-S"));
  }

  SECTION("WickTheorem: projected indices are fresh in their term") {
    // {a†_u1 a_i1}{a†_p1 a_A1}: the active projection of the surviving a†_p1
    // must not reuse the name the engine's reduction gave the dummy of
    // γ·δ(A_1), however far the global tmp-index counter has advanced
    auto in = ex<FNOperator>(cre({L"u_1"}), ann({L"i_1"})) *
              ex<FNOperator>(cre({L"p_1"}), ann({L"A_1"}));
    Index::reset_tmp_index();
    ExprPtr first;
    REQUIRE_NOTHROW(first = wick_mp(in, {.full_contractions = false}));
    const auto& u =
        get_default_context().index_basis_registry()->retrieve(L"u");
    for (int i = 0; i != 1000; ++i) Index::make_tmp_index(u);
    const auto second = wick_mp(in, {.full_contractions = false});
    REQUIRE(simplify(first - second) == ex<Constant>(0));
  }

  SECTION("WickTheorem: connectivity") {
    // partial contractions: with full ones, connecting 0 to 1 but not to 2
    // leaves 2 isolated, and no term survives
    auto in = ex<FNOperator>(cre({L"u_1"}), ann({L"u_2"})) *
              ex<FNOperator>(cre({L"u_3"}), ann({L"u_4"})) *
              ex<FNOperator>(cre({L"u_5"}), ann({L"u_6"}));
    auto all = wick_mp(in, {.full_contractions = false});
    auto filtered = wick_mp(in, {.full_contractions = false,
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

  SECTION("WickTheorem: connectivity through a cumulant") {
    // partial contractions: a 0-2 connection can be realized by a cumulant
    // block alone, which no pair of operators expresses
    auto in = ex<FNOperator>(cre({L"u_1"}), ann({L"u_2"})) *
              ex<FNOperator>(cre({L"u_3"}), ann({L"u_4"})) *
              ex<FNOperator>(cre({L"u_5"}), ann({L"u_6"}));
    auto all = wick_mp(in, {.full_contractions = false});
    auto filtered =
        wick_mp(in, {.full_contractions = false, .nop_connections = {{0, 2}}});
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

  SECTION("WickTheorem: a coefficient tensor is not a connection") {
    // h{u_1;u_6} spans operators 0 and 2, but only the factors the theorem
    // produces connect operators
    auto h = ex<Tensor>(L"h", bra{L"u_1"}, ket{L"u_6"}, Symmetry::Nonsymm,
                        BraKetSymmetry::Nonsymm, ColumnSymmetry::Symm);
    auto ops = ex<FNOperator>(cre({L"u_1"}), ann({L"u_2"})) *
               ex<FNOperator>(cre({L"u_3"}), ann({L"u_4"})) *
               ex<FNOperator>(cre({L"u_5"}), ann({L"u_6"}));
    for (const auto& opts :
         {detail::ExtendedWickOptions{.full_contractions = false,
         .nop_avoided_connections = {{0, 2}}},
          detail::ExtendedWickOptions{.full_contractions = false,
          .nop_connections = {{0, 2}}}}) {
      auto without_h = wick_mp(ops, opts);
      REQUIRE(without_h->size() > 0);
      auto lhs = wick_mp(h * ops, opts);
      REQUIRE(simplify(lhs - simplify(h * without_h)) == ex<Constant>(0));
    }
  }

  SECTION("WickTheorem: an input δ is an identification, not a connection") {
    // δ{u_2;u_5} (or s{u_2;u_5}) makes operators 0 and 2 share an index,
    // which neither realizes nor avoids a 0-2 connection
    const Index u2(L"u_2"), u5(L"u_5");
    auto ops = [](const Index& i5) {
      return ex<FNOperator>(cre({Index(L"u_1")}), ann({Index(L"u_2")})) *
             ex<FNOperator>(cre({Index(L"u_3")}), ann({Index(L"u_4")})) *
             ex<FNOperator>(cre({i5}), ann({Index(L"u_6")}));
    };
    for (const auto& opts :
         {detail::ExtendedWickOptions{.full_contractions = false,
         .nop_avoided_connections = {{0, 2}}},
          detail::ExtendedWickOptions{.full_contractions = false,
          .nop_connections = {{0, 2}}}}) {
      auto shared = wick_mp(ops(u2), opts);
      REQUIRE(shared->size() > 0);
      for (const auto& id : {make_kronecker(u2, u5), make_overlap(u2, u5)}) {
        auto result = wick_mp(id * ops(u5), opts);
        INFO("δ: " << toUtf8(to_latex(id))
                   << "\nresult: " << toUtf8(to_latex(result))
                   << "\nexpected: " << toUtf8(to_latex(shared)));
        // compared as is: simplify does not cancel a remainder carrying an
        // index as both a creator and an annihilator against itself
        REQUIRE(*result == *shared);
      }
    }
    // a δ binding an op index to a coefficient's is a rename
    auto h = [](const Index& b) {
      return ex<Tensor>(L"h", bra{b}, ket{L"p_2"}, Symmetry::Nonsymm,
                        BraKetSymmetry::Nonsymm, ColumnSymmetry::Symm);
    };
    const Index p1(L"p_1"), p7(L"p_7");
    auto pops = ex<FNOperator>(cre({p1}), ann({Index(L"p_2")})) *
                ex<FNOperator>(cre({Index(L"u_3")}), ann({Index(L"u_4")}));
    for (const bool full : {true, false}) {
      const detail::ExtendedWickOptions opts{.full_contractions = full};
      auto renamed = wick_mp(h(p1) * pops, opts);
      auto result = wick_mp(make_kronecker(p1, p7) * h(p7) * pops, opts);
      REQUIRE(simplify(result - renamed) == ex<Constant>(0));
    }
    // a δ between external indices multiplies the result
    const auto ext = make_kronecker(Index(L"p_8"), Index(L"p_9"));
    for (const auto& opts :
         {detail::ExtendedWickOptions{.full_contractions = false,
         .nop_avoided_connections = {{0, 2}}},
          detail::ExtendedWickOptions{.full_contractions = false,
          .nop_connections = {{0, 2}}}}) {
      auto without = wick_mp(ops(u5), opts);
      auto result = wick_mp(ext * ops(u5), opts);
      REQUIRE(simplify(result - simplify(ext * without)) == ex<Constant>(0));
    }
  }

  SECTION("WickTheorem: use_topology with a context-named dummy") {
    // the context may name an index that appears twice; it is then external,
    // so the pruning may not treat the operators carrying it as equivalent
    auto coeff = [](std::wstring_view label, IndexList b, IndexList k) {
      return ex<Tensor>(label, bra(b), ket(k), Symmetry::Antisymm,
                        BraKetSymmetry::Nonsymm, ColumnSymmetry::Symm);
    };
    auto in = coeff(L"g", {L"u_1", L"u_2"}, {L"u_3", L"u_4"}) *
              ex<FNOperator>(cre({L"u_1", L"u_2"}), ann({L"u_3", L"u_4"})) *
              coeff(L"t", {L"u_5", L"u_6"}, {L"u_7", L"u_8"}) *
              ex<FNOperator>(cre({L"u_5", L"u_6"}), ann({L"u_7", L"u_8"}));
    auto scope = scoped_canonicalize_options(
        CanonicalizeOptions::default_options().copy_and_set(
            container::set<Index>{Index(L"u_1"), Index(L"u_5")}));
    for (const bool full : {true, false}) {
      INFO("full=" << full);
      std::array<ExprPtr, 2> result;
      for (const bool top : {false, true}) {
        FWickTheorem wick{in->clone()};
        result[top] = wick.full_contractions(full).use_topology(top).compute();
      }
      REQUIRE(get_used_indices_with_counts(result[1]).contains(Index(L"u_1")));
      REQUIRE(simplify(result[1] - result[0]) == ex<Constant>(0));
    }
  }

  SECTION("WickTheorem: use_topology") {
    // WickTheorem under the MultiProduct vacuum with use_topology(@p top);
    // returns the result and the number of attempted contractions
    auto run = [](const ExprPtr& in, bool full, bool top,
                  const container::svector<std::pair<std::size_t, std::size_t>>&
                      connections = {}) {
      FWickTheorem wick{in};
      wick.full_contractions(full).use_topology(top).set_nop_connections(
          connections);
      auto result = wick.compute();
      return std::pair{result, wick.stats().num_attempted_contractions.load()};
    };
    auto coeff = [](std::wstring_view label, IndexList b, IndexList k) {
      return ex<Tensor>(label, bra(b), ket(k), Symmetry::Antisymm,
                        BraKetSymmetry::Nonsymm, ColumnSymmetry::Symm);
    };
    // antisymmetric 2-body coefficients summed against all-active operators:
    // the creators (annihilators) of each operator are equivalent
    auto g_op = coeff(L"g", {L"u_1", L"u_2"}, {L"u_3", L"u_4"}) *
                ex<FNOperator>(cre({L"u_1", L"u_2"}), ann({L"u_3", L"u_4"}));
    auto t_op = coeff(L"t", {L"u_5", L"u_6"}, {L"u_7", L"u_8"}) *
                ex<FNOperator>(cre({L"u_5", L"u_6"}), ann({L"u_7", L"u_8"}));
    for (const bool full : {true, false}) {
      INFO("full=" << full);
      const auto [on, attempted_on] = run(g_op * t_op, full, true);
      const auto [off, attempted_off] = run(g_op * t_op, full, false);
      REQUIRE(simplify(on - off) == ex<Constant>(0));
      REQUIRE(attempted_on > 0);
      REQUIRE(attempted_on < attempted_off);
    }

    // an index-free factor does not take part in the topology analysis
    {
      const auto c_g_t = ex<Variable>(L"c") * g_op * t_op;
      ExprPtr on, off;
      REQUIRE_NOTHROW(on = run(c_g_t, true, true).first);
      off = run(c_g_t, true, false).first;
      REQUIRE(simplify(on - off) == ex<Constant>(0));
    }

    // operators 1 and 2 are equivalent, but only 1 must connect to 0
    auto t1_op = [&](std::wstring_view b, std::wstring_view k) {
      return ex<Tensor>(L"t", bra{b}, ket{k}, Symmetry::Nonsymm,
                        BraKetSymmetry::Nonsymm, ColumnSymmetry::Symm) *
             ex<FNOperator>(cre({b}), ann({k}));
    };
    auto three = g_op * t1_op(L"u_5", L"u_6") * t1_op(L"u_7", L"u_8");
    for (const bool full : {true, false}) {
      INFO("full=" << full);
      const auto [on, attempted_on] = run(three, full, true, {{0, 1}});
      const auto [off, attempted_off] = run(three, full, false, {{0, 1}});
      REQUIRE(simplify(on - off) == ex<Constant>(0));
      REQUIRE(attempted_on < attempted_off);
    }

    // an input that vanishes by symmetry (a symmetric tensor contracted with
    // the antisymmetric creators of an operator) is returned as zero with no
    // contraction attempted, as under the other vacua; the exhaustive
    // enumeration produces terms that each vanish by their own symmetry,
    // which the canonicalizer does not yet recognize
    {
      auto zero = ex<Tensor>(L"X", bra{L"u_1", L"u_2"}, ket{L"u_3", L"u_4"},
                             Symmetry::Symm, BraKetSymmetry::Nonsymm,
                             ColumnSymmetry::Symm) *
                  ex<FNOperator>(cre({L"u_1", L"u_2"}), ann({L"u_3", L"u_4"}));
      for (const bool full : {true, false}) {
        INFO("full=" << full);
        const auto [on, attempted_on] = run(zero * t_op, full, true);
        REQUIRE(on == ex<Constant>(0));
        REQUIRE(attempted_on == 0);
      }
    }

    // a†_u1 and a†_u2 are only equivalent together with the tensors they
    // are attached to, so nothing is pruned
    auto X = [](std::wstring_view b, std::wstring_view k) {
      return ex<Tensor>(L"X", bra{b}, ket{k}, Symmetry::Nonsymm,
                        BraKetSymmetry::Nonsymm, ColumnSymmetry::Symm);
    };
    auto xx = X(L"u_1", L"u_3") * X(L"u_2", L"u_4") *
              ex<FNOperator>(cre({}), ann({L"u_3"})) *
              ex<FNOperator>(cre({}), ann({L"u_4"})) *
              ex<FNOperator>(cre({L"u_1", L"u_2"}), ann({}));
    for (const bool full : {true, false}) {
      INFO("full=" << full);
      REQUIRE(simplify(run(xx, full, true).first -
                       run(xx, full, false).first) == ex<Constant>(0));
    }
  }

  SECTION("WickTheorem: an avoided pair is rejected during contraction") {
    // a pair contraction between an avoided pair of operators is fatal the
    // moment it is attempted (no later cumulant can undo it), so the engine
    // rejects it early and attempts fewer contractions; the result must still
    // be exactly what filtering the unconstrained result would give
    auto in = ex<FNOperator>(cre({L"u_1"}), ann({L"u_2"})) *
              ex<FNOperator>(cre({L"u_3"}), ann({L"u_4"})) *
              ex<FNOperator>(cre({L"u_5"}), ann({L"u_6"}));
    auto run = [](const ExprPtr& in,
                  const container::svector<std::pair<std::size_t, std::size_t>>&
                      avoided) {
      FWickTheorem wick{in};
      wick.full_contractions(false).set_nop_avoided_connections(avoided);
      auto result = wick.compute();
      return std::pair{result, wick.stats().num_attempted_contractions.load()};
    };
    auto [all, attempted_all] = run(in, {});
    auto [avoided, attempted_avoided] = run(in, {{0, 2}});
    REQUIRE(attempted_avoided < attempted_all);

    // the terms of the unconstrained result without a 0-2 edge
    auto ord = [](const Index& i) -> int {
      const auto l = i.label();
      if (l == L"u_1" || l == L"u_2") return 0;
      if (l == L"u_3" || l == L"u_4") return 1;
      return 2;
    };
    auto expected = std::make_shared<Sum>();
    for (const auto& term : *all) {
      bool e02 = false;
      auto check = [&](const ExprPtr& f) {
        if (!f->is<Tensor>()) return;
        container::set<int> ords;
        for (const auto& idx : f->as<Tensor>().const_braket())
          ords.insert(ord(idx));
        if (ords.contains(0) && ords.contains(2)) e02 = true;
      };
      if (term->is<Product>())
        for (const auto& f : *term) check(f);
      else
        check(term);
      if (!e02) expected->append(term);
    }
    REQUIRE(expected->size() > 0);
    REQUIRE(simplify(avoided - ExprPtr(expected)) == ex<Constant>(0));
  }

  SECTION("WickTheorem: connection ordinals must name input operators") {
    auto in = ex<FNOperator>(cre({L"u_1"}), ann({L"u_2"})) *
              ex<FNOperator>(cre({L"u_3"}), ann({L"u_4"}));
    REQUIRE_THROWS_AS(wick_mp(in, {.nop_connections = {{0, 2}}}), Exception);
    REQUIRE_THROWS_AS(wick_mp(in, {.nop_avoided_connections = {{2, 1}}}),
                      Exception);
    // a term without operators is kept, as under the other vacua
    const auto c = ex<Constant>(2) *
                   ex<Tensor>(L"h", bra{L"u_5"}, ket{L"u_6"}, Symmetry::Nonsymm,
                              BraKetSymmetry::Conjugate, ColumnSymmetry::Symm);
    const detail::ExtendedWickOptions connected{.nop_connections = {{0, 1}}};
    REQUIRE(simplify(wick_mp(in + c, connected) - wick_mp(in, connected) - c) ==
            ex<Constant>(0));
  }

  SECTION("WickTheorem: Sum input") {
    auto in = ex<FNOperator>(cre({L"u_1"}), ann({L"u_2"})) *
                  ex<FNOperator>(cre({L"u_3"}), ann({L"u_4"})) +
              ex<Constant>(2) * ex<FNOperator>(cre({L"u_1"}), ann({L"u_2"})) *
                  ex<FNOperator>(cre({L"u_3"}), ann({L"u_4"}));
    auto result = wick_mp(in);
    REQUIRE_THAT(result, EquivalentTo(L"3 γ{u_4;u_1} * η{u_2;u_3} "
                                      L"+ 3 κ{u_2,u_4;u_1,u_3}"));
  }

  SECTION("WickTheorem: η between two CSV bases is their overlap") {
    // the virtual block of η is the identity on the virtuals; between
    // cluster-specific virtuals of different pairs it is their overlap, as
    // under a single-product vacuum
    const Index i1(L"i_1"), i2(L"i_2");
    const Index a1(L"a_1", {i1}), a2(L"a_2", {i2}), a3(L"a_3", {i1});
    auto in =
        ex<FNOperator>(cre({}), ann({a1})) * ex<FNOperator>(cre({a2}), ann({}));
    REQUIRE_THAT(wick_mp(in), EquivalentTo(L"s{a_1<i_1>;a_2<i_2>}"));
    in =
        ex<FNOperator>(cre({}), ann({a1})) * ex<FNOperator>(cre({a3}), ann({}));
    REQUIRE_THAT(wick_mp(in), EquivalentTo(L"δ{a_1<i_1>;a_3<i_1>}"));
  }

  SECTION("WickTheorem: η = δ - γ") {
    auto in = ex<FNOperator>(cre({L"u_1"}), ann({L"u_2"})) *
              ex<FNOperator>(cre({L"u_3"}), ann({L"u_4"}));
    auto result = wick_mp(in, {.eta_as_delta_minus_gamma = true});
    REQUIRE_THAT(result, EquivalentTo(L"γ{u_4;u_1} * δ{u_2;u_3} "
                                      L"- γ{u_4;u_1} * γ{u_2;u_3} "
                                      L"+ κ{u_2,u_4;u_1,u_3}"));
    bool has_eta = false;
    result->visit(
        [&](const ExprPtr& e) {
          if (e->is<Tensor>() && e->as<Tensor>().label() == L"η")
            has_eta = true;
        },
        /*atoms_only=*/true);
    REQUIRE(!has_eta);
  }

  SECTION("WickTheorem: η = δ - γ leaves a multi-body η alone") {
    // only the one-body η the theorem produces is δ - γ; a multi-body η in
    // the input (e.g. from an earlier result) is kept as is
    const auto eta2 = density::make_density(
        density::hole_rdm_label(), bra{L"u_5", L"u_6"}, ket{L"u_7", L"u_8"});
    auto in = eta2 * ex<FNOperator>(cre({L"u_1"}), ann({L"u_2"})) *
              ex<FNOperator>(cre({L"u_3"}), ann({L"u_4"}));
    ExprPtr result;
    REQUIRE_NOTHROW(result = wick_mp(in, {.eta_as_delta_minus_gamma = true}));
    std::size_t n_eta1 = 0, n_eta2 = 0;
    result->visit(
        [&](const ExprPtr& e) {
          if (e->is<Tensor>() && e->as<Tensor>().label() == L"η")
            (e->as<Tensor>().rank() == 1 ? n_eta1 : n_eta2) += 1;
        },
        /*atoms_only=*/true);
    REQUIRE(n_eta1 == 0);
    REQUIRE(n_eta2 == result->size());
  }

  SECTION("WickTheorem: single-reference limit") {
    // without an active space (reference occupancy == vacuum occupancy)
    // MultiProduct reduces to SingleProduct. The SR registry declares both
    // occupancies, which MultiProduct requires.
    auto sr_isr = mbpt::make_sr_spaces();
    REQUIRE(sr_isr->reference_occupied_space() ==
            sr_isr->vacuum_occupied_space());

    // inputs carry their context's vacuum, so each is built per context
    auto make_input = [](std::size_t k) {
      return k == 0 ? ex<FNOperator>(cre({L"i_1"}), ann({L"a_1"})) *
                          ex<FNOperator>(cre({L"a_2"}), ann({L"i_2"}))
                    : ex<FNOperator>(cre({L"p_1"}), ann({L"p_2"})) *
                          ex<FNOperator>(cre({L"p_3"}), ann({L"p_4"}));
    };
    auto sr_ctx = get_default_context();
    sr_ctx.set(sr_isr);
    auto with_vacuum = [&](Vacuum v) {
      auto c = sr_ctx;
      c.set(v);
      return c;
    };

    for (bool full : {true, false}) {
      for (std::size_t k = 0; k != 2; ++k) {
        // the extended theorem projects the survivors of general indices onto
        // the core and virtual parts (see "projected survivors"), which the
        // standard theorem does not: only full contractions compare
        if (!full && k == 1) continue;
        ExprPtr mp, sp;
        {
          auto r =
              set_scoped_default_context(with_vacuum(Vacuum::MultiProduct));
          mp = wick_mp(make_input(k), {.full_contractions = full});
        }
        {
          auto r =
              set_scoped_default_context(with_vacuum(Vacuum::SingleProduct));
          FWickTheorem wick{make_input(k)};
          sp = wick.full_contractions(full).compute();
        }
        mp = simplify(mp);
        sp = simplify(sp);
        // no active space: no γ, η or κ
        mp->visit(
            [&](const ExprPtr& e) {
              if (e->is<Tensor>()) {
                const auto l = e->as<Tensor>().label();
                REQUIRE(l != L"γ");
                REQUIRE(l != L"η");
                REQUIRE(l != L"κ");
              }
            },
            /*atoms_only=*/true);
        INFO("full=" << full << " k=" << k << "\nMP: " << toUtf8(to_latex(mp))
                     << "\nSP: " << toUtf8(to_latex(sp)));
        // The standard theorem spells a contraction as the overlap s, the
        // extended one as δ; they coincide for same-space indices.
        sp->visit(
            [](const ExprPtr& e) {
              if (e->is<Tensor>() &&
                  e->as<Tensor>().label() == reserved::overlap_label())
                e->as<Tensor>().set_label(reserved::kronecker_label());
            },
            /*atoms_only=*/true);
        // Surviving ã carry their context's vacuum tag, which keeps otherwise
        // equal terms apart in mp - sp: rebuild them with a common vacuum.
        for (auto* e : {&mp, &sp})
          (*e)->visit(
              [](const ExprPtr& x) {
                if (x->is<FNOperator>()) {
                  auto& nop = x->as<FNOperator>();
                  x->as<FNOperator>() =
                      FNOperator(cre(nop.creators()), ann(nop.annihilators()),
                                 Vacuum::SingleProduct);
                }
              },
              /*atoms_only=*/true);
        REQUIRE(simplify(mp - sp) == ex<Constant>(0));
      }
    }
  }

  SECTION("WickTheorem: GNO strings from elementary operators") {
    // {a†_p a_q} = {a†_p}{a_q} - ⟨{a†_p}{a_q}⟩, so a product of 1-body GNO
    // strings equals the product of these differences, each single-operator
    // string its own input operator
    auto single = [](const Index& p, const Index& q) {
      return ex<FNOperator>(cre({p}), ann({})) *
             ex<FNOperator>(cre({}), ann({q}));
    };
    // ⟨{a†_p}{a_q}⟩ with its dummies renamed apart from every other's
    auto contraction = [&](const Index& p, const Index& q) {
      auto e = wick_mp(single(p, q));
      container::map<Index, Index> fresh;
      e->visit(
          [&](const ExprPtr& x) {
            if (x->is<Tensor>())
              for (const auto& idx : x->as<Tensor>().const_braket())
                if (idx != p && idx != q && !fresh.contains(idx))
                  fresh.emplace(idx, Index::make_tmp_index(idx.space()));
          },
          /*atoms_only=*/true);
      e->visit(
          [&](const ExprPtr& x) {
            if (x->is<Tensor>()) {
              x->as<Tensor>().transform_indices(fresh);
              x->as<Tensor>().reset_tags();
            }
          },
          /*atoms_only=*/true);
      return e;
    };
    using Pairs = container::svector<std::pair<std::wstring, std::wstring>>;
    struct Case {
      Pairs pairs;
      bool full;
      std::size_t nterms;
    };
    for (const auto& [pairs, full, nterms] :
         {Case{Pairs{{L"u_1", L"u_2"}, {L"u_3", L"u_4"}, {L"u_5", L"u_6"}},
               true, 9},
          Case{Pairs{{L"p_1", L"p_2"}, {L"p_3", L"p_4"}}, true, 7},
          Case{Pairs{{L"u_1", L"u_2"}, {L"u_3", L"u_4"}}, false, 5}}) {
      ExprPtr gno = ex<Constant>(1), elementary = ex<Constant>(1);
      for (const auto& [p, q] : pairs) {
        gno = gno * ex<FNOperator>(cre({Index(p)}), ann({Index(q)}));
        elementary = elementary * (single(Index(p), Index(q)) -
                                   contraction(Index(p), Index(q)));
      }
      const detail::ExtendedWickOptions opts{.full_contractions = full};
      auto lhs = wick_mp(gno, opts);
      auto rhs = wick_mp(elementary, opts);
      INFO("lhs: " << toUtf8(to_latex(lhs))
                   << "\nrhs: " << toUtf8(to_latex(rhs)));
      REQUIRE(lhs->size() == nterms);
      REQUIRE(simplify(lhs - rhs) == ex<Constant>(0));
    }
    // the three all-active strings have both a κ2 and a κ3
    container::set<std::size_t> kappa_ranks;
    wick_mp(ex<FNOperator>(cre({L"u_1"}), ann({L"u_2"})) *
            ex<FNOperator>(cre({L"u_3"}), ann({L"u_4"})) *
            ex<FNOperator>(cre({L"u_5"}), ann({L"u_6"})))
        ->visit(
            [&](const ExprPtr& e) {
              if (e->is<Tensor>() && e->as<Tensor>().label() == L"κ")
                kappa_ranks.insert(e->as<Tensor>().bra_rank());
            },
            /*atoms_only=*/true);
    REQUIRE(kappa_ranks == container::set<std::size_t>{2, 3});
  }
}
