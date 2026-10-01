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
}
