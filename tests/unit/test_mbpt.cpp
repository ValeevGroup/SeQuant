//
// Created by Eduard Valeyev on 2019-02-19.
//

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/density.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/core/utility/expr.hpp>
#include <SeQuant/core/utility/timer.hpp>
#include <SeQuant/core/wick.hpp>
#include <SeQuant/domain/mbpt/context.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/op.hpp>
#include <SeQuant/domain/mbpt/op_registry.hpp>
#include <SeQuant/domain/mbpt/rdm.hpp>
#include <SeQuant/domain/mbpt/rules/df.hpp>
#include <SeQuant/domain/mbpt/rules/thc.hpp>
#include <SeQuant/domain/mbpt/utils.hpp>
#include <SeQuant/domain/mbpt/vac_av.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "catch2_sequant.hpp"

#include <algorithm>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

#include "SeQuant/core/utility/debug.hpp"

namespace {
/// Returns an RAII guard that, while alive, pins the (de)excitation amplitude
/// operators (t, λ, R, L) Hermitian in the default MBPT context. The
/// field/hermiticity model makes these operators bra<->ket *nonsymmetric* by
/// default (which is correct); pinning them Hermitian here lets the reference
/// equations below -- generated under the legacy conjugate-symmetric assumption
/// -- stay valid without regenerating them.
[[nodiscard]] auto scoped_hermitian_amplitudes() {
  auto reg = sequant::mbpt::get_default_mbpt_context().op_registry()->clone();
  for (std::wstring lbl : {L"t", L"λ", L"R", L"L"})
    if (reg.contains(lbl))
      reg.set_hermiticity(lbl, sequant::Hermiticity::Hermitian);
  return sequant::mbpt::set_scoped_default_mbpt_context(
      {.csv = sequant::mbpt::get_default_mbpt_context().csv(),
       .op_registry = std::move(reg)});
}

/// the first ordinal of the indices external_in_base() names: above those the
/// tests spell, below those of temporary indices
constexpr std::size_t external_in_base_ordinal = 50;

/// @return the index that stands for @p idx, one of @p externals, in the base
/// space @p base
sequant::Index external_in_base(
    const sequant::Index& idx, const sequant::IndexSpace& base,
    const sequant::container::set<sequant::Index>& externals) {
  const auto k = std::distance(externals.begin(), externals.find(idx));
  REQUIRE(external_in_base_ordinal + k < sequant::Index::min_tmp_index());
  return sequant::Index(
      base.base_key() + L"_" + std::to_wstring(external_in_base_ordinal + k),
      base);
}

/// @return @p expr spelled so that two reference expectation values that
/// differ only in how they write the same sums compare equal: every η is
/// δ - γ, every δ over a dummy is applied, and every index in a non-base
/// space (e.g. E, or O) is split into a sum over the base spaces it spans,
/// one of @p externals into indices named by external_in_base()
sequant::ExprPtr in_base_spaces(
    sequant::ExprPtr expr,
    const sequant::container::set<sequant::Index>& externals = {}) {
  using namespace sequant;
  const auto isr = get_default_context().index_space_registry();
  auto terms_of = [](const ExprPtr& e) {
    return e->is<Sum>() ? e->as<Sum>().summands() | ranges::to_vector
                        : std::vector<ExprPtr>{e};
  };
  auto first_index = [](const ExprPtr& term, auto&& pred) {
    std::optional<Index> found;
    auto look = [&](const ExprPtr& f) {
      if (found || !f->is<Tensor>()) return;
      for (const auto& idx : f->as<Tensor>().const_braket())
        if (pred(idx)) {
          found = idx;
          return;
        }
    };
    if (term->is<Tensor>())
      look(term);
    else
      term->visit(look, /*atoms_only=*/true);
    return found;
  };

  expr = expr->clone();
  expand(expr);
  for (bool split = true; split;) {
    split = false;
    auto result = std::make_shared<Sum>();
    for (const auto& term : terms_of(expr)) {
      const auto idx = first_index(
          term, [&](const Index& i) { return !isr->is_base(i.space()); });
      if (!idx) {
        result->append(term);
        continue;
      }
      split = true;
      for (const auto& base : isr->base_spaces())
        if (base.qns() == idx->space().qns() &&
            idx->space().type().includes(base.type()))
          result->append(transform_expr(
              term, {{*idx, externals.contains(*idx)
                                ? external_in_base(*idx, base, externals)
                                : Index::make_tmp_index(base)}}));
    }
    expr = result;
    expand(expr);
  }

  auto eta_to_delta_minus_gamma = [](ExprPtr& f) {
    if (f->is<Tensor>() &&
        f->as<Tensor>().label() == density::hole_rdm_label()) {
      const auto& t = f->as<Tensor>();
      f = make_kronecker(t.bra()[0], t.ket()[0]) -
          density::make_rdm(t.bra()[0], t.ket()[0]);
    }
  };
  if (expr->is_atom())
    eta_to_delta_minus_gamma(expr);
  else
    expr->visit(eta_to_delta_minus_gamma, /*atoms_only=*/true);
  expand(expr);
  auto result = std::make_shared<Sum>();
  for (auto term : terms_of(expr)) {
    if (term->is<Product>()) {
      FWickTheorem reducer{term};
      reducer.reduce(term);
    }
    result->append(canonicalize(term));
  }
  ExprPtr out = result;
  return simplify(out);
}
}  // namespace

TEST_CASE("mbpt_operator_type_id", "[mbpt]") {
  using namespace sequant;
  using mbpt::qns_t;
  const auto fid = Expr::get_type_id<mbpt::FOperator<qns_t>>();
  const auto bid = Expr::get_type_id<mbpt::BOperator<qns_t>>();
  REQUIRE(fid != bid);
  REQUIRE(Expr::get_type_id<mbpt::FOperatorBase>() != fid);
  REQUIRE(Expr::get_type_id<mbpt::BOperatorBase>() != bid);
  REQUIRE(Expr::type_rank_of(fid) == Expr::default_type_rank);
  REQUIRE(Expr::type_rank_of(bid) == Expr::default_type_rank);
  REQUIRE(mbpt::FOperator<qns_t>::static_type_name() ==
          "sequant::mbpt::Operator<sequant::mbpt::QuantumNumberChange<int64,"
          "sequant::mbpt::default_qns_tag>,FermiDirac>");
}

TEST_CASE("mbpt multireference regressions", "[mbpt]") {
  using namespace sequant;
  using namespace sequant::mbpt;
  auto ctx = get_default_context();
  ctx.set(make_mr_spaces());
  auto ctx_resetter = set_scoped_default_context(ctx);

  SECTION("commutators retain active-space contractions") {
    const auto f = op::h(1);
    const auto t = op::t(1);
    CHECK_FALSE(f->commutes_with(*t));
    const auto commutator = simplify(f * t - t * f);
    const auto expected = tensor::ref_av(tensor::h(1) * tensor::t(1) -
                                         tensor::t(1) * tensor::h(1));
    REQUIRE(expected != ex<Constant>(0));
    CHECK_THAT(op::ref_av(commutator), EquivalentTo(expected));
  }

  SECTION("single-reference excitation operators still commute") {
    auto sr_ctx = get_default_context();
    sr_ctx.set(make_sr_spaces());
    auto sr_scope = set_scoped_default_context(sr_ctx);
    CHECK(op::t(1)->commutes_with(*op::t(2)));
  }

  SECTION("RDM replacement preserves external indices") {
    for (const bool named : {false, true}) {
      auto opts = CanonicalizeOptions::default_options();
      if (named)
        opts.named_indices = container::set<Index>{
            Index(L"u_1"), Index(L"u_2"), Index(L"u_3"), Index(L"u_4")};
      auto scope = set_scoped_modified_default_context(
          [&opts](sequant::Context& ctx) { ctx.set(opts); });
      CAPTURE(named);
      CHECK_THAT(tensor::ref_av(deserialize(L"ã{u_1;u_2} * ã{u_3;u_4}")),
                 EquivalentTo(L"γ{u_1,u_3;u_2,u_4}:A-C-S + "
                              L"s{u_1;u_4} * γ{u_3;u_2}"));
      CHECK_THAT(tensor::ref_av(deserialize(L"x{u_1;u_2} * ã{u_2;u_1}")),
                 EquivalentTo(L"x{u_1;u_2} * γ{u_2;u_1}"));
    }
  }

  SECTION("external RDM indices beyond the active space are restricted") {
    for (const bool named : {false, true}) {
      auto opts = CanonicalizeOptions::default_options();
      if (named)
        opts.named_indices =
            container::set<Index>{Index(L"I_1"), Index(L"I_2"), Index(L"I_4"),
                                  Index(L"I_6"), Index(L"p_1"), Index(L"p_2")};
      auto scope = set_scoped_modified_default_context(
          [&opts](sequant::Context& ctx) { ctx.set(opts); });
      CAPTURE(named);
      CHECK_THAT(tensor::ref_av(deserialize(L"ã{I_1;I_2}")),
                 EquivalentTo(L"γ{u_1;u_2} * δ{I_1;u_1} * δ{u_2;I_2}"));
      CHECK_THAT(tensor::ref_av(deserialize(L"ã{p_1;p_2}")),
                 EquivalentTo(L"γ{u_1;u_2} * δ{p_1;u_1} * δ{u_2;p_2}"));
      CHECK_THAT(tensor::ref_av(deserialize(L"f{I_4;I_5} * ã{I_5;I_6}")),
                 EquivalentTo(L"f{I_4;u_1} * γ{u_1;u_2} * δ{u_2;I_6}"));
      // only the RDM slot of an external index is restricted
      if (named)
        CHECK_THAT(
            tensor::ref_av(deserialize(L"x{;;I_2} * ã{I_1;I_2}")),
            EquivalentTo(L"x{;;I_2} * γ{u_1;u_2} * δ{I_1;u_1} * δ{u_2;I_2}"));
    }
  }

  SECTION("a single normal operator becomes an RDM") {
    CHECK_THAT(tensor::ref_av(deserialize(L"ã{u_1;u_2}")),
               EquivalentTo(L"γ{u_1;u_2}"));
    CHECK_THAT(tensor::ref_av(deserialize(L"ã{u_1,u_3;u_2,u_4}")),
               EquivalentTo(L"γ{u_1,u_3;u_2,u_4}:A-C-S"));
  }

  SECTION("bare abstract operators retain their reference average") {
    const auto expected = tensor::ref_av(tensor::H(1));
    REQUIRE(expected != ex<Constant>(0));
    CHECK_THAT(op::ref_av(op::H(1)), EquivalentTo(expected));
  }
}

TEST_CASE("mbpt", "[mbpt][valgrind_skip]") {
  SECTION("cardinal tensor labels") {
    // every reference density label sorts as a cardinal label
    for (const auto& label : sequant::reserved::density_labels())
      REQUIRE(ranges::contains(sequant::mbpt::cardinal_tensor_labels(), label));
  }

  SECTION("registry") {
    using namespace sequant::mbpt;

    SECTION("empty-registry") {
      OpRegistry registry;
      REQUIRE(registry.ops().empty());
      REQUIRE_FALSE(registry.contains(L"T"));
    }

    SECTION("add-operators") {
      OpRegistry registry;
      registry.add(L"T", OpClass::Ex);
      registry.add(L"L", OpClass::Deex);
      registry.add(L"F", OpClass::Gen);

      REQUIRE(registry.contains(L"T"));
      REQUIRE(registry.contains(L"L"));
      REQUIRE(registry.contains(L"F"));
      REQUIRE(registry.ops().size() == 3);

      REQUIRE(registry.to_class(L"T") == OpClass::Ex);
      REQUIRE(registry.to_class(L"L") == OpClass::Deex);
      REQUIRE(registry.to_class(L"F") == OpClass::Gen);
    }

    SECTION("add-operators-with-hermiticity") {
      using sequant::Hermiticity;

      OpRegistry registry;
      // no explicit Hermiticity => default_hermiticity(class)
      registry.add(L"T", OpClass::Ex);
      registry.add(L"F", OpClass::Gen);
      REQUIRE(registry.hermiticity(L"T") == Hermiticity::NonHermitian);
      REQUIRE(registry.hermiticity(L"F") == Hermiticity::Hermitian);

      // explicit Hermiticity overrides the class default
      registry.add(L"Z", OpClass::Ex, Hermiticity::Hermitian);
      registry.add(L"G", OpClass::Gen, Hermiticity::NonHermitian);
      REQUIRE(registry.to_class(L"Z") == OpClass::Ex);
      REQUIRE(registry.to_class(L"G") == OpClass::Gen);
      REQUIRE(registry.hermiticity(L"Z") == Hermiticity::Hermitian);
      REQUIRE(registry.hermiticity(L"G") == Hermiticity::NonHermitian);

      // set_hermiticity overrides after the fact, in both directions
      registry.set_hermiticity(L"T", Hermiticity::Hermitian);
      registry.set_hermiticity(L"Z", Hermiticity::NonHermitian);
      REQUIRE(registry.hermiticity(L"T") == Hermiticity::Hermitian);
      REQUIRE(registry.hermiticity(L"Z") == Hermiticity::NonHermitian);

      // overrides survive cloning
      auto cloned = registry.clone();
      REQUIRE(cloned.hermiticity(L"Z") == Hermiticity::NonHermitian);
      REQUIRE(cloned.hermiticity(L"G") == Hermiticity::NonHermitian);
      REQUIRE(cloned.hermiticity(L"F") == Hermiticity::Hermitian);
    }

    SECTION("remove-operators") {
      OpRegistry registry;
      registry.add(L"T", OpClass::Ex);
      registry.add(L"L", OpClass::Deex);

      REQUIRE(registry.contains(L"T"));
      registry.remove(L"T");
      REQUIRE_FALSE(registry.contains(L"T"));
      REQUIRE(registry.contains(L"L"));
    }

    SECTION("clone-registry") {
      OpRegistry registry;
      registry.add(L"T", OpClass::Ex);
      registry.add(L"L", OpClass::Deex);

      auto cloned = registry.clone();
      REQUIRE(cloned.contains(L"T"));
      REQUIRE(cloned.contains(L"L"));
      REQUIRE(cloned.ops().size() == registry.ops().size());
    }

    SECTION("purge-registry") {
      OpRegistry registry;
      registry.add(L"T", OpClass::Ex);
      registry.add(L"L", OpClass::Deex);

      REQUIRE_FALSE(registry.ops().empty());
      registry.purge();
      REQUIRE(registry.ops().empty());
    }

    if (sequant::assert_behavior() == sequant::AssertBehavior::Throw) {
      SECTION("reserved-labels") {
        OpRegistry registry;
        // should not be able to add reserved labels
        REQUIRE_THROWS(
            registry.add(sequant::reserved::antisymm_label(), OpClass::Gen));
        REQUIRE_THROWS(
            registry.add(sequant::reserved::symm_label(), OpClass::Gen));
      }

      SECTION("duplicate-operators") {
        OpRegistry registry;
        registry.add(L"T", OpClass::Ex);
        // should not be able to add duplicate
        REQUIRE_THROWS(registry.add(L"T", OpClass::Deex));
      }
    }
  }  // SECTION("registry")

  SECTION("context") {
    using namespace sequant::mbpt;

    SECTION("default-context") {
      auto default_ctx = get_default_mbpt_context();
      // check default values
      REQUIRE(default_ctx.csv() == CSV::No);
      REQUIRE_NOTHROW(default_ctx.csv());
      REQUIRE_NOTHROW(default_ctx.op_registry());
    }

    SECTION("context-with-CSV") {
      auto ctx = Context({.csv = CSV::Yes});
      REQUIRE(ctx.csv() == CSV::Yes);

      ctx.set(CSV::No);
      REQUIRE(ctx.csv() == CSV::No);
    }

    SECTION("context-with-registry") {
      auto reg = make_minimal_registry();
      auto ctx = Context({.op_registry_ptr = reg});

      REQUIRE(ctx.csv() == CSV::No);
      REQUIRE(ctx.op_registry() != nullptr);
      REQUIRE(ctx.op_registry() == reg);
    }

    SECTION("set-operator-registry") {
      auto ctx = Context();
      auto reg = make_legacy_registry();

      ctx.set(reg);
      REQUIRE(ctx.op_registry() != nullptr);
      REQUIRE(ctx.op_registry() == reg);
    }

    SECTION("clone-context") {
      auto reg = std::make_shared<OpRegistry>();
      reg->add(L"T", OpClass::Ex);

      auto ctx = Context({.csv = CSV::Yes, .op_registry_ptr = reg});
      auto cloned = ctx.clone();

      REQUIRE(cloned.csv() == ctx.csv());
      REQUIRE(cloned.op_registry() != ctx.op_registry());  // Different pointer
      REQUIRE(cloned.op_registry()->contains(L"T"));
    }

    SECTION("context-equality") {
      auto reg1 = make_legacy_registry();
      auto reg2 = make_minimal_registry();
      auto ctx1 = Context({.csv = CSV::Yes, .op_registry_ptr = reg1});
      auto ctx2 = Context({.csv = CSV::Yes, .op_registry_ptr = reg1});
      auto ctx3 = Context({.csv = CSV::No, .op_registry_ptr = reg1});
      auto ctx4 = Context({.csv = CSV::Yes, .op_registry_ptr = reg2});

      REQUIRE(ctx1 == ctx2);
      REQUIRE(ctx1 != ctx3);
      REQUIRE(ctx1 != ctx4);
    }

    SECTION("scoped-context") {
      auto original_ctx = get_default_mbpt_context();
      auto original_csv = original_ctx.csv();

      {
        auto reg1 = make_legacy_registry();
        auto ctx = Context({.csv = CSV::Yes, .op_registry_ptr = reg1});
        auto _ = set_scoped_default_mbpt_context(ctx);

        auto current_ctx = get_default_mbpt_context();
        REQUIRE(current_ctx.csv() == CSV::Yes);
        REQUIRE(current_ctx.op_registry() == reg1);
        REQUIRE(current_ctx.op_registry()->contains(L"t"));
        REQUIRE(current_ctx.op_registry()->contains(L"λ"));
      }

      // after scope, original context should be restored
      auto restored_ctx = get_default_mbpt_context();
      REQUIRE(restored_ctx.csv() == original_csv);
    }
  }  // SECTION("context")

  SECTION("nbody_operators") {
    using namespace sequant;

    // reference equations below predate the hermiticity model; keep amplitudes
    // Hermitian so they remain valid (see scoped_hermitian_amplitudes)
    auto herm_amplitudes_guard = scoped_hermitian_amplitudes();

    SECTION("constructor") {
      // tests 1-space quantum number case
      {
        using namespace sequant::mbpt;

        op_t f1([]() -> std::wstring_view { return L"f"; },
                []() -> ExprPtr {
                  return ex<Tensor>(L"f", bra{L"p_1"}, ket{L"p_2"}) *
                         ex<FNOperator>(cre({L"p_1"}), ann({L"p_2"}));
                },
                [](qns_t& qns) { qns += general_type_qns(1); });

        REQUIRE(f1.label() == L"f");

        {
          // exact compare of intervals
          using namespace boost::numeric::interval_lib::compare::possible;
          REQUIRE(operator==(
              f1()[0], general_type_qns(1)[0]));  // produces single replacement
          REQUIRE(operator!=(
              f1()[0],
              general_type_qns(2)[0]));  // cannot produce double replacement
          /// TODO clearly this test does not make sense for context implicit
          /// size of qns. Need help to reimagine this test.
          // REQUIRE(operator==(f1(qns_t{5, 0}), qns_t{{5, 6}, {0, 1}})); //
        }
      }

      // tests 2-space quantum number case
      {
        using namespace sequant::mbpt;

        // this is fock operator in terms of general spaces
        op_t f_gg([]() -> std::wstring_view { return L"f"; },
                  []() -> ExprPtr {
                    return ex<Tensor>(L"f", bra{L"p_1"}, ket{L"p_2"}) *
                           ex<FNOperator>(cre({L"p_1"}), ann({L"p_2"}));
                  },
                  [](qns_t& qns) { qns += mbpt::general_type_qns(1); });
        // excitation part of the Fock operator
        op_t f_uo([]() -> std::wstring_view { return L"f"; },
                  []() -> ExprPtr {
                    return ex<Tensor>(L"f", bra{L"a_2"}, ket{L"i_2"}) *
                           ex<FNOperator>(cre({L"a_1"}), ann({L"i_2"}));
                  },
                  [](qns_t& qns) { qns += mbpt::excitation_type_qns(1); });

        REQUIRE(f_gg.label() == L"f");
        REQUIRE(f_uo.label() == L"f");

        {
          // comparison

          // exact
          REQUIRE((f_uo() == excitation_type_qns(
                                 1)));  // f_uo produces single excitations
          REQUIRE((f_gg() !=
                   excitation_type_qns(
                       1)));  // f_gg does not produce just single excitations
          /* REQUIRE(f_gg().in(excitation_type_qns(1)));  // f_gg can produce
           single excitations REQUIRE(f_gg().in(deexcitation_type_qns(1)));  //
           f_gg can also produce single de-excitations REQUIRE(f_gg().in( {1, 1,
           0, 0}));  // f_gg can produce replacements within occupieds
           REQUIRE(f_gg().in(
               {0, 0, 1, 1}));  // f_gg can produce replacements within virtuals
           REQUIRE(f_gg().in(
               {1, 1, 1, 1}));  // f_gg cannot produce this double replacements,
                                // but this returns true TODO introduce
           constraints
                                // on the total number of creators/annihilators,
                                // the interval logic does not constrain it
           REQUIRE(f_gg().in(
               {0, 0, 0, 0}));  // f_gg cannot produce a null replacement, but
           this
                                // returns true TODO introduce constraints on
           the
                                // total number of creators/annihilators, the
                                // interval logic does not constrain it
                                */
          // most of these seem like artifacts of fixed interval logic. we can
          // add them back if needed

          /*REQUIRE(
              f_uo().in(excitation_type_qns(1)));  // f_uo can produce single
          excitations REQUIRE(!f_uo().in( deexcitation_type_qns(1)));  // f_uo
          cannot produce single de-excitations REQUIRE(!f_uo().in( {1, 1, 0,
          0}));
          // f_uo can produce replacements withing occupieds REQUIRE(!f_uo().in(
              {0, 0, 1, 1}));  // f_uo can produce replacements withing virtuals
          REQUIRE(!f_uo().in(
              {1, 1, 1, 1}));  // f_uo cannot produce double replacements
          REQUIRE(
              !f_uo().in({0, 0, 0, 0}));  // f_uo cannot produce null
          replacements

          REQUIRE(f_gg({0, 1, 1, 0})
                      .in({0, 0, 0, 0}));  // f_gg can produce reference when
                                           // acting on singly-excited
          determinant REQUIRE(f_gg({0, 1, 1, 0}) .in({0, 1, 1, 0}));  // f_gg
          can produce singly-excited determinant
                                  // when acting on singly-excited determinant
          REQUIRE(
              !f_uo({0, 1, 1, 0})
                   .in({0, 0, 0, 0}));  // f_uo can't produce reference when
                                        // acting on singl-y-excited determinant
          REQUIRE(f_uo({0, 1, 1, 0})
                      .in({0, 2, 2,
                           0}));  // f_uo can produce doubly-excited determinant
                                  // when acting on singl-y-excited determinant

          //        REQUIRE(!f1(qns_t{2, 2}).in(0));  // can't produce reference
          //        when
          //                                          // acting on
          doubly-excited
           */
        }
        {
          // equal compare
          // using namespace
          // boost::numeric::interval_lib::compare::lexicographic;
          // REQUIRE(f1(qns_t{0, 0}) == qns_t{-1, 1}); // not same as below due
          // to interaction with Catch could do REQUIRE(operator==(f1(qns_t{0,
          // 0}), qns_t{-1, 1})); but equal is shorter
          //        REQUIRE(equal(f1(qns_t{0, 0}), qns_t{-1, 1}));
          //        REQUIRE(equal(f1(qns_t{-1, 1}), qns_t{-2, 2}));
        }
      }

      // tests OpParams constructor with perturbation order
      {
        using namespace sequant::mbpt;

        op_t op_order0([]() -> std::wstring_view { return L"t"; },
                       []() -> ExprPtr {
                         return ex<Tensor>(L"t", bra{L"a_1"}, ket{L"i_1"}) *
                                ex<FNOperator>(cre({L"a_1"}), ann({L"i_1"}));
                       },
                       [](qns_t& qns) { qns += excitation_type_qns(1); },
                       OpParams{.order = 0});
        REQUIRE(op_order0.label() == L"t");
        REQUIRE(op_order0.order() == 0);

        op_t op_order1([]() -> std::wstring_view { return L"t"; },
                       []() -> ExprPtr {
                         return ex<Tensor>(L"t", bra{L"a_1"}, ket{L"i_1"}) *
                                ex<FNOperator>(cre({L"a_1"}), ann({L"i_1"}));
                       },
                       [](qns_t& qns) { qns += excitation_type_qns(1); },
                       OpParams{.order = 1});
        REQUIRE(op_order1.label() == L"t¹");
        REQUIRE(op_order1.order() == 1);

        op_t op_order2([]() -> std::wstring_view { return L"t"; },
                       []() -> ExprPtr {
                         return ex<Tensor>(L"t", bra{L"a_1", L"a_2"},
                                           ket{L"i_1", L"i_2"}) *
                                ex<FNOperator>(cre({L"a_1", L"a_2"}),
                                               ann({L"i_1", L"i_2"}));
                       },
                       [](qns_t& qns) { qns += excitation_type_qns(2); },
                       OpParams{.order = 2});
        REQUIRE(op_order2.label() == L"t²");
        REQUIRE(op_order2.order() == 2);
      }
    }  // SECTION("constructor")

    SECTION("to_latex") {
      using qns_t [[maybe_unused]] = mbpt::qns_t;
      using namespace sequant::mbpt;

      auto f = F();
      auto t1 = T(1);
      auto t2 = t(2);
      auto lambda1 = λ(1);
      auto lambda2 = λ(2);
      auto r_2_1 = r(nₚ(1), nₕ(2));
      auto r_1_2 = r(nₚ(2), nₕ(1));
      auto theta2 = θ(2);
      REQUIRE(to_latex(theta2) == L"{\\hat{\\theta}_{2}}");
      REQUIRE(to_latex(f) == L"{\\hat{f}}");
      REQUIRE(to_latex(t1) == L"{\\hat{t}_{1}}");
      REQUIRE(to_latex(t2) == L"{\\hat{t}_{2}}");
      REQUIRE(to_latex(lambda1) == L"{\\hat{\\lambda}_{1}}");
      REQUIRE(to_latex(lambda2) == L"{\\hat{\\lambda}_{2}}");
      REQUIRE(to_latex(r_2_1) == L"{\\hat{R}_{2,1}}");
      REQUIRE(to_latex(r_1_2) == L"{\\hat{R}_{1,2}}");

      // projectors
      auto P2 = P(2);
      auto P_neg2 = P(-2);
      auto P_2_1 = P(nₚ(1), nₕ(2));
      auto P_neg_2_1 = P(nₚ(-1), nₕ(-2));
      REQUIRE(to_latex(P2) == L"{\\hat{A}_{-2}}");
      REQUIRE(to_latex(P_neg2) == L"{\\hat{A}_{2}}");
      REQUIRE(to_latex(P_2_1) == L"{\\hat{A}_{-2,-1}}");
      REQUIRE(to_latex(P_neg_2_1) == L"{\\hat{A}_{2,1}}");
    }  // SECTION("to_latex")

    SECTION("canonicalize") {
      using qns_t [[maybe_unused]] = mbpt::qns_t;
      using namespace sequant::mbpt;
      auto f = F();
      auto t1 = t(1);
      auto l1 = λ(1);
      auto t2 = t(2);
      auto l2 = λ(2);
      auto h_pt = Hʼ(1, {.order = 1});
      REQUIRE(to_latex(f * t1 * t2) == to_latex(canonicalize(f * t2 * t1)));
      REQUIRE(to_latex(canonicalize(f * t1 * t2)) ==
              to_latex(canonicalize(f * t2 * t1)));
      REQUIRE(to_latex(t1 * t2 * f * t1 * t2) ==
              to_latex(canonicalize(t2 * t1 * f * t2 * t1)));

      REQUIRE(to_latex(ex<Constant>(3) * f * t1 * t2) ==
              to_latex(simplify(ex<Constant>(2) * f * t2 * t1 + f * t1 * t2)));

      REQUIRE(to_latex(simplify(t1 * l1)) != to_latex(simplify(l1 * t1)));
      REQUIRE(to_latex(simplify(t1 * l2)) != to_latex(simplify(l2 * t1)));

      REQUIRE(to_latex(simplify(l2 * t1)) ==
              L"{{\\hat{\\lambda}_{2}}{\\hat{t}_{1}}}");
      REQUIRE(to_latex(simplify(t1 * l2)) ==
              L"{{\\hat{t}_{1}}{\\hat{\\lambda}_{2}}}");

      REQUIRE(to_latex(simplify(t1 + t1)) ==
              to_latex(simplify(ex<Constant>(2) * t1)));

      REQUIRE(to_latex(simplify(t1 + t1 + t2)) ==
              to_latex(simplify(ex<Constant>(2) * t1 + t2)));

      REQUIRE(to_latex(simplify(f + f)) ==
              to_latex(simplify(ex<Constant>(2) * f)));

      REQUIRE(to_latex(simplify(h_pt + h_pt + f)) ==
              to_latex(simplify(ex<Constant>(2) * h_pt + f)));

      auto t = t1 + t2;

      {
        //      std::wcout << "to_latex(simplify(f * t * t)): "
        //                 << to_latex(simplify(f * t * t)) << std::endl;
        REQUIRE_THAT(simplify(f * t * t),
                     EquivalentTo(f * t2 * t2 + ex<Constant>(2) * f * t1 * t2 +
                                  f * t1 * t1));
      }

      {
        //      std::wcout << "to_latex(simplify(f * t * t * t): "
        //                 << to_latex(simplify(f * t * t * t)) << std::endl;
        REQUIRE_THAT(simplify(f * t * t * t),
                     EquivalentTo(f * t1 * t1 * t1 + f * t2 * t2 * t2 +
                                  ex<Constant>(3) * f * t1 * t2 * t2 +
                                  ex<Constant>(3) * f * t1 * t1 * t2));
      }
    }  // SECTION("canonicalize")

    SECTION("adjoint") {
      using qns_t = mbpt::qns_t;
      using op_t = mbpt::Operator<qns_t>;
      using namespace mbpt;
      op_t f = F()->as<op_t>();
      op_t t1 = t(1)->as<op_t>();
      op_t lambda2 = λ(2)->as<op_t>();
      op_t r_1_2 = r(nₚ(2), nₕ(1))->as<op_t>();

      REQUIRE_NOTHROW(adjoint(f));
      REQUIRE_NOTHROW(adjoint(t1));
      REQUIRE_NOTHROW(adjoint(lambda2));
      REQUIRE_NOTHROW(adjoint(r_1_2));

      REQUIRE(adjoint(f)() == mbpt::general_type_qns(1));
      REQUIRE(adjoint(t1)() == mbpt::deexcitation_type_qns(1));
      REQUIRE(adjoint(lambda2)() == mbpt::excitation_type_qns(2));
      REQUIRE(adjoint(r_1_2)() == l(nₚ(2), nₕ(1))->as<op_t>()());

      // adjoint(adjoint(Op)) = Op
      REQUIRE(adjoint(adjoint(t1))() == t1());
      REQUIRE(adjoint(adjoint(r_1_2))() == r_1_2());

      // tensor_form()
      REQUIRE((simplify(adjoint(t1).tensor_form())) ==
              (simplify(adjoint(t1.tensor_form()))));

      REQUIRE(simplify(adjoint(lambda2).tensor_form()) ==
              simplify(adjoint(lambda2.tensor_form())));
      REQUIRE(simplify(adjoint(r_1_2).tensor_form()) ==
              simplify(adjoint(r_1_2.tensor_form())));

      // to_latex()
      REQUIRE(to_latex(adjoint(f).as<Expr>()) == L"{\\hat{f⁺}}");
      REQUIRE(to_latex(adjoint(t1).as<Expr>()) == L"{\\hat{t⁺}^{1}}");
      REQUIRE(to_latex(adjoint(lambda2).as<Expr>()) ==
              L"{\\hat{\\lambda⁺}^{2}}");
      REQUIRE(to_latex(adjoint(r_1_2).as<Expr>()) == L"{\\hat{R⁺}^{1,2}}");

      // adjoint(adjoint(op)) == op
      auto t1_adj = adjoint(t1);
      auto r_1_2_adj = adjoint(r_1_2);
      auto lambda2_adj = adjoint(lambda2);
      REQUIRE(to_latex(adjoint(t1_adj).as<Expr>()) == L"{\\hat{t}_{1}}");
      REQUIRE(to_latex(adjoint(r_1_2_adj).as<Expr>()) == L"{\\hat{R}_{1,2}}");
      REQUIRE(to_latex(adjoint(lambda2_adj).as<Expr>()) ==
              L"{\\hat{\\lambda}_{2}}");

      // adjoint should preserve perturbation order
      auto t1_order2 = tʼ(1, {.order = 2});
      REQUIRE(t1_order2.as<op_t>().order() == 2);
      auto t1_order2_adj = adjoint(t1_order2);
      REQUIRE(t1_order2_adj.as<op_t>().order() == 2);

      auto λ1_order3 = λʼ(1, {.order = 3});
      REQUIRE(λ1_order3.as<op_t>().order() == 3);
      auto λ1_order3_adj = adjoint(λ1_order3);
      REQUIRE(λ1_order3_adj.as<op_t>().order() == 3);

      // each call to tensor_form() must generate fresh dummy indices, otherwise
      // using an operator more than once (e.g. squaring it in a similarity
      // transform) yields tensors that share indices. Non-adjoint
      // operators satisfy this:
      op_t t2op = t(2)->as<op_t>();
      REQUIRE(to_latex(t2op.tensor_form()) != to_latex(t2op.tensor_form()));
      // adjoint operators must satisfy it too:
      op_t λ2adj = adjoint(λ(2))->as<op_t>();
      REQUIRE(to_latex(λ2adj.tensor_form()) != to_latex(λ2adj.tensor_form()));

    }  // SECTION("adjoint")

    SECTION("screen") {
      using namespace sequant::mbpt;
      auto g_t2_t2 = h(2) * t(2) * t(2);
      REQUIRE(raises_vacuum_to_rank(g_t2_t2, 2));
      REQUIRE(raises_vacuum_up_to_rank(g_t2_t2, 2));

      auto g_t2 = h(2) * t(2);
      REQUIRE(raises_vacuum_to_rank(g_t2, 3));

      auto lambda2_f = λ(2) * h(1);
      REQUIRE(lowers_rank_to_vacuum(lambda2_f, 2));

      auto expr1 = P(nₚ(0), nₕ(1)) * H() * R(nₚ(0), nₕ(1));
      auto expr1_tnsr = lower_to_tensor_form(expr1);
      auto vev1_op = op::vac_av(expr1);
      auto vev1_t = tensor::vac_av(expr1_tnsr);  // no operator level screening
      REQUIRE(to_latex(vev1_op) == to_latex(vev1_t));

      auto expr2 = P(nₚ(2), nₕ(1)) * H() * R(nₚ(1), nₕ(0));
      auto expr2_tnsr = lower_to_tensor_form(expr2);
      auto vev2_op = op::vac_av(expr2);
      auto vev2_t = tensor::vac_av(expr2_tnsr);  // no operator level screening
      REQUIRE(to_latex(vev2_op) == to_latex(vev2_t));

      // Test screen_vac_av
      // CCD Hbar in connected-product form
      auto hbar = mbpt::lst(H(), t(2), 4, {.use_connected_form = true});
      auto screened_hbar = screen_vac_av(hbar);
      auto expected = h(2) * t(2);
      REQUIRE(simplify(screened_hbar - expected) == ex<Constant>(0));

      auto expr3 = P(2) * hbar * r(nₚ(2), nₕ(2));
      auto screened_expr3 = screen_vac_av(expr3);
      auto expected3 = op::P(2) * (h(2) * t(2) + h(2) + h(1)) * r(nₚ(2), nₕ(2));
      REQUIRE(simplify(screened_expr3 - expected3) == ex<Constant>(0));

      auto expr4 = P(nₚ(2), nₕ(1)) * hbar * R(nₚ(1), nₕ(0));
      auto screened_expr4 = screen_vac_av(expr4);
      auto expected4 = P(nₚ(2), nₕ(1)) *
                       (h(1) + h(1) * t(2) + h(2) * t(2) + h(2)) *
                       R(nₚ(1), nₕ(0));
      REQUIRE(simplify(screened_expr4 - expected4) == ex<Constant>(0));

    }  // SECTION("screen")

    SECTION("lst") {
      using namespace sequant::mbpt;

      auto commutator = [](const ExprPtr& A, const ExprPtr& B) {
        return A * B - B * A;
      };

      // non-unitary, rank 3
      auto expr1 =
          lst(H(), t(2), 3, {.unitary = false, .use_connected_form = true});
      auto expected1 =
          H() *
          (ex<Constant>(1) + t(2) + ex<Constant>(rational{1, 2}) * t(2) * t(2) +
           ex<Constant>(rational{1, 6}) * t(2) * t(2) * t(2));
      REQUIRE(simplify(expr1 - expected1) == ex<Constant>(0));

      auto expr2 =
          lst(H(), t(2), 3, {.unitary = false, .use_connected_form = false});
      auto expected2 =
          H() + commutator(H(), t(2)) +
          ex<Constant>(rational{1, 2}) *
              commutator(commutator(H(), t(2)), t(2)) +
          ex<Constant>(rational{1, 6}) *
              commutator(commutator(commutator(H(), t(2)), t(2)), t(2));
      REQUIRE(simplify(expr2 - expected2) == ex<Constant>(0));

      // the default LSTOptions must produce the explicit commutator form, i.e.
      // lst() must never assume the caller supplies operator connectivity
      REQUIRE(simplify(lst(H(), t(2), 3) - expected2) == ex<Constant>(0));
      REQUIRE(simplify(lst(H(), t(2), 3, {}) - expected2) == ex<Constant>(0));

      // unitary, rank 2
      using sequant::adjoint;
      auto expr3 =
          lst(H(), t(2), 2, {.unitary = true, .use_connected_form = true});
      auto expected3 =
          H() + H() * t(2) + adjoint(t(2)) * H() + adjoint(t(2)) * H() * t(2) +
          H() * ex<Constant>(rational{1, 2}) * t(2) * t(2) +
          ex<Constant>(rational{1, 2}) * adjoint(t(2)) * adjoint(t(2)) * H();
      REQUIRE(simplify(expr3 - expected3) == ex<Constant>(0));

      auto expr4 =
          lst(H(), t(2), 2, {.unitary = true, .use_connected_form = false});
      auto generator = commutator(H(), t(2)) - commutator(H(), adjoint(t(2)));
      auto expected4 =
          H() + generator +
          ex<Constant>(rational{1, 2}) * (commutator(generator, t(2)) -
                                          commutator(generator, adjoint(t(2))));
      REQUIRE(simplify(expr4 - expected4) == ex<Constant>(0));

      // `unitary` must no longer select the commutator representation: it only
      // swaps B for B - B^+, leaving use_connected_form at its default
      REQUIRE(simplify(lst(H(), t(2), 2, {.unitary = true}) - expected4) ==
              ex<Constant>(0));
    }  // SECTION("lst")

    SECTION("predefined") {
      // P.S. ref outputs produced with complete canonicalization
      auto ctx = get_default_context();
      ctx.set(CanonicalizeOptions{.method = CanonicalizationMethod::Complete});
      auto _ = set_scoped_default_context(ctx);
      using namespace sequant::mbpt;

      auto theta1 = θ(1)->as<op_t>();
      // std::wcout << "theta1: " << to_latex(simplify(theta1.tensor_form()));
      REQUIRE(to_latex(simplify(theta1.tensor_form())) ==
              L"{{\\theta^{{p_2}}_{{p_1}}}{\\tilde{a}^{{p_1}}_{{p_2}}}}");

      {
        // replacement operator: label "ã", lowering yields a bare
        // FNOperator with `rank` cre and `rank` ann over the complete space
        auto ã_2 = ã(2)->as<op_t>();
        REQUIRE(ã_2.label() == L"ã");

        auto ã_2_tform = ã_2.tensor_form();
        REQUIRE(ã_2_tform->is<FNOperator>());
        const auto& ã_2_fnop = ã_2_tform->as<FNOperator>();
        REQUIRE(ã_2_fnop.ncreators() == 2);
        REQUIRE(ã_2_fnop.nannihilators() == 2);

        const auto complete_space = get_complete_space(Spin::any);
        for (const auto& o : ã_2_fnop.creators())
          REQUIRE(o.index().space() == complete_space);
        for (const auto& o : ã_2_fnop.annihilators())
          REQUIRE(o.index().space() == complete_space);

        // callers depend on this exact free-index mapping
        REQUIRE(to_latex(ã_2_tform) ==
                L"{\\tilde{a}^{{p_1}{p_2}}_{{p_3}{p_4}}}");

        // regression: to_latex on the operator form throws if the label is
        // missing from the registry
        REQUIRE(to_latex(ã(1)) == L"{\\hat{\\tilde{a}}}");

        // rank 0 is an assert, not a throw
        if (sequant::assert_behavior() == sequant::AssertBehavior::Throw)
          REQUIRE_THROWS_AS(ã(0), Exception);
      }

      auto R_2 = r(2)->as<op_t>();
      //    std::wcout << "R_2: " << to_latex(simplify(R_2.tensor_form())) <<
      //    std::endl;
      REQUIRE(
          to_latex(simplify(R_2.tensor_form())) ==
          L"{{{\\frac{1}{4}}}{\\bar{R}^{{i_1}{i_2}}_{{a_1}{a_2}}}{\\tilde{a}^"
          L"{{a_1}{"
          L"a_2}}_{{i_1}{i_2}}}}");

      auto L_3 = l(3)->as<op_t>();
      //    std::wcout << "L_3: " << to_latex(simplify(L_3.tensor_form())) <<
      //    std::endl;
      REQUIRE(
          to_latex(simplify(L_3.tensor_form())) ==
          L"{{{\\frac{1}{36}}}{\\bar{L}^{{a_1}{a_2}{a_3}}_{{i_1}{i_2}{i_3}}}{"
          L"\\tilde{a}^{{i_1}{i_2}{i_3}}_{{a_1}{a_2}{a_3}}}}");

      auto R_2_3 = r(nₚ(3), nₕ(2))->as<op_t>();
      //    std::wcout << "R_2_3: " << to_latex(simplify(R_2_3.tensor_form()))
      //    << std::endl;
      REQUIRE(to_latex(simplify(R_2_3.tensor_form())) ==
              L"{{{\\frac{1}{12}}}{\\bar{R}^{{i_1}{i_2}}_{{a_1}{a_2}{a_3}}}{"
              L"\\tilde{a}^{"
              L"{a_1}{a_2}{a_3}}_{\\textvisiblespace\\,{i_1}{i_2}}}}");

      auto L_1_2 = l(nₚ(1), nₕ(2))->as<op_t>();
      // std::wcout << "l(1,2): " << to_latex(simplify(L_1_2.tensor_form())) <<
      // std::endl;
      REQUIRE(
          to_latex(simplify(L_1_2.tensor_form())) ==
          L"{{{\\frac{1}{2}}}{\\bar{L}^{{a_1}}_{{i_1}{i_2}}}{\\tilde{a}^{{i_"
          L"1}{i_2}}"
          L"_{\\textvisiblespace\\,{a_1}}}}");

      auto A_2_1 = A(nₚ(2), nₕ(1))->as<op_t>();
      //    std::wcout << "A_2_1: " << to_latex(simplify(A_2_1.tensor_form()))
      //               << std::endl;
      REQUIRE(to_latex(simplify(A_2_1.tensor_form())) ==
              L"{{\\hat{A}^{{i_1}}_{{a_1}{a_2}}}{\\tilde{a}^{{a_"
              L"1}{a_2}}"
              L"_{\\textvisiblespace\\,{i_1}}}}");

      auto P_0_1 = P(nₚ(0), nₕ(1))->as<op_t>();
      //    std::wcout << "P_0_1: " << to_latex(simplify(P_0_1.tensor_form()))
      //               << std::endl;
      REQUIRE(to_latex(simplify(P_0_1.tensor_form())) ==
              L"{{\\hat{A}^{}_{{i_1}}}{\\tilde{a}^{{i_1}}}}");

      auto P_2_1 = P(nₚ(2), nₕ(1))->as<op_t>();
      //    std::wcout << "P_2_1: " << to_latex(simplify(P_2_1.tensor_form()))
      //               << std::endl;
      REQUIRE(to_latex(simplify(P_2_1.tensor_form())) ==
              L"{{\\hat{A}^{{a_1}{a_2}}_{{i_1}}}{\\tilde{a}^{"
              L"\\textvisiblespace\\,{i_1}}_{{a_1}{a_2}}}}");

      auto P_2_3 = P(nₚ(2), nₕ(3))->as<op_t>();
      //    std::wcout << "P_2_3: " << to_latex(simplify(P_3_2.tensor_form()))
      //               << std::endl;
      REQUIRE(to_latex(simplify(P_2_3.tensor_form())) ==
              L"{{\\hat{A}^{{a_1}{a_2}}_{{i_1}{i_2}{i_3}}}{"
              L"\\tilde{a}^{"
              L"{i_1}{i_2}{i_3}}_{\\textvisiblespace\\,{a_1}{a_2}}}}");

      auto R33 = R(3);
      lower_to_tensor_form(R33);
      simplify(R33);
      //    std::wcout << "R33: " << to_latex(R33) << std::endl;
      REQUIRE(to_latex(R33) ==
              L"{ "
              L"\\bigl({{R^{{i_1}}_{{a_1}}}{\\tilde{a}^{{a_1}}_{{i_1}}}} + "
              L"{{{\\frac{1}{36}}}{\\bar{R}^{{i_1}{i_2}{i_3}}_{{a_1}{a_2}{a_3}}"
              L"}{\\tilde{a}^{{a_1}{a_2}{a_3}}_{{i_1}{i_2}{i_3}}}} + "
              L"{{{\\frac{1}{4}}}{\\bar{R}^{{i_1}{i_2}}_{{a_1}{a_2}}}{\\tilde{"
              L"a}^{{a_1}{a_2}}_{{i_1}{i_2}}}}\\bigr) }");

      auto R12 = R(nₚ(2), nₕ(1));
      lower_to_tensor_form(R12);
      simplify(R12);
      //    std::wcout << "R12: " << to_latex(R12) << std::endl;
      REQUIRE(to_latex(R12) ==
              L"{ \\bigl({{R^{}_{{a_1}}}{\\tilde{a}^{{a_1}}}} + "
              L"{{{\\frac{1}{2}}}{\\bar{R}^{{i_1}}_{{a_1}{a_2}}}{\\tilde{a}^{{"
              L"a_1}{a_"
              L"2}}_{\\textvisiblespace\\,{i_1}}}}\\bigr) }");

      auto R21 = R(nₚ(1), nₕ(2));
      lower_to_tensor_form(R21);
      simplify(R21);
      //    std::wcout << "R21: " << to_latex(R21) << std::endl;
      REQUIRE(to_latex(R21) ==
              L"{ "
              L"\\bigl({{{\\frac{1}{2}}}{\\bar{R}^{{i_1}{i_2}}_{{a_1}}}{"
              L"\\tilde{a}^{"
              L"\\textvisiblespace\\,{a_1}}_{{i_1}{i_2}}}} + "
              L"{{R^{{i_1}}_{}}{\\tilde{a}_{{i_1}}}}\\bigr) }");

      auto L23 = L(nₚ(2), nₕ(3));
      lower_to_tensor_form(L23);
      simplify(L23);
      // std::wcout << "L23: " << to_latex(L23) << std::endl;
      REQUIRE(to_latex(L23) ==
              L"{ "
              L"\\bigl({{L^{}_{{i_1}}}{\\tilde{a}^{{i_1}}}} + "
              L"{{{\\frac{1}{2}}}{\\bar{L}^{{a_1}}_{{i_1}{i_2}}}{\\tilde{a}^{{"
              L"i_1}{i_2}}_{\\textvisiblespace\\,{a_1}}}} + "
              L"{{{\\frac{1}{12}}}{\\bar{L}^{{a_1}{a_2}}_{{i_1}{i_2}{i_"
              L"3}}}{\\tilde{a}^{{i_1}{i_2}{i_3}}_{\\textvisiblespace\\,{a_1}{"
              L"a_2}}}}\\bigr) }");

      // perturbation ops
      REQUIRE_NOTHROW(Hʼ(1, {.order = 1}));
      REQUIRE_NOTHROW(Hʼ(2, {.order = 2}));
      REQUIRE_NOTHROW(Λʼ(3, {.order = 5, .skip1 = true}));
      REQUIRE_NOTHROW(Tʼ(2, {.order = 9}));
      if (sequant::assert_behavior() == sequant::AssertBehavior::Throw) {
        REQUIRE_THROWS(Hʼ(1, {.order = 10}));  // invalid order
      }

      auto h0 = Hʼ(1);
      REQUIRE(to_latex(h0) == L"{\\hat{h¹}}");
      auto h1 = Hʼ(1, {.order = 1});
      auto h2 = Hʼ(1, {.order = 2});
      auto t2 = tʼ(2, {.order = 2});

      REQUIRE(h1 != h2);
      REQUIRE(h0 == h1);

      REQUIRE(to_latex(simplify(h0 + h1)) == L"{{{2}}{\\hat{h¹}}}");
      REQUIRE(to_latex(simplify(h1 + h2)) ==
              L"{ \\bigl({\\hat{h¹}} + {\\hat{h²}}\\bigr) }");

      REQUIRE(to_latex(simplify(h1 * t2)) == L"{{\\hat{h¹}}{\\hat{t²}_{2}}}");
      REQUIRE(to_latex(simplify(h2 * t2)) == L"{{\\hat{h²}}{\\hat{t²}_{2}}}");

      // δl δr ops
      {  // Spinor basis
        auto dl2 = tensor::δl(2);
        REQUIRE(simplify(dl2 - rational{1, 2} * tensor::P(2)) ==
                ex<Constant>(0));
        auto dr2 = tensor::δr(2);
        REQUIRE(simplify(dr2 - rational{1, 2} * tensor::P(-2)) ==
                ex<Constant>(0));

        auto dl23 = tensor::δl(nₚ(2), nₕ(3));
        REQUIRE(simplify(dl23 - ex<Power>(rational{1, 12}, rational{1, 2}) *
                                    tensor::P(nₚ(2), nₕ(3))) ==
                ex<Constant>(0));
      }

      {  // spinfree basis
        auto ctx = get_default_context();
        auto ctx_resetter =
            set_scoped_default_context(ctx.set(SPBasis::Spinfree));
        auto dl2 = tensor::δl(2);
        REQUIRE(simplify(dl2 - ex<Power>(rational{1, 2}, rational{1, 2}) *
                                   tensor::P(2)) == ex<Constant>(0));
        auto dr3 = tensor::δr(3);
        REQUIRE(simplify(dr3 - ex<Power>(rational{1, 6}, rational{1, 2}) *
                                   tensor::P(-3)) == ex<Constant>(0));
      }

    }  // SECTION("predefined")

    SECTION("batching") {
      // update context to use batching index
      auto isr = sequant::mbpt::make_legacy_spaces();
      mbpt::add_batching_spaces(isr);
      auto ctx_resetter =
          set_scoped_default_context({.index_space_registry_shared_ptr = isr,
                                      .vacuum = Vacuum::SingleProduct});
      REQUIRE_NOTHROW(
          get_default_context().index_space_registry()->retrieve(L"z"));

      using namespace mbpt;
      REQUIRE_NOTHROW(op::Hʼ(1, {.order = 1, .nbatch = 1}));
      REQUIRE_NOTHROW(op::Hʼ(2, {.order = 1, .nbatch = 2}));
      REQUIRE_NOTHROW(
          op::Λʼ(3, {.order = 1, .batch_ordinals = {3, 4, 5}, .skip1 = true}));
      REQUIRE_NOTHROW(op::Tʼ(1, {.order = 1, .nbatch = 20}));

      // invalid usages
      if (sequant::assert_behavior() == sequant::AssertBehavior::Throw) {
        // cannot set both nbatch and batch_ordinals
        REQUIRE_THROWS_AS(
            op::Hʼ(2, {.order = 1, .nbatch = 2, .batch_ordinals = {1, 2}}),
            sequant::Exception);
        // all ordinals must be unique
        REQUIRE_THROWS_AS(op::Hʼ(2, {.order = 1, .batch_ordinals = {1, 2, 2}}),
                          sequant::Exception);
        // ordinals must be sorted
        REQUIRE_THROWS_AS(op::Hʼ(1, {.order = 1, .batch_ordinals = {3, 2}}),
                          sequant::Exception);
      }

      // operations
      auto h0 = op::Hʼ(1);
      REQUIRE(to_latex(h0) == L"{\\hat{h¹}}");

      auto h1 = op::Hʼ(1, {.order = 1, .nbatch = 1});
      auto h1_2 = op::Hʼ(1, {.order = 1, .batch_ordinals = {1, 2}});
      auto pt1 = op::Tʼ(2, {.order = 1, .batch_ordinals = {1}});

      auto sum0 = h0 + h1;
      simplify(sum0);
      REQUIRE(to_latex(sum0) ==
              L"{ \\bigl({\\hat{h¹}}{[{z}_{1}]} + {\\hat{h¹}}\\bigr) }");

      auto sum1 = h1 + h1;
      simplify(sum1);
      REQUIRE(to_latex(sum1) == L"{{{2}}{\\hat{h¹}}{[{z}_{1}]}}");
      auto sum2 = h1 + pt1;
      simplify(sum2);
      // std::wcout << "sum2:  " << to_latex(sum2) << std::endl;
      REQUIRE(to_latex(sum2) ==
              L"{ \\bigl({\\hat{h¹}}{[{z}_{1}]} + {\\hat{t¹}_{2}}{[{z}_{1}]} + "
              L"{\\hat{t¹}_{1}}{[{z}_{1}]}\\bigr) }");

      auto sum3 = h1 + h1_2;
      simplify(sum3);
      // std::wcout << "sum3:  " << to_latex(sum3) << std::endl;
      REQUIRE(to_latex(sum3) ==
              L"{ \\bigl({\\hat{h¹}}{[{z}_{1},{z}_{2}]} + "
              L"{\\hat{h¹}}{[{z}_{1}]}\\bigr) }");

      auto pdt1 = h1 * h1_2;
      simplify(pdt1);
      // std::wcout << "pdt1: " << to_latex(pdt1) << std::endl;
      REQUIRE(to_latex(pdt1) ==
              L"{{\\hat{h¹}}{[{z}_{1}]}{\\hat{h¹}}{[{z}_{1},{z}_{2}]}}");

      auto pdt2 = h1 * pt1;
      simplify(pdt2);
      // std::wcout << "pdt1: " << to_latex(pdt2) << std::endl;
      REQUIRE(to_latex(pdt2) ==
              L"{ \\bigl({{\\hat{h¹}}{[{z}_{1}]}{\\hat{t¹}_{1}}{[{z}_{1}]}} + "
              L"{{\\hat{h¹}}{[{z}_{1}]}{\\hat{t¹}_{2}}{[{z}_{1}]}}\\bigr) }");

      // lowering to tensor form
      auto sum1_t = simplify(lower_to_tensor_form(sum1));
      // std::wcout << "sum1_t: " << to_latex(simplify(sum1_t)) << std::endl;
      REQUIRE_THAT(sum1_t, EquivalentTo(L"2 * h¹{κ1;κ2;z1}:A-C-S * ã{κ2;κ1}"));

      auto sum2_t = simplify(lower_to_tensor_form(sum2));
      // std::wcout << "sum2_t: " << to_latex(sum2_t) << std::endl;
      REQUIRE_THAT(
          sum2_t,
          EquivalentTo(L"h¹{κ2;κ1;z1}:A-C-S * ã{κ1;κ2} + t¹{a1;i1;z1}:A-C-S * "
                       L"ã{i1;a1} + "
                       "(1/4) * t¹{a1,a2;i1,i2;z1}:A-C-S * ã{i1,i2;a1,a2}"));

      auto expr3_t = simplify(lower_to_tensor_form(sum2 * h1_2));
      // std::wcout << "expr3_t: " << to_latex(expr3_t) << std::endl;
      REQUIRE_THAT(
          expr3_t,
          EquivalentTo(
              L"1/4 ã{i1,i2;a1,a2} * t¹{a1,a2;i1,i2;z1}:A-C-S * "
              L"h¹{κ2;κ1;z1,z2}:A-C-S * ã{κ1;κ2} + h¹{κ2;κ1;z1,z2}:A-C-S * "
              L"ã{i1;a1} * ã{κ1;κ2} * t¹{a1;i1;z1}:A-C-S + h¹{κ4;κ3;z1}:A-C-S "
              L"* h¹{κ2;κ1;z1,z2}:A-C-S * ã{κ3;κ4} * ã{κ1;κ2}"));

      // batching + adjoint: verify batch_ordinals are preserved
      auto h1_adj = adjoint(h1);
      REQUIRE(h1_adj.as<op_t>().batch_ordinals().has_value());
      REQUIRE(h1_adj.as<op_t>().batch_ordinals().value() ==
              h1.as<op_t>().batch_ordinals().value());

      auto h1_2_adj = adjoint(h1_2);
      REQUIRE(h1_2_adj.as<op_t>().batch_ordinals().value() ==
              h1_2.as<op_t>().batch_ordinals().value());

      // custom batch ordinals
      auto t = op::tʼ(2, {.batch_ordinals = {5, 10, 15}});
      auto t_adj = adjoint(t);
      REQUIRE(t_adj.as<op_t>().batch_ordinals().value() ==
              t.as<op_t>().batch_ordinals().value());
    }  // SECTION("batching")
  }

  SECTION("wick") {
    using namespace sequant;
    using namespace sequant::mbpt;
    namespace o = sequant::mbpt::op;
    namespace t = sequant::mbpt::tensor;

    SECTION("expectation values default to empty connectivity") {
      // Requiring every h to connect to each later t
      // removes this nonzero product.
      const auto expr = o::h(1) * o::t(1) * o::h(1) * o::t(1);
      const auto unconstrained = o::vac_av(expr, {});
      REQUIRE(unconstrained != ex<Constant>(0));
      CHECK_THAT(o::vac_av(expr), EquivalentTo(unconstrained));
      CHECK_THAT(o::ref_av(expr), EquivalentTo(unconstrained));
      CHECK(o::vac_av(expr, {.connect = {{L"f", L"t"}}}) == ex<Constant>(0));
    }

    SECTION("ref_av retains connectivity when reference equals vacuum") {
      REQUIRE(t::ref_av(t::h(1) * t::t(1) * t::h(1) * t::t(1),
                        {.connect = {{0, 3}}}) == ex<Constant>(0));
    }

    SECTION("SRSO"){
        // H**T12**T12 -> R2
        SECTION("wick(H**T12**T12 -> R2)"){
            auto result = t::vac_av(t::A(nₚ(-2)) * t::H(2) * t::T(2) * t::T(2),
                                    {.connect = {{1, 2}, {1, 3}}});

    //      std::wcout << "H*T12*T12 -> R2 = " << to_latex_align(result, 20)
    //                 << std::endl;
    REQUIRE(result->size() == 15);

    {
      // check against op
      auto result_op = o::vac_av(o::P(nₚ(2)) * o::H() * o ::T(2) * o::T(2),
                                 {.connect = default_op_connections()});
      REQUIRE(result_op->size() == result->size());  // as compact as result ..
      REQUIRE(simplify(result_op - result) ==
              ex<Constant>(0));  // .. and equivalent to it
    }
  }

  // H2**T3**T3 -> R4
  SECTION("wick(H2**T3**T3 -> R4)") {
    auto result = t::vac_av(t::A(nₚ(-4)) * t::h(2) * t::t(3) * t::t(3),
                            {.connect = {{1, 2}, {1, 3}}});

    // std::wcout << "H2**T3**T3 -> R4 = " << to_latex_align(result, 20)
    //            << std::endl;
    REQUIRE(result->size() == 4);
  }

#ifndef SEQUANT_SKIP_LONG_TESTS
  // the longest term in CCSDTQP
  // H2**T2**T2**T3 -> R5
  {
    ExprPtr ref_result;
    SECTION("wick(H2**T2**T2**T3 -> R5)") {
      ref_result = t::vac_av(t::A(-5) * t::h(2) * t::t(2) * t::t(2) * t::t(3),
                             {.connect = {{1, 2}, {1, 3}, {1, 4}}});
      REQUIRE(ref_result->size() == 7);
    }
  }
#endif  // !defined(SEQUANT_SKIP_LONG_TESTS)
}  // SECTION ("SRSO")

SECTION("SRSO Fock") {
  // <2p1h|H2|1p> ->
  SECTION("wick(<2p1h|H2|1p>)") {
    auto input = t::l(nₚ(2), nₕ(1)) * t::h(2) * t::r(nₚ(1), nₕ(0));
    auto result = t::vac_av(input);

    REQUIRE(result->is<Product>());  // product ...
    REQUIRE(result->size() == 3);    // ... of 3 factors
  }

  // <2p1h|H2|2p1h(c)> ->
  SECTION("wick(<2p1h|H2|2p1h(c)>)") {
    auto input = t::l(nₚ(2), nₕ(1)) * t::H() * t::r(nₚ(2), nₕ(1));
    auto result = t::vac_av(input);

    // std::wcout << "<2p1h|H|2p1h(c)> = " << to_latex(result)
    //            << std::endl;
    REQUIRE(result->is<Sum>());    // sub ...
    REQUIRE(result->size() == 4);  // ... of 4 factors
  }
}  // SECTION("SRSO Fock")

SECTION("vac_av with Power") {
  // keep amplitudes Hermitian so the reference equation below stays valid
  // (see scoped_hermitian_amplitudes)
  auto herm_amplitudes_guard = scoped_hermitian_amplitudes();
  auto expr1 = ex<Power>(ex<Variable>(L"α"), rational{2}) * t::h(2) * t::t(2);
  auto vev1 = t::vac_av(expr1, {.connect = {{0, 1}}});
  simplify(vev1);
  REQUIRE_THAT(
      vev1, EquivalentTo(
                L"1/4 * α^(2) * g{i1,i2;a1,a2}:A-C-S * t{a1,a2;i1,i2}:A-C-S"));

  auto expr2 = ex<Power>(ex<Constant>(2), rational{2}) * t::h(2) * t::t(2);
  auto vev2 = op::vac_av(expr2);
  simplify(vev2);
  // should not have Power object, 1/4 * 2^{2}
  REQUIRE_THAT(vev2,
               EquivalentTo(L"g{i1,i2;a1,a2}:A-C-S * t{a1,a2;i1,i2}:A-C-S"));
}

SECTION("SRSO-PNO") {
  using sequant::mbpt::Context;
  auto mbpt_ctx = sequant::mbpt::set_scoped_default_mbpt_context(
      Context({.csv = CSV::Yes, .op_registry_ptr = make_minimal_registry()}));

  // H2**T2 -> E
  SECTION("wick(H2**T2 -> E)") {
    // the contracted occupied indices of h and t are identified, including
    // where they are protoindices of t's virtuals; their overlap must not stand
    REQUIRE_THAT(t::vac_av(t::h(2) * t::t(2)),
                 EquivalentTo(L"1/4 g{i_1,i_2;a_1<i_1,i_2>,a_2<i_1,i_2>}:A-C-S "
                              L"* t{a_1<i_1,i_2>,a_2<i_1,i_2>;i_1,i_2}:A-N-S"));
  }

  // H2**T2**T2 -> R2
  SECTION("wick(H2**T2**T2 -> R2)") {
    auto result = t::vac_av(t::A(nₚ(-2)) * t::h(2) * t::t(2) * t::t(2),
                            {.connect = {{1, 2}, {1, 3}}});

    REQUIRE(result->size() == 4);
  }
}  // SECTION("SRSO-PNO")

SECTION("SRSF") {
  auto ctx = get_default_context();
  ctx.set(SPBasis::Spinfree);
  auto ctx_resetter = set_scoped_default_context(ctx);

  // H2 -> R2
  SECTION("wick(H2 -> R2)") {
    auto result = t::vac_av(t::S(-2) * t::h(2));

    {
      // std::wcout << "H2 -> R2 = " << to_latex_align(result, 0, 1)
      //            << std::endl;
      REQUIRE(result->is<Sum>());
      REQUIRE(result->size() == 2);
    }
  }

  // H2**T2 -> R2
  SECTION("wick(H2**T2 -> R2)") {
    auto result =
        t::vac_av(t::S(-2) * t::h(2) * t::t(2), {.connect = {{1, 2}}});

    {
      // std::wcout << "H2**T2 -> R2 = " << to_latex_align(result, 0, 1)
      //            << std::endl;
      REQUIRE(result->is<Sum>());
      REQUIRE(result->size() == 12);
    }
  }
  SECTION("S operator action") {
    auto expr = tensor::S(2) * deserialize(L"f{a1,a2;i1,i2}");
    simplify(expr);
    auto scalar = expr->as<Product>().scalar();
    REQUIRE(scalar == 1);
    REQUIRE_THAT(expr,
                 EquivalentTo("Ŝ{a_3,a_4;i_3,i_4} * "
                              "f{a_1,a_2;i_1,i_2}:N-N-S * ã{i_3,i_4;a_3,a_4}"));
  }
}  // SECTION("SRSF")

SECTION("MRSO") {
  auto ctx = get_default_context();
  ctx.set(mbpt::make_mr_spaces());
  auto ctx_resetter = set_scoped_default_context(ctx);

  SECTION("ref_av rejects connectivity when reference differs from vacuum") {
    const auto expr = t::h(1) * t::t(1);
    const auto unconstrained = t::ref_av(expr);
    REQUIRE(unconstrained != ex<Constant>(0));
    if (assert_behavior() != AssertBehavior::Abort) {
      // match the message: under THROW an unrelated SEQUANT_ASSERT throws
      // the same type
      const auto rejects_connections =
          Catch::Matchers::MessageMatches(Catch::Matchers::ContainsSubstring(
              "connect and do_not_connect must be empty"));
      REQUIRE_THROWS_MATCHES(t::ref_av(expr, {.connect = {{0, 1}}}), Exception,
                             rejects_connections);
      REQUIRE_THROWS_MATCHES(t::ref_av(expr, {.do_not_connect = {{0, 1}}}),
                             Exception, rejects_connections);
      // Operator requests must be rejected before screening or label lowering.
      REQUIRE_THROWS_MATCHES(
          o::ref_av(ex<Constant>(0), {.connect = {{L"f", L"t"}}}), Exception,
          rejects_connections);
      REQUIRE_THROWS_MATCHES(
          o::ref_av(o::h(1) * o::t(1), {.do_not_connect = {{L"f", L"t"}}}),
          Exception, rejects_connections);
    }
  }

  SECTION("wick(H2**T2 -> 0)") {
    {
      auto result = t::ref_av(t::h(2) * t::t(2));

      auto result_wo_top =
          t::ref_av(t::h(2) * t::t(2), {.use_topology = false});
      REQUIRE(simplify(result - result_wo_top) == ex<Constant>(0));
    }

    // now compute using physical vacuum
    {
      auto ctx = get_default_context();
      ctx.set(mbpt::make_mr_spaces());
      ctx.set(Vacuum::Physical);
      auto ctx_resetter = set_scoped_default_context(ctx);
      auto result_phys = t::ref_av(t::h(2) * t::t(2));
    }
  }

  // H2 ** T2 ** T2 -> 0
#ifndef SEQUANT_SKIP_LONG_TESTS
  // Cross-checks that use_topology's contraction-symmetry exploitation
  // agrees with the brute-force (use_topology=false) path; the brute-force
  // side dominates this section's runtime.
  SECTION("wick(H2**T2**T2 -> 0)") {
    // first without use of topology
    auto result =
        t::ref_av(t::h(2) * t::t(2) * t::t(2), {.use_topology = false});
    // now with topology use
    auto result_top =
        t::ref_av(t::h(2) * t::t(2) * t ::t(2), {.use_topology = true});

    REQUIRE(simplify(result - result_top) == ex<Constant>(0));
  }
#endif  // !defined(SEQUANT_SKIP_LONG_TESTS)

  // non-normal-ordered one-body operator in product form: h * a†_p * a_q
  SECTION("ref_av of non-normal-ordered one-body product") {
    const Index p{L"p_1"};  // complete space (spans core + active + virtual)
    const Index q{L"p_2"};
    // creator index (p) -> tensor bra, annihilator index (q) -> tensor ket
    // column symmetry is spelled out to match how mbpt::OpMaker builds h
    // everywhere else; the programmatic default is the conservative Nonsymm
    auto H1 = ex<Tensor>(L"h", bra{p}, ket{q}, Symmetry::Nonsymm,
                         BraKetSymmetry::Conjugate, ColumnSymmetry::Symm) *
              fcrex(p) * fannx(q);
    ExprPtr result;
    const auto* index_comparer = &get_default_context().index_comparer();
    REQUIRE_NOTHROW(result = t::ref_av(H1));
    // the active-first comparer is scoped to ref_av
    CHECK(&get_default_context().index_comparer() == index_comparer);
    REQUIRE_THAT(result, SimplifiesTo(L"h{O_1;O_1}:N-C-S + "
                                      L"h{u_2;u_1}:N-C-S * γ{u_1;u_2}"));
  }

#if 0
    // H**T12 -> R2
    SECTION("wick(H**T2 -> R2)") {
      auto result = t::ref_av(t::A(-2) * t::H() * t::t(2), {.connect = {{1, 2}}});

      {
        std::wcout << "H*T2 -> R2 = " << to_latex_align(result, 0, 1)
                   << std::endl;
      }
    }
#endif
}  // SECTION("MRSO")

SECTION("MRSO-MultiProduct") {
  auto ctx = get_default_context();
  ctx.set(mbpt::make_mr_spaces());
  ctx.set(Vacuum::MultiProduct);
  auto ctx_resetter = set_scoped_default_context(ctx);

  // one-body: same expectation as the core-vacuum path at MRSO
  SECTION("ref_av of non-normal-ordered one-body product") {
    const Index p{L"p_1"};
    const Index q{L"p_2"};
    auto H1 = ex<Tensor>(L"h", bra{p}, ket{q}, Symmetry::Nonsymm,
                         BraKetSymmetry::Conjugate, ColumnSymmetry::Symm) *
              fcrex(p) * fannx(q);
    ExprPtr result;
    REQUIRE_NOTHROW(result = t::ref_av(H1));
    REQUIRE_THAT(result, SimplifiesTo(L"h{O_1;O_1}:N-C-S + "
                                      L"h{u_2;u_1}:N-C-S * γ{u_1;u_2}"));
  }

  // the mbpt operators (ã) are normal-ordered relative to the context vacuum,
  // i.e. to the reference here and to the core under SingleProduct, so
  // t::h(k)·t::t(k) is a different operator on the two paths; only products
  // of elementary operators, whose normal order is immaterial, are compared.
  // connect is not used: it means different things on the two paths (a core
  // or virtual δ between the operators vs any density linking them)
  SECTION("elementary operators match the core-vacuum path") {
    const Index p1{L"p_1"}, p2{L"p_2"}, p3{L"p_3"}, p4{L"p_4"}, p5{L"p_5"},
        p6{L"p_6"};
    auto coeff = [](IndexList b, IndexList k) {
      return ex<Tensor>(L"h", bra(b), ket(k), Symmetry::Nonsymm,
                        BraKetSymmetry::Nonsymm, ColumnSymmetry::Nonsymm);
    };
    // the operators carry the vacuum of the context they are built under,
    // so the input is built under each
    auto check = [](auto&& make_x) {
      const auto mp =
          mbpt::decompositions::cumulants_to_densities(t::ref_av(make_x()));
      ExprPtr sp;
      {
        auto sp_ctx = get_default_context();
        sp_ctx.set(Vacuum::SingleProduct);
        auto sp_resetter = set_scoped_default_context(sp_ctx);
        sp = t::ref_av(make_x());
      }
      // the two paths spell the same sums differently, e.g. h{E;O} vs
      // h{a;O} + h{g;O} + h{u;O}, and η vs δ - γ
      REQUIRE(simplify(in_base_spaces(mp) - in_base_spaces(sp)) ==
              ex<Constant>(0));
    };
    // one-body
    check([&] { return coeff({p1}, {p2}) * fcrex(p1) * fannx(p2); });
    // two-body: up to κ₂
    check([&] {
      return coeff({p1, p2}, {p3, p4}) * fcrex(p1) * fcrex(p2) * fannx(p4) *
             fannx(p3);
    });
    // a two-body times a one-body string: up to κ₃
    check([&] {
      return coeff({p1, p2, p5}, {p3, p4, p6}) * fcrex(p1) * fcrex(p2) *
             fannx(p4) * fannx(p3) * fcrex(p5) * fannx(p6);
    });
#ifndef SEQUANT_SKIP_LONG_TESTS
    // two two-body strings: up to κ₄
    const Index p7{L"p_7"}, p8{L"p_8"};
    check([&] {
      return coeff({p1, p2, p5, p6}, {p3, p4, p7, p8}) * fcrex(p1) * fcrex(p2) *
             fannx(p4) * fannx(p3) * fcrex(p5) * fcrex(p6) * fannx(p8) *
             fannx(p7);
    });
#endif  // !defined(SEQUANT_SKIP_LONG_TESTS)
  }

  // relative to a fixed-N reference, a GNO string {x} is defined by
  // x = Σ ± Π κ(B) {x minus the legs of every B} over all sets of disjoint
  // blocks B of the ops of x, each with as many creators as annihilators
  // (a pair is a γ or an η), the sign being that of moving the blocks' legs,
  // in order, to the front; κ(B) is the ordered cumulant of the reference
  // averages of B's substrings. Built so from core-vacuum averages, {x} is an
  // elementary-operator expression, which makes products of GNO strings,
  // number-conserving or not, comparable on the two paths
  SECTION("GNO products match the core-vacuum path") {
    struct Leg {
      Index index;
      Action action;
    };
    using Legs = container::svector<Leg>;
    // a GNO string from its creators and its annihilators in storage order
    auto gno_string = [](IndexList cre_idxs, IndexList ann_idxs) {
      Legs legs;
      for (const auto& i : cre_idxs) legs.push_back({i, Action::Create});
      for (const auto& i : ann_idxs) legs.push_back({i, Action::Annihilate});
      return legs;
    };
    auto balanced = [](const Legs& legs, const container::svector<int>& pos) {
      long net = 0;
      for (auto p : pos) net += legs[p].action == Action::Create ? 1 : -1;
      return net == 0;
    };
    // the parity of the permutation that lists @p order
    auto parity = [](const container::svector<int>& order) {
      int sign = 1;
      for (std::size_t i = 0; i != order.size(); ++i)
        for (std::size_t j = i + 1; j != order.size(); ++j)
          if (order[i] > order[j]) sign = -sign;
      return sign;
    };
    // calls f(sign, blocks, rest) for every partition of @p pos into blocks
    // of ≥ 2 legs with as many creators as annihilators and, if
    // @p with_rest, a rest of single legs
    using Blocks = container::svector<container::svector<int>>;
    auto for_each_partition = [&](const Legs& legs,
                                  const container::svector<int>& pos,
                                  bool with_rest, auto&& f) {
      Blocks blocks;
      container::svector<int> rest;
      container::svector<bool> used(legs.size(), false);
      auto recurse = [&](auto&& self) -> void {
        auto first = std::find_if(pos.begin(), pos.end(),
                                  [&](int p) { return !used[p]; });
        if (first == pos.end()) {
          container::svector<int> order;
          for (const auto& b : blocks)
            order.insert(order.end(), b.begin(), b.end());
          order.insert(order.end(), rest.begin(), rest.end());
          f(parity(order), blocks, rest);
          return;
        }
        const int p0 = *first;
        used[p0] = true;
        if (with_rest) {
          rest.push_back(p0);
          self(self);
          rest.pop_back();
        }
        container::svector<int> others;
        for (auto p : pos)
          if (!used[p]) others.push_back(p);
        for (std::size_t m = 1; m < (std::size_t{1} << others.size()); ++m) {
          container::svector<int> block{p0};
          for (std::size_t k = 0; k != others.size(); ++k)
            if (m & (std::size_t{1} << k)) block.push_back(others[k]);
          if (!balanced(legs, block)) continue;
          for (auto p : block) used[p] = true;
          blocks.push_back(block);
          self(self);
          blocks.pop_back();
          for (auto p : block) used[p] = p == p0;
        }
        used[p0] = false;
      };
      recurse(recurse);
    };
    // an elementary-operator expression: Σ c·(string of elementary operators)
    using Expansion = container::svector<std::pair<ExprPtr, Legs>>;
    auto subset = [](const Legs& legs, const container::svector<int>& pos) {
      Legs result;
      for (auto p : pos) result.push_back(legs[p]);
      return result;
    };
    // the reference average of an elementary string via the standard theorem
    // under the core vacuum: a surviving all-active string is a density, any
    // other survivor averages to 0, and so does an unbalanced string, the
    // reference having a fixed N; call under the SingleProduct vacuum
    auto core_vacuum_average = [&](const Legs& legs) -> ExprPtr {
      container::svector<int> all(legs.size());
      std::iota(all.begin(), all.end(), 0);
      if (!balanced(legs, all)) return ex<Constant>(0);
      if (legs.empty()) return ex<Constant>(1);
      ExprPtr string = ex<Constant>(1);
      for (const auto& l : legs)
        string = string *
                 (l.action == Action::Create ? fcrex(l.index) : fannx(l.index));
      auto contracted = FWickTheorem{string}.full_contractions(false).compute();
      const auto isr = get_default_context().index_space_registry();
      const auto& active = isr->retrieve(L"u");
      auto average = [&](ExprPtr& f) {
        if (f->is<FNOperator>()) {
          const auto& nop = f->as<FNOperator>();
          const bool all_active = ranges::all_of(nop, [&](const auto& op) {
            return op.index().space() == active;
          });
          f = all_active && nop.ncreators() == nop.nannihilators()
                  ? density::rdm_from_nop(nop, density::rdm_label())
                  : ex<Constant>(0);
        } else if (f->is<Tensor>() &&
                   f->as<Tensor>().label() == reserved::overlap_label()) {
          f = make_kronecker(f->as<Tensor>().bra()[0],
                             f->as<Tensor>().ket()[0]);
        }
      };
      ExprPtr result = std::make_shared<Sum>();
      for (const auto& term : contracted->is<Sum>()
                                  ? contracted->as<Sum>().summands() |
                                        ranges::to<container::svector<ExprPtr>>
                                  : container::svector<ExprPtr>{contracted}) {
        ExprPtr t = term->is<Product>()
                        ? term->clone()
                        : ex<Product>(ExprPtrList{term->clone()});
        t->visit(average, /*atoms_only=*/true);
        result = result + t;
      }
      return simplify(result);
    };
    // {x} as an elementary-operator expression with core-vacuum averages;
    // call under the SingleProduct vacuum
    auto gno_from_core_vacuum = [&](const Legs& legs) {
      std::map<container::svector<int>, ExprPtr> cumulants;
      std::map<container::svector<int>, Expansion> gnos;
      auto cumulant = [&](auto&& self,
                          const container::svector<int>& pos) -> ExprPtr {
        if (auto it = cumulants.find(pos); it != cumulants.end())
          return it->second->clone();
        ExprPtr result = core_vacuum_average(subset(legs, pos));
        for_each_partition(legs, pos, /*with_rest=*/false,
                           [&](int sign, const Blocks& blocks, const auto&) {
                             if (blocks.size() < 2) return;
                             ExprPtr term = ex<Constant>(-sign);
                             for (const auto& b : blocks)
                               term = term * self(self, b);
                             result = result + term;
                           });
        expand(result);
        cumulants.emplace(pos, result);
        return result->clone();
      };
      auto gno = [&](auto&& self,
                     const container::svector<int>& pos) -> Expansion {
        if (auto it = gnos.find(pos); it != gnos.end()) return it->second;
        Expansion result{{ex<Constant>(1), subset(legs, pos)}};
        for_each_partition(
            legs, pos, /*with_rest=*/true,
            [&](int sign, const Blocks& blocks,
                const container::svector<int>& rest) {
              if (blocks.empty()) return;
              ExprPtr c = ex<Constant>(-sign);
              for (const auto& b : blocks) c = c * cumulant(cumulant, b);
              for (const auto& [c_rest, string] : self(self, rest))
                result.emplace_back(c * c_rest->clone(), string);
            });
        gnos.emplace(pos, result);
        return result;
      };
      container::svector<int> all(legs.size());
      std::iota(all.begin(), all.end(), 0);
      return gno(gno, all);
    };
    // a leg over a union of base spaces (e.g. I = i ∪ u) is split into them on
    // the core-vacuum side, whose averages take base-space legs; every density
    // leg of the MultiProduct result must be active
    const auto isr = get_default_context().index_space_registry();
    auto check = [&](std::initializer_list<Legs> strings_il) {
      container::svector<Legs> strings(strings_il);
      container::set<Index> externals;
      for (const auto& legs : strings)
        for (const auto& l : legs)
          if (!isr->is_base(l.index.space())) externals.insert(l.index);
      ExprPtr mp_input = ex<Constant>(1);
      for (const auto& legs : strings) {
        container::svector<Index> cre_idxs, ann_idxs;
        for (const auto& l : legs)
          (l.action == Action::Create ? cre_idxs : ann_idxs).push_back(l.index);
        // the ctor takes annihilators in particle order, the reverse of
        // storage
        std::reverse(ann_idxs.begin(), ann_idxs.end());
        mp_input = mp_input * ex<FNOperator>(cre(cre_idxs), ann(ann_idxs));
      }
      const auto mp = mbpt::decompositions::cumulants_to_densities(
          FWickTheorem{mp_input}.compute());
      // all base-space assignments of the legs; each union leg multiplies
      // them by the number of its base spaces
      container::svector<container::svector<Legs>> assignments{strings};
      for (std::size_t s_i = 0; s_i != strings.size(); ++s_i)
        for (std::size_t l_i = 0; l_i != strings[s_i].size(); ++l_i) {
          const auto& idx = strings[s_i][l_i].index;
          if (isr->is_base(idx.space())) continue;
          decltype(assignments) next;
          for (const auto& asg : assignments)
            for (const auto& base : isr->base_spaces())
              if (base.qns() == idx.space().qns() &&
                  idx.space().type().includes(base.type())) {
                auto a2 = asg;
                a2[s_i][l_i].index = external_in_base(idx, base, externals);
                next.push_back(a2);
              }
          assignments = std::move(next);
        }
      ExprPtr sp = ex<Constant>(0);
      {
        auto sp_ctx = get_default_context();
        sp_ctx.set(Vacuum::SingleProduct);
        auto sp_resetter = set_scoped_default_context(sp_ctx);
        for (const auto& asg : assignments) {
          Expansion product{{ex<Constant>(1), Legs{}}};
          for (const auto& legs : asg) {
            Expansion next;
            for (const auto& [c1, s1] : product)
              for (const auto& [c2, s2] : gno_from_core_vacuum(legs)) {
                Legs st = s1;
                st.insert(st.end(), s2.begin(), s2.end());
                next.emplace_back(c1->clone() * c2->clone(), std::move(st));
              }
            product = std::move(next);
          }
          for (const auto& [c, st] : product)
            sp = sp + c->clone() * core_vacuum_average(st);
        }
      }
      // every density leg of the MultiProduct result is active
      bool active_only = true;
      mp->visit(
          [&](const ExprPtr& e) {
            if (!e->is<Tensor>()) return;
            const auto& t = e->as<Tensor>();
            if (t.label() == L"γ" || t.label() == L"η" || t.label() == L"κ")
              for (const auto& idx : t.const_braket())
                active_only = active_only && idx.space() == isr->retrieve(L"u");
          },
          /*atoms_only=*/true);
      INFO("mp: " << toUtf8(to_latex(mp)));
      CHECK(active_only);
      REQUIRE(simplify(in_base_spaces(mp, externals) -
                       in_base_spaces(sp, externals)) == ex<Constant>(0));
    };
    const Index u1{L"u_1"}, u2{L"u_2"}, u3{L"u_3"}, u4{L"u_4"}, u5{L"u_5"},
        u6{L"u_6"}, i1{L"i_1"}, i2{L"i_2"}, a1{L"a_1"}, a2{L"a_2"};
    // number-conserving
    check({gno_string({u1}, {u2}), gno_string({u3}, {u4})});
    check({gno_string({u1}, {i1}), gno_string({a1}, {u2})});
    // non-conserving
    check({gno_string({u1, u2}, {u3}), gno_string({}, {u4})});
    check({gno_string({}, {u1}), gno_string({u2}, {})});
    check({gno_string({u1}, {u3, u2}), gno_string({u4, u5}, {})});
    check({gno_string({}, {u1}), gno_string({u2}, {u3}), gno_string({u4}, {})});
    check({gno_string({u1}, {}), gno_string({u2}, {})});
    check({gno_string({}, {i1}), gno_string({i2}, {u1}), gno_string({u2}, {})});
    check({gno_string({a1}, {i1, u1}), gno_string({u2}, {a2}),
           gno_string({i2}, {})});
    // up to κ₃; the 2-body string's own κ₂ enters its definition
    check({gno_string({u1, u2}, {u4, u3}), gno_string({u5}, {u6})});
    check({gno_string({u1, u2}, {u3}), gno_string({u4}, {u6, u5})});
    // legs over unions of the active space with others
    const Index I1{L"I_1"}, I2{L"I_2"}, A1{L"A_1"}, A2{L"A_2"}, M1{L"M_1"},
        M2{L"M_2"}, E1{L"E_1"}, E2{L"E_2"}, p1{L"p_1"}, p2{L"p_2"}, p3{L"p_3"},
        p4{L"p_4"};
    check({gno_string({I1}, {A1}), gno_string({A2}, {I2})});
    check({gno_string({M1}, {E1}), gno_string({E2}, {M2})});
    check({gno_string({M1}, {M2}), gno_string({E1}, {E2})});
    check({gno_string({p1}, {p2}), gno_string({p3}, {p4})});
    // with a κ₃
    check({gno_string({I1, A1}, {I2, A2}), gno_string({u5}, {u6})});
    // non-conserving
    check({gno_string({M1, E1}, {I1}), gno_string({}, {A1})});

    // partial contractions are operator-valued, so the two paths are compared
    // in a common normal form: the GNO remainders of the MultiProduct result
    // are replaced by their elementary-operator expansions and both sides are
    // normal-ordered relative to the core vacuum by the standard theorem
    // a lone term is returned as is: a Product whose factor is a one-term Sum
    // is not flattened by expand, and the theorem would take it for a scalar
    auto to_expr = [](const Expansion& expansion) {
      auto result = std::make_shared<Sum>();
      for (const auto& [c, legs] : expansion) {
        ExprPtr string = c->clone();
        for (const auto& l : legs)
          string = string * (l.action == Action::Create ? fcrex(l.index)
                                                        : fannx(l.index));
        result->append(string);
      }
      return result->size() == 1 ? result->summand(0) : ExprPtr(result);
    };
    // call under the SingleProduct vacuum; the standard theorem spells a
    // contraction as the overlap s, the extended one as δ
    auto core_vacuum_normal_form = [&](ExprPtr x) {
      expand(x);
      auto result = FWickTheorem{x}.full_contractions(false).compute();
      auto overlap_to_kronecker = [](ExprPtr& f) {
        if (f->is<Tensor>() &&
            f->as<Tensor>().label() == reserved::overlap_label())
          f = make_kronecker(f->as<Tensor>().bra()[0],
                             f->as<Tensor>().ket()[0]);
      };
      if (result->is_atom())
        overlap_to_kronecker(result);
      else
        result->visit(overlap_to_kronecker, /*atoms_only=*/true);
      return in_base_spaces(result);
    };
    auto check_partial = [&](std::initializer_list<Legs> strings) {
      ExprPtr mp_input = ex<Constant>(1);
      for (const auto& legs : strings) {
        container::svector<Index> cre_idxs, ann_idxs;
        for (const auto& l : legs)
          (l.action == Action::Create ? cre_idxs : ann_idxs).push_back(l.index);
        std::reverse(ann_idxs.begin(), ann_idxs.end());
        mp_input = mp_input * ex<FNOperator>(cre(cre_idxs), ann(ann_idxs));
      }
      const auto mp = mbpt::decompositions::cumulants_to_densities(
          FWickTheorem{mp_input}.full_contractions(false).compute());
      ExprPtr lhs, rhs;
      {
        auto sp_ctx = get_default_context();
        sp_ctx.set(Vacuum::SingleProduct);
        auto sp_resetter = set_scoped_default_context(sp_ctx);
        auto expand_gno = [&](const ExprPtr& f) {
          Legs legs;
          for (const auto& op : f->as<FNOperator>())
            legs.push_back({op.index(), op.action()});
          return to_expr(gno_from_core_vacuum(legs));
        };
        // each term's GNO remainder, if any, is replaced by its expansion
        ExprPtr elementary = ex<Constant>(0);
        for (const auto& term :
             mp->is<Sum>() ? mp->as<Sum>().summands() |
                                 ranges::to<container::svector<ExprPtr>>
                           : container::svector<ExprPtr>{mp}) {
          if (term->is<FNOperator>()) {
            elementary = elementary + expand_gno(term);
          } else if (term->is<Product>()) {
            ExprPtr coefficient = ex<Constant>(term->as<Product>().scalar());
            ExprPtr remainder;
            for (const auto& factor : term->as<Product>().factors())
              if (factor->is<FNOperator>())
                remainder = expand_gno(factor);
              else
                coefficient = coefficient * factor->clone();
            elementary = elementary +
                         (remainder ? coefficient * remainder : coefficient);
          } else {
            elementary = elementary + term->clone();
          }
        }
        lhs = core_vacuum_normal_form(elementary);
        ExprPtr product = ex<Constant>(1);
        for (const auto& legs : strings)
          product = product * to_expr(gno_from_core_vacuum(legs));
        rhs = core_vacuum_normal_form(product);
      }
      INFO("lhs: " << toUtf8(to_latex(lhs))
                   << "\nrhs: " << toUtf8(to_latex(rhs)));
      REQUIRE(!lhs->is<Constant>());
      REQUIRE(simplify(lhs - rhs) == ex<Constant>(0));
    };
    // number-conserving
    check_partial({gno_string({u1}, {u2}), gno_string({u3}, {u4})});
    // non-conserving
    check_partial({gno_string({u1, u2}, {u3}), gno_string({}, {u4})});
    check_partial({gno_string({}, {u1}), gno_string({u2}, {})});
    check_partial({gno_string({u1}, {u3, u2}), gno_string({u4, u5}, {})});
    check_partial(
        {gno_string({}, {i1}), gno_string({i2}, {u1}), gno_string({u2}, {})});
    check_partial({gno_string({a1}, {i1, u1}), gno_string({u2}, {a2}),
                   gno_string({i2}, {})});
    check_partial({gno_string({u1, u2}, {u3}), gno_string({u4}, {u6, u5})});
  }

  SECTION("wick(H2**T2) runs in generalized normal order") {
    ExprPtr result;
    REQUIRE_NOTHROW(result =
                        t::ref_av(t::h(2) * t::t(2), {.connect = {{0, 1}}}));
    REQUIRE(!result->is<Constant>());
    // the product reaches κ₄
    ExprPtr densities;
    REQUIRE_NOTHROW(densities =
                        mbpt::decompositions::cumulants_to_densities(result));
    REQUIRE(!densities->is<Constant>());
  }

  SECTION("topology on/off agree") {
    auto a = t::ref_av(t::h(2) * t::t(2), {.connect = {{0, 1}}});
    auto b = t::ref_av(t::h(2) * t::t(2),
                       {.connect = {{0, 1}}, .use_topology = false});
    REQUIRE(simplify(a - b) == ex<Constant>(0));
  }

  SECTION("topology prunes contractions") {
    auto attempted = [](bool top) {
      FWickTheorem wick{simplify(t::h(2) * t::t(2))};
      wick.use_topology(top).set_nop_connections({{0, 1}});
      wick.compute();
      return wick.stats().num_attempted_contractions.load();
    };
    const auto attempted_on = attempted(true);
    REQUIRE(attempted_on > 0);
    REQUIRE(attempted_on < attempted(false));
  }

#ifndef SEQUANT_SKIP_LONG_TESTS
  SECTION("wick(H2**T2**T2) topology on/off agree") {
    auto a = t::ref_av(t::h(2) * t::t(2) * t::t(2), {.connect = {{0, 1}}});
    auto b = t::ref_av(t::h(2) * t::t(2) * t::t(2),
                       {.connect = {{0, 1}}, .use_topology = false});
    REQUIRE(simplify(a - b) == ex<Constant>(0));
  }
#endif  // !defined(SEQUANT_SKIP_LONG_TESTS)

  SECTION("operator-level ref_av agrees with tensor-level") {
    auto result_op = o::ref_av(o::h(2) * o::t(2));
    auto result_t = t::ref_av(t::h(2) * t::t(2), {.connect = {{0, 1}}});
    REQUIRE(simplify(result_op - result_t) == ex<Constant>(0));
  }

  SECTION("ref_av honors connectivity through densities and cumulants") {
    // ⟨{a†_u1 a_u2}{a†_u3 a_u4}⟩ = γ η + κ: every term connects the two
    // operators, through a density or a cumulant only
    const auto x = ex<FNOperator>(cre({L"u_1"}), ann({L"u_2"})) *
                   ex<FNOperator>(cre({L"u_3"}), ann({L"u_4"}));
    const auto all = t::ref_av(x);
    REQUIRE(all != ex<Constant>(0));
    REQUIRE(simplify(t::ref_av(x, {.connect = {{0, 1}}}) - all) ==
            ex<Constant>(0));
    REQUIRE(t::ref_av(x, {.do_not_connect = {{0, 1}}}) == ex<Constant>(0));
  }

  SECTION("ref_av connectivity partitions the terms") {
    // h·t·t has terms with the two t connected and terms without; each list
    // keeps one part, at the operator and the tensor level, and vac_av
    // agrees with ref_av
    const auto x = t::h(2) * t::t(1) * t::t(1);
    const auto all = t::ref_av(x);
    const auto connected = t::ref_av(x, {.connect = {{1, 2}}});
    const auto disconnected = t::ref_av(x, {.do_not_connect = {{1, 2}}});
    REQUIRE(connected != ex<Constant>(0));
    REQUIRE(disconnected != ex<Constant>(0));
    REQUIRE(simplify(connected + disconnected - all) == ex<Constant>(0));
    REQUIRE(simplify(t::vac_av(x, {.connect = {{1, 2}}}) - connected) ==
            ex<Constant>(0));
    const auto x_op = o::h(2) * o::t(1) * o::t(1);
    REQUIRE(simplify(o::ref_av(x_op, {.connect = {{L"t", L"t"}}}) -
                     connected) == ex<Constant>(0));
    REQUIRE(simplify(o::ref_av(x_op, {.do_not_connect = {{L"t", L"t"}}}) -
                     disconnected) == ex<Constant>(0));
  }
}  // SECTION("MRSO-MultiProduct")

SECTION("MRSF") {
  // now compute using (closed) Fermi vacuum + spinfree basis
  auto ctx = get_default_context();
  ctx.set(mbpt::make_mr_spaces());
  ctx.set(SPBasis::Spinfree);
  auto ctx_resetter = set_scoped_default_context(ctx);

  SECTION("wick(H2**T2 -> 0)") {
    auto result = t::ref_av(t::h(2) * t::t(2));

    {
      // make sure get same result without use of topology
      auto result_wo_top =
          t::ref_av(t::h(2) * t::t(2), {.use_topology = false});

      REQUIRE(simplify(result - result_wo_top) == ex<Constant>(0));
    }

    {
      // make sure get same result using operators
      auto result_op = o::ref_av(o::h(2) * o::t(2));

      REQUIRE(result_op->size() == result->size());
      REQUIRE(simplify(result - result_op) == ex<Constant>(0));
    }
  }
}  // SECTION("MRSF")
}

SECTION("rules") {
  using namespace sequant;

  SECTION("density-fit") {
    const std::vector<std::wstring> inputs = {
        L"t{a1,a2;i1,i2} t{a3;i3}",
        L"t{a1,a2;i1,i2} g{a3;i3}",
        L"t{a1,a2;i1,i2} g{i1,i2;a1,a2}",
        L"t{a1,a2;i1,i2} g{i1,i2;a1,a2}:A",
    };
    const std::vector<std::wstring> expected = {
        L"t{a1,a2;i1,i2} t{a3;i3}",
        L"t{a1,a2;i1,i2} g{a3;i3}",
        L"t{a1,a2;i1,i2} B{i1;a1;x_1}:N-C-S B{i2;a2;x_1}:N-C-S",
        L"t{a1,a2;i1,i2} (B{i1;a1;x_1}:N-C-S B{i2;a2;x_1}:N-C-S "
        "- B{i2;a1;x_1}:N-C-S B{i1;a2;x_1}:N-C-S)",
    };

    REQUIRE(inputs.size() == expected.size());

    for (std::size_t i = 0; i < inputs.size(); ++i) {
      CAPTURE(inputs.at(i));

      ExprPtr input_expr = deserialize(inputs.at(i));

      const IndexSpace aux_space =
          get_default_context().index_space_registry()->retrieve(L"x");

      ExprPtr actual = mbpt::density_fit(input_expr, aux_space, L"g", L"B");

      REQUIRE_THAT(actual, EquivalentTo(expected.at(i)));
    }
  }

  SECTION("tensor-hypercontract") {
    const std::vector<std::wstring> inputs = {
        L"t{a1,a2;i1,i2} t{a3;i3}",
        L"t{a1,a2;i1,i2} g{i1,i2;a1,a2}",
        L"t{a1,a2;i1,i2} g{i1,i2;a1,a2}:A",
    };
    const std::vector<std::wstring> expected = {
        L"t{a1,a2;i1,i2} t{a3;i3}",
        L"t{a1,a2;i1,i2} B{i1;;x_1} B{;a1;x_1} C{;;x_1,x_2} B{i2;;x_2} "
        L"B{;a2;x_2}",
        L"t{a1,a2;i1,i2} (B{i1;;x_1} B{;a1;x_1} C{;;x_1,x_2} B{i2;;x_2} "
        L"B{;a2;x_2}"
        " - B{i2;;x_1} B{;a1;x_1} C{;;x_1,x_2} B{i1;;x_2} B{;a2;x_2})",
    };

    REQUIRE(inputs.size() == expected.size());

    for (std::size_t i = 0; i < inputs.size(); ++i) {
      CAPTURE(inputs.at(i));

      ExprPtr input_expr = deserialize(inputs.at(i));

      const IndexSpace aux_space =
          get_default_context().index_space_registry()->retrieve(L"x");

      ExprPtr actual =
          mbpt::tensor_hypercontract(input_expr, aux_space, L"g", L"B", L"C");

      REQUIRE_THAT(actual, EquivalentTo(expected.at(i)));
    }
  }
}  // SECTION("rules")

SECTION("manuscript-examples") {
  using namespace sequant::mbpt;

  using sequant::reserved::antisymm_label;

  /// When order == 0, returns sim. transformed of zeroth order Hamiltonian.
  /// When order > 0, returns sim. transformed perturbation operator or given
  /// order.
  auto H̅ = [](size_t order = 0) {
    auto hbar0 = lst(op::H(), T(2), 4, {.use_connected_form = true});
    if (order == 0) return hbar0;
    // only one-body perturbation operator
    auto hbar_pt = lst(op::Hʼ(/*rank*/ 1, {.order = order}), T(2), 2,
                       {.use_connected_form = true});
    return hbar_pt;
  };

  SECTION("CCD Term") {
    auto expr = ref_av(P(2) * H() * t(2) * t(2),
                       {.connect = {{L"f", L"t"}, {L"g", L"t"}}});
    REQUIRE(expr.size() == 4);
  }

  SECTION("Ground State Amplitudes") {
    // connectivity info for t and λ amplitude equations
    const auto t_connect = default_op_connections();
    const auto l_connect =
        concat(default_op_connections(),
               OpConnections<std::wstring>{{L"h", antisymm_label()},
                                           {L"f", antisymm_label()},
                                           {L"g", antisymm_label()}});

    auto t = ref_av(P(2) * H̅(), {.connect = t_connect});
    auto λ = ref_av((1 + Λ(2)) * H̅() * P(-2), {.connect = l_connect});

    // numer of terms are verified against srcc results
    REQUIRE(t.size() == 31);
    REQUIRE(λ.size() == 32);
  }

  SECTION("CC LR Function") {
    const int N = 2;  // CC rank

    auto θ̅ = lst(θ(1), T(N), 2, {.use_connected_form = true});
    auto expr = (1 + Λ(N)) * θ̅ * Tʼ(N) + Λʼ(N) * θ̅;
    auto result = ref_av(expr, {.connect = {{L"θ", L"t"}, {L"θ", L"t¹"}}});
    // number of terms is verified against MPQC4 implementation
    REQUIRE(result.size() == 21);
  }

#ifndef SEQUANT_SKIP_LONG_TESTS
  SECTION("EOM-CC Equations") {
    // connectivity info for right and left amplitude equations
    const auto r_connect = concat(
        default_op_connections(),
        OpConnections<std::wstring>{{L"h", L"R"}, {L"f", L"R"}, {L"g", L"R"}});
    const auto l_connect =
        concat(default_op_connections(),
               OpConnections<std::wstring>{{L"h", antisymm_label()},
                                           {L"f", antisymm_label()},
                                           {L"g", antisymm_label()}});

    // EE
    auto r_EE = ref_av(P(2) * H̅() * R(2), {.connect = r_connect});
    auto l_EE = ref_av(L(2) * H̅() * P(-2), {.connect = l_connect});
    // EA
    auto r_EA =
        ref_av(P(nₚ(2), nₕ(1)) * H̅() * R(nₚ(2), nₕ(1)), {.connect = r_connect});
    auto l_EA = ref_av(L(nₚ(2), nₕ(1)) * H̅() * P(nₚ(-2), nₕ(-1)),
                       {.connect = l_connect});
    // IP
    auto r_IP =
        ref_av(P(nₚ(1), nₕ(2)) * H̅() * R(nₚ(1), nₕ(2)), {.connect = r_connect});
    auto l_IP = ref_av(L(nₚ(1), nₕ(2)) * H̅() * P(nₚ(-1), nₕ(-2)),
                       {.connect = l_connect});

    // number of terms are verified against eomcc results
    REQUIRE(r_EE.size() == 53);
    REQUIRE(l_EE.size() == 31);
    REQUIRE(r_EA.size() == 32);
    REQUIRE(l_EA.size() == 24);
    REQUIRE(r_IP.size() == 32);
    REQUIRE(l_IP.size() == 24);
  }
#endif  // !defined(SEQUANT_SKIP_LONG_TESTS)

#ifndef SEQUANT_SKIP_LONG_TESTS
  SECTION("CC Perturbed Amplitudes") {
    // connectivity info for perturbed t and λ amplitude equations
    const auto t_connect =
        concat(default_op_connections(),
               OpConnections<std::wstring>{
                   {L"h", L"t¹"}, {L"f", L"t¹"}, {L"g", L"t¹"}, {L"h¹", L"t"}});

    const auto l_connect =
        concat(default_op_connections(),
               OpConnections<std::wstring>{{L"h", L"t¹"},
                                           {L"f", L"t¹"},
                                           {L"g", L"t¹"},
                                           {L"h¹", L"t"},
                                           {L"h", antisymm_label()},
                                           {L"f", antisymm_label()},
                                           {L"g", antisymm_label()},
                                           {L"h¹", antisymm_label()}});

    // perturbed t amplitudes (Eq 18 in SQ Manuscript #2)
    auto t = ref_av(P(2) * (H̅(1) + H̅() * Tʼ(2) - "ω" * tʼ(2)),
                    {.connect = t_connect});
    // perturbed λ amplitudes (Eq 19 in SQ Manuscript #2)
    auto λ = ref_av(
        ((1 + Λ(2)) * (H̅(1) + H̅() * Tʼ(2)) + Λʼ(2) * H̅() + "ω" * λʼ(2)) * P(-2),
        {.connect = l_connect});

    // number of terms are verified against MPQC4 implementation
    REQUIRE(t.size() == 58);
    REQUIRE(λ.size() == 63);
  }
#endif  // !defined(SEQUANT_SKIP_LONG_TESTS)

  SECTION("Custom MBPT Operators") {
    using namespace sequant;
    using namespace sequant::mbpt;
    // get two IndexSpaces for custom excitation operator, could be any spaces
    const auto& cre_space = get_particle_space(Spin::any);
    const auto& ann_space = get_hole_space(Spin::any);

    OpRegistry registry;
    registry.add(L"f", OpClass::Gen, Hermiticity::Hermitian)
        .add(L"g", OpClass::Gen, Hermiticity::Hermitian)
        .add(L"t", OpClass::Ex, Hermiticity::NonHermitian)
        .add(L"x", OpClass::Ex, Hermiticity::NonHermitian)
        .add(L"y", OpClass::Ex, Hermiticity::NonHermitian);

    // set MBPT context
    auto ctx_resetter = set_scoped_default_mbpt_context(
        {.csv = CSV::No, .op_registry = registry});

    // use OpMaker to define a custom excitation operator of rank 2
    auto x = OpMaker<Statistics::FermiDirac>(L"x", 2)();

    // particle non-conserving excitation operator with custom IndexSpaces
    auto y = OpMaker<Statistics::FermiDirac>(L"y", ncre(2), nann(1),
                                             cre(cre_space), ann(ann_space))();

    auto expr = tensor::H() * x * y;
    REQUIRE_THAT(
        simplify(expr),
        EquivalentTo(
            L"1/8 f{p1;p2}:A-C-S x{a1,a2;i1,i2}:A-N-S y{a3,a4;i3}:A-N-S "
            L"ã{p2;p1} ã{i1,i2;a1,a2} ã{i3;a3,a4} + 1/32 g{p3,p4;p1,p2}:A-C-S "
            L"x{a1,a2;i1,i2}:A-N-S y{a3,a4;i3}:A-N-S ã{p1,p2;p3,p4} "
            L"ã{i1,i2;a1,a2} ã{i3;a3,a4}"));
  }
}  // SECTION("manuscript-examples")

SECTION("avoided-connections") {
  using namespace sequant::mbpt;
  using sequant::reserved::antisymm_label;

  // keep amplitudes Hermitian so the reference equation below stays valid
  // (see scoped_hermitian_amplitudes)
  auto herm_amplitudes_guard = scoped_hermitian_amplitudes();

  // h(1) * t(1): two operators, avoid the only possible contraction
  auto expr1 = tensor::h(1) * tensor::t(1);
  auto res1 = tensor::vac_av(expr1, {.do_not_connect = {{0, 1}}});
  REQUIRE(res1 == sequant::ex<sequant::Constant>(0));  // result should be zero

  // P(1) * H() * T(2): avoid connections between projector and Hamiltonian
  auto expr2 = tensor::P(1) * tensor::H() * tensor::T(2);
  auto res2_full = tensor::vac_av(expr2, {.connect = {{1, 2}}});
  auto res2 =
      tensor::vac_av(expr2, {.connect = {{1, 2}}, .do_not_connect = {{0, 1}}});
  REQUIRE(res2_full.size() == 6);
  // only one term with no A-{f,g} connection
  REQUIRE(res2.is<sequant::Product>());
  const std::wstring expected2 =
      L"-1 Â{i_1;a_1} t{a_1,a_2;i_2,i_1}:A-C-S f{i_2;a_2}:A-C-S";
  REQUIRE_THAT(sequant::simplify(res2), EquivalentTo(expected2));

  // same test as above but from Operator level and labels for connectivity
  using namespace sequant::mbpt::op;
  auto expr3 = op::P(1) * op::H(2) * op::T(2);
  auto res3 = op::vac_av(expr3, {.connect = {{L"f", L"t"}, {L"g", L"t"}},
                                 .do_not_connect = {{antisymm_label(), L"f"},
                                                    {antisymm_label(), L"g"}}});
  REQUIRE_THAT(sequant::simplify(res3), EquivalentTo(expected2));

  // projectors are never connected
  auto expr4 = op::P(1) * op::H() * op::t(2) * op::P(-1);
  auto res4_full = op::vac_av(expr4);
  auto res4 = op::vac_av(
      expr4, {.connect = op::default_op_connections(),
              .do_not_connect = {{antisymm_label(), antisymm_label()}}});
  REQUIRE(res4_full.size() == 4);
  REQUIRE(res4.is<sequant::Product>());  // only single term survives
  const std::wstring expected4 =
      L"Â{i_1;a_2} Â{a_1;i_2} g{i_3,i_2;a_3,a_1}:A-C-S "
      L"t{a_3,a_2;i_3,i_1}:A-C-S";
  REQUIRE_THAT(simplify(res4), EquivalentTo(expected4));
}

SECTION("rdm-decomposition symmetries") {
  using namespace sequant;

  // an RDM is Hermitian and particle (column) symmetric by definition: the γ
  // that the decompositions in mbpt/rdm.cpp build has to equal the γ that
  // expectation_value_impl() (mbpt/op.cpp) and the parser build, or
  // otherwise-equal terms stop merging. The symmetries take part in the tensor
  // hash, so a mismatch in either attribute is enough to break it.
  auto ctx_resetter = set_scoped_default_context(
      Context({.index_space_registry_shared_ptr = mbpt::make_sr_spaces(),
               .vacuum = Vacuum::SingleProduct}));

  const auto kappa =
      density::make_cumulant(bra{Index(L"i_1")}, ket{Index(L"i_2")});
  const auto gamma = mbpt::decompositions::cumulant_to_density(kappa);
  REQUIRE(gamma->is<Tensor>());
  REQUIRE(gamma->as<Tensor>().label() == L"γ");
  REQUIRE(gamma->as<Tensor>().hermiticity() == Hermiticity::Hermitian);
  REQUIRE(gamma->as<Tensor>().column_symmetry() == ColumnSymmetry::Symm);

  // the decomposition, the factory in SeQuant/core/density.hpp and the parser
  // agree on the spelling
  const auto gamma_factory = density::make_rdm(Index(L"i_1"), Index(L"i_2"));
  REQUIRE(*gamma == *gamma_factory);
  REQUIRE(*deserialize(L"γ{i_1;i_2}") == *gamma_factory);
  REQUIRE(*deserialize(L"κ{i_1;i_2}") == *kappa);

  const auto eta = density::make_hole_rdm(Index(L"i_1"), Index(L"i_2"));
  REQUIRE(eta->as<Tensor>().label() == L"η");
  REQUIRE(eta->as<Tensor>().hermiticity() == Hermiticity::Hermitian);
  REQUIRE(eta->as<Tensor>().column_symmetry() == ColumnSymmetry::Symm);

  // κ from a 2-body normal operator: bra = annihilators, ket = creators
  const FNOperator nop2(cre({L"i_1", L"i_2"}), ann({L"i_3", L"i_4"}));
  const auto kappa2 = density::make_cumulant(nop2);
  REQUIRE(kappa2->as<Tensor>().label() == L"κ");
  REQUIRE(kappa2->as<Tensor>().symmetry() == Symmetry::Antisymm);
  REQUIRE(kappa2->as<Tensor>().bra()[0] == Index(L"i_3"));
  REQUIRE(kappa2->as<Tensor>().bra()[1] == Index(L"i_4"));
  REQUIRE(kappa2->as<Tensor>().ket()[0] == Index(L"i_1"));
  REQUIRE(kappa2->as<Tensor>().ket()[1] == Index(L"i_2"));
  REQUIRE(*deserialize(L"κ{i_3,i_4;i_1,i_2}") == *kappa2);
  // antisymmetrize() generates each distinct pairing once, also under the
  // topological canonicalization that the test suite defaults to
  {
    auto ctx = get_default_context();
    ctx.set(CanonicalizeOptions::default_options().copy_and_set(
        CanonicalizationMethod::Topological));
    auto topological = set_scoped_default_context(ctx);
    const Index i1(L"i_1"), i2(L"i_2"), i3(L"i_3"), i4(L"i_4");
    const auto kappa = density::make_cumulant(bra{i1, i3}, ket{i2, i4});
    const auto expected =
        density::make_rdm(bra{i1, i3}, ket{i2, i4}) -
        density::make_rdm(i1, i2) * density::make_rdm(i3, i4) +
        density::make_rdm(i1, i4) * density::make_rdm(i3, i2);
    REQUIRE(simplify(mbpt::decompositions::cumulant2_to_density(kappa) -
                     expected) == ex<Constant>(0));
  }
  // every density a decomposition builds has these symmetries, including
  // those that antisymmetrize() builds by permuting indices
  const auto densities2 =
      simplify(mbpt::decompositions::cumulant2_to_density(kappa2));
  densities2->visit(
      [](const ExprPtr& e) {
        if (!e->is<Tensor>()) return;
        REQUIRE(e->as<Tensor>().label() == L"γ");
        REQUIRE(e->as<Tensor>().hermiticity() == Hermiticity::Hermitian);
        REQUIRE(e->as<Tensor>().column_symmetry() == ColumnSymmetry::Symm);
      },
      /*atoms_only=*/true);
  // a spin-orbital multi-body density is antisymmetric, like the one
  // expectation_value_impl() builds from a leftover normal operator
  const auto gamma2 = density::rdm_from_nop(nop2, density::rdm_label());
  REQUIRE(simplify(densities2 - gamma2)->size() == densities2->size() - 1);
  // η is a reserved label, hence an mbpt operator label
  REQUIRE(mbpt::to_op_class(L"η") == mbpt::OpClass::Gen);

  // every multi-body κ that a decomposition builds is the κ that the extended
  // Wick theorem builds
  {
    const FNOperator nop3(cre({L"i_1", L"i_2", L"i_3"}),
                          ann({L"i_4", L"i_5", L"i_6"}));
    const auto decomp =
        simplify(mbpt::decompositions::three_body_decomp(ex<FNOperator>(nop3),
                                                         /*approx=*/false)
                     .first);
    std::size_t n_multibody_kappa = 0;
    decomp->visit(
        [&n_multibody_kappa](const ExprPtr& e) {
          if (!e->is<Tensor>() || e->as<Tensor>().label() != L"κ" ||
              e->as<Tensor>().rank() < 2)
            return;
          ++n_multibody_kappa;
          REQUIRE(e->as<Tensor>().symmetry() == Symmetry::Antisymm);
          REQUIRE(e->as<Tensor>().hermiticity() == Hermiticity::Hermitian);
          REQUIRE(e->as<Tensor>().column_symmetry() == ColumnSymmetry::Symm);
        },
        /*atoms_only=*/true);
    REQUIRE(n_multibody_kappa > 0);
    REQUIRE(simplify(decomp - density::make_cumulant(nop3))->size() ==
            decomp->size() - 1);
  }
}

SECTION("cumulant-to-density decompositions") {
  using namespace sequant;
  auto ctx_resetter = set_scoped_default_context(
      Context({.index_space_registry_shared_ptr = mbpt::make_sr_spaces(),
               .vacuum = Vacuum::SingleProduct}));
  auto gamma = [](std::vector<Index> b, std::vector<Index> k) {
    return b.size() == 1
               ? density::make_rdm(b[0], k[0])
               : density::make_rdm(bra(std::move(b)), ket(std::move(k)));
  };

  // κ₃ = γ₃ - Σ γ₁γ₂ (9 terms) + 2 Σ γ₁γ₁γ₁ (6 terms)
  const FNOperator nop3(cre({L"i_1", L"i_2", L"i_3"}),
                        ann({L"i_4", L"i_5", L"i_6"}));
  const auto densities3 = simplify(
      mbpt::decompositions::cumulant3_to_density(density::make_cumulant(nop3)));
  REQUIRE(densities3->is<Sum>());
  std::size_t n_gamma3 = 0, n_gamma1_gamma2 = 0, n_gamma1_cubed = 0;
  for (const auto& term : *densities3) {
    if (term->is<Tensor>()) {
      REQUIRE(term->as<Tensor>().label() == L"γ");
      REQUIRE(term->as<Tensor>().rank() == 3);
      REQUIRE(term->as<Tensor>().symmetry() == Symmetry::Antisymm);
      ++n_gamma3;
      continue;
    }
    const auto& product = term->as<Product>();
    for (const auto& f : product) REQUIRE(f->as<Tensor>().label() == L"γ");
    if (product.size() == 2) {
      REQUIRE(abs(product.scalar()) == 1);
      ++n_gamma1_gamma2;
    } else {
      REQUIRE(product.size() == 3);
      REQUIRE(abs(product.scalar()) == 2);
      ++n_gamma1_cubed;
    }
  }
  REQUIRE(n_gamma3 == 1);
  REQUIRE(n_gamma1_gamma2 == 9);
  REQUIRE(n_gamma1_cubed == 6);
  // the identity pairing of γ₁γ₁γ₁ comes with +2
  const auto identity = ex<Constant>(2) *
                        density::make_rdm(Index(L"i_4"), Index(L"i_1")) *
                        density::make_rdm(Index(L"i_5"), Index(L"i_2")) *
                        density::make_rdm(Index(L"i_6"), Index(L"i_3"));
  REQUIRE(simplify(densities3 - identity)->size() == densities3->size() - 1);

  // the general cumulant_to_density reproduces the hand-derived κ₂ and κ₃
  {
    using mbpt::antisymmetrize;
    using mbpt::decompositions::cumulant_to_density;
    const Index i1(L"i_1"), i2(L"i_2"), i3(L"i_3"), i4(L"i_4"), i5(L"i_5"),
        i6(L"i_6");
    const auto kappa2_ref =
        gamma({i3, i4}, {i1, i2}) -
        antisymmetrize(gamma({i3}, {i1}) * gamma({i4}, {i2})).result;
    REQUIRE(simplify(cumulant_to_density(density::make_cumulant(
                         FNOperator(cre({i1, i2}), ann({i3, i4})))) -
                     kappa2_ref) == ex<Constant>(0));
    const auto kappa3_ref =
        gamma({i4, i5, i6}, {i1, i2, i3}) -
        antisymmetrize(gamma({i4}, {i1}) * gamma({i5, i6}, {i2, i3})).result +
        ex<Constant>(2) * antisymmetrize(gamma({i4}, {i1}) * gamma({i5}, {i2}) *
                                         gamma({i6}, {i3}))
                              .result;
    REQUIRE(simplify(cumulant_to_density(density::make_cumulant(nop3)) -
                     kappa3_ref) == ex<Constant>(0));
  }

  // antisymmetrize() keeps the bra-ket symmetry of the tensors it permutes
  // as spelled, not only their Hermiticity: over a complex field a tensor is
  // bra-ket Symm without being Conjugate
  {
    auto cidx = [](std::wstring_view label) {
      Index i{label};
      IndexSpace space = i.space();
      space.field(Field::Complex);
      return Index(space, i.ordinal());
    };
    auto t = [&](std::wstring_view b, std::wstring_view k) {
      return ex<Tensor>(L"t", bra{cidx(b)}, ket{cidx(k)}, Symmetry::Nonsymm,
                        BraKetSymmetry::Symm, ColumnSymmetry::Symm);
    };
    const auto permuted =
        mbpt::antisymmetrize(t(L"i_3", L"i_1") * t(L"i_4", L"i_2")).result;
    REQUIRE(permuted->is<Sum>());
    REQUIRE(permuted->size() == 2);
    permuted->visit(
        [](const ExprPtr& e) {
          if (e->is<Tensor>())
            REQUIRE(e->as<Tensor>().braket_symmetry() == BraKetSymmetry::Symm);
        },
        /*atoms_only=*/true);
  }

  // cumulants_to_densities rewrites every κ in an expression
  using mbpt::decompositions::cumulants_to_densities;
  const FNOperator nop2(cre({L"i_1", L"i_2"}), ann({L"i_3", L"i_4"}));
  const auto kappa2 = density::make_cumulant(nop2);
  const auto kappa1 =
      density::make_cumulant(bra{Index(L"i_5")}, ket{Index(L"i_6")});
  const auto g = ex<Tensor>(L"g", bra{L"i_1", L"i_2"}, ket{L"i_3", L"i_4"},
                            Symmetry::Antisymm);
  const auto h = ex<Tensor>(L"h", bra{L"i_6"}, ket{L"i_5"});
  const auto expr = g * kappa2 + h * kappa1;
  const auto expected = g * mbpt::decompositions::cumulant2_to_density(kappa2) +
                        h * mbpt::decompositions::cumulant_to_density(kappa1);
  const auto result = cumulants_to_densities(expr);
  result->visit(
      [](const ExprPtr& e) {
        if (e->is<Tensor>()) REQUIRE(e->as<Tensor>().label() != L"κ");
      },
      /*atoms_only=*/true);
  REQUIRE(simplify(result - expected) == ex<Constant>(0));
  // ... including a κ that is the whole expression
  REQUIRE(simplify(cumulants_to_densities(density::make_cumulant(nop3)) -
                   densities3) == ex<Constant>(0));

  // κ₄ = γ₄ - Σ γ₁γ₃ (4·4 terms) - Σ γ₂γ₂ (6·6/2 terms)
  //      + 2 Σ γ₁γ₁γ₂ (12·12/2 terms) - 6 Σ γ₁γ₁γ₁γ₁ (4!·4!/4! terms)
  const FNOperator nop4(cre({L"i_1", L"i_2", L"i_3", L"i_4"}),
                        ann({L"i_5", L"i_6", L"i_7", L"i_8"}));
  const auto densities4 = simplify(
      mbpt::decompositions::cumulant_to_density(density::make_cumulant(nop4)));
  REQUIRE(densities4->is<Sum>());
  std::size_t n_gamma4 = 0, n_gamma1_gamma3 = 0, n_gamma2_gamma2 = 0,
              n_gamma1_gamma1_gamma2 = 0, n_gamma1_fourth = 0;
  for (const auto& term : *densities4) {
    if (term->is<Tensor>()) {
      REQUIRE(term->as<Tensor>().label() == L"γ");
      REQUIRE(term->as<Tensor>().rank() == 4);
      REQUIRE(term->as<Tensor>().symmetry() == Symmetry::Antisymm);
      ++n_gamma4;
      continue;
    }
    const auto& product = term->as<Product>();
    for (const auto& f : product) REQUIRE(f->as<Tensor>().label() == L"γ");
    switch (product.size()) {
      case 2:
        REQUIRE(abs(product.scalar()) == 1);
        if (product.factors()[0]->as<Tensor>().rank() == 2)
          ++n_gamma2_gamma2;
        else
          ++n_gamma1_gamma3;
        break;
      case 3:
        REQUIRE(abs(product.scalar()) == 2);
        ++n_gamma1_gamma1_gamma2;
        break;
      default:
        REQUIRE(product.size() == 4);
        REQUIRE(abs(product.scalar()) == 6);
        ++n_gamma1_fourth;
    }
  }
  REQUIRE(n_gamma4 == 1);
  REQUIRE(n_gamma1_gamma3 == 16);
  REQUIRE(n_gamma2_gamma2 == 18);
  REQUIRE(n_gamma1_gamma1_gamma2 == 72);
  REQUIRE(n_gamma1_fourth == 24);
  // the identity pairing of each product comes with (-1)^(m-1) (m-1)!
  {
    const Index i1(L"i_1"), i2(L"i_2"), i3(L"i_3"), i4(L"i_4"), i5(L"i_5"),
        i6(L"i_6"), i7(L"i_7"), i8(L"i_8");
    for (const auto& identity :
         {ex<Constant>(-1) * gamma({i5}, {i1}) *
              gamma({i6, i7, i8}, {i2, i3, i4}),
          ex<Constant>(-1) * gamma({i5, i6}, {i1, i2}) *
              gamma({i7, i8}, {i3, i4}),
          ex<Constant>(2) * gamma({i5}, {i1}) * gamma({i6}, {i2}) *
              gamma({i7, i8}, {i3, i4}),
          ex<Constant>(-6) * gamma({i5}, {i1}) * gamma({i6}, {i2}) *
              gamma({i7}, {i3}) * gamma({i8}, {i4})}) {
      REQUIRE(simplify(densities4 - identity)->size() ==
              densities4->size() - 1);
    }
    // ... and the expansions of κ₂, κ₃ and κ₄ invert the expansion of γ₄ in
    // cumulants, γ₄ = κ₄ + A[γ₁κ₃] + A[κ₂κ₂] + A[γ₁γ₁κ₂] + A[γ₁γ₁γ₁γ₁]
    using mbpt::antisymmetrize;
    auto kappa = [](std::vector<Index> b, std::vector<Index> k) {
      return density::make_cumulant(bra(std::move(b)), ket(std::move(k)));
    };
    const auto moments =
        kappa({i5, i6, i7, i8}, {i1, i2, i3, i4}) +
        antisymmetrize(gamma({i5}, {i1}) * kappa({i6, i7, i8}, {i2, i3, i4}))
            .result +
        antisymmetrize(kappa({i5, i6}, {i1, i2}) * kappa({i7, i8}, {i3, i4}))
            .result +
        antisymmetrize(gamma({i5}, {i1}) * gamma({i6}, {i2}) *
                       kappa({i7, i8}, {i3, i4}))
            .result +
        antisymmetrize(gamma({i5}, {i1}) * gamma({i6}, {i2}) *
                       gamma({i7}, {i3}) * gamma({i8}, {i4}))
            .result;
    REQUIRE(simplify(cumulants_to_densities(moments) -
                     gamma({i5, i6, i7, i8}, {i1, i2, i3, i4})) ==
            ex<Constant>(0));
  }
  REQUIRE(simplify(cumulants_to_densities(density::make_cumulant(nop4)) -
                   densities4) == ex<Constant>(0));
}
}
