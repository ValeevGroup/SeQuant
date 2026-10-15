//
// Wick's theorem on operator products whose legs carry a basis instance.
//

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include "catch2_sequant.hpp"
#include "csv_test_utils.hpp"

#include <SeQuant/core/basis.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/wick.hpp>
#include <SeQuant/domain/mbpt/bernoulli.hpp>
#include <SeQuant/domain/mbpt/context.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/op.hpp>
#include <SeQuant/domain/mbpt/vac_av.hpp>

#include <algorithm>
#include <cstddef>
#include <map>
#include <optional>

TEST_CASE("basis-wick-contractions", "[algorithms][wick][basis]") {
  using namespace sequant;
  using namespace sequant::mbpt;
  using sequant::tests::csv::message_contains;
  using sequant::tests::csv::standing_metrics;
  using sequant::tests::csv::tensors_labelled;
  using sequant::tests::csv::with_leg_instance;
  namespace t = op::tensor;

  const auto csv = GENERATE(CSV::Yes, CSV::No);
  const bool yes = csv == CSV::Yes;
  auto ctx = set_scoped_default_context(
      sequant::Context({.index_basis_registry_shared_ptr = make_min_sr_spaces(),
                        .vacuum = Vacuum::SingleProduct,
                        .spbasis = SPBasis::Spinor}));
  auto mbpt_ctx = set_scoped_default_mbpt_context(
      mbpt::Context({.csv = csv, .op_registry_ptr = make_minimal_registry()}));
  // one product: Wick runs on the calling thread, so a throw reaches the test
  auto vac_av = [](ExprPtr const& product) {
    Index::reset_tmp_index();
    return t::vac_av(product);
  };

  SECTION("K1: an integral leg meets a granted amplitude leg") {
    CHECK(serialize(vac_av(t::h(2) * with_leg_instance(t::t(2), 1))) ==
          (yes ? L"1/4 g{i_1,i_2;a_1<i_1,i_2;1>,a_2<i_1,i_2;1>}:A-C-S * "
                 L"t{a_1<i_1,i_2;1>,a_2<i_1,i_2;1>;i_1,i_2}:A-N-S"
               : L"1/4 g{i_1,i_2;a_1<;1>,a_2<;1>}:A-C-S * "
                 L"t{a_1<;1>,a_2<;1>;i_1,i_2}:A-N-S"));
  }

  SECTION("K2: projector and amplitude legs of different pairs, one family") {
    CHECK(serialize(vac_av(with_leg_instance(t::P(nₚ(1)), 1) * t::h(1) *
                           with_leg_instance(t::t(2), 1))) ==
          (yes ? L"-1 Â{i_1;a_1<i_1;1>}:A-N-S * f{i_2;a_2<i_1,i_2;1>}:A-C-S * "
                 L"s{a_1<i_1;1>;a_3<i_1,i_2;1>}:N-C-S * "
                 L"t{a_2<i_1,i_2;1>,a_3<i_1,i_2;1>;i_1,i_2}:A-N-S"
               : L"Â{i_1;a_1<;1>}:A-N-S * f{i_2;a_2<;1>}:A-C-S * "
                 L"t{a_1<;1>,a_2<;1>;i_1,i_2}:A-N-S"));
  }

  SECTION("K3: λ and t legs of one family merge") {
    CHECK(serialize(vac_av(with_leg_instance(t::λ(2), 1) *
                           with_leg_instance(t::t(2), 1))) ==
          (yes ? L"1/4 t{a_1<i_1,i_2;1>,a_2<i_1,i_2;1>;i_1,i_2}:A-N-S * "
                 L"λ{i_1,i_2;a_1<i_1,i_2;1>,a_2<i_1,i_2;1>}:A-N-S"
               : L"1/4 t{a_1<;1>,a_2<;1>;i_1,i_2}:A-N-S * "
                 L"λ{i_1,i_2;a_1<;1>,a_2<;1>}:A-N-S"));
  }

  SECTION("K4: λ and t legs of two families leave the overlaps standing") {
    const auto r =
        vac_av(with_leg_instance(t::λ(2), 2) * with_leg_instance(t::t(2), 1));
    CHECK(serialize(r) ==
          (yes ? L"1/4 s{a_1<i_1,i_2;2>;a_2<i_1,i_2;1>}:N-C-S * "
                 L"s{a_3<i_1,i_2;2>;a_4<i_1,i_2;1>}:N-C-S * "
                 L"t{a_2<i_1,i_2;1>,a_4<i_1,i_2;1>;i_1,i_2}:A-N-S * "
                 L"λ{i_1,i_2;a_1<i_1,i_2;2>,a_3<i_1,i_2;2>}:A-N-S"
               : L"1/4 s{a_1<;2>;a_2<;1>}:N-C-S * s{a_3<;2>;a_4<;1>}:N-C-S * "
                 L"t{a_2<;1>,a_4<;1>;i_1,i_2}:A-N-S * "
                 L"λ{i_1,i_2;a_1<;2>,a_3<;2>}:A-N-S"));
    CHECK(standing_metrics(r) == 2);
  }

  SECTION("K5: ungranted legs next to granted ones take their instance") {
    const auto products = {t::λ(2) * with_leg_instance(t::t(2), 1),
                           with_leg_instance(t::λ(2), 1) * t::t(2)};
    for (auto const& product : products) {
      CHECK(serialize(vac_av(product)) ==
            (yes ? L"1/4 t{a_1<i_1,i_2;1>,a_2<i_1,i_2;1>;i_1,i_2}:A-N-S * "
                   L"λ{i_1,i_2;a_1<i_1,i_2;1>,a_2<i_1,i_2;1>}:A-N-S"
                 : L"1/4 t{a_1<;1>,a_2<;1>;i_1,i_2}:A-N-S * "
                   L"λ{i_1,i_2;a_1<;1>,a_2<;1>}:A-N-S"));
    }
  }

  SECTION("K6: general-space legs of two families contract to one overlap") {
    // neither leg is pure, so the contraction is δ{p_1<;1>;q<;1>}
    // s{q<;1>;q'<;2>} δ{q'<;2>;p_2<;2>}: a projection onto the quasiparticle
    // space keeps the instance of the leg it projects, and the overlap
    // between the two families stands
    const Index i1(L"i_1"), i2(L"i_2");
    const auto p1 = Index(L"p_1").replace_basis_instance(1);
    const auto p2 = Index(L"p_2").replace_basis_instance(2);
    {
      const auto c = FWickTheorem::contract(FOp(p1, Action::Annihilate),
                                            FOp(p2, Action::Create));
      for (auto const* t : tensors_labelled(c, reserved::kronecker_label()))
        CHECK(t->_bra()[0].basis().basis_instance() ==
              t->_ket()[0].basis().basis_instance());
      const auto metrics = tensors_labelled(c, reserved::overlap_label());
      REQUIRE(metrics.size() == 1);
      CHECK(metrics[0]->_bra()[0].basis().basis_instance() == 1);
      CHECK(metrics[0]->_ket()[0].basis().basis_instance() == 2);
    }
    // legs in one instance (or none) contract to a Kronecker delta, and so
    // do a generic leg and a specific one: the identity narrows the generic
    // one (an index without an instance runs over every instance)
    for (auto const& [l, r] :
         {std::pair{p1, Index(L"p_2").replace_basis_instance(1)},
          std::pair{Index(L"p_1"), Index(L"p_2")},
          std::pair{Index(L"p_1"), Index(L"p_2").replace_basis_instance(1)}}) {
      const auto c = FWickTheorem::contract(FOp(l, Action::Annihilate),
                                            FOp(r, Action::Create));
      CHECK(tensors_labelled(c, reserved::overlap_label()).empty());
      CHECK(tensors_labelled(c, reserved::kronecker_label()).size() == 3);
    }
    const auto expr = ex<Tensor>(L"t", bra{i1}, ket{p1}, Symmetry::Nonsymm) *
                      ex<FNOperator>(cre{i1}, ann{p1}) *
                      ex<Tensor>(L"u", bra{p2}, ket{i2}, Symmetry::Nonsymm) *
                      ex<FNOperator>(cre{p2}, ann{i2});
    Index::reset_tmp_index();
    ExprPtr r;
    REQUIRE_NOTHROW(r = FWickTheorem{expr}.full_contractions(true).compute());
    CHECK(tensors_labelled(r, reserved::kronecker_label()).empty());
    const auto metrics = tensors_labelled(r, reserved::overlap_label());
    REQUIRE(metrics.size() == 1);
    container::svector<IndexBasis::optional_instance> instances;
    for (auto const& idx : metrics[0]->_slots())
      instances.push_back(idx.basis().basis_instance());
    CHECK(
        (instances == container::svector<IndexBasis::optional_instance>{1, 2} ||
         instances == container::svector<IndexBasis::optional_instance>{2, 1}));
  }

  SECTION("a projection resolves its instance to the registered name") {
    // a name belongs to one (space, instance) pair: the projection of a leg
    // in p̃ (p in instance 1) onto the particle space is in instance 1 of a,
    // i.e. in ã if that is registered, not in p̃
    auto isr = make_min_sr_spaces();
    isr->add(L"p̃", IndexBasis{isr->retrieve_space(L"p"), 1});
    isr->add(L"ã", IndexBasis{isr->retrieve_space(L"a"), 1});
    auto nctx = get_default_context();
    nctx.set(std::move(isr));
    auto nctx_resetter = set_scoped_default_context(nctx);
    const Index x(L"p̃_1");
    REQUIRE(x.basis().name() == L"p̃");
    const auto c = FWickTheorem::contract(
        FOp(x, Action::Annihilate),
        FOp(Index(L"p_2").replace_basis_instance(2), Action::Create));
    const auto deltas = tensors_labelled(c, reserved::kronecker_label());
    REQUIRE(deltas.size() == 2);
    std::size_t projected = 0;
    for (auto const* t : deltas)
      for (auto const& idx : t->_slots())
        if (idx.basis().basis_instance() == 1 && idx != x) {
          CHECK(idx.basis().name() == L"ã");
          ++projected;
        }
    CHECK(projected == 1);
  }

  SECTION("K0: no instance anywhere") {
    CHECK(serialize(vac_av(t::h(2) * t::t(2))) ==
          (yes ? L"1/4 g{i_1,i_2;a_1<i_1,i_2>,a_2<i_1,i_2>}:A-C-S * "
                 L"t{a_1<i_1,i_2>,a_2<i_1,i_2>;i_1,i_2}:A-N-S"
               : L"1/4 g{i_1,i_2;a_1,a_2}:A-C-S * t{a_1,a_2;i_1,i_2}:A-N-S"));
    CHECK(serialize(vac_av(t::P(nₚ(1)) * t::h(1) * t::t(2))) ==
          (yes ? L"-1 Â{i_1;a_1<i_1>}:A-N-S * f{i_2;a_2<i_1,i_2>}:A-C-S * "
                 L"s{a_1<i_1>;a_3<i_1,i_2>}:N-C-S * "
                 L"t{a_2<i_1,i_2>,a_3<i_1,i_2>;i_1,i_2}:A-N-S"
               : L"Â{i_1;a_1}:A-N-S * f{i_2;a_2}:A-C-S * "
                 L"t{a_1,a_2;i_1,i_2}:A-N-S"));
    CHECK(serialize(vac_av(t::λ(2) * t::t(2))) ==
          (yes ? L"1/4 t{a_1<i_1,i_2>,a_2<i_1,i_2>;i_1,i_2}:A-N-S * "
                 L"λ{i_1,i_2;a_1<i_1,i_2>,a_2<i_1,i_2>}:A-N-S"
               : L"1/4 t{a_1,a_2;i_1,i_2}:A-N-S * λ{i_1,i_2;a_1,a_2}:A-N-S"));
  }

  SECTION("K14: CSV mode without instances, legacy registry") {
    if (yes) {
      auto legacy_ctx = set_scoped_default_mbpt_context(mbpt::Context(
          {.csv = CSV::Yes, .op_registry_ptr = make_legacy_registry()}));
      Index::reset_tmp_index();
      CHECK(serialize(op::vac_av(op::H(1) * op::T(1))) ==
            L"f{i_1;a_1<i_1>}:A-C-S * t{a_1<i_1>;i_1}:A-N-S");
    }
  }
}

TEST_CASE("basis-wick-reduce", "[algorithms][wick][basis]") {
  using namespace sequant;
  using sequant::tests::csv::message_contains;
  using sequant::tests::csv::standing_metrics;
  using sequant::tests::csv::tensors_labelled;

  auto ctx = set_scoped_default_context(sequant::Context(
      {.index_basis_registry_shared_ptr = mbpt::make_min_sr_spaces(),
       .vacuum = Vacuum::SingleProduct,
       .spbasis = SPBasis::Spinor}));

  const Index p1(L"p_1"), q1(L"p_2"), a1(L"a_1"), a2(L"a_2"), i1(L"i_1"),
      i2(L"i_2");
  const auto a1_1 = a1.replace_basis_instance(1);
  const auto a2_1 = a2.replace_basis_instance(1);
  const auto a2_2 = a2.replace_basis_instance(2);
  // indices repeated in the product are dummy, the others external
  auto reduce = [](ExprPtr expr) {
    Index::reset_tmp_index();
    FWickTheorem(expr).reduce(expr);
    return expr;
  };
  auto g = [&](Index const& x, Index const& y) {
    return ex<Tensor>(L"g", bra{x, y}, ket{i1, i2}, Symmetry::Nonsymm);
  };
  auto instances = [](AbstractTensor const& t) {
    container::svector<IndexBasis::optional_instance> result;
    for (auto const& idx : t._slots())
      result.push_back(idx.basis().basis_instance());
    return result;
  };
  auto slots = [](AbstractTensor const& t) {
    container::svector<Index> result;
    for (auto const& idx : t._slots()) result.push_back(idx);
    return result;
  };

  SECTION("a generic index meeting a specific one takes its instance") {
    // s{p_1;a_1} s{p_1;a_2<;1>} g{a_1,a_2<;1>;i_1,i_2}: all three virtual
    // indices become one, in instance 1, whichever overlap comes first
    const auto s1 = make_overlap(p1, a1), s2 = make_overlap(p1, a2_1);
    for (auto const& product : {s1 * s2 * g(a1, a2_1), s2 * s1 * g(a1, a2_1)}) {
      const auto r = reduce(product);
      REQUIRE(tensors_labelled(r, reserved::overlap_label()).empty());
      const auto gs = tensors_labelled(r, L"g");
      REQUIRE(gs.size() == 1);
      const auto g_slots = slots(*gs[0]);
      CHECK(g_slots[0] == g_slots[1]);
      CHECK(instances(*gs[0]) ==
            container::svector<IndexBasis::optional_instance>{1, 1, {}, {}});
    }
  }

  SECTION("an overlap between two instances stands, even through a delta") {
    // s{p_1;a_1<;1>} δ{p_1;p_2} s{p_2;a_2<;2>} g{a_1<;1>,a_2<;2>;i_1,i_2}
    // = <a_1<;1>|a_2<;2>> g: one overlap between the two families stands
    const auto s1 = make_overlap(p1, a1_1), d = make_kronecker(p1, q1),
               s2 = make_overlap(q1, a2_2);
    for (auto const& product :
         {s1 * d * s2 * g(a1_1, a2_2), s2 * d * s1 * g(a1_1, a2_2)}) {
      const auto r = reduce(product);
      CHECK(tensors_labelled(r, reserved::kronecker_label()).empty());
      const auto metrics = tensors_labelled(r, reserved::overlap_label());
      REQUIRE(metrics.size() == 1);
      CHECK(standing_metrics(r) == 1);
      const auto gs = tensors_labelled(r, L"g");
      REQUIRE(gs.size() == 1);
      CHECK(instances(*gs[0]) ==
            container::svector<IndexBasis::optional_instance>{1, 2, {}, {}});
      // the overlap connects g's two virtual slots
      const auto g_slots = slots(*gs[0]);
      const auto m_slots = slots(*metrics[0]);
      CHECK(((m_slots[0] == g_slots[0] && m_slots[1] == g_slots[1]) ||
             (m_slots[0] == g_slots[1] && m_slots[1] == g_slots[0])));
    }
  }

  SECTION("an overlap between two generic indices is resolved last") {
    // δ{a_1<;1>;p_1} s{p_1;p_2} δ{p_2;a_2<;2>} g{a_1<;1>,a_2<;2>;i_1,i_2}: the
    // deltas identify p_1 and p_2 with indices of two families, so the overlap
    // between the generic indices is <a_1<;1>|a_2<;2>> and stands, whichever
    // factor comes first (Wick itself emits the overlap first)
    const auto d1 = make_kronecker(a1_1, p1), s = make_overlap(p1, q1),
               d2 = make_kronecker(q1, a2_2);
    for (auto const& product :
         {d1 * d2 * s * g(a1_1, a2_2), d1 * s * d2 * g(a1_1, a2_2),
          s * d1 * d2 * g(a1_1, a2_2)}) {
      INFO(toUtf8(serialize(product)));
      ExprPtr r;
      REQUIRE_NOTHROW(r = reduce(product));
      CHECK(tensors_labelled(r, reserved::kronecker_label()).empty());
      const auto metrics = tensors_labelled(r, reserved::overlap_label());
      REQUIRE(metrics.size() == 1);
      const auto gs = tensors_labelled(r, L"g");
      REQUIRE(gs.size() == 1);
      CHECK(instances(*gs[0]) ==
            container::svector<IndexBasis::optional_instance>{1, 2, {}, {}});
      const auto g_slots = slots(*gs[0]);
      const auto m_slots = slots(*metrics[0]);
      CHECK(((m_slots[0] == g_slots[0] && m_slots[1] == g_slots[1]) ||
             (m_slots[0] == g_slots[1] && m_slots[1] == g_slots[0])));
    }
  }

  SECTION("an overlap between two instances is not a Kronecker delta") {
    CHECK_FALSE(is_kronecker_equivalent(a1_1, a2_2));
    CHECK(is_kronecker_equivalent(a1_1, a2_1));
    CHECK(is_kronecker_equivalent(a1, a2_2));
    // f{i_1;a_1<;1>} s{a_1<;1>;a_2<;2>} t{a_2<;2>;i_1}: the overlap stands
    const auto r =
        reduce(ex<Tensor>(L"f", bra{i1}, ket{a1_1}) * make_overlap(a1_1, a2_2) *
               ex<Tensor>(L"t", bra{a2_2}, ket{i1}));
    const auto metrics = tensors_labelled(r, reserved::overlap_label());
    REQUIRE(metrics.size() == 1);
    CHECK(slots(*metrics[0]) == container::svector<Index>{a1_1, a2_2});
  }

  SECTION("an overlap inside a non-orthonormal basis stands") {
    // ã: a basis of a with a general metric; an overlap between two of its
    // functions is the metric, a Kronecker delta is still an identity
    auto isr = mbpt::make_min_sr_spaces();
    isr->add(L"ã", IndexBasis{isr->retrieve_space(L"a"), 1},
             IndexBasisMetric::General);
    auto gctx = get_default_context();
    gctx.set(std::move(isr));
    auto gctx_resetter = set_scoped_default_context(gctx);
    const Index x1(L"ã_1"), x2(L"ã_2"), x3(L"ã_3");
    REQUIRE(x1.basis().metric() == IndexBasisMetric::General);
    auto t = [&](Index const& x) { return ex<Tensor>(L"t", bra{i1}, ket{x}); };
    // int + ext and int + int: the overlap stands
    auto r = reduce(t(x1) * make_overlap(x1, x2));
    CHECK(tensors_labelled(r, reserved::overlap_label()).size() == 1);
    CHECK(slots(*tensors_labelled(r, L"t")[0])[1] == x1);
    r = reduce(t(x1) * make_overlap(x1, x3) *
               ex<Tensor>(L"g", bra{x3}, ket{i2}));
    CHECK(tensors_labelled(r, reserved::overlap_label()).size() == 1);
    // a Kronecker delta between its functions reduces
    r = reduce(t(x1) * make_kronecker(x1, x2));
    CHECK(tensors_labelled(r, reserved::kronecker_label()).empty());
    CHECK(slots(*tensors_labelled(r, L"t")[0])[1] == x2);
    // the same overlap in an orthonormal instance of a is an identity
    r = reduce(t(a1_1) * make_overlap(a1_1, a2_1));
    CHECK(tensors_labelled(r, reserved::overlap_label()).empty());
    CHECK(slots(*tensors_labelled(r, L"t")[0])[1] == a2_1);
  }

  SECTION("a delta or overlap between disjoint spaces is zero in any bases") {
    const auto i1_1 = i1.replace_basis_instance(1);
    const auto a1_2 = a1.replace_basis_instance(2);
    auto t = [&] { return ex<Tensor>(L"t", bra{a1_2}, ket{i1_1}); };
    CHECK(reduce(make_overlap(i1_1, a1_2) * t()) == ex<Constant>(0));
    CHECK(reduce(make_kronecker(i1_1, a1_2) * t()) == ex<Constant>(0));
    // also when the spaces become disjoint through a rule
    CHECK(reduce(make_kronecker(p1, i1_1) * make_overlap(p1, a1_2) * t()) ==
          ex<Constant>(0));
    CHECK(reduce(make_overlap(p1, a1_2) * make_kronecker(p1, i1_1) * t()) ==
          ex<Constant>(0));
  }

  SECTION("a Kronecker delta between two instances is an error") {
    CHECK_THROWS_MATCHES(
        reduce(make_kronecker(a1_1, a2_2) * g(a1_1, a2_2)), Exception,
        message_contains("Kronecker delta between a_1<;1> and a_2<;2>") &&
            message_contains("different basis instances (1 and 2)"));
    // also when the instances meet through a generic index
    CHECK_THROWS_MATCHES(
        reduce(make_kronecker(p1, a1_1) * make_kronecker(p1, a2_2) *
               g(a1_1, a2_2)),
        Exception, message_contains("different basis instances (1 and 2)"));
  }

  SECTION("an external generic index never replaces a specific one") {
    // s{a_1;a_2<;1>} t{a_2<;1>;i_1}, a_1 external: a_2<;1> cannot become a_1,
    // the overlap stands
    {
      const auto r =
          reduce(make_overlap(a1, a2_1) * ex<Tensor>(L"t", bra{a2_1}, ket{i1}));
      const auto metrics = tensors_labelled(r, reserved::overlap_label());
      REQUIRE(metrics.size() == 1);
      const auto ts = tensors_labelled(r, L"t");
      REQUIRE(ts.size() == 1);
      const auto t_slots = slots(*ts[0]);
      CHECK(t_slots[0] != a1);
      CHECK(t_slots[0].basis().basis_instance() == 1);
      CHECK(slots(*metrics[0]) == container::svector<Index>{a1, t_slots[0]});
    }
    // s{a_1<;1>;a_2} t{a_2;i_1}, a_1<;1> external: a_2 becomes a_1<;1>
    {
      const auto r =
          reduce(make_overlap(a1_1, a2) * ex<Tensor>(L"t", bra{a2}, ket{i1}));
      CHECK(tensors_labelled(r, reserved::overlap_label()).empty());
      const auto ts = tensors_labelled(r, L"t");
      REQUIRE(ts.size() == 1);
      CHECK(slots(*ts[0]) == container::svector<Index>{a1_1, i1});
    }
  }
}

TEST_CASE("basis-wick-bernoulli-commutator", "[algorithms][wick][basis]") {
  using namespace sequant;
  using namespace sequant::mbpt;
  using sequant::tests::csv::instance_histogram;
  using sequant::tests::csv::tensors_labelled;
  using sequant::tests::csv::with_leg_instance;

  auto ctx = set_scoped_default_context(
      sequant::Context({.index_basis_registry_shared_ptr = make_min_sr_spaces(),
                        .vacuum = Vacuum::SingleProduct,
                        .spbasis = SPBasis::Spinor}));
  auto mbpt_ctx = set_scoped_default_mbpt_context(mbpt::Context(
      {.csv = CSV::No, .op_registry_ptr = make_minimal_registry()}));
  const auto& isr = get_default_context().index_basis_registry();

  // one commutator [V, σ] of the Bernoulli UCC derivation
  const auto ab = bernoulli::detail::wick_commutator(
      op::tensor::h(2), with_leg_instance(op::tensor::t(2), 3));

  const auto ts = tensors_labelled(ab, L"t");
  REQUIRE(!ts.empty());
  std::size_t misses = 0;
  for (auto const* t : ts)
    for (Index const& idx : t->_slots())
      if (isr->is_pure_unoccupied(idx.space()) &&
          idx.basis().basis_instance() != 3)
        ++misses;
  CHECK(misses == 0);
  CHECK(instance_histogram(ab) ==
        std::map<IndexBasis::optional_instance, std::size_t>{{std::nullopt, 60},
                                                             {3, 32}});
  REQUIRE(ab->is<Sum>());
  CHECK(ab->size() == 8);
}

TEST_CASE("basis-wick-extended", "[algorithms][wick][basis]") {
  using namespace sequant;
  using sequant::tests::csv::standing_metrics;
  using sequant::tests::csv::tensors_labelled;

  auto ctx = get_default_context();
  ctx.set(mbpt::make_mr_spaces());
  ctx.set(Vacuum::MultiProduct);
  auto ctx_resetter = set_scoped_default_context(ctx);

  // a Kronecker delta between two indices in different basis instances
  // identifies functions of different bases, which is meaningless
  auto has_delta_across_instances = [](ExprPtr const& expr) {
    return std::ranges::any_of(
        tensors_labelled(expr, reserved::kronecker_label()), [](auto const* t) {
          return different_instances(t->_bra()[0].basis(),
                                     t->_ket()[0].basis());
        });
  };

  SECTION("a contraction between two instances leaves an overlap standing") {
    for (IndexBasis::instance_type ket_instance : {1, 2}) {
      CAPTURE(ket_instance);
      const auto p1 = Index(L"p_1").replace_basis_instance(1);
      const auto p2 = Index(L"p_2").replace_basis_instance(ket_instance);
      const auto expr = ex<Tensor>(L"h", bra{p1}, ket{p2}, Symmetry::Nonsymm) *
                        ex<FNOperator>(cre{p1}, ann{}) *
                        ex<FNOperator>(cre{}, ann{p2});
      ExprPtr r;
      REQUIRE_NOTHROW(r = FWickTheorem{expr}.compute());
      CHECK(!has_delta_across_instances(r));
      // the core term keeps one overlap between the two bases
      CHECK(standing_metrics(r) == (ket_instance == 1 ? 0 : 1));
    }
  }

  SECTION("η = δ - γ between two instances is an overlap minus γ") {
    for (IndexBasis::instance_type bra_instance : {1, 2}) {
      CAPTURE(bra_instance);
      const auto u1 = Index(L"u_1").replace_basis_instance(1);
      const auto u2 = Index(L"u_2").replace_basis_instance(bra_instance);
      const auto expr =
          ex<FNOperator>(cre{}, ann{u2}) * ex<FNOperator>(cre{u1}, ann{});
      FWickTheorem wick{expr};
      wick.eta_as_delta_minus_gamma(true);
      ExprPtr r;
      REQUIRE_NOTHROW(r = wick.compute());
      CHECK(!has_delta_across_instances(r));
      CHECK(standing_metrics(r) == (bra_instance == 1 ? 0 : 1));
      CHECK(tensors_labelled(r, reserved::rdm_label()).size() == 1);
    }
  }

  SECTION("η = δ - γ in a non-orthonormal basis is an overlap minus γ") {
    // ũ: a basis of the active space with a general metric
    auto isr = mbpt::make_mr_spaces();
    isr->add(L"ũ", IndexBasis{isr->retrieve_space(L"u"), 1},
             IndexBasisMetric::General);
    auto metric_ctx = get_default_context();
    metric_ctx.set(std::move(isr));
    auto metric_resetter = set_scoped_default_context(metric_ctx);
    const Index u1(L"ũ_1"), u2(L"ũ_2");
    REQUIRE(u1.basis().metric() == IndexBasisMetric::General);
    const auto expr =
        ex<FNOperator>(cre{}, ann{u2}) * ex<FNOperator>(cre{u1}, ann{});
    FWickTheorem wick{expr};
    wick.eta_as_delta_minus_gamma(true);
    ExprPtr r;
    REQUIRE_NOTHROW(r = wick.compute());
    CHECK(tensors_labelled(r, reserved::kronecker_label()).empty());
    CHECK(tensors_labelled(r, reserved::overlap_label()).size() == 1);
    CHECK(tensors_labelled(r, reserved::rdm_label()).size() == 1);
  }
}
