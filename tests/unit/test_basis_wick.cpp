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
      sequant::Context({.index_space_registry_shared_ptr = make_min_sr_spaces(),
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
    // neither leg is pure, so the contraction is s{q;q'} δ{p_1<;1>;q}
    // δ{q';p_2<;2>} with fresh generic q, q': one overlap between the two
    // families must stand, and nothing may throw
    const Index i1(L"i_1"), i2(L"i_2");
    const auto p1 = Index(L"p_1").replace_basis_instance(1);
    const auto p2 = Index(L"p_2").replace_basis_instance(2);
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
      {.index_space_registry_shared_ptr = mbpt::make_min_sr_spaces(),
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
    const auto unit = IndexSpaceMetric::Unit;
    CHECK_FALSE(is_kronecker_equivalent(a1_1, a2_2, unit));
    CHECK(is_kronecker_equivalent(a1_1, a2_1, unit));
    CHECK(is_kronecker_equivalent(a1, a2_2, unit));
    // f{i_1;a_1<;1>} s{a_1<;1>;a_2<;2>} t{a_2<;2>;i_1}: the overlap stands
    const auto r =
        reduce(ex<Tensor>(L"f", bra{i1}, ket{a1_1}) * make_overlap(a1_1, a2_2) *
               ex<Tensor>(L"t", bra{a2_2}, ket{i1}));
    const auto metrics = tensors_labelled(r, reserved::overlap_label());
    REQUIRE(metrics.size() == 1);
    CHECK(slots(*metrics[0]) == container::svector<Index>{a1_1, a2_2});
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
      sequant::Context({.index_space_registry_shared_ptr = make_min_sr_spaces(),
                        .vacuum = Vacuum::SingleProduct,
                        .spbasis = SPBasis::Spinor}));
  auto mbpt_ctx = set_scoped_default_mbpt_context(mbpt::Context(
      {.csv = CSV::No, .op_registry_ptr = make_minimal_registry()}));
  const auto& isr = get_default_context().index_space_registry();

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

  SECTION("η = δ - γ under a non-unit metric is an overlap minus γ") {
    auto metric_ctx = get_default_context();
    metric_ctx.set(IndexSpaceMetric::General);
    auto metric_resetter = set_scoped_default_context(metric_ctx);
    const Index u1(L"u_1"), u2(L"u_2");
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
