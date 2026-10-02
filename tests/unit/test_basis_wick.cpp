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

#include <cstddef>
#include <map>
#include <optional>

TEST_CASE("basis-wick-contractions", "[algorithms][wick][basis]") {
  using namespace sequant;
  using namespace sequant::mbpt;
  using sequant::tests::csv::message_contains;
  using sequant::tests::csv::ScopedDefaultCardinalLabels;
  using sequant::tests::csv::standing_metrics;
  using sequant::tests::csv::with_leg_instance;
  namespace t = op::tensor;

  const auto csv = GENERATE(CSV::Yes, CSV::No);
  const ScopedDefaultCardinalLabels cardinal_labels;
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
