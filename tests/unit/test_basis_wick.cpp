//
// Wick's theorem on operator products whose legs carry a basis instance.
//

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include "catch2_sequant.hpp"
#include "csv_test_utils.hpp"

#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/index_basis.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/utility/exception.hpp>
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

  SECTION("K5: an ungranted family next to a granted one") {
    const auto products = {t::λ(2) * with_leg_instance(t::t(2), 1),
                           with_leg_instance(t::λ(2), 1) * t::t(2)};
    for (auto const& product : products) {
      if (yes)
        CHECK_THROWS_MATCHES(
            vac_av(product), sequant::Exception,
            message_contains("carries proto indices and no basis instance") &&
                message_contains("a_100<i_100, i_101> carries") &&
                message_contains("a_102<i_100, i_101;1>"));
      else
        CHECK(serialize(vac_av(product)) ==
              L"1/4 t{a_1<;1>,a_2<;1>;i_1,i_2}:A-N-S * "
              L"λ{i_1,i_2;a_1<;1>,a_2<;1>}:A-N-S");
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
