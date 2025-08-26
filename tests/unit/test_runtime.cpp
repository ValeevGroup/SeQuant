//
// Created by Eduard Valeyev on 10/12/22.
//

#include <catch2/catch_test_macros.hpp>

#include "catch2_sequant.hpp"

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>

#include <iostream>

TEST_CASE("context", "[runtime]") {
  using namespace sequant;

  SECTION("constructors") {
    REQUIRE_NOTHROW(Context{});

    REQUIRE_NOTHROW(Context({.vacuum = Vacuum::Invalid}));
    REQUIRE(Context({.vacuum = Vacuum::Invalid}).vacuum() == Vacuum::Invalid);
  }

  SECTION("default context") {
    CHECK_NOTHROW(get_default_context());
    auto initial_ctx = get_default_context();

    // basic set_default_context test
    CHECK_NOTHROW(set_default_context(
        {.index_space_registry_shared_ptr = mbpt::make_sr_spaces(),
         .vacuum = Vacuum::SingleProduct,
         .metric = IndexSpaceMetric::Unit,
         .braket_symmetry = BraKetSymmetry::symm,
         .spbasis = SPBasis::spinfree}));
    CHECK(get_default_context().vacuum() == Vacuum::SingleProduct);
    CHECK(get_default_context().metric() == IndexSpaceMetric::Unit);
    CHECK(get_default_context().braket_symmetry() == BraKetSymmetry::symm);
    CHECK(get_default_context().spbasis() == SPBasis::spinfree);

    // set distinct contexts for fermi and bose statistics
    auto [fermi_isr, bose_isr] = mbpt::make_fermi_and_bose_spaces();
    CHECK(fermi_isr->spaces() ==
          bose_isr->spaces());  // fermi_isr and bose_isr share the space set
    CHECK_NOTHROW(set_default_context(
        {{Statistics::FermiDirac,
          Context({.index_space_registry_shared_ptr = fermi_isr,
                   .vacuum = Vacuum::SingleProduct})},
         {Statistics::BoseEinstein,
          Context({.index_space_registry_shared_ptr = bose_isr,
                   .vacuum = Vacuum::Physical})}}));
    CHECK(get_default_context(Statistics::Arbitrary).vacuum() ==
          Vacuum::SingleProduct);
    CHECK(get_default_context(Statistics::FermiDirac).vacuum() ==
          Vacuum::SingleProduct);
    CHECK(get_default_context(Statistics::BoseEinstein).vacuum() ==
          Vacuum::Physical);

    // reset back to default
    CHECK_NOTHROW(reset_default_context());
    CHECK(get_default_context().vacuum() == Vacuum::Physical);
    CHECK(get_default_context().metric() == IndexSpaceMetric::Unit);
    CHECK(get_default_context().braket_symmetry() == BraKetSymmetry::conjugate);
    CHECK(get_default_context().spbasis() == SPBasis::spinor);

    // reset back to initial context
    CHECK_NOTHROW(set_default_context(initial_ctx));
    CHECK(get_default_context() == initial_ctx);
    CHECK(!(get_default_context() != initial_ctx));

    // scoped changes to default context
    {
      // if we do not save the resetter context is reset back immediately
      CHECK_NOTHROW(set_scoped_default_context(
          {.index_space_registry_shared_ptr = mbpt::make_sr_spaces(),
           .vacuum = Vacuum::SingleProduct,
           .metric = IndexSpaceMetric::Unit,
           .braket_symmetry = BraKetSymmetry::symm,
           .spbasis = SPBasis::spinfree}));
      CHECK(get_default_context() == initial_ctx);

      auto ctx = get_default_context();
      ctx.set(mbpt::make_sr_spaces());
      ctx.set(BraKetSymmetry::symm);
      ctx.set(SPBasis::spinfree);
      const auto ctx_copy = ctx;
      auto resetter = set_scoped_default_context(ctx);
      CHECK(get_default_context() == ctx_copy);
    }
    // leaving scope resets the context back
    CHECK(get_default_context() == initial_ctx);
  }
}
