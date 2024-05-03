//
// Created by Eduard Valeyev on 10/12/22.
//

#include <catch2/catch_test_macros.hpp>

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>

#include <iostream>

TEST_CASE("Context", "[runtime]") {
  using namespace sequant;

  CHECK_NOTHROW(get_default_context());
  auto initial_ctx = get_default_context();

  // basic set_default_context test
  CHECK_NOTHROW(set_default_context(Context(
      Vacuum::SingleProduct, mbpt::make_sr_subspaces(), IndexSpaceMetric::Unit,
      BraKetSymmetry::symm, SPBasis::spinfree)));
  CHECK(get_default_context().vacuum() == Vacuum::SingleProduct);
  CHECK(get_default_context().metric() == IndexSpaceMetric::Unit);
  CHECK(get_default_context().braket_symmetry() == BraKetSymmetry::symm);
  CHECK(get_default_context().spbasis() == SPBasis::spinfree);

  // default default context
  CHECK_NOTHROW(reset_default_context());
  CHECK(get_default_context().vacuum() == Vacuum::Physical);
  CHECK(get_default_context().metric() == IndexSpaceMetric::Unit);
  CHECK(get_default_context().braket_symmetry() == BraKetSymmetry::conjugate);
  CHECK(get_default_context().spbasis() == SPBasis::spinorbital);

  // reset back to initial context
  CHECK_NOTHROW(set_default_context(initial_ctx));
  CHECK(get_default_context() == initial_ctx);
  CHECK(!(get_default_context() != initial_ctx));

  // scoped changes to default context
  {
    // if we do not save the resetter context is reset back immediately
    CHECK_NOTHROW(set_scoped_default_context(Context(
        Vacuum::SingleProduct, mbpt::make_sr_subspaces(),
        IndexSpaceMetric::Unit, BraKetSymmetry::symm, SPBasis::spinfree)));
    CHECK(get_default_context() == initial_ctx);

    auto resetter = set_scoped_default_context(Context(
        Vacuum::SingleProduct, mbpt::make_sr_subspaces(),
        IndexSpaceMetric::Unit, BraKetSymmetry::symm, SPBasis::spinfree));
    CHECK(get_default_context() ==
          Context(Vacuum::SingleProduct, mbpt::make_sr_subspaces(),
                  IndexSpaceMetric::Unit, BraKetSymmetry::symm,
                  SPBasis::spinfree));
  }
  // leaving scope resets the context back
  CHECK(get_default_context() == initial_ctx);

}  // TEST_CASE("Index")
