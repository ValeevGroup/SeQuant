//
// Created by Eduard Valeyev on 3/20/18.
//

#define CATCH_CONFIG_RUNNER
#include <catch2/catch_session.hpp>
#include <catch2/catch_test_macros.hpp>

#include "catch2_sequant.hpp"

#include <SeQuant/core/logger.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/runtime.hpp>
#include <SeQuant/core/space.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/domain/mbpt/context.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/op.hpp>

#ifdef SEQUANT_HAS_TILEDARRAY
#include <tiledarray.h>
#endif

#include <clocale>

int main(int argc, char* argv[]) {
  using namespace std;
  using namespace sequant;

  Catch::Session session;

  // global setup...
  std::wcout.precision(std::numeric_limits<double>::max_digits10);
  std::wcerr.precision(std::numeric_limits<double>::max_digits10);
  sequant::set_locale();
  sequant::detail::OpIdRegistrar op_id_registrar;
  sequant::set_default_context(
      {.index_space_registry_shared_ptr = sequant::mbpt::make_sr_spaces(),
       .vacuum = Vacuum::SingleProduct,
       .metric = IndexSpaceMetric::Unit,
       .braket_symmetry = BraKetSymmetry::conjugate,
       .spbasis = SPBasis::spinor,
       .first_dummy_index_ordinal = 100,
       .braket_typesetting = BraKetTypesetting::ContraSub,
       // to_latex() reference outputs predominantly assume the original
       // (naive) convention
       .braket_slot_typesetting = BraKetSlotTypesetting::Naive});
  TensorCanonicalizer::set_cardinal_tensor_labels(
      sequant::mbpt::cardinal_tensor_labels());
  // uncomment to enable verbose output ...
  // Logger::set_instance(1);
  // ... or can instead selectively set/unset particular logging flags
  // Logger::instance().wick_contract = true;

#ifdef SEQUANT_HAS_TILEDARRAY
  auto& world = TA::initialize(argc, argv);
  TA::set_default_world(world);
#endif

  int result = session.run(argc, argv);

#ifdef SEQUANT_HAS_TILEDARRAY
  TA::finalize();
#endif

  return result;
}
