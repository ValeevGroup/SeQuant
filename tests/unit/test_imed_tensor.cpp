#include "catch.hpp"

#include <SeQuant/core/tensor.hpp>
#include <SeQuant/core/tensor_network.cpp>
#include <SeQuant/domain/factorize/imed_tensor.hpp>

using namespace sequant;
using namespace sequant::factorize;

TEST_CASE("ImedTensor Constructors") {
  SECTION("Default constructor") {
    REQUIRE_NOTHROW(ImedTensor());
    REQUIRE_NOTHROW(ex<ImedTensor>());

    ExprPtr imed1 = ex<ImedTensor>();
    REQUIRE(imed1->is<Tensor>());
    REQUIRE_NOTHROW(imed1->as<Tensor>());
    REQUIRE(imed1->is<ImedTensor>());
    // REQUIRE_NOTHROW(imed1->as<ImedTensor>());
  }
}
