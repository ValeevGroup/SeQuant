//
// Created by Eduard Valeyev on 3/23/18.
//

#include "catch.hpp"

#include <iostream>
#include "../../src/SeQuant2/tensor_network.hpp"

TEST_CASE("TensorNetwork", "[elements]") {

  using namespace sequant2;
  IndexSpace::register_standard_instances();

  SECTION("constructors") {

    auto t1 = ex<Tensor>(L"F", WstrList{L"i_1"}, WstrList{L"i_2"});
    auto t2 = ex<Tensor>(L"t", WstrList{L"i_2"}, WstrList{L"i_1"});
    auto t1_x_t2 = t1 * t2;
    REQUIRE_NOTHROW(TensorNetwork(*t1_x_t2));

    auto t1_x_t2_p_t2 = t1 * (t2 + t2); // can only use a flat tensor product
    REQUIRE_THROWS_AS(TensorNetwork(*t1_x_t2_p_t2), std::logic_error);
  }  // SECTION("constructors")

}  // TEST_CASE("Tensor")
