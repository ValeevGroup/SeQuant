//
// Created by Eduard Valeyev on 3/20/18.
//

#include "catch.hpp"

#include <SeQuant/core/space.hpp>

TEST_CASE("IndexSpace", "[elements]") {
  using namespace sequant;

  SECTION("register_instance") {
    REQUIRE_THROWS(IndexSpace::register_instance(
        L"p", IndexSpace::all));  // already registered standard instances
    REQUIRE_THROWS(
        IndexSpace::register_instance(L"p_1", IndexSpace::all));  // ditto
    // since convention does not include ALL standard spaces
    // (see https://github.com/ValeevGroup/SeQuant/issues/97 )
    // user can register additional instances
    REQUIRE_NOTHROW(
        IndexSpace::register_instance(L"m'", IndexSpace::frozen_occupied));
    REQUIRE(IndexSpace::instance_exists(L"m'"));
  }

  // these are loaded in test_main.cpp
  SECTION("register_standard_instances") {
    REQUIRE(IndexSpace::instance_exists(L"i"));
    REQUIRE(IndexSpace::instance_exists(L"m"));
    REQUIRE(IndexSpace::instance_exists(L"p"));
    REQUIRE(IndexSpace::instance_exists(L"a_17"));
    REQUIRE(IndexSpace::instance_exists(L"e_19"));
    REQUIRE(IndexSpace::instance_exists(L"α_21"));
    REQUIRE(IndexSpace::instance_exists(L"α'_32"));
    REQUIRE(IndexSpace::instance_exists(L"κ_48"));
  }

  SECTION("register_key") {
    REQUIRE_NOTHROW(IndexSpace::register_key(
        L"w",
        IndexSpace::all));  // can assign additional key to a space already
                            // registered, this does not redefine base key
    REQUIRE(IndexSpace::instance(L"w") == IndexSpace::instance(L"p"));
  }

  SECTION("equality") {
    REQUIRE(IndexSpace::instance(L"i") == IndexSpace::instance(L"i"));
    REQUIRE(IndexSpace::instance(L"i") != IndexSpace::instance(L"p"));

    REQUIRE(IndexSpace::null_instance() ==
            IndexSpace::instance(IndexSpace::null_key()));

    REQUIRE(IndexSpace::instance(L"i").type() == IndexSpace::active_occupied);
    REQUIRE(IndexSpace::instance(L"i") == IndexSpace::active_occupied);
    REQUIRE(IndexSpace::active_occupied == IndexSpace::instance(L"i").type());
    REQUIRE(IndexSpace::active_occupied == IndexSpace::instance(L"i"));
    REQUIRE(IndexSpace::instance(L"i").qns() == IndexSpace::nullqns);
    REQUIRE(IndexSpace::instance(L"i") == IndexSpace::nullqns);
    REQUIRE(IndexSpace::nullqns == IndexSpace::instance(L"i").qns());
    REQUIRE(IndexSpace::nullqns == IndexSpace::instance(L"i"));

    // use nondefault mask
    TypeAttr::used_bits = 0b100;
    REQUIRE(IndexSpace::active_occupied == IndexSpace::occupied);
    REQUIRE(IndexSpace::active_occupied != IndexSpace::inactive_occupied);
    REQUIRE(IndexSpace::active_occupied != IndexSpace::frozen_occupied);
    REQUIRE(IndexSpace::active_occupied == IndexSpace::all);
    REQUIRE(IndexSpace::active_occupied == IndexSpace::all_active);
    TypeAttr::used_bits = 0xffff;
  }

  SECTION("ordering") {
    REQUIRE(!(IndexSpace::instance(L"i") < IndexSpace::instance(L"i")));
    REQUIRE(IndexSpace::instance(L"i") < IndexSpace::instance(L"a"));
    REQUIRE(!(IndexSpace::instance(L"a") < IndexSpace::instance(L"i")));
    REQUIRE(IndexSpace::instance(L"i") < IndexSpace::instance(L"m"));
    REQUIRE(!(IndexSpace::instance(L"m") < IndexSpace::instance(L"m")));
    REQUIRE(IndexSpace::instance(L"m") < IndexSpace::instance(L"a"));
    REQUIRE(IndexSpace::instance(L"m") < IndexSpace::instance(L"p"));
    REQUIRE(IndexSpace::instance(L"a") < IndexSpace::instance(L"p"));
    REQUIRE(IndexSpace::instance(L"p") < IndexSpace::instance(L"α"));

    // test ordering with quantum numbers
    {
      auto i = IndexSpace::instance(L"i");
      [[maybe_unused]] auto a = IndexSpace::instance(L"a");
      auto iA =
          IndexSpace::instance(IndexSpace::active_occupied, IndexSpace::alpha);
      auto iB =
          IndexSpace::instance(IndexSpace::active_occupied, IndexSpace::beta);
      auto aA = IndexSpace::instance(IndexSpace::active_unoccupied,
                                     IndexSpace::alpha);
      auto aB =
          IndexSpace::instance(IndexSpace::active_unoccupied, IndexSpace::beta);

      REQUIRE(iA < aA);
      REQUIRE(iB < aB);
      REQUIRE(iA < aB);
      REQUIRE(iA < aB);
      REQUIRE(iA < iB);
      REQUIRE(aA < aB);
      REQUIRE(!(iA < iA));
      REQUIRE(i < iA);
      REQUIRE(i < iB);
    }

    // use nondefault mask
    TypeAttr::used_bits = 0b100;
    REQUIRE(!(IndexSpace::active_occupied < IndexSpace::occupied));
    REQUIRE(!(IndexSpace::inactive_occupied < IndexSpace::frozen_occupied));
    REQUIRE(IndexSpace::inactive_occupied < IndexSpace::active_occupied);
    REQUIRE(!(IndexSpace::active_occupied < IndexSpace::all));
    REQUIRE(!(IndexSpace::active_occupied < IndexSpace::all_active));
    TypeAttr::used_bits = 0xffff;
  }

  SECTION("set operations") {
    REQUIRE(
        IndexSpace::instance(L"i") ==
        intersection(IndexSpace::instance(L"i"), IndexSpace::instance(L"p")));
    REQUIRE(
        IndexSpace::null_instance() ==
        intersection(IndexSpace::instance(L"a"), IndexSpace::instance(L"i")));
    REQUIRE(
        IndexSpace::null_instance() ==
        intersection(IndexSpace::instance(L"a"), IndexSpace::instance(L"α'")));

    REQUIRE(IndexSpace::instance(L"κ") ==
            unIon(IndexSpace::instance(L"p"), IndexSpace::instance(L"α'")));

    REQUIRE(includes(IndexSpace::instance(L"κ"), IndexSpace::instance(L"m")));
    REQUIRE(!includes(IndexSpace::instance(L"m"), IndexSpace::instance(L"κ")));
    REQUIRE(includes(IndexSpace::instance(L"α"), IndexSpace::instance(L"a")));
  }

  SECTION("occupancy_class") {
    REQUIRE(occupancy_class(IndexSpace::instance(L"i")) == -1);
    REQUIRE(occupancy_class(IndexSpace::instance(L"a")) == +1);
    REQUIRE(occupancy_class(IndexSpace::instance(L"p")) == 0);
    REQUIRE(occupancy_class(IndexSpace::instance(L"u")) == +1);
  }
}
