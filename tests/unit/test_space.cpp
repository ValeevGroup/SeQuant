//
// Created by Eduard Valeyev on 3/20/18.
//

#include "catch.hpp"

#include "../../src/SeQuant2/space.hpp"

TEST_CASE("IndexSpace", "[elements]") {

  using namespace sequant2;

  SECTION("register_instance") {
    REQUIRE_NOTHROW(IndexSpace::register_instance(L"p", IndexSpace::all));
    REQUIRE_THROWS(IndexSpace::register_instance(L"p_1", IndexSpace::all));
  }

  SECTION("register_standard_instances") {
    REQUIRE_NOTHROW(IndexSpace::register_standard_instances());
    REQUIRE_NOTHROW(IndexSpace::register_standard_instances());  // can do twice
    REQUIRE(IndexSpace::instance_exists(L"i"));
    REQUIRE(IndexSpace::instance_exists(L"m"));
    REQUIRE(IndexSpace::instance_exists(L"p"));
    REQUIRE(IndexSpace::instance_exists(L"a_17"));
    REQUIRE(IndexSpace::instance_exists(L"e_19"));
    REQUIRE(IndexSpace::instance_exists(L"⍺_21"));
    REQUIRE(IndexSpace::instance_exists(L"⍺'_32"));
    REQUIRE(IndexSpace::instance_exists(L"κ_48"));
  }

  SECTION("equality") {
    REQUIRE(IndexSpace::instance(L"i") == IndexSpace::instance(L"i"));
    REQUIRE(IndexSpace::instance(L"i") != IndexSpace::instance(L"p"));

    REQUIRE(IndexSpace::null_instance() == IndexSpace::instance(IndexSpace::null_key()));

    REQUIRE(IndexSpace::instance(L"i").type() == IndexSpace::active_occupied);
    REQUIRE(IndexSpace::instance(L"i") == IndexSpace::active_occupied);
    REQUIRE(IndexSpace::active_occupied == IndexSpace::instance(L"i").type());
    REQUIRE(IndexSpace::active_occupied == IndexSpace::instance(L"i"));
    REQUIRE(IndexSpace::instance(L"i").qns() == IndexSpace::nullqns);
    REQUIRE(IndexSpace::instance(L"i") == IndexSpace::nullqns);
    REQUIRE(IndexSpace::nullqns == IndexSpace::instance(L"i").qns());
    REQUIRE(IndexSpace::nullqns == IndexSpace::instance(L"i"));
  }

  SECTION("ordering") {
    REQUIRE(IndexSpace::instance(L"i") < IndexSpace::instance(L"a"));
    REQUIRE(!(IndexSpace::instance(L"a") < IndexSpace::instance(L"i")));
    REQUIRE(IndexSpace::instance(L"i") < IndexSpace::instance(L"m"));
    REQUIRE(!(IndexSpace::instance(L"m") < IndexSpace::instance(L"m")));
    REQUIRE(IndexSpace::instance(L"m") < IndexSpace::instance(L"a"));
    REQUIRE(IndexSpace::instance(L"m") < IndexSpace::instance(L"p"));
    REQUIRE(IndexSpace::instance(L"a") < IndexSpace::instance(L"p"));
    REQUIRE(IndexSpace::instance(L"p") < IndexSpace::instance(L"⍺"));
  }

  SECTION("set operations") {
    REQUIRE(IndexSpace::instance(L"i") == intersection(IndexSpace::instance(L"i"), IndexSpace::instance(L"p")));
    REQUIRE(IndexSpace::null_instance() == intersection(IndexSpace::instance(L"a"), IndexSpace::instance(L"i")));
    REQUIRE(IndexSpace::null_instance() == intersection(IndexSpace::instance(L"a"), IndexSpace::instance(L"⍺'")));

    REQUIRE(IndexSpace::instance(L"κ") == unIon(IndexSpace::instance(L"m"), IndexSpace::instance(L"⍺")));

    REQUIRE(includes(IndexSpace::instance(L"κ"), IndexSpace::instance(L"m")));
    REQUIRE(!includes(IndexSpace::instance(L"m"), IndexSpace::instance(L"κ")));
    REQUIRE(includes(IndexSpace::instance(L"⍺"), IndexSpace::instance(L"a")));
  }

  SECTION("occupancy_class") {
    REQUIRE(occupancy_class(IndexSpace::instance(L"i")) == -1);
    REQUIRE(occupancy_class(IndexSpace::instance(L"a")) == +1);
    REQUIRE(occupancy_class(IndexSpace::instance(L"p")) == 0);
  }

}
