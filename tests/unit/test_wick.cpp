//
// Created by Eduard Valeyev on 3/23/18.
//

#include "catch.hpp"

#include <iostream>
#include "../../src/SeQuant2/wick.hpp"

TEST_CASE("WickTheorem", "[algorithms]") {

  using namespace sequant2;

  SECTION("constructors") {

    REQUIRE_NOTHROW(FWickTheorem{FNOperatorSeq{}});

    auto opseq1 =
        FNOperatorSeq({FNOperator({L"i_1"}, {L"i_2"}), FNOperator({L"i_3"}, {L"i_4"}), FNOperator({L"i_5"}, {L"i_6"})});
    REQUIRE_NOTHROW(FWickTheorem{opseq1});
    auto wick1 = FWickTheorem{opseq1};
    REQUIRE_THROWS(wick1.full_contractions(false).compute());
    REQUIRE_THROWS(wick1.spinfree(true).compute());
    REQUIRE_NOTHROW(wick1.full_contractions(true).spinfree(false).compute());

  }  // SECTION("constructors")

}  // TEST_CASE("WickTheorem")