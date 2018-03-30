//
// Created by Eduard Valeyev on 3/23/18.
//

#include <iostream>

#include "catch.hpp"
#include "timer.hpp"
#include "../../src/SeQuant2/wick.hpp"

#define SEQUANT2_SKIP_LONG_TESTS 1

#if 1
TEST_CASE("WickTheorem", "[algorithms]") {

  using namespace sequant2;
  IndexSpace::register_standard_instances();

  SECTION("Op contractions") {

    REQUIRE(  FWickTheorem::can_contract(fann(L"i_1"),fcre(L"i_2"), Vacuum::Physical) );
    REQUIRE( !FWickTheorem::can_contract(fcre(L"i_1"),fcre(L"i_2"), Vacuum::Physical) );
    REQUIRE( !FWickTheorem::can_contract(fcre(L"i_1"),fann(L"i_2"), Vacuum::Physical) );
    REQUIRE( !FWickTheorem::can_contract(fann(L"i_1"),fann(L"i_2"), Vacuum::Physical) );

    REQUIRE( !FWickTheorem::can_contract(fann(L"i_1"),fcre(L"i_2"), Vacuum::SingleProduct) );
    REQUIRE( !FWickTheorem::can_contract(fcre(L"i_1"),fcre(L"i_2"), Vacuum::SingleProduct) );
    REQUIRE(  FWickTheorem::can_contract(fcre(L"i_1"),fann(L"i_2"), Vacuum::SingleProduct) );
    REQUIRE( !FWickTheorem::can_contract(fann(L"i_1"),fann(L"i_2"), Vacuum::SingleProduct) );

    REQUIRE( !FWickTheorem::can_contract(fann(L"i_1"),fcre(L"a_2"), Vacuum::SingleProduct) );
    REQUIRE( !FWickTheorem::can_contract(fcre(L"i_1"),fcre(L"a_2"), Vacuum::SingleProduct) );
    REQUIRE( !FWickTheorem::can_contract(fcre(L"i_1"),fann(L"a_2"), Vacuum::SingleProduct) );
    REQUIRE( !FWickTheorem::can_contract(fann(L"i_1"),fann(L"a_2"), Vacuum::SingleProduct) );

    REQUIRE( !FWickTheorem::can_contract(fann(L"a_1"),fcre(L"i_2"), Vacuum::SingleProduct) );
    REQUIRE( !FWickTheorem::can_contract(fcre(L"a_1"),fcre(L"i_2"), Vacuum::SingleProduct) );
    REQUIRE( !FWickTheorem::can_contract(fcre(L"a_1"),fann(L"i_2"), Vacuum::SingleProduct) );
    REQUIRE( !FWickTheorem::can_contract(fann(L"a_1"),fann(L"i_2"), Vacuum::SingleProduct) );

    REQUIRE(  FWickTheorem::can_contract(fann(L"a_1"),fcre(L"a_2"), Vacuum::SingleProduct) );
    REQUIRE( !FWickTheorem::can_contract(fcre(L"a_1"),fcre(L"a_2"), Vacuum::SingleProduct) );
    REQUIRE( !FWickTheorem::can_contract(fcre(L"a_1"),fann(L"a_2"), Vacuum::SingleProduct) );
    REQUIRE( !FWickTheorem::can_contract(fann(L"a_1"),fann(L"a_2"), Vacuum::SingleProduct) );

    REQUIRE( !FWickTheorem::can_contract(fann(L"p_1"),fcre(L"i_2"), Vacuum::SingleProduct) );
    REQUIRE( !FWickTheorem::can_contract(fcre(L"p_1"),fcre(L"i_2"), Vacuum::SingleProduct) );
    REQUIRE(  FWickTheorem::can_contract(fcre(L"p_1"),fann(L"i_2"), Vacuum::SingleProduct) );
    REQUIRE( !FWickTheorem::can_contract(fann(L"p_1"),fann(L"i_2"), Vacuum::SingleProduct) );
  }

  SECTION("constructors") {

    REQUIRE_NOTHROW(FWickTheorem{FNOperatorSeq{}});

    {
      auto opseq1 =
          FNOperatorSeq({FNOperator({L"i_1"}, {L"i_2"}), FNOperator({L"i_3"}, {L"i_4"}),
                         FNOperator({L"i_5"}, {L"i_6"})});
      REQUIRE_NOTHROW(FWickTheorem{opseq1});
      auto wick1 = FWickTheorem{opseq1};
      REQUIRE_THROWS(wick1.full_contractions(false).compute());
      REQUIRE_THROWS(wick1.spinfree(true).compute());
    }

  }  // SECTION("constructors")

  SECTION("physical vacuum") {

    constexpr Vacuum V = Vacuum::Physical;

    {
      auto opseq1 =
          FNOperatorSeq({FNOperator({L"i_1"}, {L"i_2"}, V), FNOperator({L"i_3"}, {L"i_4"}, V),
                         FNOperator({L"i_5"}, {L"i_6"}, V)});
      auto wick1 = FWickTheorem{opseq1};
      REQUIRE_NOTHROW(wick1.full_contractions(true).spinfree(false).compute());
      REQUIRE(FWickTheorem{opseq1}.full_contractions(true).spinfree(false).compute().size() == 0);
    }

    // two 2-body operators
    {
      auto opseq =
          FNOperatorSeq({FNOperator({}, {L"i_1", L"i_2"}, V), FNOperator({L"i_3", L"i_4"}, {}, V)});
      auto wick = FWickTheorem{opseq};
      REQUIRE_NOTHROW(wick.full_contractions(true).spinfree(false).compute());
      auto result = wick.full_contractions(true).spinfree(false).compute();
      REQUIRE(result.size() == 2);
    }

    // two 3-body operators
    {
      auto opseq =
          FNOperatorSeq({FNOperator({}, {L"i_1", L"i_2", L"i_3"}, V), FNOperator({L"i_4", L"i_5", L"i_6"}, {}, V)});
      auto wick = FWickTheorem{opseq};
      REQUIRE_NOTHROW(wick.full_contractions(true).spinfree(false).compute());
      auto result = wick.full_contractions(true).spinfree(false).compute();
      REQUIRE(result.size() == 6);
    }

    // two 4-body operators
    {
      auto opseq =
          FNOperatorSeq({FNOperator({}, {L"i_1", L"i_2", L"i_3", L"i_4"}, V), FNOperator({L"i_5", L"i_6", L"i_7", L"i_8"}, {}, V)});
      auto wick = FWickTheorem{opseq};
      REQUIRE_NOTHROW(wick.full_contractions(true).spinfree(false).compute());
      auto result = wick.full_contractions(true).spinfree(false).compute();
      REQUIRE(result.size() == 24);
    }

    // three 2-body operators
    {
      auto opseq =
          FNOperatorSeq({FNOperator({}, {L"i_1"}, V),
                         FNOperator({L"i_2"}, {L"i_3"}, V),
                         FNOperator({L"i_4"}, {}, V)
                        });
      auto wick = FWickTheorem{opseq};
      REQUIRE_NOTHROW(wick.full_contractions(true).spinfree(false).compute());
      auto result = wick.full_contractions(true).spinfree(false).compute();
      REQUIRE(result.size() == 1);
    }


  }  // SECTION("physical vacuum")

  SECTION("fermi vacuum") {

    constexpr Vacuum V = Vacuum::SingleProduct;

    // two (pure qp) 1-body operators
    {
      auto opseq =
          FNOperatorSeq({FNOperator({L"i_1"}, {L"a_1"}, V), FNOperator({L"a_2"}, {L"i_2"}, V)});
      auto wick = FWickTheorem{opseq};
      REQUIRE_NOTHROW(wick.full_contractions(true).spinfree(false).compute());
      auto result = wick.full_contractions(true).spinfree(false).compute();
      REQUIRE(result.size() == 1);
    }

    // two general 1-body operators
    {
      auto opseq =
          FNOperatorSeq({FNOperator({L"p_1"}, {L"p_2"}, V), FNOperator({L"p_3"}, {L"p_4"}, V)});
      auto wick = FWickTheorem{opseq};
      REQUIRE_NOTHROW(wick.full_contractions(true).spinfree(false).compute());
      auto result = wick.full_contractions(true).spinfree(false).compute();
      REQUIRE(result.size() == 1);
    }

    // two (pure qp) 2-body operators
    {
      auto opseq =
          FNOperatorSeq({FNOperator({L"i_1", L"i_2"}, {L"a_1", L"a_2"}, V), FNOperator({L"a_3", L"a_4"}, {L"i_3", L"i_4"}, V)});
      auto wick = FWickTheorem{opseq};
      REQUIRE_NOTHROW(wick.full_contractions(true).spinfree(false).compute());
      auto result = wick.full_contractions(true).spinfree(false).compute();
      REQUIRE(result.size() == 4);
    }

    // two general 2-body operators
    {
      auto opseq =
          FNOperatorSeq({FNOperator({L"p_1", L"p_2"}, {L"p_3", L"p_4"}, V),
                         FNOperator({L"p_5", L"p_6"}, {L"p_7", L"p_8"}, V)});
      auto wick = FWickTheorem{opseq};
      REQUIRE_NOTHROW(wick.full_contractions(true).spinfree(false).compute());
      auto result = wick.full_contractions(true).spinfree(false).compute();
      REQUIRE(result.size() == 4);
    }

    // two general 3-body operators
    {
      auto opseq =
          FNOperatorSeq({FNOperator({L"p_1", L"p_2", L"p_3"}, {L"p_4", L"p_5", L"p_6"}, V),
                         FNOperator({L"p_7", L"p_8", L"p_9"}, {L"p_10", L"p_11", L"p_12"}, V)});
      auto wick = FWickTheorem{opseq};
      REQUIRE_NOTHROW(wick.full_contractions(true).spinfree(false).compute());
      auto result = wick.full_contractions(true).spinfree(false).compute();
      REQUIRE(result.size() == 36);
    }

    // two N-nonconserving operators
    {
      auto opseq =
          FNOperatorSeq({FNOperator({L"p_1", L"p_2", L"p_3"}, {L"p_4", L"p_5"}, V),
                         FNOperator({L"p_7", L"p_8"}, {L"p_10", L"p_11", L"p_12"}, V)});
      auto wick = FWickTheorem{opseq};
      REQUIRE_NOTHROW(wick.full_contractions(true).spinfree(false).compute());
      auto result = wick.full_contractions(true).spinfree(false).compute();
      REQUIRE(result.size() == 12);
    }

    // odd number of ops -> full contraction is 0
    {
      auto opseq =
          FNOperatorSeq({FNOperator({L"p_1", L"p_2"}, {L"p_4", L"p_5"}, V),
                         FNOperator({L"p_7", L"p_8"}, {L"p_10", L"p_11", L"p_12"}, V)});
      auto wick = FWickTheorem{opseq};
      REQUIRE_NOTHROW(wick.full_contractions(true).spinfree(false).compute());
      auto result = wick.full_contractions(true).spinfree(false).compute();
      REQUIRE(result.size() == 0);
    }

    // 4-body ^ 4-body
    SEQUANT2_PROFILE_SINGLE("wick(4^4)", {
      auto opseq =
          FNOperatorSeq({FNOperator({L"p_1", L"p_2", L"p_3", L"p_4"}, {L"p_5", L"p_6", L"p_7", L"p_8"}, V),
                         FNOperator({L"p_21", L"p_22", L"p_23", L"p_24"}, {L"p_25", L"p_26", L"p_27", L"p_28"}, V)
                        });
      auto wick = FWickTheorem{opseq};
      auto result = wick.full_contractions(true).spinfree(false).compute(true);
      REQUIRE(result.size() == 576);
    }
    )

    // three general 1-body operators
    {
      auto opseq =
          FNOperatorSeq({FNOperator({L"p_1"}, {L"p_2"}, V),
                         FNOperator({L"p_3"}, {L"p_4"}, V),
                         FNOperator({L"p_5"}, {L"p_6"}, V)
                        });
      auto wick = FWickTheorem{opseq};
      REQUIRE_NOTHROW(wick.full_contractions(true).spinfree(false).compute());
      auto result = wick.full_contractions(true).spinfree(false).compute();
      REQUIRE(result.size() == 2);
    }

    // 4-body ^ 2-body ^ 2-body
    {
      auto opseq =
          FNOperatorSeq({FNOperator({L"p_1", L"p_2", L"p_3", L"p_4"}, {L"p_5", L"p_6", L"p_7", L"p_8"}, V),
                         FNOperator({L"p_9", L"p_10"}, {L"p_11", L"p_12"}, V),
                         FNOperator({L"p_13", L"p_14"}, {L"p_15", L"p_16"}, V)
                        });
      auto wick = FWickTheorem{opseq};
      auto result = wick.full_contractions(true).spinfree(false).compute();
      REQUIRE(result.size() == 576);
    }

    // 2-body ^ 2-body ^ 2-body
    {
      auto opseq =
          FNOperatorSeq({FNOperator({L"p_1", L"p_2"}, {L"p_5", L"p_6"}, V),
                         FNOperator({L"p_9", L"p_10"}, {L"p_11", L"p_12"}, V),
                         FNOperator({L"p_17", L"p_18"}, {L"p_19", L"p_20"}, V)
                        });
      auto wick = FWickTheorem{opseq};
      auto result = wick.full_contractions(true).spinfree(false).compute();
      REQUIRE(result.size() == 80);
    }

    // 2-body ^ 2-body ^ 2-body ^ 2-body
    SEQUANT2_PROFILE_SINGLE("wick(2^2^2^2)", {
      auto opseq =
          FNOperatorSeq({FNOperator({L"p_1", L"p_2"}, {L"p_5", L"p_6"}, V),
                         FNOperator({L"p_9", L"p_10"}, {L"p_11", L"p_12"}, V),
                         FNOperator({L"p_13", L"p_14"}, {L"p_15", L"p_16"}, V),
                         FNOperator({L"p_17", L"p_18"}, {L"p_19", L"p_20"}, V)
                        });
      auto wick = FWickTheorem{opseq};
      auto result = wick.full_contractions(true).spinfree(false).compute(true);
      REQUIRE(result.size() == 4752);
    }
    )

    // 4-body ^ 2-body ^ 2-body ^ 2-body
    SEQUANT2_PROFILE_SINGLE("wick(4^2^2^2)", {
      auto opseq =
          FNOperatorSeq({FNOperator({L"p_1", L"p_2", L"p_3", L"p_4"}, {L"p_5", L"p_6", L"p_7", L"p_8"}, V),
                         FNOperator({L"p_9", L"p_10"}, {L"p_11", L"p_12"}, V),
                         FNOperator({L"p_13", L"p_14"}, {L"p_15", L"p_16"}, V),
                         FNOperator({L"p_17", L"p_18"}, {L"p_19", L"p_20"}, V)
                        });
      auto wick = FWickTheorem{opseq};
      auto result = wick.full_contractions(true).spinfree(false).compute(true);
      REQUIRE(result.size() == 117504);
    }
    )

    // 3-body ^ 2-body ^ 2-body ^ 3-body
    SEQUANT2_PROFILE_SINGLE("wick(3^2^2^3)", {
      auto opseq =
          FNOperatorSeq({FNOperator({L"p_1", L"p_2", L"p_3"}, {L"p_5", L"p_6", L"p_7"}, V),
                         FNOperator({L"p_9", L"p_10"}, {L"p_11", L"p_12"}, V),
                         FNOperator({L"p_13", L"p_14"}, {L"p_15", L"p_16"}, V),
                         FNOperator({L"p_17", L"p_18", L"p_19"}, {L"p_20", L"p_21", L"p_22"}, V)
                        });
      auto wick = FWickTheorem{opseq};
      auto result = wick.full_contractions(true).spinfree(false).compute(true);
      REQUIRE(result.size() == 202320);
    }
    )

    // 4-body ^ 2-body ^ 4-body
    SEQUANT2_PROFILE_SINGLE("wick(4^2^4)", {
      auto opseq =
          FNOperatorSeq({FNOperator({L"p_1", L"p_2", L"p_3", L"p_4"}, {L"p_5", L"p_6", L"p_7", L"p_8"}, V),
                         FNOperator({L"p_9", L"p_10"}, {L"p_11", L"p_12"}, V),
                         FNOperator({L"p_21", L"p_22", L"p_23", L"p_24"}, {L"p_25", L"p_26", L"p_27", L"p_28"}, V)
                        });
      auto wick = FWickTheorem{opseq};
      auto result = wick.full_contractions(true).spinfree(false).compute(true);
      REQUIRE(result.size() == 50688);
    }
    )

    // 4-body ^ 4-body ^ 4-body
#ifndef SEQUANT2_SKIP_LONG_TESTS
    SEQUANT2_PROFILE_SINGLE("wick(4^4^4)", {
      auto opseq =
          FNOperatorSeq({FNOperator({L"p_1", L"p_2", L"p_3", L"p_4"}, {L"p_5", L"p_6", L"p_7", L"p_8"}, V),
                         FNOperator({L"p_11", L"p_12", L"p_13", L"p_14"}, {L"p_15", L"p_16", L"p_17", L"p_18"}, V),
                         FNOperator({L"p_21", L"p_22", L"p_23", L"p_24"}, {L"p_25", L"p_26", L"p_27", L"p_28"}, V)
                        });
      auto wick = FWickTheorem{opseq};
      auto result = wick.full_contractions(true).spinfree(false).compute(true);
      REQUIRE(result.size() == 4783104);
    }
    )
#endif

#if 0
    // impossible: 4-body ^ 4-body ^ 4-body ^ 4-body ^ 4-body ^ 4-body
    {
      auto opseq =
          FNOperatorSeq({FNOperator({L"p_1", L"p_2", L"p_3", L"p_4"}, {L"p_5", L"p_6", L"p_7", L"p_8"}, V),
                         FNOperator({L"p_11", L"p_12", L"p_13", L"p_14"}, {L"p_15", L"p_16", L"p_17", L"p_18"}, V),
                         FNOperator({L"p_21", L"p_22", L"p_23", L"p_24"}, {L"p_25", L"p_26", L"p_27", L"p_28"}, V),
                         FNOperator({L"p_31", L"p_32", L"p_33", L"p_34"}, {L"p_35", L"p_36", L"p_37", L"p_38"}, V),
                         FNOperator({L"p_41", L"p_42", L"p_43", L"p_44"}, {L"p_45", L"p_46", L"p_47", L"p_48"}, V),
                         FNOperator({L"p_51", L"p_52", L"p_53", L"p_54"}, {L"p_55", L"p_56", L"p_57", L"p_58"}, V)
                        });
      auto wick = FWickTheorem{opseq};
      auto result = wick.full_contractions(true).spinfree(false).compute(true);
      REQUIRE(result.size() == 0);
    }
#endif

  }

  }  // TEST_CASE("WickTheorem")
#endif
