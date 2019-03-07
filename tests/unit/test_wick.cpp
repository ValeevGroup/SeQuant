//
// Created by Eduard Valeyev on 3/23/18.
//

#include <iostream>

#include "catch.hpp"
#include "timer.hpp"
#include "../../src/SeQuant/wick.hpp"

// comment out if want to run long tests
//#define SEQUANT_SKIP_LONG_TESTS 1

namespace sequant {
struct WickAccessor {};

template <>
template <>
struct WickTheorem<Statistics::FermiDirac>::template access_by<WickAccessor> {
  auto compute_nontensor_wick(WickTheorem<Statistics::FermiDirac>& wick) {
    return wick.compute_nontensor_wick(false);
  }
};

auto compute_nontensor_wick(WickTheorem<Statistics::FermiDirac>& wick) {
  return WickTheorem<Statistics::FermiDirac>::access_by<WickAccessor>{}
      .compute_nontensor_wick(wick);
}

}  // namespace sequant

#if 1
TEST_CASE("WickTheorem", "[algorithms]") {

  using namespace sequant;

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
      REQUIRE(! FWickTheorem{opseq1}.full_contractions(true).spinfree(false).compute() );
    }

    // two 2-body operators
    {
      auto opseq =
          FNOperatorSeq({FNOperator({}, {L"i_1", L"i_2"}, V), FNOperator({L"i_3", L"i_4"}, {}, V)});
      auto wick = FWickTheorem{opseq};
      REQUIRE_NOTHROW(wick.full_contractions(true).spinfree(false).compute());
      auto result = wick.full_contractions(true).spinfree(false).compute();
      REQUIRE(result->size() == 2);
    }

    // two 3-body operators
    {
      auto opseq =
          FNOperatorSeq({FNOperator({}, {L"i_1", L"i_2", L"i_3"}, V), FNOperator({L"i_4", L"i_5", L"i_6"}, {}, V)});
      auto wick = FWickTheorem{opseq};
      REQUIRE_NOTHROW(wick.full_contractions(true).spinfree(false).compute());
      auto result = wick.full_contractions(true).spinfree(false).compute();
      REQUIRE(result->size() == 6);
    }

    // two 4-body operators
    {
      auto opseq =
          FNOperatorSeq({FNOperator({}, {L"i_1", L"i_2", L"i_3", L"i_4"}, V), FNOperator({L"i_5", L"i_6", L"i_7", L"i_8"}, {}, V)});
      auto wick = FWickTheorem{opseq};
      REQUIRE_NOTHROW(wick.full_contractions(true).spinfree(false).compute());
      auto result = wick.full_contractions(true).spinfree(false).compute();
      REQUIRE(result->size() == 24);
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
      REQUIRE(result->size() == 2);  // 1 term = product of 2 terms
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
      REQUIRE(result->size() == 2);  // 1 term = product of 2 terms
    }

    // two general 1-body operators
    {
      auto opseq =
          FNOperatorSeq({FNOperator({L"p_1"}, {L"p_2"}, V), FNOperator({L"p_3"}, {L"p_4"}, V)});
      auto wick = FWickTheorem{opseq};
      REQUIRE_NOTHROW(wick.full_contractions(true).spinfree(false).compute());
      auto result = wick.full_contractions(true).spinfree(false).compute();
      REQUIRE(result->size() == 2
          * 2);  // 1 product of 4 terms (since each contraction of 2 *general* indices produces 2 overlaps)
      REQUIRE(to_latex(result)
                  == L"{{S^{{p_1}}_{{m_{102}}}}{S^{{m_{102}}}_{{p_4}}}{S^{{e_{103}}}_{{p_2}}}{S^{{p_3}}_{{e_{103}}}}}");
    }
    // two (pure qp) 2-body operators
    {
      auto opseq =
          FNOperatorSeq({FNOperator({L"i_1", L"i_2"}, {L"a_1", L"a_2"}, V), FNOperator({L"a_3", L"a_4"}, {L"i_3", L"i_4"}, V)});
      auto wick = FWickTheorem{opseq};
      REQUIRE_NOTHROW(wick.full_contractions(true).spinfree(false).compute());
      auto result = wick.full_contractions(true).spinfree(false).compute();
      REQUIRE(result->size() == 4);
    }
    // two (pure qp) 3-body operators
    {
      auto opseq =
          FNOperatorSeq({FNOperator({L"i_1", L"i_2", L"i_3"}, {L"a_1", L"a_2", L"a_3"}, V),
                         FNOperator({L"a_4", L"a_5", L"a_6"}, {L"i_4", L"i_5", L"i_6"}, V)});
      auto wick = FWickTheorem{opseq};
      REQUIRE_NOTHROW(wick.full_contractions(true).spinfree(false).compute());
      auto result = wick.full_contractions(true).spinfree(false).compute();
      REQUIRE(result->size() == 36);
    }

    // two general 2-body operators
    {
      auto opseq =
          FNOperatorSeq({FNOperator({L"p_1", L"p_2"}, {L"p_3", L"p_4"}, V),
                         FNOperator({L"p_5", L"p_6"}, {L"p_7", L"p_8"}, V)});
      auto wick = FWickTheorem{opseq};
      REQUIRE_NOTHROW(wick.full_contractions(true).spinfree(false).compute());
      auto result = wick.full_contractions(true).spinfree(false).compute();
      REQUIRE(result->size() == 4);
    }

    // two general 3-body operators
    {
      auto opseq =
          FNOperatorSeq({FNOperator({L"p_1", L"p_2", L"p_3"}, {L"p_4", L"p_5", L"p_6"}, V),
                         FNOperator({L"p_7", L"p_8", L"p_9"}, {L"p_10", L"p_11", L"p_12"}, V)});
      auto wick = FWickTheorem{opseq};
      REQUIRE_NOTHROW(wick.full_contractions(true).spinfree(false).compute());
      auto result = wick.full_contractions(true).spinfree(false).compute();
      REQUIRE(result->size() == 36);
    }

    // two N-nonconserving operators
    {
      auto opseq =
          FNOperatorSeq({FNOperator({L"p_1", L"p_2", L"p_3"}, {L"p_4", L"p_5"}, V),
                         FNOperator({L"p_7", L"p_8"}, {L"p_10", L"p_11", L"p_12"}, V)});
      auto wick = FWickTheorem{opseq};
      REQUIRE_NOTHROW(wick.full_contractions(true).spinfree(false).compute());
      auto result = wick.full_contractions(true).spinfree(false).compute();
      REQUIRE(result->size() == 12);
    }

    // odd number of ops -> full contraction is 0
    {
      auto opseq =
          FNOperatorSeq({FNOperator({L"p_1", L"p_2"}, {L"p_4", L"p_5"}, V),
                         FNOperator({L"p_7", L"p_8"}, {L"p_10", L"p_11", L"p_12"}, V)});
      auto wick = FWickTheorem{opseq};
      REQUIRE_NOTHROW(wick.full_contractions(true).spinfree(false).compute());
      auto result = wick.full_contractions(true).spinfree(false).compute();
      REQUIRE(!result);
    }

    // 4-body ^ 4-body
    SEQUANT_PROFILE_SINGLE("wick(4^4)", {
      auto opseq =
          FNOperatorSeq({FNOperator({L"p_1", L"p_2", L"p_3", L"p_4"}, {L"p_5", L"p_6", L"p_7", L"p_8"}, V),
                         FNOperator({L"p_21", L"p_22", L"p_23", L"p_24"}, {L"p_25", L"p_26", L"p_27", L"p_28"}, V)
                        });
      auto wick = FWickTheorem{opseq};
      auto result = wick.full_contractions(true).spinfree(false).compute(true);
      REQUIRE(result->size() == 576);
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
      REQUIRE(result->size() == 2);
    }

    // 4 general 1-body operators
    {
      auto opseq =
          FNOperatorSeq({FNOperator({L"p_1"}, {L"p_2"}, V),
                         FNOperator({L"p_3"}, {L"p_4"}, V),
                         FNOperator({L"p_5"}, {L"p_6"}, V),
                         FNOperator({L"p_7"}, {L"p_8"}, V)
                        });
      auto ext_indices =
          make_indices<std::vector<Index>>(WstrList{L"p_1", L"p_2", L"p_3", L"p_4", L"p_5", L"p_6", L"p_7", L"p_8"});
      auto wick1 = FWickTheorem{opseq};
      auto result1 = wick1.full_contractions(true).set_external_indices(ext_indices).spinfree(false).compute();
//      std::wcout << "G1*G1*G1*G1 = " << to_latex_align(result1) << std::endl;
      REQUIRE(result1->size() == 9);
      auto wick2 = FWickTheorem{opseq};
      auto result2 =
          wick2.full_contractions(true).set_external_indices(ext_indices).set_op_connections({
              {1, 2}, {1, 3}}).spinfree(false).compute();
//      std::wcout << "G1*G1*G1*G1(1->2, 1->3) = " << to_latex_align(result2) << std::endl;
      REQUIRE(result1->size() == 9);
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
      REQUIRE(result->size() == 576);
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
      REQUIRE(result->size() == 80);
    }

    // 2-body ^ 2-body ^ 2-body ^ 2-body
    SEQUANT_PROFILE_SINGLE("wick(2^2^2^2)", {
      auto opseq =
          FNOperatorSeq({FNOperator({L"p_1", L"p_2"}, {L"p_5", L"p_6"}, V),
                         FNOperator({L"p_9", L"p_10"}, {L"p_11", L"p_12"}, V),
                         FNOperator({L"p_13", L"p_14"}, {L"p_15", L"p_16"}, V),
                         FNOperator({L"p_17", L"p_18"}, {L"p_19", L"p_20"}, V)
                        });
      auto wick = FWickTheorem{opseq};
      auto result = wick.full_contractions(true).spinfree(false).compute(true);
      REQUIRE(result->size() == 4752);
    }
    )

#ifndef SEQUANT_SKIP_LONG_TESTS
    // 4-body ^ 2-body ^ 2-body ^ 2-body
    SEQUANT_PROFILE_SINGLE("wick(4^2^2^2)", {
      auto opseq =
          FNOperatorSeq({FNOperator({L"p_1", L"p_2", L"p_3", L"p_4"}, {L"p_5", L"p_6", L"p_7", L"p_8"}, V),
                         FNOperator({L"p_9", L"p_10"}, {L"p_11", L"p_12"}, V),
                         FNOperator({L"p_13", L"p_14"}, {L"p_15", L"p_16"}, V),
                         FNOperator({L"p_17", L"p_18"}, {L"p_19", L"p_20"}, V)
                        });
      auto wick = FWickTheorem{opseq};
      auto result = wick.full_contractions(true)
                        .spinfree(false)
                        .use_topology(true)
                        .compute(true);
      REQUIRE(result->size() == 2088);
    })

    // 3-body ^ 2-body ^ 2-body ^ 3-body
    SEQUANT_PROFILE_SINGLE("wick(3^2^2^3)", {
      auto opseq = FNOperatorSeq(
          {FNOperator({L"p_1", L"p_2", L"p_3"}, {L"p_5", L"p_6", L"p_7"}, V),
           FNOperator({L"p_9", L"p_10"}, {L"p_11", L"p_12"}, V),
           FNOperator({L"p_13", L"p_14"}, {L"p_15", L"p_16"}, V),
           FNOperator({L"p_17", L"p_18", L"p_19"}, {L"p_20", L"p_21", L"p_22"},
                      V)});
      auto wick = FWickTheorem{opseq};
      auto result = wick.full_contractions(true)
                        .spinfree(false)
                        .use_topology(true)
                        .compute(true);
      REQUIRE(result->size() == 694);
    })

    // 4-body ^ 2-body ^ 4-body
    SEQUANT_PROFILE_SINGLE("wick(4^2^4)", {
      auto opseq =
          FNOperatorSeq({FNOperator({L"p_1", L"p_2", L"p_3", L"p_4"},
                                    {L"p_5", L"p_6", L"p_7", L"p_8"}, V),
                         FNOperator({L"p_9", L"p_10"}, {L"p_11", L"p_12"}, V),
                         FNOperator({L"p_21", L"p_22", L"p_23", L"p_24"},
                                    {L"p_25", L"p_26", L"p_27", L"p_28"}, V)});
      auto wick = FWickTheorem{opseq};
      auto result = wick.full_contractions(true)
                        .spinfree(false)
                        .use_topology(true)
                        .compute(true);
      REQUIRE(result->size() == 28);
    })

    // 4-body ^ 4-body ^ 4-body
    SEQUANT_PROFILE_SINGLE("wick(4^4^4)", {
      auto opseq =
          FNOperatorSeq({FNOperator({L"p_1", L"p_2", L"p_3", L"p_4"},
                                    {L"p_5", L"p_6", L"p_7", L"p_8"}, V),
                         FNOperator({L"p_11", L"p_12", L"p_13", L"p_14"},
                                    {L"p_15", L"p_16", L"p_17", L"p_18"}, V),
                         FNOperator({L"p_21", L"p_22", L"p_23", L"p_24"},
                                    {L"p_25", L"p_26", L"p_27", L"p_28"}, V)});
      auto wick = FWickTheorem{opseq};
      auto result = wick.full_contractions(true)
                        .spinfree(false)
                        .use_topology(true)
                        .compute(true);
      REQUIRE(result->size() == 70);
    })
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
      auto result = wick.full_contractions(true).spinfree(false).use_topology(true).compute(true);
      REQUIRE(result->size() == 0);
    }
#endif
  }  // SECTION("fermi vacuum")

  auto print = [](const auto &lead, const auto &expr) {
    std::wcout << lead << to_latex(expr) << std::endl;
  };

  SECTION("Expression Reduction") {
    constexpr Vacuum V = Vacuum::SingleProduct;

    // make these Tensor labels special ...
    // these are already lexicographically ordered, but tensors that do not
    // appear on this list (e.g. "S") will appear after Tensors with these
    // labels
    TensorCanonicalizer::set_cardinal_tensor_labels({L"A", L"f", L"g", L"t"});

    // 2-body ^ 2-body
    SEQUANT_PROFILE_SINGLE("wick(H2*T2)", {
      auto opseq =
          FNOperatorSeq({FNOperator({L"p_1", L"p_2"}, {L"p_3", L"p_4"}, V),
                         FNOperator({L"a_4", L"a_5"}, {L"i_4", L"i_5"}, V)});
      auto wick = FWickTheorem{opseq};
      auto wick_result = wick.full_contractions(true).spinfree(false).compute();
      REQUIRE(wick_result->size() == 4);

      // multiply tensor factors and expand
      auto wick_result_2 =
          ex<Tensor>(L"g", WstrList{L"p_1", L"p_2"}, WstrList{L"p_3", L"p_4"},
                     Symmetry::antisymm) *
          ex<Tensor>(L"t", WstrList{L"a_4", L"a_5"}, WstrList{L"i_4", L"i_5"},
                     Symmetry::antisymm) *
          wick_result;
      expand(wick_result_2);
      REQUIRE(to_latex(wick_result_2) ==
              L"{ "
              L"\\left({{g^{{p_3}{p_4}}_{{p_1}{p_2}}}{t^{{i_4}{i_5}}_{{a_4}{a_"
              L"5}}}{S^{{p_1}}_{{i_5}}}{S^{{p_2}}_{{i_4}}}{S^{{a_4}}_{{p_4}}}{"
              L"S^{{a_5}}_{{p_3}}}} + {{{-1}} \\times "
              L"{g^{{p_3}{p_4}}_{{p_1}{p_2}}}{t^{{i_4}{i_5}}_{{a_4}{a_5}}}{S^{{"
              L"p_1}}_{{i_5}}}{S^{{p_2}}_{{i_4}}}{S^{{a_5}}_{{p_4}}}{S^{{a_4}}_"
              L"{{p_3}}}} + {{{-1}} \\times "
              L"{g^{{p_3}{p_4}}_{{p_1}{p_2}}}{t^{{i_4}{i_5}}_{{a_4}{a_5}}}{S^{{"
              L"p_1}}_{{i_4}}}{S^{{p_2}}_{{i_5}}}{S^{{a_4}}_{{p_4}}}{S^{{a_5}}_"
              L"{{p_3}}}} + "
              L"{{g^{{p_3}{p_4}}_{{p_1}{p_2}}}{t^{{i_4}{i_5}}_{{a_4}{a_5}}}{S^{"
              L"{p_1}}_{{i_4}}}{S^{{p_2}}_{{i_5}}}{S^{{a_5}}_{{p_4}}}{S^{{a_4}}"
              L"_{{p_3}}}}\\right) }");
      wick.reduce(wick_result_2);
      rapid_simplify(wick_result_2);
      TensorCanonicalizer::register_instance(
          std::make_shared<DefaultTensorCanonicalizer>());
      canonicalize(wick_result_2);
      rapid_simplify(wick_result_2);

      std::wcout << L"H2*T2 = " << to_latex(wick_result_2) << std::endl;
      std::wcout << L"H2*T2 = " << to_wolfram(wick_result_2) << std::endl;
      REQUIRE(to_latex(wick_result_2) ==
              L"{{{4}} \\times "
              L"{g^{{a_1}{a_2}}_{{i_1}{i_2}}}{t^{{i_1}{i_2}}_{{a_1}{a_2}}}}");
    });

    // 2-body ^ 1-body ^ 1-body
    SEQUANT_PROFILE_SINGLE("wick(H2*T1*T1)", {
      constexpr bool use_op_partitions = true;
      auto opseq =
          FNOperatorSeq({FNOperator({L"p_1", L"p_2"}, {L"p_3", L"p_4"}, V),
                         FNOperator({L"a_4"}, {L"i_4"}, V),
                         FNOperator({L"a_5"}, {L"i_5"}, V)});
      auto wick = FWickTheorem{opseq};
      wick.full_contractions(true).spinfree(false).use_topology(true);
      if (use_op_partitions)
        wick.set_op_partitions({{1,2}});
      auto wick_result = wick.compute();
      REQUIRE(wick_result->size() == (use_op_partitions ? 2 : 4));

      // multiply tensor factors and expand
      auto wick_result_2 =
          ex<Tensor>(L"g", WstrList{L"p_1", L"p_2"}, WstrList{L"p_3", L"p_4"},
                     Symmetry::antisymm) *
          ex<Tensor>(L"t", WstrList{L"a_4"}, WstrList{L"i_4"},
                     Symmetry::antisymm) *
          ex<Tensor>(L"t", WstrList{L"a_5"}, WstrList{L"i_5"},
                     Symmetry::antisymm) *
          wick_result;
      expand(wick_result_2);
      wick.reduce(wick_result_2);
      rapid_simplify(wick_result_2);
      TensorCanonicalizer::register_instance(
          std::make_shared<DefaultTensorCanonicalizer>(std::vector<Index>{}));
      canonicalize(wick_result_2);
      rapid_simplify(wick_result_2);

      print("H2*T1*T1 = ", wick_result_2);
      REQUIRE(to_latex(wick_result_2) ==
              L"{{{4}} \\times "
              L"{g^{{a_1}{a_2}}_{{i_1}{i_2}}}{t^{{i_1}}_{{a_1}}}{t^{{i_2}}_{{a_"
              L"2}}}}");
    });

    // 2=body ^ 1-body ^ 2-body with dependent (PNO) indices
    SEQUANT_PROFILE_SINGLE("wick(P2*H1*T2)", {
      auto opseq = FNOperatorSeq({FNOperator(IndexList{L"i_1", L"i_2"},
                                             {Index(L"a_1", {L"i_1", L"i_2"}),
                                              Index(L"a_2", {L"i_1", L"i_2"})},
                                             V),
                                  FNOperator({L"p_1"}, {L"p_2"}, V),
                                  FNOperator({Index(L"a_3", {L"i_3", L"i_4"}),
                                              Index(L"a_4", {L"i_3", L"i_4"})},
                                             IndexList{L"i_3", L"i_4"}, V)});
      auto wick = FWickTheorem{opseq};
      auto wick_result = wick.full_contractions(true).spinfree(false).compute();
      REQUIRE(wick_result->size() == 16);

      // multiply tensor factors and expand
      auto wick_result_2 =
          ex<Tensor>(
              L"A", IndexList{L"i_1", L"i_2"},
              IndexList{{L"a_1", {L"i_1", L"i_2"}}, {L"a_2", {L"i_1", L"i_2"}}},
              Symmetry::antisymm) *
          ex<Tensor>(L"f", WstrList{L"p_1"}, WstrList{L"p_2"},
                     Symmetry::antisymm) *
          ex<Tensor>(
              L"t",
              IndexList{{L"a_3", {L"i_3", L"i_4"}}, {L"a_4", {L"i_3", L"i_4"}}},
              IndexList{L"i_3", L"i_4"}, Symmetry::antisymm) *
          wick_result;
      expand(wick_result_2);
      wick.reduce(wick_result_2);
      rapid_simplify(wick_result_2);
      TensorCanonicalizer::register_instance(
          std::make_shared<DefaultTensorCanonicalizer>());
      canonicalize(wick_result_2);
      rapid_simplify(wick_result_2);

      std::wcout << L"P2*H1*T2(PNO) = " << to_latex_align(wick_result_2)
                 << std::endl;
      REQUIRE(to_latex(wick_result_2) ==
              L"{ \\left({{{-8}} \\times "
              L"{A^{{a_1^{{i_1}{i_2}}}{a_2^{{i_1}{i_2}}}}_{{i_1}{i_2}}}{f^{{a_"
              L"3^{{i_1}{i_2}}}}_{{a_1^{{i_1}{i_2}}}}}{t^{{i_1}{i_2}}_{{a_2^{{"
              L"i_1}{i_2}}}{a_3^{{i_1}{i_2}}}}}} + {{{8}} \\times "
              L"{A^{{a_1^{{i_1}{i_2}}}{a_2^{{i_1}{i_2}}}}_{{i_1}{i_2}}}{f^{{i_"
              L"1}}_{{i_3}}}{t^{{i_2}{i_3}}_{{a_3^{{i_2}{i_3}}}{a_4^{{i_2}{i_3}"
              L"}}}}{S^{{a_3^{{i_2}{i_3}}}}_{{a_1^{{i_1}{i_2}}}}}{S^{{a_4^{{i_"
              L"2}{i_3}}}}_{{a_2^{{i_1}{i_2}}}}}}\\right) }");
    });

    // 2=body ^ 2-body ^ 2-body ^ 2-body with dependent (PNO) indices
    SEQUANT_PROFILE_SINGLE("wick(P2*H2*T2*T2)", {
      constexpr bool use_op_partitions = true;
      auto opseq =
          FNOperatorSeq({FNOperator(IndexList{L"i_1", L"i_2"},
                                    {Index(L"a_1", {L"i_1", L"i_2"}),
                                     Index(L"a_2", {L"i_1", L"i_2"})},
                                    V),
                         FNOperator({L"p_1", L"p_2"}, {L"p_3", L"p_4"}, V),
                         FNOperator({Index(L"a_3", {L"i_3", L"i_4"}),
                                     Index(L"a_4", {L"i_3", L"i_4"})},
                                    IndexList{L"i_3", L"i_4"}, V),
                         FNOperator({Index(L"a_5", {L"i_5", L"i_6"}),
                                     Index(L"a_6", {L"i_5", L"i_6"})},
                                    IndexList{L"i_5", L"i_6"}, V)});
      auto wick = FWickTheorem{opseq};
      wick.full_contractions(true)
              .spinfree(false)
              .set_op_connections({{1, 2}, {1, 3}})
              .use_topology(true);

      if (use_op_partitions)
        wick.set_op_partitions({{2,3}});
      auto wick_result =
          wick.compute();
      REQUIRE(wick_result->size() == (use_op_partitions ? 17 : 34));

      // multiply tensor factors and expand
      auto wick_result_2 =
          ex<Constant>(1. / 256) *
          ex<Tensor>(
              L"A", IndexList{L"i_1", L"i_2"},
              IndexList{{L"a_1", {L"i_1", L"i_2"}}, {L"a_2", {L"i_1", L"i_2"}}},
              Symmetry::antisymm) *
          ex<Tensor>(L"g", WstrList{L"p_1", L"p_2"}, WstrList{L"p_3", L"p_4"},
                     Symmetry::antisymm) *
          ex<Tensor>(
              L"t",
              IndexList{{L"a_3", {L"i_3", L"i_4"}}, {L"a_4", {L"i_3", L"i_4"}}},
              IndexList{L"i_3", L"i_4"}, Symmetry::antisymm) *
          ex<Tensor>(
              L"t",
              IndexList{{L"a_5", {L"i_5", L"i_6"}}, {L"a_6", {L"i_5", L"i_6"}}},
              IndexList{L"i_5", L"i_6"}, Symmetry::antisymm) *
          wick_result;
      expand(wick_result_2);
      wick.reduce(wick_result_2);
      rapid_simplify(wick_result_2);
      TensorCanonicalizer::register_instance(
          std::make_shared<DefaultTensorCanonicalizer>());
      canonicalize(wick_result_2);
      canonicalize(wick_result_2);
      canonicalize(wick_result_2);
      rapid_simplify(wick_result_2);

      std::wcout << L"P2*H2*T2*T2(PNO) = " << to_latex_align(wick_result_2, 20)
                 << std::endl;
      REQUIRE(wick_result_2->size() == 4);
    });

#if 1
    // 3-body ^ 2-body ^ 2-body ^ 3-body
    SEQUANT_PROFILE_SINGLE("wick(P3*H2*T2*T3)", {
      constexpr bool connected_only = true;
      constexpr bool topology = true;
      auto P3 =
          ex<Constant>(1. / 36) *
          ex<Tensor>(L"A", WstrList{L"i_1", L"i_2", L"i_3"},
                     WstrList{L"a_1", L"a_2", L"a_3"}, Symmetry::antisymm) *
          ex<FNOperator>(WstrList{L"i_1", L"i_2", L"i_3"},
                         WstrList{L"a_1", L"a_2", L"a_3"}, V);
      auto H2 =
          ex<Constant>(1. / 4) *
          ex<Tensor>(L"g", WstrList{L"p_1", L"p_2"}, WstrList{L"p_3", L"p_4"},
                     Symmetry::antisymm) *
          ex<FNOperator>(WstrList{L"p_1", L"p_2"}, WstrList{L"p_3", L"p_4"}, V);
      auto T2 =
          ex<Constant>(1. / 4) *
          ex<Tensor>(L"t", WstrList{L"a_4", L"a_5"}, WstrList{L"i_4", L"i_5"},
                     Symmetry::antisymm) *
          ex<FNOperator>(WstrList{L"a_4", L"a_5"}, WstrList{L"i_4", L"i_5"}, V);
      auto T3 =
          ex<Constant>(1. / 36) *
          ex<Tensor>(L"t", WstrList{L"a_6", L"a_7", L"a_8"},
                     WstrList{L"i_6", L"i_7", L"i_8"}, Symmetry::antisymm) *
          ex<FNOperator>(WstrList{L"a_6", L"a_7", L"a_8"},
                         WstrList{L"i_6", L"i_7", L"i_8"}, V);
      auto ext_indices = make_indices<std::vector<Index>>({});
      auto wick = FWickTheorem{P3 * H2 * T2 * T3}
                      .full_contractions(true)
                      .set_external_indices(ext_indices)
                      .spinfree(false)
                      .use_topology(topology);
      if (connected_only)
        wick.set_op_connections({{1, 2}, {1, 3}});
      auto wick_result = wick.compute();

      std::wcout << "P3*H2*T2*T3 = " << to_latex_align(wick_result, 20)
                 << std::endl;
      REQUIRE(
          wick_result->size() ==
          (connected_only ? 7 : 9));  // 9 = 2 disconnected + 7 connected terms
    });
#endif
  }

  }  // TEST_CASE("WickTheorem")
#endif
