//
// Created by Eduard Valeyev on 3/23/18.
//

#include <iostream>

#include "SeQuant/core/timer.hpp"
#include "SeQuant/core/wick.hpp"
#include "catch.hpp"
#include "test_config.hpp"

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
TEST_CASE("WickTheorem", "[algorithms][wick]") {
  using namespace sequant;

  TensorCanonicalizer::register_instance(
      std::make_shared<DefaultTensorCanonicalizer>());
  Index::reset_tmp_index();

  SECTION("Op contractions") {
    REQUIRE(FWickTheorem::can_contract(fann(L"i_1"), fcre(L"i_2"),
                                       Vacuum::Physical));
    REQUIRE(!FWickTheorem::can_contract(fcre(L"i_1"), fcre(L"i_2"),
                                        Vacuum::Physical));
    REQUIRE(!FWickTheorem::can_contract(fcre(L"i_1"), fann(L"i_2"),
                                        Vacuum::Physical));
    REQUIRE(!FWickTheorem::can_contract(fann(L"i_1"), fann(L"i_2"),
                                        Vacuum::Physical));

    REQUIRE(!FWickTheorem::can_contract(fann(L"i_1"), fcre(L"i_2"),
                                        Vacuum::SingleProduct));
    REQUIRE(!FWickTheorem::can_contract(fcre(L"i_1"), fcre(L"i_2"),
                                        Vacuum::SingleProduct));
    REQUIRE(FWickTheorem::can_contract(fcre(L"i_1"), fann(L"i_2"),
                                       Vacuum::SingleProduct));
    REQUIRE(!FWickTheorem::can_contract(fann(L"i_1"), fann(L"i_2"),
                                        Vacuum::SingleProduct));

    REQUIRE(!FWickTheorem::can_contract(fann(L"i_1"), fcre(L"a_2"),
                                        Vacuum::SingleProduct));
    REQUIRE(!FWickTheorem::can_contract(fcre(L"i_1"), fcre(L"a_2"),
                                        Vacuum::SingleProduct));
    REQUIRE(!FWickTheorem::can_contract(fcre(L"i_1"), fann(L"a_2"),
                                        Vacuum::SingleProduct));
    REQUIRE(!FWickTheorem::can_contract(fann(L"i_1"), fann(L"a_2"),
                                        Vacuum::SingleProduct));

    REQUIRE(!FWickTheorem::can_contract(fann(L"a_1"), fcre(L"i_2"),
                                        Vacuum::SingleProduct));
    REQUIRE(!FWickTheorem::can_contract(fcre(L"a_1"), fcre(L"i_2"),
                                        Vacuum::SingleProduct));
    REQUIRE(!FWickTheorem::can_contract(fcre(L"a_1"), fann(L"i_2"),
                                        Vacuum::SingleProduct));
    REQUIRE(!FWickTheorem::can_contract(fann(L"a_1"), fann(L"i_2"),
                                        Vacuum::SingleProduct));

    REQUIRE(FWickTheorem::can_contract(fann(L"a_1"), fcre(L"a_2"),
                                       Vacuum::SingleProduct));
    REQUIRE(!FWickTheorem::can_contract(fcre(L"a_1"), fcre(L"a_2"),
                                        Vacuum::SingleProduct));
    REQUIRE(!FWickTheorem::can_contract(fcre(L"a_1"), fann(L"a_2"),
                                        Vacuum::SingleProduct));
    REQUIRE(!FWickTheorem::can_contract(fann(L"a_1"), fann(L"a_2"),
                                        Vacuum::SingleProduct));

    REQUIRE(!FWickTheorem::can_contract(fann(L"p_1"), fcre(L"i_2"),
                                        Vacuum::SingleProduct));
    REQUIRE(!FWickTheorem::can_contract(fcre(L"p_1"), fcre(L"i_2"),
                                        Vacuum::SingleProduct));
    REQUIRE(FWickTheorem::can_contract(fcre(L"p_1"), fann(L"i_2"),
                                       Vacuum::SingleProduct));
    REQUIRE(!FWickTheorem::can_contract(fann(L"p_1"), fann(L"i_2"),
                                        Vacuum::SingleProduct));

    REQUIRE(BWickTheorem::can_contract(bann(L"i_1"), bcre(L"i_2"),
                                       Vacuum::Physical));
    REQUIRE(!BWickTheorem::can_contract(bcre(L"i_1"), bcre(L"i_2"),
                                        Vacuum::Physical));
    REQUIRE(!BWickTheorem::can_contract(bcre(L"i_1"), bann(L"i_2"),
                                        Vacuum::Physical));
    REQUIRE(!BWickTheorem::can_contract(bann(L"i_1"), bann(L"i_2"),
                                        Vacuum::Physical));
  }

  SECTION("constructors") {
    REQUIRE_NOTHROW(FWickTheorem{FNOperatorSeq{}});
    REQUIRE_NOTHROW(BWickTheorem{BNOperatorSeq{}});

    {
      auto opseq1 = FNOperatorSeq({FNOperator({L"i_1"}, {L"i_2"}),
                                   FNOperator({L"i_3"}, {L"i_4"}),
                                   FNOperator({L"i_5"}, {L"i_6"})});
      REQUIRE_NOTHROW(FWickTheorem{opseq1});
      auto wick1 = FWickTheorem{opseq1};
      REQUIRE_THROWS(wick1.spinfree(true).compute());
    }

  }  // SECTION("constructors")

  SECTION("physical vacuum") {
    constexpr Vacuum V = Vacuum::Physical;

    // number operator
    {{auto opseq1 = FNOperatorSeq(
          {FNOperator({L"i_1"}, {}, V), FNOperator({}, {L"i_2"}, V)});
    auto wick1 = FWickTheorem{opseq1};
    REQUIRE_NOTHROW(wick1.spinfree(false).compute());
    // full contractions = null (N is already in normal form)
    auto full_contractions = FWickTheorem{opseq1}.spinfree(false).compute();
    REQUIRE(full_contractions->is<Constant>());
    REQUIRE(full_contractions->as<Constant>().value<int>() == 0);
    // partial contractions = N
    auto partial_contractions =
        FWickTheorem{opseq1}.full_contractions(false).spinfree(false).compute();
    // std::wcout << "partial_contractions=" << to_latex(partial_contractions)
    // << std::endl;
    REQUIRE(partial_contractions->is<Product>());
    REQUIRE(partial_contractions->as<Product>().size() == 1);
  }
  {
    auto opseq1 = BNOperatorSeq(
        {BNOperator({L"i_1"}, {}, V), BNOperator({}, {L"i_2"}, V)});
    auto wick1 = BWickTheorem{opseq1};
    REQUIRE_NOTHROW(wick1.spinfree(false).compute());
    // full contractions = null
    auto full_contractions = BWickTheorem{opseq1}.spinfree(false).compute();
    REQUIRE(full_contractions->is<Constant>());
    REQUIRE(full_contractions->as<Constant>().value<int>() == 0);
    // partial contractions = N
    auto partial_contractions =
        BWickTheorem{opseq1}.full_contractions(false).spinfree(false).compute();
    // std::wcout << "partial_contractions=" << to_latex(partial_contractions)
    // << std::endl;
    REQUIRE(partial_contractions->is<Product>());
    REQUIRE(partial_contractions->as<Product>().size() == 1);
  }
}

// hole number operator
{{auto opseq1 =
      FNOperatorSeq({FNOperator({}, {L"i_1"}, V), FNOperator({L"i_2"}, {}, V)});
auto wick1 = FWickTheorem{opseq1};
REQUIRE_NOTHROW(wick1.spinfree(false).compute());
// full contractions = delta
auto full_contractions = FWickTheorem{opseq1}.spinfree(false).compute();
REQUIRE(full_contractions->is<Product>());
REQUIRE(full_contractions->as<Product>().size() == 1);
// partial contractions = delta - N
auto partial_contractions =
    FWickTheorem{opseq1}.full_contractions(false).spinfree(false).compute();
// std::wcout << "partial_contractions=" << to_latex(partial_contractions) <<
// std::endl;
REQUIRE(partial_contractions->is<Sum>());
REQUIRE(partial_contractions->as<Sum>().size() == 2);
REQUIRE(to_latex(partial_contractions) ==
        L"{ \\bigl({{s^{{i_2}}_{{i_1}}}} - {{a^{{i_2}}_{{i_1}}}}\\bigr) }");
}
{
  auto opseq1 =
      BNOperatorSeq({BNOperator({}, {L"i_1"}, V), BNOperator({L"i_2"}, {}, V)});
  auto wick1 = BWickTheorem{opseq1};
  REQUIRE_NOTHROW(wick1.spinfree(false).compute());
  // full contractions = delta
  auto full_contractions = BWickTheorem{opseq1}.spinfree(false).compute();
  REQUIRE(full_contractions->is<Product>());
  REQUIRE(full_contractions->as<Product>().size() == 1);
  // partial contractions = delta + N
  auto partial_contractions =
      BWickTheorem{opseq1}.full_contractions(false).spinfree(false).compute();
  // std::wcout << "partial_contractions=" << to_latex(partial_contractions) <<
  // std::endl;
  REQUIRE(partial_contractions->is<Sum>());
  REQUIRE(partial_contractions->as<Sum>().size() == 2);
  REQUIRE(to_latex(partial_contractions) ==
          L"{ \\bigl({{s^{{i_2}}_{{i_1}}}} + {{b^{{i_2}}_{{i_1}}}}\\bigr) }");
}
}

// three 1-body operators
{
  auto opseq1 = FNOperatorSeq({FNOperator({L"i_1"}, {L"i_2"}, V),
                               FNOperator({L"i_3"}, {L"i_4"}, V),
                               FNOperator({L"i_5"}, {L"i_6"}, V)});
  auto wick1 = FWickTheorem{opseq1};
  REQUIRE_NOTHROW(wick1.spinfree(false).compute());
  auto result = FWickTheorem{opseq1}.spinfree(false).compute();
  REQUIRE(result->is<Constant>());
  REQUIRE(result->as<Constant>().value<int>() == 0);
}

// two 2-body operators
{
  auto opseq = FNOperatorSeq({FNOperator({}, {L"i_1", L"i_2"}, V),
                              FNOperator({L"i_3", L"i_4"}, {}, V)});
  auto wick = FWickTheorem{opseq};
  REQUIRE_NOTHROW(wick.spinfree(false).compute());
  auto result = wick.spinfree(false).compute();
  REQUIRE(result->is<Sum>());
  REQUIRE(result->size() == 2);
}

// two 3-body operators
{
  auto opseq = FNOperatorSeq({FNOperator({}, {L"i_1", L"i_2", L"i_3"}, V),
                              FNOperator({L"i_4", L"i_5", L"i_6"}, {}, V)});
  auto wick = FWickTheorem{opseq};
  REQUIRE_NOTHROW(wick.spinfree(false).compute());
  auto result = wick.spinfree(false).compute();
  REQUIRE(result->is<Sum>());
  REQUIRE(result->size() == 6);
}

// two 4-body operators
{
  auto opseq =
      FNOperatorSeq({FNOperator({}, {L"i_1", L"i_2", L"i_3", L"i_4"}, V),
                     FNOperator({L"i_5", L"i_6", L"i_7", L"i_8"}, {}, V)});
  auto wick = FWickTheorem{opseq};
  REQUIRE_NOTHROW(wick.spinfree(false).compute());
  auto result = wick.spinfree(false).compute();
  REQUIRE(result->is<Sum>());
  REQUIRE(result->size() == 24);
}

// 1/2 * 1 * 1/2 body ops, full contraction
{
  auto opseq = FNOperatorSeq({FNOperator({}, {L"i_1"}, V),
                              FNOperator({L"i_2"}, {L"i_3"}, V),
                              FNOperator({L"i_4"}, {}, V)});
  auto wick = FWickTheorem{opseq};
  REQUIRE_NOTHROW(wick.spinfree(false).compute());
  auto result = wick.spinfree(false).compute();
  REQUIRE(result->is<Product>());
  REQUIRE(result->size() == 2);  // product of 2 terms
}

// 1/2 * 1 * 1/2 body ops, partial contraction
{
  auto opseq = FNOperatorSeq({FNOperator({}, {L"i_1"}, V),
                              FNOperator({L"i_2"}, {L"i_3"}, V),
                              FNOperator({L"i_4"}, {}, V)});
  auto wick = FWickTheorem{opseq};
  REQUIRE_NOTHROW(wick.full_contractions(false).spinfree(false).compute());
  auto result = wick.full_contractions(false).spinfree(false).compute();
  REQUIRE(result->is<Sum>());
  REQUIRE(result->size() == 5);  // sum of 4 terms
  REQUIRE(to_latex(result) ==
          L"{ \\bigl( - {{s^{{i_2}}_{{i_1}}}{a^{{i_4}}_{{i_3}}}} + "
          L"{{s^{{i_2}}_{{i_1}}}{s^{{i_4}}_{{i_3}}}} + "
          L"{{s^{{i_4}}_{{i_1}}}{a^{{i_2}}_{{i_3}}}} - "
          L"{{s^{{i_4}}_{{i_3}}}{a^{{i_2}}_{{i_1}}}} - "
          L"{{a^{{i_2}{i_4}}_{{i_3}{i_1}}}}\\bigr) }");
}

// three 1-body operators, partial contraction
{
  auto opseq1 = FNOperatorSeq({FNOperator({L"i_1"}, {L"i_2"}, V),
                               FNOperator({L"i_3"}, {L"i_4"}, V),
                               FNOperator({L"i_5"}, {L"i_6"}, V)});
  auto wick1 = FWickTheorem{opseq1};
  REQUIRE_NOTHROW(wick1.full_contractions(false).spinfree(false).compute());
  auto result =
      FWickTheorem{opseq1}.full_contractions(false).spinfree(false).compute();
  REQUIRE(result->is<Sum>());
  REQUIRE(result->size() == 5);
  REQUIRE(to_latex(result) ==
          L"{ \\bigl({{s^{{i_3}}_{{i_2}}}{a^{{i_1}{i_5}}_{{i_4}{i_6}}}} + "
          L"{{s^{{i_3}}_{{i_2}}}{s^{{i_5}}_{{i_4}}}{a^{{i_1}}_{{i_6}}}} + "
          L"{{s^{{i_5}}_{{i_2}}}{a^{{i_1}{i_3}}_{{i_6}{i_4}}}} + "
          L"{{s^{{i_5}}_{{i_4}}}{a^{{i_1}{i_3}}_{{i_2}{i_6}}}} + "
          L"{{a^{{i_1}{i_3}{i_5}}_{{i_2}{i_4}{i_6}}}}\\bigr) }");
}

// two 2-body operators, partial contraction: Eq. 9b of DOI 10.1063/1.474405
{
  auto opseq =
      FNOperatorSeq({FNOperator({L"i_1", L"i_2"}, {L"i_3", L"i_4"}, V),
                     FNOperator({L"i_5", L"i_6"}, {L"i_7", L"i_8"}, V)});
  auto wick = FWickTheorem{opseq};
  REQUIRE_NOTHROW(wick.full_contractions(false).spinfree(false).compute());
  auto result = wick.full_contractions(false).spinfree(false).compute();
  REQUIRE(result->is<Sum>());
  REQUIRE(result->size() == 7);
  REQUIRE(to_latex(result) ==
          L"{ "
          L"\\bigl({{s^{{i_5}}_{{i_4}}}{a^{{i_1}{i_2}{i_6}}_{{i_3}{i_7}{i_8}}}}"
          L" + "
          L"{{s^{{i_5}}_{{i_4}}}{s^{{i_6}}_{{i_3}}}{a^{{i_1}{i_2}}_{{i_8}{i_7}}"
          L"}} + {{s^{{i_6}}_{{i_4}}}{a^{{i_1}{i_2}{i_5}}_{{i_3}{i_8}{i_7}}}} "
          L"+ "
          L"{{s^{{i_6}}_{{i_4}}}{s^{{i_5}}_{{i_3}}}{a^{{i_1}{i_2}}_{{i_7}{i_8}}"
          L"}} + {{s^{{i_5}}_{{i_3}}}{a^{{i_1}{i_2}{i_6}}_{{i_7}{i_4}{i_8}}}} "
          L"+ {{s^{{i_6}}_{{i_3}}}{a^{{i_1}{i_2}{i_5}}_{{i_8}{i_4}{i_7}}}} + "
          L"{{a^{{i_1}{i_2}{i_5}{i_6}}_{{i_3}{i_4}{i_7}{i_8}}}}\\bigr) }");
}

}  // SECTION("physical vacuum")

SECTION("fermi vacuum") {
  constexpr Vacuum V = Vacuum::SingleProduct;

  // two (pure qp) 1-body operators
  {
    auto opseq = FNOperatorSeq(
        {FNOperator({L"i_1"}, {L"a_1"}, V), FNOperator({L"a_2"}, {L"i_2"}, V)});
    auto wick = FWickTheorem{opseq};
    REQUIRE_NOTHROW(wick.spinfree(false).compute());
    auto result = wick.spinfree(false).compute();
    REQUIRE(result->is<Product>());
    REQUIRE(result->size() == 2);  // product of 2 terms
  }

  // two (pure qp) N-nonconserving 2-body operators
  {
    auto opseq = FNOperatorSeq({FNOperator({L"i_1", L"i_2"}, {L"a_1"}, V),
                                FNOperator({L"a_2"}, {L"i_3", L"i_4"}, V)});
    auto wick = FWickTheorem{opseq};
    REQUIRE_NOTHROW(wick.spinfree(false).compute());
    auto result = wick.spinfree(false).compute();
    REQUIRE(result->is<Sum>());
    REQUIRE(result->size() == 2);
  }

  // two general 1-body operators
  {
    auto opseq = FNOperatorSeq(
        {FNOperator({L"p_1"}, {L"p_2"}, V), FNOperator({L"p_3"}, {L"p_4"}, V)});
    auto wick = FWickTheorem{opseq};
    REQUIRE_NOTHROW(wick.spinfree(false).compute());
    auto result = wick.spinfree(false).compute();
    REQUIRE(result->is<Product>());
    REQUIRE(result->size() ==
            2 * 2);  // product of 4 terms (since each contraction of 2
                     // *general* indices produces 2 overlaps)
    REQUIRE(to_latex(result) ==
            L"{{s^{{p_1}}_{{m_{102}}}}{s^{{m_{102}}}_{{p_4}}}{s^{{e_{103}}}_{"
            L"{p_2}}}{s^{{p_3}}_{{e_{103}}}}}");
  }
  // two general 1-body operators, partial contractions: Eq. 21a of
  // DOI 10.1063/1.474405
  {
    auto opseq = FNOperatorSeq(
        {FNOperator({L"p_1"}, {L"p_2"}, V), FNOperator({L"p_3"}, {L"p_4"}, V)});
    auto wick = FWickTheorem{opseq};
    REQUIRE_NOTHROW(wick.full_contractions(false).spinfree(false).compute());
    auto result = wick.full_contractions(false).spinfree(false).compute();
    REQUIRE(result->is<Sum>());
    REQUIRE(result->size() == 4);
    REQUIRE(
        to_latex(result) ==
        L"{ \\bigl( - "
        L"{{s^{{p_1}}_{{m_{107}}}}{s^{{m_{107}}}_{{p_4}}}{\\tilde{a}^{{p_3}}_"
        L"{{p_2}}}} + "
        L"{{s^{{p_1}}_{{m_{107}}}}{s^{{m_{107}}}_{{p_4}}}{s^{{e_{108}}}_{{p_"
        L"2}}}{s^{{p_3}}_{{e_{108}}}}} + "
        L"{{s^{{e_{109}}}_{{p_2}}}{s^{{p_3}}_{{e_{109}}}}{\\tilde{a}^{{p_1}}_"
        L"{{p_4}}}} + {{\\tilde{a}^{{p_1}{p_3}}_{{p_2}{p_4}}}}\\bigr) }");
  }

  // two (pure qp) 2-body operators
  {
    auto opseq =
        FNOperatorSeq({FNOperator({L"i_1", L"i_2"}, {L"a_1", L"a_2"}, V),
                       FNOperator({L"a_3", L"a_4"}, {L"i_3", L"i_4"}, V)});
    auto wick = FWickTheorem{opseq};
    REQUIRE_NOTHROW(wick.spinfree(false).compute());
    auto result = wick.spinfree(false).compute();
    REQUIRE(result->is<Sum>());
    REQUIRE(result->size() == 4);
  }
  // two (pure qp) 3-body operators
  {
    auto opseq = FNOperatorSeq(
        {FNOperator({L"i_1", L"i_2", L"i_3"}, {L"a_1", L"a_2", L"a_3"}, V),
         FNOperator({L"a_4", L"a_5", L"a_6"}, {L"i_4", L"i_5", L"i_6"}, V)});
    auto wick = FWickTheorem{opseq};
    REQUIRE_NOTHROW(wick.spinfree(false).compute());
    auto result = wick.spinfree(false).compute();
    REQUIRE(result->is<Sum>());
    REQUIRE(result->size() == 36);
  }

  // one general 1-body operator + one general 2-body operator, partial
  // contraction: Eq. 9 of DOI 10.1063/1.474405
  {
    auto opseq =
        FNOperatorSeq({FNOperator({L"p_1"}, {L"p_2"}, V),
                       FNOperator({L"p_3", L"p_4"}, {L"p_5", L"p_6"}, V)});
    auto wick = FWickTheorem{opseq};
    REQUIRE_NOTHROW(wick.full_contractions(false).spinfree(false).compute());
    auto result = wick.full_contractions(false).spinfree(false).compute();
    REQUIRE(result->is<Sum>());
    REQUIRE(result->size() == 9);
  }

  // two general 2-body operators
  {
    auto opseq =
        FNOperatorSeq({FNOperator({L"p_1", L"p_2"}, {L"p_3", L"p_4"}, V),
                       FNOperator({L"p_5", L"p_6"}, {L"p_7", L"p_8"}, V)});
    auto wick = FWickTheorem{opseq};
    REQUIRE_NOTHROW(wick.spinfree(false).compute());
    auto result = wick.spinfree(false).compute();
    REQUIRE(result->is<Sum>());
    REQUIRE(result->size() == 4);
  }
  // two general 2-body operators, partial contractions: Eqs. 22 of
  // DOI 10.1063/1.474405
  {
    auto opseq =
        FNOperatorSeq({FNOperator({L"p_1", L"p_2"}, {L"p_3", L"p_4"}, V),
                       FNOperator({L"p_5", L"p_6"}, {L"p_7", L"p_8"}, V)});
    auto wick = FWickTheorem{opseq};
    REQUIRE_NOTHROW(wick.full_contractions(false).spinfree(false).compute());
    auto result = wick.full_contractions(false).spinfree(false).compute();
    REQUIRE(result->is<Sum>());
    REQUIRE(result->size() == 49);  // the MK paper only gives 47 terms,
                                    // misses the 2 double-hole contractions
  }

  // two general 3-body operators
  {
    auto opseq = FNOperatorSeq(
        {FNOperator({L"p_1", L"p_2", L"p_3"}, {L"p_4", L"p_5", L"p_6"}, V),
         FNOperator({L"p_7", L"p_8", L"p_9"}, {L"p_10", L"p_11", L"p_12"}, V)});
    auto wick = FWickTheorem{opseq};
    REQUIRE_NOTHROW(wick.spinfree(false).compute());
    auto result = wick.spinfree(false).compute();
    REQUIRE(result->is<Sum>());
    REQUIRE(result->size() == 36);
  }

  // two N-nonconserving operators
  {
    auto opseq = FNOperatorSeq(
        {FNOperator({L"p_1", L"p_2", L"p_3"}, {L"p_4", L"p_5"}, V),
         FNOperator({L"p_7", L"p_8"}, {L"p_10", L"p_11", L"p_12"}, V)});
    auto wick = FWickTheorem{opseq};
    REQUIRE_NOTHROW(wick.spinfree(false).compute());
    auto result = wick.spinfree(false).compute();
    REQUIRE(result->is<Sum>());
    REQUIRE(result->size() == 12);
  }

  // more N-nonconserving operators
  {
    auto input = ex<FNOperator>(WstrList{L"i_1"}, WstrList{L"a_3", L"a_4"}, V) *
                 (ex<Constant>(rational{1, 4}) *
                  ex<Tensor>(L"g", WstrList{L"p_1", L"p_2"},
                             WstrList{L"p_3", L"p_4"}, Symmetry::antisymm) *
                  ex<FNOperator>(WstrList{L"p_1", L"p_2"},
                                 WstrList{L"p_3", L"p_4"}, V)) *
                 ex<FNOperator>(WstrList{L"a_2"}, WstrList{}, V);
    auto wick = FWickTheorem{input};
    wick.spinfree(false)
        .set_external_indices(IndexList{L"i_1", L"a_3", L"a_4", L"a_2"})
        .use_topology(true);
    ExprPtr result;
    REQUIRE_NOTHROW(result = wick.compute());
    // std::wcout << "result = " << to_latex(result) << std::endl;
    REQUIRE(to_latex(result) ==
            L"{{{-1}}{\\bar{g}^{{i_1}{a_2}}_{{a_3}{a_4}}}}");
  }

  // odd number of ops -> full contraction is 0
  {
    auto opseq = FNOperatorSeq(
        {FNOperator({L"p_1", L"p_2"}, {L"p_4", L"p_5"}, V),
         FNOperator({L"p_7", L"p_8"}, {L"p_10", L"p_11", L"p_12"}, V)});
    auto wick = FWickTheorem{opseq};
    REQUIRE_NOTHROW(wick.spinfree(false).compute());
    auto result = wick.spinfree(false).compute();
    REQUIRE(result->is<Constant>());
    REQUIRE(result->as<Constant>().value<int>() == 0);
  }

  // 4-body ^ 4-body
  SEQUANT_PROFILE_SINGLE(
      "wick(4^4)",
      {
        auto opseq = FNOperatorSeq(
            {FNOperator({L"p_1", L"p_2", L"p_3", L"p_4"},
                        {L"p_5", L"p_6", L"p_7", L"p_8"}, V),
             FNOperator({L"p_21", L"p_22", L"p_23", L"p_24"},
                        {L"p_25", L"p_26", L"p_27", L"p_28"}, V)});
        auto wick = FWickTheorem{opseq};
        auto result = wick.spinfree(false).compute(true);
        REQUIRE(result->is<Constant>());
        REQUIRE(result->as<Constant>().value<int>() == 576);
      })

  // three general 1-body operators
  {
    auto opseq = FNOperatorSeq({FNOperator({L"p_1"}, {L"p_2"}, V),
                                FNOperator({L"p_3"}, {L"p_4"}, V),
                                FNOperator({L"p_5"}, {L"p_6"}, V)});
    auto wick = FWickTheorem{opseq};
    REQUIRE_NOTHROW(wick.spinfree(false).compute());
    auto result = wick.spinfree(false).compute();
    REQUIRE(result->is<Sum>());
    REQUIRE(result->size() == 2);
  }

  // 4 general 1-body operators
  {
    auto opseq = FNOperatorSeq(
        {FNOperator({L"p_1"}, {L"p_2"}, V), FNOperator({L"p_3"}, {L"p_4"}, V),
         FNOperator({L"p_5"}, {L"p_6"}, V), FNOperator({L"p_7"}, {L"p_8"}, V)});
    auto ext_indices = make_indices<std::vector<Index>>(WstrList{
        L"p_1", L"p_2", L"p_3", L"p_4", L"p_5", L"p_6", L"p_7", L"p_8"});
    auto wick1 = FWickTheorem{opseq};
    auto result1 =
        wick1.set_external_indices(ext_indices).spinfree(false).compute();
    REQUIRE(result1->is<Sum>());
    REQUIRE(result1->size() == 9);
    auto wick2 = FWickTheorem{opseq};
    auto result2 = wick2.set_external_indices(ext_indices)
                       .set_nop_connections({{1, 2}, {1, 3}})
                       .spinfree(false)
                       .compute();
    REQUIRE(result2->is<Sum>());
    REQUIRE(result2->size() == 2);
  }

  // 4-body ^ 2-body ^ 2-body
  {
    auto opseq =
        FNOperatorSeq({FNOperator({L"p_1", L"p_2", L"p_3", L"p_4"},
                                  {L"p_5", L"p_6", L"p_7", L"p_8"}, V),
                       FNOperator({L"p_9", L"p_10"}, {L"p_11", L"p_12"}, V),
                       FNOperator({L"p_13", L"p_14"}, {L"p_15", L"p_16"}, V)});
    auto wick = FWickTheorem{opseq};
    auto result = wick.spinfree(false).compute();
    REQUIRE(result->is<Sum>());
    REQUIRE(result->size() == 576);
  }

  // 2-body ^ 2-body ^ 2-body
  {
    auto opseq =
        FNOperatorSeq({FNOperator({L"p_1", L"p_2"}, {L"p_5", L"p_6"}, V),
                       FNOperator({L"p_9", L"p_10"}, {L"p_11", L"p_12"}, V),
                       FNOperator({L"p_17", L"p_18"}, {L"p_19", L"p_20"}, V)});
    auto wick = FWickTheorem{opseq};
    auto result = wick.spinfree(false).compute();
    REQUIRE(result->is<Sum>());
    REQUIRE(result->size() == 80);
  }

  // 2-body ^ 2-body ^ 2-body ^ 2-body
  SEQUANT_PROFILE_SINGLE("wick(2^2^2^2)", {
    auto opseq =
        FNOperatorSeq({FNOperator({L"p_1", L"p_2"}, {L"p_5", L"p_6"}, V),
                       FNOperator({L"p_9", L"p_10"}, {L"p_11", L"p_12"}, V),
                       FNOperator({L"p_13", L"p_14"}, {L"p_15", L"p_16"}, V),
                       FNOperator({L"p_17", L"p_18"}, {L"p_19", L"p_20"}, V)});
    auto wick = FWickTheorem{opseq};
    auto result = wick.spinfree(false).compute(true);
    REQUIRE(result->is<Constant>());
    REQUIRE(result->as<Constant>().value<int>() == 4752);
  })

#ifndef SEQUANT_SKIP_LONG_TESTS
  // 4-body ^ 2-body ^ 2-body ^ 2-body
  SEQUANT_PROFILE_SINGLE("wick(4^2^2^2)", {
    auto opseq =
        FNOperatorSeq({FNOperator({L"p_1", L"p_2", L"p_3", L"p_4"},
                                  {L"p_5", L"p_6", L"p_7", L"p_8"}, V),
                       FNOperator({L"p_9", L"p_10"}, {L"p_11", L"p_12"}, V),
                       FNOperator({L"p_13", L"p_14"}, {L"p_15", L"p_16"}, V),
                       FNOperator({L"p_17", L"p_18"}, {L"p_19", L"p_20"}, V)});
    auto wick = FWickTheorem{opseq};
    auto result = wick.spinfree(false).use_topology(true).compute(true);
    REQUIRE(result->is<Constant>());
    REQUIRE(result->as<Constant>().value<int>() == 2088);
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
    auto result = wick.spinfree(false).use_topology(true).compute(true);
    REQUIRE(result->is<Constant>());
    REQUIRE(result->as<Constant>().value<int>() == 694);
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
    auto result = wick.spinfree(false).use_topology(true).compute(true);
    REQUIRE(result->is<Constant>());
    REQUIRE(result->as<Constant>().value<int>() == 28);
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
    auto result = wick.spinfree(false).use_topology(true).compute(true);
    REQUIRE(result->is<Constant>());
    REQUIRE(result->as<Constant>().value<int>() == 70);
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
      auto result = wick.spinfree(false).use_topology(true).compute(true);
    }
#endif
}  // SECTION("fermi vacuum")

auto print = [](const auto& lead, const auto& expr) {
  std::wcout << lead << to_latex(expr) << std::endl;
};

SECTION("Expression Reduction") {
  constexpr Vacuum V = Vacuum::SingleProduct;

  // 2-body ^ 2-body
  SEQUANT_PROFILE_SINGLE("wick(H2*T2)", {
    auto opseq =
        FNOperatorSeq({FNOperator({L"p_1", L"p_2"}, {L"p_3", L"p_4"}, V),
                       FNOperator({L"a_4", L"a_5"}, {L"i_4", L"i_5"}, V)});
    auto wick = FWickTheorem{opseq};
    auto wick_result = wick.spinfree(false).compute();
    REQUIRE(wick_result->is<Sum>());
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
            L"\\bigl({{\\bar{g}^{{p_3}{p_4}}_{{p_1}{p_2}}}{\\bar{t}^{{i_4}{i_"
            L"5}}_{{a_4}{a_"
            L"5}}}{s^{{p_1}}_{{i_5}}}{s^{{p_2}}_{{i_4}}}{s^{{a_4}}_{{p_4}}}{"
            L"s^{{a_5}}_{{p_3}}}} - {"
            L"{\\bar{g}^{{p_3}{p_4}}_{{p_1}{p_2}}}{\\bar{t}^{{i_4}{i_5}}_{{a_"
            L"4}{a_5}}}{s^{{"
            L"p_1}}_{{i_5}}}{s^{{p_2}}_{{i_4}}}{s^{{a_5}}_{{p_4}}}{s^{{a_4}}_"
            L"{{p_3}}}} - {"
            L"{\\bar{g}^{{p_3}{p_4}}_{{p_1}{p_2}}}{\\bar{t}^{{i_4}{i_5}}_{{a_"
            L"4}{a_5}}}{s^{{"
            L"p_1}}_{{i_4}}}{s^{{p_2}}_{{i_5}}}{s^{{a_4}}_{{p_4}}}{s^{{a_5}}_"
            L"{{p_3}}}} + "
            L"{{\\bar{g}^{{p_3}{p_4}}_{{p_1}{p_2}}}{\\bar{t}^{{i_4}{i_5}}_{{"
            L"a_4}{a_5}}}{s^{"
            L"{p_1}}_{{i_4}}}{s^{{p_2}}_{{i_5}}}{s^{{a_5}}_{{p_4}}}{s^{{a_4}}"
            L"_{{p_3}}}}\\bigr) }");
    wick.reduce(wick_result_2);
    rapid_simplify(wick_result_2);
    TensorCanonicalizer::register_instance(
        std::make_shared<DefaultTensorCanonicalizer>());
    canonicalize(wick_result_2);
    rapid_simplify(wick_result_2);

    std::wcout << L"H2*T2 = " << to_latex(wick_result_2) << std::endl;
    std::wcout << L"H2*T2 = " << to_wolfram(wick_result_2) << std::endl;
    REQUIRE(to_latex(wick_result_2) ==
            L"{{{4}}"
            L"{\\bar{g}^{{a_1}{a_2}}_{{i_1}{i_2}}}{\\bar{t}^{{i_1}{i_2}}_{{a_"
            L"1}{a_2}}}}");
  });

  // 2-body ^ 1-body ^ 1-body, with/without using topology
  SEQUANT_PROFILE_SINGLE("wick(H2*T1*T1)", {
    for (auto&& use_nop_partitions : {false}) {
      for (auto&& use_op_partitions : {true, false}) {
        std::wostringstream oss;
        oss << "use_op_partitions=" << use_op_partitions << "}: H2*T1*T1 = ";

        auto opseq =
            FNOperatorSeq({FNOperator({L"p_1", L"p_2"}, {L"p_3", L"p_4"}, V),
                           FNOperator({L"a_4"}, {L"i_4"}, V),
                           FNOperator({L"a_5"}, {L"i_5"}, V)});
        auto wick = FWickTheorem{opseq};
        wick.spinfree(false).use_topology(use_nop_partitions ||
                                          use_op_partitions);
        // if (use_nop_partitions) wick.set_nop_partitions({{1, 2}});
        if (use_op_partitions) wick.set_op_partitions({{0, 1}, {2, 3}});
        auto wick_result = wick.compute();
        // print(oss.str() + L" (nopseq only) ", wick_result);
        if (use_op_partitions) {
          REQUIRE(wick_result->is<Product>());
          REQUIRE(wick_result->size() == 4 /* factors */);
        } else {
          REQUIRE(wick_result->is<Sum>());
          REQUIRE(wick_result->size() == 4 /* summands */);
        }

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

        // print(oss.str(), wick_result_2);
        REQUIRE(wick_result_2->size() == 3 /* factors */);
        REQUIRE(to_latex(wick_result_2) ==
                L"{{{4}}"
                L"{\\bar{g}^{{a_1}{a_2}}_{{i_1}{i_2}}}{t^{{i_1}}_{{a_1}}}{"
                L"t^{{i_"
                L"2}}_{{a_"
                L"2}}}}");
      }  // use_op_partitions
    }    // use_nop_partitions
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
    auto wick_result = wick.spinfree(false).compute();
    REQUIRE(wick_result->is<Sum>());
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

    //    std::wcout << L"P2*H1*T2(PNO) = " << to_latex_align(wick_result_2)
    //               << std::endl;
    // it appears that the two terms are swapped when using gcc 8 on linux
    // TODO investigate why sum canonicalization seems to produce
    // platform-dependent results.
    //      REQUIRE(to_latex(wick_result_2) ==
    //              L"{ \\bigl( - {{{8}}"
    //              L"{A^{{a_1^{{i_1}{i_2}}}{a_2^{{i_1}{i_2}}}}_{{i_1}{i_2}}}{f^{{a_"
    //              L"3^{{i_1}{i_2}}}}_{{a_1^{{i_1}{i_2}}}}}{t^{{i_1}{i_2}}_{{a_2^{{"
    //              L"i_1}{i_2}}}{a_3^{{i_1}{i_2}}}}}} + {{{8}}"
    //              L"{A^{{a_1^{{i_1}{i_2}}}{a_2^{{i_1}{i_2}}}}_{{i_1}{i_2}}}{f^{{i_"
    //              L"1}}_{{i_3}}}{t^{{i_2}{i_3}}_{{a_3^{{i_2}{i_3}}}{a_4^{{i_2}{i_3}"
    //              L"}}}}{s^{{a_3^{{i_2}{i_3}}}}_{{a_1^{{i_1}{i_2}}}}}{s^{{a_4^{{i_"
    //              L"2}{i_3}}}}_{{a_2^{{i_1}{i_2}}}}}}\\bigr) }");
  });

  // 2=body ^ 2-body ^ 2-body ^ 2-body with dependent (PNO) indices
  SEQUANT_PROFILE_SINGLE("wick(P2*H2*T2*T2)", {
    for (auto&& use_nop_partitions : {false}) {
      for (auto&& use_op_partitions : {true, false}) {
        std::wostringstream oss;
        oss << "use_{nop,op}_partitions={" << use_nop_partitions << ","
            << use_op_partitions << "}: P2*H2*T2*T2(PNO) = ";

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
        wick.spinfree(false)
            .set_nop_connections({{1, 2}, {1, 3}})
            .use_topology(true);

        if (use_nop_partitions) wick.set_nop_partitions({{2, 3}});
        if (use_op_partitions)
          wick.set_op_partitions({{0, 1},
                                  {2, 3},
                                  {4, 5},
                                  {6, 7},
                                  {8, 9},
                                  {10, 11},
                                  {12, 13},
                                  {14, 15}});
        auto wick_result = wick.compute();
        REQUIRE(wick_result->is<Sum>());
        if (use_op_partitions) {
          REQUIRE(wick_result->size() == 7);
        } else {
          REQUIRE(wick_result->size() == 544);
        }
        std::wcout << oss.str() << L" (nop-only) "
                   << to_latex_align(wick_result, 0, 1) << std::endl;

        // multiply tensor factors and expand
        auto wick_result_2 =
            ex<Constant>(rational{1, 256}) *
            ex<Tensor>(L"A", IndexList{L"i_1", L"i_2"},
                       IndexList{{L"a_1", {L"i_1", L"i_2"}},
                                 {L"a_2", {L"i_1", L"i_2"}}},
                       Symmetry::antisymm) *
            ex<Tensor>(L"g", WstrList{L"p_1", L"p_2"}, WstrList{L"p_3", L"p_4"},
                       Symmetry::antisymm) *
            ex<Tensor>(L"t",
                       IndexList{{L"a_3", {L"i_3", L"i_4"}},
                                 {L"a_4", {L"i_3", L"i_4"}}},
                       IndexList{L"i_3", L"i_4"}, Symmetry::antisymm) *
            ex<Tensor>(L"t",
                       IndexList{{L"a_5", {L"i_5", L"i_6"}},
                                 {L"a_6", {L"i_5", L"i_6"}}},
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

        //        std::wcout << oss.str() << to_latex_align(wick_result_2, 20)
        //                   << std::endl;
        REQUIRE(wick_result_2->is<Sum>());
        REQUIRE(wick_result_2->size() == 4);
      }  // use_op_partitions
    }    // use_nop_partitions
  });

#if 1
  // 3-body ^ 2-body ^ 2-body ^ 3-body
  SEQUANT_PROFILE_SINGLE("wick(P3*H2*T2*T3)", {
    constexpr bool connected_only = true;
    constexpr bool topology = true;
    auto P3 = ex<Constant>(rational{1, 36}) *
              ex<Tensor>(L"A", WstrList{L"i_1", L"i_2", L"i_3"},
                         WstrList{L"a_1", L"a_2", L"a_3"}, Symmetry::antisymm) *
              ex<FNOperator>(WstrList{L"i_1", L"i_2", L"i_3"},
                             WstrList{L"a_1", L"a_2", L"a_3"}, V);
    auto H2 =
        ex<Constant>(rational{1, 4}) *
        ex<Tensor>(L"g", WstrList{L"p_1", L"p_2"}, WstrList{L"p_3", L"p_4"},
                   Symmetry::antisymm) *
        ex<FNOperator>(WstrList{L"p_1", L"p_2"}, WstrList{L"p_3", L"p_4"}, V);
    auto T2 =
        ex<Constant>(rational{1, 4}) *
        ex<Tensor>(L"t", WstrList{L"a_4", L"a_5"}, WstrList{L"i_4", L"i_5"},
                   Symmetry::antisymm) *
        ex<FNOperator>(WstrList{L"a_4", L"a_5"}, WstrList{L"i_4", L"i_5"}, V);
    auto T3 = ex<Constant>(rational{1, 36}) *
              ex<Tensor>(L"t", WstrList{L"a_6", L"a_7", L"a_8"},
                         WstrList{L"i_6", L"i_7", L"i_8"}, Symmetry::antisymm) *
              ex<FNOperator>(WstrList{L"a_6", L"a_7", L"a_8"},
                             WstrList{L"i_6", L"i_7", L"i_8"}, V);
    FWickTheorem wick{P3 * H2 * T2 * T3};
    wick.spinfree(false).use_topology(topology);
    if (connected_only) wick.set_nop_connections({{1, 2}, {1, 3}});
    auto wick_result = wick.compute();

    std::wcout << "P3*H2*T2*T3 = " << to_latex_align(wick_result, 20)
               << std::endl;
    REQUIRE(wick_result->is<Sum>());
    REQUIRE(
        wick_result->size() ==
        (connected_only ? 7 : 9));  // 9 = 2 disconnected + 7 connected terms
  });
#endif
}

}  // TEST_CASE("WickTheorem")
#endif
