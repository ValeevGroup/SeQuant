//
// Created by Eduard Valeyev on 3/20/18.
//

#include "catch.hpp"

#include <iostream>
#include "../../src/SeQuant2/op.hpp"

TEST_CASE("Op", "[elements]") {

  using namespace sequant2;

  SECTION("constructors") {
    REQUIRE_NOTHROW(FOp{});

    REQUIRE_NOTHROW(FOp(Index(L"i_1"), Action::create));
    FOp o1(Index(L"i_1"), Action::create);
    REQUIRE(o1.statistics == Statistics::FermiDirac);
    REQUIRE(o1.index() == Index(L"i_1"));
    REQUIRE(o1.action() == Action::create);

    REQUIRE_NOTHROW(bann(Index(L"i_2")));
    auto o2 = bann(Index(L"i_2"));
    REQUIRE(o2.statistics == Statistics::BoseEinstein);
    REQUIRE(o2.index() == Index(L"i_2"));
    REQUIRE(o2.action() == Action::annihilate);

    REQUIRE_NOTHROW(bcre(L"i_3"));
    auto o3 = bcre(L"i_3");
    REQUIRE(o3.statistics == Statistics::BoseEinstein);
    REQUIRE(o3.index() == Index(L"i_3"));
    REQUIRE(o3.action() == Action::create);

    REQUIRE_NOTHROW(fann(L"i_4", {L"a_1", L"a_2"}));
    auto o4 = fann(L"i_4", {L"a_1", L"a_2"});
    REQUIRE(o4.statistics == Statistics::FermiDirac);
    REQUIRE(o4.index() == Index(L"i_4", {L"a_1", L"a_2"}));
    REQUIRE(o4.action() == Action::annihilate);

    REQUIRE_NOTHROW(FOperator{fcre(L"i_1"), fann(L"i_1")});
    auto oper1 = FOperator{fcre(L"i_1"), fann(L"i_1")};
    REQUIRE(oper1.statistics == Statistics::FermiDirac);
    REQUIRE(oper1.size() == 2);
    REQUIRE(oper1[0] == fcre(L"i_1"));
    REQUIRE(oper1[1] == fann(L"i_1"));
  }

  SECTION("latex") {
    FOp o1(Index(L"i_1"), Action::create);
    REQUIRE(to_latex(o1) == L"{a^{\\dagger}_{i_1}}");

    BOp o2(Index(L"a_1"), Action::create);
    REQUIRE(to_latex(o2) == L"{b^{\\dagger}_{a_1}}");

    auto oper1 = FOperator{fcre(L"i_1"), fann(L"i_1")};
    REQUIRE(to_latex(oper1) == L"{{a^{\\dagger}_{i_1}}{a_{i_1}}}");
  }

}  // TEST_CASE("Op")