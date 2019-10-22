//
// Created by Eduard Valeyev on 3/23/18.
//

#include "catch.hpp"

#include <iostream>
#include "SeQuant/core/wick.hpp"

TEST_CASE("Tensor", "[elements]") {

  using namespace sequant;

  SECTION("constructors") {

    REQUIRE_NOTHROW(Tensor{});
    auto t1 = Tensor{};
    REQUIRE(t1.bra_rank() == 0);
    REQUIRE(t1.ket_rank() == 0);
    REQUIRE(t1.rank() == 0);
    REQUIRE(t1.symmetry() == Symmetry::invalid);
    REQUIRE(t1.braket_symmetry() == BraKetSymmetry::invalid);
    REQUIRE(t1.label() == L"");

    REQUIRE_NOTHROW(Tensor(L"F", {L"i_1"}, {L"i_1"}));
    auto t2 = Tensor(L"F", {L"i_1"}, {L"i_1"});
    REQUIRE(t2.bra_rank() == 1);
    REQUIRE(t2.ket_rank() == 1);
    REQUIRE(t2.rank() == 1);
    REQUIRE(t2.symmetry() == Symmetry::nonsymm);
    REQUIRE(t2.braket_symmetry() == BraKetSymmetry::conjugate);
    REQUIRE(t2.label() == L"F");

    REQUIRE_NOTHROW(Tensor(L"N", {L"i_1"}, {}));
    auto t3 = Tensor(L"N", {L"i_1"}, {});
    REQUIRE(t3.bra_rank() == 1);
    REQUIRE(t3.ket_rank() == 0);
    REQUIRE_THROWS(t3.rank());
    REQUIRE(t3.symmetry() == Symmetry::nonsymm);
    REQUIRE(t3.braket_symmetry() == BraKetSymmetry::conjugate);
    REQUIRE(t3.label() == L"N");

    REQUIRE_NOTHROW(Tensor(L"g", {Index{L"i_1"}, Index{L"i_2"}}, {Index{L"i_3"}, Index{L"i_4"}}, Symmetry::antisymm, BraKetSymmetry::symm));
    auto t4 = Tensor(L"g", {Index{L"i_1"}, Index{L"i_2"}}, {Index{L"i_3"}, Index{L"i_4"}}, Symmetry::antisymm, BraKetSymmetry::symm);
    REQUIRE(t4.bra_rank() == 2);
    REQUIRE(t4.ket_rank() == 2);
    REQUIRE(t4.rank() == 2);
    REQUIRE(t4.symmetry() == Symmetry::antisymm);
    REQUIRE(t4.braket_symmetry() == BraKetSymmetry::symm);
    REQUIRE(t4.label() == L"g");
  }  // SECTION("constructors")

  SECTION("index transformation") {
    auto t = Tensor(L"g", {Index{L"i_1"}, Index{L"i_2"}},
                    {Index{L"i_3"}, Index{L"i_4"}}, Symmetry::antisymm);
    std::map<Index, Index> idxmap = {{Index{L"i_1"}, Index{L"i_2"}},
                                     {Index{L"i_2"}, Index{L"i_1"}}};
    REQUIRE(t.transform_indices(idxmap, true));
    REQUIRE(t.bra()[0].tag().has_value());
    REQUIRE(t.bra()[1].tag().has_value());
    REQUIRE(!t.ket()[0].tag().has_value());
    REQUIRE(!t.ket()[1].tag().has_value());
    REQUIRE(t == Tensor(L"g", {Index{L"i_2"}, Index{L"i_1"}},
                        {Index{L"i_3"}, Index{L"i_4"}}, Symmetry::antisymm));
    // tagged indices are protected, so no replacements the second goaround
    REQUIRE(!t.transform_indices(idxmap, true));
    t.reset_tags();
    REQUIRE(t.transform_indices(idxmap, true));
    REQUIRE(t == Tensor(L"g", {Index{L"i_1"}, Index{L"i_2"}},
                        {Index{L"i_3"}, Index{L"i_4"}}, Symmetry::antisymm));
    t.reset_tags();
    REQUIRE(!t.bra()[0].tag().has_value());
    REQUIRE(!t.bra()[1].tag().has_value());
    REQUIRE(!t.ket()[0].tag().has_value());
    REQUIRE(!t.ket()[1].tag().has_value());
  }  // SECTION("index transformation")

  SECTION("latex") {
    auto t1 = Tensor(L"F", {L"i_1"}, {L"i_2"});
    REQUIRE(to_latex(t1) == L"{F^{{i_2}}_{{i_1}}}");

    auto h1 = ex<Tensor>(L"F", WstrList{L"i_1"}, WstrList{L"i_2"}) * ex<FNOperator>(WstrList{L"i_1"}, WstrList{L"i_2"});
    REQUIRE(to_latex(h1) == L"{{F^{{i_2}}_{{i_1}}}{\\tilde{a}^{{i_1}}_{{i_2}}}}");

  }  // SECTION("latex")

}  // TEST_CASE("Tensor")
