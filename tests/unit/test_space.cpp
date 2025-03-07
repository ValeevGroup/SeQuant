//
// Created by Eduard Valeyev on 3/20/18.
//

#include <catch2/catch_test_macros.hpp>

#include "catch2_sequant.hpp"

#include <SeQuant/core/space.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>
#include "SeQuant/domain/mbpt/lcao.hpp"

TEST_CASE("IndexSpace", "[elements]") {
  using namespace sequant;

  SECTION("registry synopsis") {
    auto sr_isr = sequant::mbpt::make_sr_spaces();
    REQUIRE_NOTHROW(sr_isr->retrieve(L"i"));
    REQUIRE_NOTHROW(sr_isr->remove(L"i"));
    IndexSpace active_occupied(L"i", 0b0010);
    REQUIRE_NOTHROW(sr_isr->add(active_occupied));
    REQUIRE_THROWS(sr_isr->add(
        active_occupied));  // cannot add a space that is already present
    IndexSpace jactive_occupied(L"j", 0b0010);
    REQUIRE_THROWS(sr_isr->add(
        jactive_occupied));  // cannot add a space with duplicate attr

    // can use bytestrings too
    REQUIRE_NOTHROW(sr_isr->retrieve("a"));  // N.B. string
    REQUIRE_NOTHROW(sr_isr->remove('a'));    // N.B. char
  }

  SECTION("registry construction") {
    auto isr = std::make_shared<IndexSpaceRegistry>();

    // similar make_sr_spaces, but no spin AND occupied and unoccupied bits are
    // NOT contiguous!
    REQUIRE_NOTHROW(
        isr->add(L"o", 0b0001, 3)       // approximate size
            .add("i", 0b0100, is_hole)  // N.B. narrow string
            .add(L'a', 0b0010, is_particle, QuantumNumbersAttr{},
                 50)  // N.B. wchar_t + explicit quantum numbers + size
            .add('g', 0b1000)
            .add_union(L"m", {L"o", L"i"}, is_vacuum_occupied,
                       is_reference_occupied)
            .add_union(L"e", {L"a", L"g"})
            .add_unIon(L"p", {L"m", L"e"}, is_complete)  // N.B. unIon
    );
    REQUIRE(isr->retrieve(L"o").approximate_size() == 3);
    REQUIRE(isr->retrieve(L"m").type().to_int32() == 0b0101);
    REQUIRE(isr->retrieve(L"m").approximate_size() ==
            isr->retrieve(L"o").approximate_size() +
                isr->retrieve(L"i").approximate_size());
    REQUIRE(isr->retrieve(L"p").type().to_int32() == 0b1111);
    REQUIRE(isr->retrieve(L"p").approximate_size() ==
            isr->retrieve(L"m").approximate_size() +
                isr->retrieve(L"e").approximate_size());
    REQUIRE(isr->retrieve(L"p").approximate_size() == 73);
    REQUIRE_NOTHROW(
        isr->add_union(L"iag", {L"i", L"a", L"g"})
            .add_union(L"oia", {"o", "i", "a"})  // N.B. narrow strings
            .add_intersection(L"x", {L"oia", L"iag"}));
    REQUIRE(isr->retrieve(L"x").type().to_int32() == 0b0110);
    REQUIRE(isr->retrieve(L"x").approximate_size() ==
            isr->retrieve(L"i").approximate_size() +
                isr->retrieve(L"a").approximate_size());
    REQUIRE(isr->vacuum_occupied_space(IndexSpace::QuantumNumbers::null) ==
            isr->retrieve(L"m"));
    REQUIRE(isr->vacuum_unoccupied_space(IndexSpace::QuantumNumbers::null) ==
            isr->retrieve(L"e"));
  }

  SECTION("equality") {
    auto sr_isr = sequant::mbpt::make_sr_spaces();
    REQUIRE(sr_isr->retrieve(L"i") == sr_isr->retrieve(L"i"));
    REQUIRE(sr_isr->retrieve(L"i") != sr_isr->retrieve(L"a"));
  }

  SECTION("ordering") {
    auto sr_isr = sequant::mbpt::make_sr_spaces();
    REQUIRE(!(sr_isr->retrieve(L"i") < sr_isr->retrieve(L"i")));
    REQUIRE(sr_isr->retrieve(L"i") < sr_isr->retrieve(L"a"));
    REQUIRE(!(sr_isr->retrieve(L"a") < sr_isr->retrieve(L"i")));
    REQUIRE(sr_isr->retrieve(L"i") < sr_isr->retrieve(L"m"));
    REQUIRE(!(sr_isr->retrieve(L"m") < sr_isr->retrieve(L"m")));
    REQUIRE(sr_isr->retrieve(L"m") < sr_isr->retrieve(L"a"));
    REQUIRE(sr_isr->retrieve(L"m") < sr_isr->retrieve(L"p"));
    REQUIRE(sr_isr->retrieve(L"a") < sr_isr->retrieve(L"p"));

    // test ordering with quantum numbers
    {
      auto i = IndexSpace(L"i", 0b01);
      auto a = IndexSpace(L"a", 0b10);
      auto iA = IndexSpace(L"i", 0b01, 0b01);
      auto iB = IndexSpace(L"i", 0b01, 0b10);
      auto aA = IndexSpace(L"a", 0b10, 0b01);
      auto aB = IndexSpace(L"a", 0b10, 0b10);

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
  }

  SECTION("set operations") {
    auto isr = sequant::mbpt::make_F12_sr_spaces();
    REQUIRE(isr->retrieve(L"i") ==
            isr->intersection(isr->retrieve(L"i"), isr->retrieve(L"p")));
    REQUIRE(!isr->intersection(isr->retrieve(L"p↑"), isr->retrieve(L"p↓")));
    REQUIRE(isr->intersection(isr->retrieve(L"i↑"), isr->retrieve(L"p")) ==
            isr->retrieve(L"i↑"));
    REQUIRE(!isr->intersection(isr->retrieve(L"a"), isr->retrieve(L"i")));
    REQUIRE(!isr->intersection(isr->retrieve(L"a"), isr->retrieve(L"α'")));

    REQUIRE(isr->retrieve(L"κ") ==
            isr->unIon(isr->retrieve(L"p"), isr->retrieve(L"α'")));

    REQUIRE(includes(isr->retrieve(L"κ"), isr->retrieve(L"m")));
    REQUIRE(!includes(isr->retrieve(L"m"), isr->retrieve(L"κ")));
    REQUIRE(includes(isr->retrieve(L"α"), isr->retrieve(L"a")));

    REQUIRE(isr->valid_intersection(isr->retrieve(L"i"), isr->retrieve(L"p")));

    REQUIRE(isr->valid_unIon(isr->retrieve(L"i"), isr->retrieve(L"a")));
    REQUIRE(isr->valid_unIon(isr->retrieve(L"i↑"), isr->retrieve(L"i↓")));
    REQUIRE(!isr->valid_unIon(isr->retrieve(L"i↑"), isr->retrieve(L"i↑")));
    REQUIRE(!isr->valid_unIon(isr->retrieve(L"p"), isr->retrieve(L"a")));
    REQUIRE(!isr->valid_unIon(isr->retrieve(L"p↑"), isr->retrieve(L"a↓")));
  }

  SECTION("occupancy_validation") {
    auto sr_isr = sequant::mbpt::make_sr_spaces();
    auto mr_isr = sequant::mbpt::make_mr_spaces();

    REQUIRE(sr_isr->is_pure_occupied(sr_isr->retrieve(L"i")));
    REQUIRE(sr_isr->is_pure_unoccupied(sr_isr->retrieve(L"a")));

    REQUIRE(mr_isr->is_pure_occupied(mr_isr->retrieve(L"i")));
    REQUIRE(mr_isr->is_pure_unoccupied(mr_isr->retrieve(L"a")));
    REQUIRE(!mr_isr->is_pure_occupied(mr_isr->retrieve(L"I")));
    // REQUIRE(!mr_isr->is_pure_unoccupied(mr_isr->retrieve(L"E")));

    REQUIRE(sr_isr->contains_occupied(sr_isr->retrieve(L"i")));
    REQUIRE(sr_isr->contains_unoccupied(sr_isr->retrieve(L"a")));

    REQUIRE(mr_isr->contains_occupied(mr_isr->retrieve(L"M")));
    REQUIRE(mr_isr->contains_unoccupied(mr_isr->retrieve(L"E")));
    REQUIRE(!mr_isr->contains_occupied(mr_isr->retrieve(L"a")));
    REQUIRE(!mr_isr->contains_unoccupied(mr_isr->retrieve(L"i")));
  }

  SECTION("base_space") {
    auto isr = sequant::mbpt::make_F12_sr_spaces();
    const auto& f12_base_space_types = isr->base_space_types();
    REQUIRE(f12_base_space_types.size() == 5);
    REQUIRE(f12_base_space_types[0] == 0b00001);
    REQUIRE(f12_base_space_types[1] == 0b00010);
    REQUIRE(f12_base_space_types[2] == 0b00100);
    REQUIRE(f12_base_space_types[3] == 0b01000);
    REQUIRE(f12_base_space_types[4] == 0b10000);
    const auto& f12_base_spaces = isr->base_spaces();
    auto f12_base_spaces_sf = f12_base_spaces |
                              ranges::views::filter([](const auto& space) {
                                return space.qns() == mbpt::Spin::any;
                              }) |
                              ranges::to_vector;
    REQUIRE(f12_base_spaces_sf[0].base_key() == L"o");
    REQUIRE(f12_base_spaces_sf[1].base_key() == L"i");
    REQUIRE(f12_base_spaces_sf[2].base_key() == L"a");
    REQUIRE(f12_base_spaces_sf[3].base_key() == L"g");
    REQUIRE(f12_base_spaces_sf[4].base_key() == L"α'");
  }

  SECTION("AO spaces") {
    auto isr = sequant::mbpt::make_min_sr_spaces();
    REQUIRE_NOTHROW(mbpt::add_ao_spaces(isr));

    // OBS AO space ...
    REQUIRE_NOTHROW(isr->retrieve(L"μ"));
    auto μ = isr->retrieve(L"μ");
    REQUIRE(bitset_t(μ.qns()) & bitset_t(mbpt::LCAOQNS::ao));
    REQUIRE(!(bitset_t(μ.qns()) & bitset_t(mbpt::LCAOQNS::lcao)));

    /// has same type as occupied space but differ by QNS
    REQUIRE_NOTHROW(isr->retrieve(L"i"));
    auto i = isr->retrieve(L"i");
    REQUIRE(μ.type() == i.type());
    REQUIRE(μ.qns() != i.qns());
    REQUIRE(μ != i);
  }
}
