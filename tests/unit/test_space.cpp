//
// Created by Eduard Valeyev on 3/20/18.
//

#include <catch2/catch_test_macros.hpp>

#include <SeQuant/core/space.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>

TEST_CASE("IndexSpace", "[elements]") {
  using namespace sequant;
  SECTION("registry_core_functionality") {
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
    REQUIRE(isr->nullspace ==
            isr->intersection(isr->retrieve(L"a"), isr->retrieve(L"i")));
    REQUIRE(isr->nullspace ==
            isr->intersection(isr->retrieve(L"a"), isr->retrieve(L"α'")));

    REQUIRE(isr->retrieve(L"κ") ==
            isr->unIon(isr->retrieve(L"p"), isr->retrieve(L"α'")));

    REQUIRE(includes(isr->retrieve(L"κ"), isr->retrieve(L"m")));
    REQUIRE(!includes(isr->retrieve(L"m"), isr->retrieve(L"κ")));
    REQUIRE(includes(isr->retrieve(L"α"), isr->retrieve(L"a")));

    REQUIRE(isr->has_non_overlapping_spaces(isr->retrieve(L"i"),
                                            isr->retrieve(L"p")));
    REQUIRE(isr->has_non_overlapping_spaces(isr->retrieve(L"i"),
                                            isr->retrieve(L"a")));
    REQUIRE(!isr->has_non_overlapping_spaces(isr->retrieve(L"i"),
                                             isr->retrieve(L"i")));

    REQUIRE(isr->non_overlapping_spaces(isr->retrieve(L"g"),
                                        isr->retrieve(L"α"))[0] ==
            isr->retrieve(L"a"));
    REQUIRE(isr->non_overlapping_spaces(isr->retrieve(L"g"),
                                        isr->retrieve(L"α"))[1] ==
            isr->retrieve(L"α'"));

    auto union_func = [](int32_t a, int32_t b) { return a + b; };
    REQUIRE(!isr->valid_bitop(isr->retrieve(L"i"), isr->retrieve(L"α'"),
                              union_func));
    REQUIRE(
        isr->valid_bitop(isr->retrieve(L"i"), isr->retrieve(L"a"), union_func));

    auto intersection_func = [](int32_t a, int32_t b) { return a & b; };
    REQUIRE(!isr->valid_bitop(isr->retrieve(L"i"), isr->retrieve(L"g"),
                              intersection_func));
    REQUIRE(isr->valid_bitop(isr->retrieve(L"i"), isr->retrieve(L"m"),
                             intersection_func));
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
    auto f12_base_spaces = isr->base_spaces();
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
}
