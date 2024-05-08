//
// Created by Eduard Valeyev on 3/20/18.
//

#include <catch2/catch_test_macros.hpp>

#include <SeQuant/core/space.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>

TEST_CASE("IndexSpace", "[elements]") {
  using namespace sequant;
  SECTION("registry_core_functionality") {
    auto standard_registry = sequant::mbpt::make_sr_subspaces();
    REQUIRE_NOTHROW(standard_registry.retrieve(L"i"));

    REQUIRE_NOTHROW(standard_registry.relabel(L"i", L"j"));
    REQUIRE_THROWS(standard_registry.relabel(
        L"i", L"l"));  // cannot relabel a removed entry

    REQUIRE_NOTHROW(standard_registry.remove(L"j"));
    IndexSpace active_occupied(L"i", {0b0010});
    REQUIRE_NOTHROW(standard_registry.add(active_occupied));
    REQUIRE_THROWS(standard_registry.add(
        active_occupied));  // cannot add a space that is already present
    IndexSpace jactive_occupied(L"j", {0b0010});
    REQUIRE_THROWS(standard_registry.add(
        jactive_occupied));  // cannot add a space with duplicate attr

    auto spinlabeled_registry = sequant::mbpt::make_min_sr_so_subspaces();
  }

  SECTION("equality") {
    auto standard_registry = sequant::mbpt::make_sr_subspaces();
    REQUIRE(standard_registry.retrieve(L"i") ==
            standard_registry.retrieve(L"i"));
    REQUIRE(standard_registry.retrieve(L"i") !=
            standard_registry.retrieve(L"a"));
  }

  SECTION("ordering") {
    auto standard_registry = sequant::mbpt::make_sr_subspaces();
    REQUIRE(
        !(standard_registry.retrieve(L"i") < standard_registry.retrieve(L"i")));
    REQUIRE(standard_registry.retrieve(L"i") <
            standard_registry.retrieve(L"a"));
    REQUIRE(
        !(standard_registry.retrieve(L"a") < standard_registry.retrieve(L"i")));
    REQUIRE(standard_registry.retrieve(L"i") <
            standard_registry.retrieve(L"m"));
    REQUIRE(
        !(standard_registry.retrieve(L"m") < standard_registry.retrieve(L"m")));
    REQUIRE(standard_registry.retrieve(L"m") <
            standard_registry.retrieve(L"a"));
    REQUIRE(standard_registry.retrieve(L"m") <
            standard_registry.retrieve(L"p"));
    REQUIRE(standard_registry.retrieve(L"a") <
            standard_registry.retrieve(L"p"));

    // test ordering with quantum numbers
    {
      auto i = IndexSpace(L"i", {0b01});
      auto a = IndexSpace(L"a", {0b10});
      auto iA = IndexSpace(L"i", {0b01}, 0b01);
      auto iB = IndexSpace(L"i", {0b01}, 0b10);
      auto aA = IndexSpace(L"a", {0b10}, 0b01);
      auto aB = IndexSpace(L"a", {0b10}, 0b10);

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
    auto f12_registry = sequant::mbpt::make_F12_sr_subspaces();
    REQUIRE(f12_registry.retrieve(L"i") ==
            f12_registry.intersection(f12_registry.retrieve(L"i"),
                                      f12_registry.retrieve(L"p")));
    REQUIRE(f12_registry.nulltype_() ==
            f12_registry.intersection(f12_registry.retrieve(L"a"),
                                      f12_registry.retrieve(L"i")));
    REQUIRE(f12_registry.nulltype_() ==
            f12_registry.intersection(f12_registry.retrieve(L"a"),
                                      f12_registry.retrieve(L"α'")));

    REQUIRE(f12_registry.retrieve(L"κ") ==
            f12_registry.unIon(f12_registry.retrieve(L"p"),
                               f12_registry.retrieve(L"α'")));

    REQUIRE(includes(f12_registry.retrieve(L"κ"), f12_registry.retrieve(L"m")));
    REQUIRE(
        !includes(f12_registry.retrieve(L"m"), f12_registry.retrieve(L"κ")));
    REQUIRE(includes(f12_registry.retrieve(L"α"), f12_registry.retrieve(L"a")));

    REQUIRE(f12_registry.has_non_overlapping_spaces(
        f12_registry.retrieve(L"i"), f12_registry.retrieve(L"p")));
    REQUIRE(f12_registry.has_non_overlapping_spaces(
        f12_registry.retrieve(L"i"), f12_registry.retrieve(L"a")));
    REQUIRE(!f12_registry.has_non_overlapping_spaces(
        f12_registry.retrieve(L"i"), f12_registry.retrieve(L"i")));

    REQUIRE(f12_registry.non_overlapping_spaces(
                f12_registry.retrieve(L"g"), f12_registry.retrieve(L"α"))[0] ==
            f12_registry.retrieve(L"a"));
    REQUIRE(f12_registry.non_overlapping_spaces(
                f12_registry.retrieve(L"g"), f12_registry.retrieve(L"α"))[1] ==
            f12_registry.retrieve(L"α'"));

    auto union_func = [](int32_t a, int32_t b) { return a + b; };
    REQUIRE(!f12_registry.valid_bitop(
        f12_registry.retrieve(L"i"), f12_registry.retrieve(L"α'"), union_func));
    REQUIRE(f12_registry.valid_bitop(f12_registry.retrieve(L"i"),
                                     f12_registry.retrieve(L"a"), union_func));

    auto intersection_func = [](int32_t a, int32_t b) { return a & b; };
    REQUIRE(!f12_registry.valid_bitop(f12_registry.retrieve(L"i"),
                                      f12_registry.retrieve(L"g"),
                                      intersection_func));
    REQUIRE(f12_registry.valid_bitop(f12_registry.retrieve(L"i"),
                                     f12_registry.retrieve(L"m"),
                                     intersection_func));
  }

  SECTION("occupancy_validation") {
    auto standard_registry = sequant::mbpt::make_sr_subspaces();
    auto multireference_registry = sequant::mbpt::make_mr_subspaces();

    REQUIRE(
        standard_registry.is_pure_occupied(standard_registry.retrieve(L"i")));
    REQUIRE(
        standard_registry.is_pure_unoccupied(standard_registry.retrieve(L"a")));

    REQUIRE(multireference_registry.is_pure_occupied(
        multireference_registry.retrieve(L"i")));
    REQUIRE(multireference_registry.is_pure_unoccupied(
        multireference_registry.retrieve(L"a")));
    REQUIRE(!multireference_registry.is_pure_occupied(
        multireference_registry.retrieve(L"I")));
    // REQUIRE(!multireference_registry.is_pure_unoccupied(multireference_registry.retrieve(L"E")));

    REQUIRE(
        standard_registry.contains_occupied(standard_registry.retrieve(L"i")));
    REQUIRE(standard_registry.contains_unoccupied(
        standard_registry.retrieve(L"a")));

    REQUIRE(multireference_registry.contains_occupied(
        multireference_registry.retrieve(L"M")));
    REQUIRE(multireference_registry.contains_unoccupied(
        multireference_registry.retrieve(L"E")));
    REQUIRE(!multireference_registry.contains_occupied(
        multireference_registry.retrieve(L"a")));
    REQUIRE(!multireference_registry.contains_unoccupied(
        multireference_registry.retrieve(L"i")));
  }

  SECTION("base_space") {
    auto f12_registry = sequant::mbpt::make_F12_sr_subspaces();
    auto f12_base_spaces = f12_registry.base_spaces_label();
    REQUIRE(f12_base_spaces[0].second == L"o");
    REQUIRE(f12_base_spaces[1].second == L"i");
    REQUIRE(f12_base_spaces[2].second == L"a");
    REQUIRE(f12_base_spaces[3].second == L"g");
    REQUIRE(f12_base_spaces[4].second == L"α'");
  }
}
