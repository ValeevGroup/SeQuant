//
// Created by Eduard Valeyev on 3/24/18.
//

#include <catch2/catch_test_macros.hpp>

#include "catch2_sequant.hpp"

#include <SeQuant/core/algorithm.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/ranges.hpp>

#include <algorithm>
#include <initializer_list>
#include <iterator>
#include <vector>

#include <range/v3/iterator/basic_iterator.hpp>

TEST_CASE("iterators", "[elements]") {
  using namespace sequant;

  SECTION("constructor") {
    {
      REQUIRE_NOTHROW(flattened_rangenest<FNOperatorSeq>{});
      [[maybe_unused]] auto rng0 = flattened_rangenest<FNOperatorSeq>{};
    }

    {
      auto opseq1 =
          FNOperatorSeq{{FNOperator(cre({"i_1", "i_3"}), ann({"i_2", "i_4"})),
                         FNOperator(cre({"i_5"}), ann({"i_6"}))}};

      REQUIRE_NOTHROW(flattened_rangenest<FNOperatorSeq>{&opseq1});
      auto rng1 = flattened_rangenest<FNOperatorSeq>{&opseq1};
      using std::begin;
      REQUIRE_NOTHROW(begin(rng1));
      auto it1 = begin(rng1);
      for (auto i = 0; i != 6; ++i) {
        REQUIRE(ranges::get_cursor(it1).ordinal() == i);
        switch (i) {
          case 0: {
            REQUIRE(*it1 == fcre("i_1"));
            break;
          }
          case 1: {
            REQUIRE(*it1 == fcre("i_3"));
            break;
          }
          case 2: {
            REQUIRE(*it1 == fann("i_4"));
            break;
          }
          case 3: {
            REQUIRE(*it1 == fann("i_2"));
            break;
          }
          case 4: {
            REQUIRE(*it1 == fcre("i_5"));
            break;
          }
          case 5: {
            REQUIRE(*it1 == fann("i_6"));
            break;
          }
        }
        REQUIRE_NOTHROW(++it1);
      }
      using std::end;
      REQUIRE_NOTHROW(end(rng1));
      REQUIRE(it1 == end(rng1));
      REQUIRE(ranges::get_cursor(end(rng1)).ordinal() == 6);
    }

    // flatten recursively
    {
      auto opseq1 = std::vector<FNOperatorSeq>{
          FNOperatorSeq{{FNOperator(cre({"i_1", "i_3"}), ann({"i_2", "i_4"})),
                         FNOperator(cre({"i_5"}), ann({"i_6"}))}},
          FNOperatorSeq{{FNOperator(cre({"i_7"}), ann({"i_8"})),
                         FNOperator(cre({"i_9"}), ann({"i_10"}))}}};

      auto rng1 = flattened_rangenest<decltype(opseq1)>{&opseq1};
      auto rng2 = flattened_rangenest<decltype(rng1)>{&rng1};
      using std::begin;
      using std::end;
      REQUIRE_NOTHROW(begin(rng1));
      REQUIRE_NOTHROW(end(rng1));
      REQUIRE_NOTHROW(std::distance(begin(rng1), end(rng1)) == 2);
      REQUIRE_NOTHROW(begin(rng2));
      REQUIRE_NOTHROW(end(rng2));
      REQUIRE_NOTHROW(std::distance(begin(rng2), end(rng2)) == 10);
      auto it2 = begin(rng2);
      for (auto i = 0; i != 10; ++i) {
        REQUIRE(ranges::get_cursor(it2).ordinal() == i);
        switch (i) {
          case 0: {
            REQUIRE(*it2 == fcre("i_1"));
            break;
          }
          case 1: {
            REQUIRE(*it2 == fcre("i_3"));
            break;
          }
          case 2: {
            REQUIRE(*it2 == fann("i_4"));
            break;
          }
          case 3: {
            REQUIRE(*it2 == fann("i_2"));
            break;
          }
          case 4: {
            REQUIRE(*it2 == fcre("i_5"));
            break;
          }
          case 5: {
            REQUIRE(*it2 == fann("i_6"));
            break;
          }
          case 6: {
            REQUIRE(*it2 == fcre("i_7"));
            break;
          }
          case 7: {
            REQUIRE(*it2 == fann("i_8"));
            break;
          }
          case 8: {
            REQUIRE(*it2 == fcre("i_9"));
            break;
          }
          case 9: {
            REQUIRE(*it2 == fann("i_10"));
            break;
          }
        }
        REQUIRE_NOTHROW(++it2);
      }
      REQUIRE(it2 == end(rng2));
      REQUIRE(ranges::get_cursor(end(rng2)).ordinal() == 10);
    }
  }
}
