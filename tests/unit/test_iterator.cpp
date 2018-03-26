//
// Created by Eduard Valeyev on 3/24/18.
//

#include "catch.hpp"

#include "../../src/SeQuant2/iterator.hpp"
#include "../../src/SeQuant2/op.hpp"

TEST_CASE("Iterators", "[elements]") {

  using namespace sequant2;

  SECTION("constructor") {
    REQUIRE_NOTHROW(flattened_rangenest<FNOperatorSeq>{});
    auto rng0 = flattened_rangenest<FNOperatorSeq>{};

    {
      auto opseq1 = FNOperatorSeq{{FNOperator({L"i_1", L"i_3"}, {L"i_2", L"i_4"}),
                                   FNOperator({L"i_5"}, {L"i_6"})}};

      REQUIRE_NOTHROW(flattened_rangenest<FNOperatorSeq>{&opseq1});
      auto rng1 = flattened_rangenest<FNOperatorSeq>{&opseq1};
      REQUIRE_NOTHROW(begin(rng1));
      auto it1 = begin(rng1);
      for (auto i = 0; i != 6; ++i) {
        REQUIRE(ranges::get_cursor(it1).index() == i);
        switch (i) {
          case 0: {
            REQUIRE(*it1 == fcre(L"i_1"));
            break;
          }
          case 1: {
            REQUIRE(*it1 == fcre(L"i_3"));
            break;
          }
          case 2: {
            REQUIRE(*it1 == fann(L"i_4"));
            break;
          }
          case 3: {
            REQUIRE(*it1 == fann(L"i_2"));
            break;
          }
          case 4: {
            REQUIRE(*it1 == fcre(L"i_5"));
            break;
          }
          case 5: {
            REQUIRE(*it1 == fann(L"i_6"));
            break;
          }
        }
        REQUIRE_NOTHROW(++it1);
      }
      REQUIRE_NOTHROW(end(rng1));
      REQUIRE(it1 == end(rng1));
      REQUIRE(ranges::get_cursor(end(rng1)).index() == 6);
    }
  }

}  // TEST_CASE("Iterators")