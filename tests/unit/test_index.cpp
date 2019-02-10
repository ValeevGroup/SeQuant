//
// Created by Eduard Valeyev on 3/20/18.
//

#include "catch.hpp"

#include "../../src/SeQuant2/index.hpp"

TEST_CASE("Index", "[elements]") {

  using namespace sequant2;
  IndexSpace::register_standard_instances();

  SECTION("constructors") {
    Index i{};

    REQUIRE_NOTHROW(
        Index(std::wstring(L"i_") + std::to_wstring(Index::min_tmp_index())));
    REQUIRE_THROWS(Index(std::wstring(L"i_") +
                         std::to_wstring(Index::min_tmp_index() + 1)));

    Index i1(L"i_1");
    REQUIRE(i1.label() == L"i_1");
    REQUIRE(i1.space() == IndexSpace::instance(IndexSpace::active_occupied));

    Index i2(L"i_2", IndexSpace::instance(IndexSpace::active_occupied));
    REQUIRE(i2.label() == L"i_2");
    REQUIRE(i2.space() == IndexSpace::instance(IndexSpace::active_occupied));

    // examples with proto indices
    {
      REQUIRE_NOTHROW(Index(L"i_3", IndexSpace::instance(IndexSpace::active_occupied), {i1, i2}));
      Index i3(L"i_3", IndexSpace::instance(IndexSpace::active_occupied), {i1, i2});
      REQUIRE(i3.label() == L"i_3");
      REQUIRE(i3.space() == IndexSpace::instance(IndexSpace::active_occupied));
      REQUIRE(i3.has_proto_indices());
      REQUIRE(i3.proto_indices().size() == 2);
      REQUIRE(i3.proto_indices()[0] == i1);
      REQUIRE(i3.proto_indices()[1] == i2);

      REQUIRE_NOTHROW(Index(L"i_4", {L"i_1", L"i_2"}));
      Index i4(L"i_4", {L"i_1", L"i_2"});
      REQUIRE(i4.label() == L"i_4");
      REQUIRE(i4.space() == IndexSpace::instance(IndexSpace::active_occupied));
      REQUIRE(i4.has_proto_indices());
      REQUIRE(i4.proto_indices().size() == 2);
      REQUIRE(i4.proto_indices()[0] == i1);
      REQUIRE(i4.proto_indices()[1] == i2);

      REQUIRE_NOTHROW(Index(L"i_5", {i2, i1}));
      Index i5(L"i_5", {i2, i1});
      REQUIRE(i5.label() == L"i_5");
      REQUIRE(i5.space() == IndexSpace::instance(IndexSpace::active_occupied));
      REQUIRE(i5.has_proto_indices());
      REQUIRE(i5.proto_indices().size() == 2);
      REQUIRE(i5.proto_indices()[0] == i2);
      REQUIRE(i5.proto_indices()[1] == i1);

      // one of the proto indices is a proto index
      REQUIRE_NOTHROW(Index(L"i_6", {i1, i5}));
      Index i6(L"i_6", {i1, i5});
      REQUIRE(i6.label() == L"i_6");
      REQUIRE(i6.space() == IndexSpace::instance(IndexSpace::active_occupied));
      REQUIRE(i6.has_proto_indices());
      REQUIRE(i6.proto_indices().size() == 2);
      REQUIRE(i6.proto_indices()[0] == i1);
      REQUIRE(i6.proto_indices()[1] == i5);

      REQUIRE_THROWS(Index(L"i_4", IndexSpace::instance(IndexSpace::active_occupied), {i1, i1}));
      REQUIRE_THROWS(Index(L"i_5", {L"i_1", L"i_1"}));
    }
  }

  SECTION("equality") {
    Index i1(L"i_1");
    Index i2(L"i_2");
    Index i3(L"i_1");
    REQUIRE(i1 == i1);
    REQUIRE(!(i1 == i2));
    REQUIRE(i1 != i2);
    REQUIRE(i1 == i3);
    REQUIRE(!(i1 != i3));
  }

  SECTION("ordering") {
    Index i1(L"i_1");
    Index i2(L"i_2");
    REQUIRE(i1 < i2);
    REQUIRE(!(i2 < i1));
    REQUIRE(!(i1 < i1));
    Index a1(L"a_2");
    REQUIRE(i1 < a1);
    REQUIRE(!(a1 < i1));
  }

  SECTION("hashing") {
    REQUIRE_NOTHROW(hash_value(Index{}));
    Index i1(L"i_1");
    Index i2(L"i_2");
    REQUIRE_NOTHROW(hash_value(i1));
    REQUIRE_NOTHROW(hash_value(i2));
    REQUIRE(hash_value(i1) != hash_value(Index{}));
    REQUIRE(hash_value(i1) != hash_value(i2));
  }

  SECTION("transform") {
    Index i0(L"i_0");
    Index i1(L"i_1");
    Index i0_13(L"i_0", {L"i_1", L"i_3"});
    Index i1_13(L"i_1", {L"i_1", L"i_3"});
    std::map<Index, Index> map{std::make_pair(Index{L"i_1"}, Index{L"i_2"})};
    REQUIRE(!i0.transform(map));
    REQUIRE(i0 == Index{L"i_0"});
    REQUIRE(i1.transform(map));
    REQUIRE(i1 == Index{L"i_2"});
    REQUIRE(i0_13.transform(map));
    REQUIRE(i0_13 == Index{L"i_0", {L"i_2", L"i_3"}});
    REQUIRE(i1_13.transform(map));
    REQUIRE(i1_13 == Index{L"i_1", {L"i_2", L"i_3"}});
  }

  SECTION("latex") {
    Index i1(L"i_1");
    std::wstring i1_str;
    REQUIRE_NOTHROW(i1_str = to_latex(i1));
    REQUIRE(i1_str == L"{i_1}");

    Index i2(L"i_2", {L"i_3", L"i_4"});
    std::wstring i2_str = to_latex(i2);
    REQUIRE(i2_str == L"{i_2^{{i_3}{i_4}}}");

    Index a1(L"a_1", {i1, i2});
    std::wstring a1_str = to_latex(a1);
    REQUIRE(a1_str == L"{a_1^{{i_1}{i_2^{{i_3}{i_4}}}}}");
  }

}  // TEST_CASE("Index")