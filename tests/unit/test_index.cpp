//
// Created by Eduard Valeyev on 3/20/18.
//
#include <catch2/catch_test_macros.hpp>

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/latex.hpp>
#include <SeQuant/core/wolfram.hpp>

TEST_CASE("Index", "[elements][index]") {
  using namespace sequant;

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
      REQUIRE_NOTHROW(Index(
          L"i_3", IndexSpace::instance(IndexSpace::active_occupied), {i1, i2}));
      Index i3(L"i_3", IndexSpace::instance(IndexSpace::active_occupied),
               {i1, i2});
      REQUIRE(i3.label() == L"i_3");
      REQUIRE(i3.to_string() == "i_3");
      REQUIRE(i3.space() == IndexSpace::instance(IndexSpace::active_occupied));
      REQUIRE(i3.has_proto_indices());
      REQUIRE(i3.proto_indices().size() == 2);
      REQUIRE(i3.proto_indices()[0] == i1);
      REQUIRE(i3.proto_indices()[1] == i2);

      REQUIRE_NOTHROW(Index(L"i_4", {L"i_1", L"i_2"}));
      Index i4(L"i_4", {L"i_1", L"i_2"});
      REQUIRE(i4.label() == L"i_4");
      REQUIRE(i4.to_string() == "i_4");
      REQUIRE(i4.space() == IndexSpace::instance(IndexSpace::active_occupied));
      REQUIRE(i4.has_proto_indices());
      REQUIRE(i4.proto_indices().size() == 2);
      REQUIRE(i4.proto_indices()[0] == i1);
      REQUIRE(i4.proto_indices()[1] == i2);

      // nonsymmetric proto indices
      REQUIRE_NOTHROW(Index(L"i_5", {i2, i1}, false));
      Index i5(L"i_5", {i2, i1}, false);
      REQUIRE(i5.label() == L"i_5");
      REQUIRE(i5.space() == IndexSpace::instance(IndexSpace::active_occupied));
      REQUIRE(i5.has_proto_indices());
      REQUIRE(i5.proto_indices().size() == 2);
      REQUIRE(i5.proto_indices()[0] == i2);
      REQUIRE(i5.proto_indices()[1] == i1);

      // one of the proto indices is a proto index
      REQUIRE_NOTHROW(Index(L"i_6", {i1, i5}, false));
      Index i6(L"i_6", {i1, i5}, false);
      REQUIRE(i6.label() == L"i_6");
      REQUIRE(i6.space() == IndexSpace::instance(IndexSpace::active_occupied));
      REQUIRE(i6.has_proto_indices());
      REQUIRE(i6.proto_indices().size() == 2);
      REQUIRE(i6.proto_indices()[0] == i1);
      REQUIRE(i6.proto_indices()[1] == i5);

      // symmetric proto indices
      REQUIRE_NOTHROW(Index(L"i_7", {i2, i1}));
      Index i7(L"i_7", {i2, i1});
      REQUIRE(i7.label() == L"i_7");
      REQUIRE(i7.space() == IndexSpace::instance(IndexSpace::active_occupied));
      REQUIRE(i7.has_proto_indices());
      REQUIRE(i7.proto_indices().size() == 2);
      REQUIRE(i7.proto_indices()[0] == i1);  // !!
      REQUIRE(i7.proto_indices()[1] == i2);  // !!

#ifndef NDEBUG
      REQUIRE_THROWS(Index(
          L"i_4", IndexSpace::instance(IndexSpace::active_occupied), {i1, i1}));
      REQUIRE_THROWS(Index(L"i_5", {L"i_1", L"i_1"}));
#endif
    }

    // 'g' is not a standard base key, but we can associate it with an existing
    // space to be able to extend the index vocabulary
#ifndef NDEBUG
    REQUIRE_THROWS(Index{L"h"});
#endif
    REQUIRE_NOTHROW(IndexSpace::register_key(
        L"h",
        IndexSpace::all));  // can assign additional key to a space already
                            // registered, this does not redefine base key
    // and now ...
    REQUIRE_NOTHROW(Index{L"h"});

    // copy constructor
    REQUIRE_NOTHROW(Index(i1));

    // ctor that replaces the space
    {
      Index j1, j2;
      REQUIRE_NOTHROW(j1 = Index(L"j_1", IndexSpace::instance(
                                             IndexSpace::active_unoccupied)));
      REQUIRE_NOTHROW(
          j2 = Index(j1, IndexSpace::instance(IndexSpace::active_occupied)));
      REQUIRE(j1.space() == IndexSpace::active_unoccupied);
      REQUIRE(j2.label() == j1.label());
      REQUIRE(j2.space() == IndexSpace::active_occupied);

      // same but using replace_space
      REQUIRE_NOTHROW(
          j1.replace_space(IndexSpace::instance(IndexSpace::active_occupied)));
      j2 = j1.replace_space(IndexSpace::instance(IndexSpace::active_occupied));
      REQUIRE(j2.label() == j1.label());
      REQUIRE(j2.space() == IndexSpace::active_occupied);

      // same but using replace_space
      REQUIRE_NOTHROW(j1.replace_qns(IndexSpace::alpha));
      j2 = j1.replace_qns(IndexSpace::alpha);
      REQUIRE(j2.label() == j1.label());
      REQUIRE(j2.space() == IndexSpace::active_unoccupied);
      REQUIRE(j2.space().type() == j1.space().type());
      REQUIRE(j2.space().qns() == IndexSpace::alpha);
      REQUIRE(j2.space().qns() != j1.space().qns());
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

    // check copy ctor
    Index i4(i2);
    REQUIRE(i2 == i4);
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

  SECTION("qns ordering") {
    auto p1A = Index(L"p↑_1",
                     IndexSpace::instance(IndexSpace::all, IndexSpace::alpha));
    auto p1B =
        Index(L"p↓_1", IndexSpace::instance(IndexSpace::all, IndexSpace::beta));
    auto p2A = Index(L"p↑_2",
                     IndexSpace::instance(IndexSpace::all, IndexSpace::alpha));
    auto p2B =
        Index(L"p↓_2", IndexSpace::instance(IndexSpace::all, IndexSpace::beta));
    REQUIRE(p1A.space().qns() == IndexSpace::alpha);
    REQUIRE(p2A.space().qns() == IndexSpace::alpha);
    REQUIRE(p1B.space().qns() == IndexSpace::beta);
    REQUIRE(p2B.space().qns() == IndexSpace::beta);
    REQUIRE(p1A < p1B);
    REQUIRE(p2A < p1B);
    REQUIRE(p1A < p2A);
    REQUIRE(p1B < p2B);
  }

  SECTION("hashing") {
    REQUIRE_NOTHROW(hash_value(Index{}));
    Index i1(L"i_1");
    Index i2(L"i_2");
    Index i3(L"i_1", {L"i_4", L"i_5"});
    REQUIRE_NOTHROW(hash_value(i1));
    REQUIRE_NOTHROW(hash_value(i2));
    REQUIRE_NOTHROW(hash_value(i3));
    REQUIRE(hash_value(i1) != hash_value(Index{}));
    REQUIRE(hash_value(i1) != hash_value(i2));
    REQUIRE(i1.label() == i3.label());
    REQUIRE(hash_value(i1) != hash_value(i3));
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

  SECTION("to_string") {
    Index alpha(L"α");
    REQUIRE(alpha.to_string() == "α");

    SEQUANT_PRAGMA_CLANG(diagnostic push)
    SEQUANT_PRAGMA_CLANG(diagnostic ignored "-Wdeprecated-declarations")
    SEQUANT_PRAGMA_GCC(diagnostic push)
    SEQUANT_PRAGMA_GCC(diagnostic ignored "-Wdeprecated-declarations")

    REQUIRE(alpha.ascii_label() == "alpha");

    SEQUANT_PRAGMA_GCC(diagnostic pop)
    SEQUANT_PRAGMA_CLANG(diagnostic pop)
  }

  SECTION("label manipulation") {
    Index alpha(L"α");
    Index alpha1(L"α_1");
    Index alpha_up(L"α↑");
    Index alpha1_up(L"α↑_1");
    REQUIRE_NOTHROW(alpha.make_label_plus_suffix(L'↑'));
    REQUIRE(alpha.make_label_plus_suffix(L'↑') == L"α↑");
    REQUIRE_NOTHROW(alpha.make_label_minus_substring(L'↑'));
    REQUIRE(alpha.make_label_minus_substring(L'↑') == L"α");
    REQUIRE_NOTHROW(alpha1.make_label_minus_substring(L'↑'));
    REQUIRE(alpha1.make_label_minus_substring(L'↑') == L"α_1");
    REQUIRE_NOTHROW(alpha_up.make_label_minus_substring(L'↑'));
    REQUIRE(alpha_up.make_label_minus_substring(L'↑') == L"α");
    REQUIRE_NOTHROW(alpha1_up.make_label_minus_substring(L'↑'));
    REQUIRE(alpha1_up.make_label_minus_substring(L'↑') == L"α_1");
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

    Index a1_up(L"a↑_1", {i1, i2});
    std::wstring a1_up_str = to_latex(a1_up);
    REQUIRE(a1_up_str == L"{a↑_1^{{i_1}{i_2^{{i_3}{i_4}}}}}");

    Index alpha1_up(L"α↑_1", {i1, i2});
    std::wstring alpha1_up_str = to_latex(alpha1_up);
    REQUIRE(alpha1_up_str == L"{\\alpha↑_1^{{i_1}{i_2^{{i_3}{i_4}}}}}");

    Index a12_up(L"a↑_12", {i1, i2});
    std::wstring a12_up_str = to_latex(a12_up);
    REQUIRE(a12_up_str == L"{a↑_{12}^{{i_1}{i_2^{{i_3}{i_4}}}}}");
  }

  SECTION("wolfram") {
    Index i1(L"i_1");
    std::wstring i1_str;
    REQUIRE_NOTHROW(i1_str = i1.to_wolfram());
    REQUIRE(i1_str ==
            L"particleIndex[\"\\!\\(\\*SubscriptBox[\\(i\\), "
            L"\\(1\\)]\\)\",particleSpace[occupied]]");
    REQUIRE(i1.to_wolfram(Action::create) ==
            L"particleIndex[\"\\!\\(\\*SubscriptBox[\\(i\\), "
            L"\\(1\\)]\\)\",particleSpace[occupied],indexType[cre]]");
    REQUIRE(i1.to_wolfram(BraKetPos::ket) ==
            L"particleIndex[\"\\!\\(\\*SubscriptBox[\\(i\\), "
            L"\\(1\\)]\\)\",particleSpace[occupied],indexType[ket]]");
    REQUIRE(i1.to_wolfram(Action::annihilate) ==
            L"particleIndex[\"\\!\\(\\*SubscriptBox[\\(i\\), "
            L"\\(1\\)]\\)\",particleSpace[occupied],indexType[ann]]");
    REQUIRE(i1.to_wolfram(BraKetPos::bra) ==
            L"particleIndex[\"\\!\\(\\*SubscriptBox[\\(i\\), "
            L"\\(1\\)]\\)\",particleSpace[occupied],indexType[bra]]");

    REQUIRE(Index(L"a_1").to_wolfram() ==
            L"particleIndex[\"\\!\\(\\*SubscriptBox[\\(a\\), "
            L"\\(1\\)]\\)\",particleSpace[virtual]]");
    REQUIRE(Index(L"p_1").to_wolfram() ==
            L"particleIndex[\"\\!\\(\\*SubscriptBox[\\(p\\), "
            L"\\(1\\)]\\)\",particleSpace[occupied,virtual]]");
    REQUIRE(Index(L"α'_1").to_wolfram() ==
            L"particleIndex[\"\\!\\(\\*SubscriptBox[\\(α'\\), "
            L"\\(1\\)]\\)\",particleSpace[othervirtual]]");
    REQUIRE(Index(L"α_1").to_wolfram() ==
            L"particleIndex[\"\\!\\(\\*SubscriptBox[\\(α\\), "
            L"\\(1\\)]\\)\",particleSpace[virtual,othervirtual]]");
    REQUIRE(Index(L"κ_1").to_wolfram() ==
            L"particleIndex[\"\\!\\(\\*SubscriptBox[\\(κ\\), "
            L"\\(1\\)]\\)\",particleSpace[occupied,virtual,othervirtual]]");
  }

}  // TEST_CASE("Index")
