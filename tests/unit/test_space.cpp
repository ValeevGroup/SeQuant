//
// Created by Eduard Valeyev on 3/20/18.
//

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "catch2_sequant.hpp"

#include <SeQuant/core/basis.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/index_basis_registry.hpp>
#include <SeQuant/core/space.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>
#include "SeQuant/domain/mbpt/space_qns.hpp"

#include <range/v3/view/filter.hpp>

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <limits>
#include <optional>
#include <string_view>
#include <type_traits>

TEST_CASE("index_space", "[elements]") {
  using namespace sequant;

  SECTION("constructor") {
    REQUIRE_NOTHROW(IndexSpace{});
    REQUIRE(IndexSpace{} == IndexSpace::null);

    REQUIRE_NOTHROW(IndexSpace(L"i", 0b0010, 20));
    IndexSpace active_occupied(L"i", 0b0010, 20);

    // move leave null in its wake
    REQUIRE_NOTHROW(IndexSpace(std::move(active_occupied)));
    REQUIRE(active_occupied == IndexSpace::null);
    active_occupied = IndexSpace(L"i", 0b0010, 20);
    REQUIRE_NOTHROW(IndexSpace{} = std::move(active_occupied));
    REQUIRE(active_occupied == IndexSpace::null);
  }

  SECTION("to_string") {
    IndexSpace active_occupied(L"i", 0b0010, 20);
    REQUIRE_NOTHROW(to_string(active_occupied));
    const auto str = to_string(active_occupied);
    REQUIRE(str == "{attr={type=0x2,qns=0x14},base_key=i,dimension=10}");
  }

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

  // the table constructor takes the entries as given, but their labels must
  // still be base keys, or indices in them could not be parsed back
  SECTION("table constructor validates labels") {
    auto isr = sequant::mbpt::make_min_sr_spaces();
    CHECK_NOTHROW(IndexBasisRegistry{isr->bases()});  // the null space's
                                                      // empty key is fine
    auto table = isr->bases();
    table.emplace(L"a1", IndexBasis{isr->retrieve(L"a"), 1});
    CHECK_THROWS_AS(IndexBasisRegistry{table}, Exception);
    auto table2 = isr->bases();
    table2.emplace(L"x_1", IndexBasis{IndexSpace{L"x_1", 0b100000}});
    CHECK_THROWS_AS(IndexBasisRegistry{table2}, Exception);
  }

  // the pre-rename spellings compile and reach the same registry
  SECTION("deprecated spellings") {
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto isr = sequant::mbpt::make_min_sr_spaces();
    Context by_shared_ptr({.index_space_registry_shared_ptr = isr});
    CHECK(*by_shared_ptr.index_space_registry() == *isr);
    CHECK(by_shared_ptr.index_space_registry() ==
          by_shared_ptr.index_basis_registry());
    Context by_value({.index_space_registry = IndexBasisRegistry(*isr)});
    CHECK(*by_value.index_space_registry() == *isr);
    auto scoped = set_scoped_default_context(by_shared_ptr);
    CHECK(*get_default_index_space_registry() == *isr);
    CHECK(get_default_index_space_registry() ==
          get_default_index_basis_registry());
    // the current spelling wins when both are given
    auto other = sequant::mbpt::make_sr_spaces();
    Context both({.index_basis_registry_shared_ptr = other,
                  .index_space_registry_shared_ptr = isr});
    CHECK(*both.index_basis_registry() == *other);
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
  }

  SECTION("named basis instances") {
    auto isr = sequant::mbpt::make_min_sr_spaces();
    const IndexSpace a = isr->retrieve(L"a");
    const IndexBasis pao{a, 2147483647};
    const auto nbase = isr->base_spaces().size();
    const auto nentries = isr->bases().size();
    static_assert(std::is_same_v<decltype(isr->bases()),
                                 const IndexBasisRegistry::table_type&>);

    // registration: label, basis, own size; the space keeps its own size
    CHECK_FALSE(isr->basis_label(pao));  // nothing named yet
    REQUIRE_NOTHROW(isr->add(L"μ̃", pao, 120ul));
    CHECK(isr->bases().size() == nentries + 1);
    CHECK(isr->contains(L"μ̃"));
    CHECK(isr->contains(pao));
    CHECK(isr->contains(IndexBasis{a}));
    CHECK_FALSE(isr->contains(IndexBasis{a, 7}));
    CHECK(std::as_const(*isr).retrieve_ptr(L"μ̃") ==
          nullptr);  // a name is not a space
    CHECK_THROWS_AS(isr->retrieve(L"μ̃"), IndexBasisRegistry::not_a_space);
    CHECK_THROWS_AS(isr->retrieve(L"μ̃"),
                    IndexSpace::bad_key);  // not_a_space is a bad_key
    REQUIRE(isr->retrieve_basis_ptr(L"μ̃") != nullptr);
    REQUIRE(isr->retrieve_basis_ptr(L"a") !=
            nullptr);  // both kinds live in the one table
    CHECK(*isr->retrieve_basis_ptr(L"a") == IndexBasis{a});
    CHECK(isr->retrieve_basis_ptr(L"ζ") == nullptr);
    const IndexBasis got =
        isr->retrieve_basis(L"μ̃_3");  // an Index label reduces to its base key
    CHECK(got == IndexBasis(a, 2147483647, L"μ̃"));  // the name is identity,
    CHECK(got != pao);                              // the unnamed pair is not
    CHECK(same_instance(got, pao));
    CHECK(got.extent() == 120);  // the entry carries the metadata
    CHECK(got.space().dimension() == a.dimension());
    CHECK(got.metric() == IndexSpaceMetric::Unit);
    CHECK(got.field() == a.field());
    CHECK(isr->retrieve(L"a").dimension() == a.dimension());
    CHECK(isr->basis_label(pao) == std::optional<std::wstring_view>{L"μ̃"});
    CHECK_FALSE(isr->basis_label(IndexBasis{a}));
    CHECK_FALSE(isr->basis_label(IndexBasis{a, 7}));
    CHECK(isr->resolve(IndexBasis{a, 2147483647}).extent() == 120);
    CHECK(isr->resolve(IndexBasis{a, 7}) == IndexBasis{a, 7});
    CHECK(isr->resolve(IndexBasis{a}) == IndexBasis{a});
    CHECK(isr->bases().find(std::wstring_view(L"μ̃"))->second == got);

    // the space views do not see the name: the space a appears once, not also
    // as the named entry's copy of it
    const auto is_a = [](const IndexSpace& s) { return s.base_key() == L"a"; };
    CHECK(std::ranges::count_if(*isr, is_a) == 1);
    CHECK(std::ranges::count_if(isr->spaces(), is_a) == 1);
    CHECK(std::distance(isr->spaces().begin(), isr->spaces().end()) ==
          static_cast<std::ptrdiff_t>(nentries));
    CHECK(isr->base_spaces().size() == nbase);  // unchanged by the registration
    // type/qns and attr lookups return the space entry, never the named
    // entry's copy
    REQUIRE(std::as_const(*isr).retrieve_ptr(L"a") != nullptr);
    CHECK(isr->retrieve_ptr(a.type(), a.qns()) ==
          std::as_const(*isr).retrieve_ptr(L"a"));
    CHECK(isr->retrieve_ptr(a.attr()) ==
          std::as_const(*isr).retrieve_ptr(L"a"));
    // ... also when the name sorts before the space's key (Ĩ < i): localized
    // occupied orbitals as basis instance 1 of i, with a size of their own
    const IndexSpace i = isr->retrieve(L"i");
    REQUIRE(i.dimension() != 50);
    REQUIRE_NOTHROW(isr->add(L"Ĩ", IndexBasis{i, 1}, 50ul));
    REQUIRE(isr->bases().find(std::wstring_view(L"Ĩ")) <
            isr->bases().find(std::wstring_view(L"i")));
    REQUIRE(std::as_const(*isr).retrieve_ptr(L"i") != nullptr);
    CHECK(isr->retrieve_ptr(i.attr()) ==
          std::as_const(*isr).retrieve_ptr(L"i"));
    CHECK(isr->retrieve_ptr(i.type(), i.qns()) ==
          std::as_const(*isr).retrieve_ptr(L"i"));
    CHECK(isr->retrieve_ptr(i.attr())->dimension() == i.dimension());
    CHECK(isr->retrieve_ptr(i.type(), i.qns())->dimension() == i.dimension());
    REQUIRE_NOTHROW(isr->remove(L"Ĩ"));

    // metadata by label, either kind
    isr->extent(L"μ̃", 77)
        .field(L"μ̃", Field::Real)
        .metric(L"μ̃", IndexSpaceMetric::General);
    CHECK(isr->retrieve_basis(L"μ̃").extent() == 77);
    CHECK(isr->retrieve_basis(L"μ̃").field() == Field::Real);
    CHECK(isr->retrieve_basis(L"μ̃").metric() == IndexSpaceMetric::General);
    // the basis's metadata is its own: the space is untouched
    CHECK(isr->retrieve_basis(L"μ̃").space().dimension() == a.dimension());
    CHECK(isr->retrieve_basis(L"μ̃").space().field() == a.field());
    CHECK(isr->retrieve(L"a").field() == a.field());
    // the own basis of a space is orthonormal
    CHECK_THROWS_AS(isr->metric(L"a", IndexSpaceMetric::General), Exception);
    REQUIRE(std::ranges::find(isr->base_spaces(), a) !=
            isr->base_spaces().end());  // memoizes the base spaces
    // ã: a named basis of a with no extent or field of its own
    REQUIRE_NOTHROW(isr->add(L"ã", IndexBasis{a, 5}));
    REQUIRE(isr->retrieve_basis(L"ã").extent() == a.dimension());
    REQUIRE(isr->retrieve_basis(L"ã").field() == a.field());
    CHECK_FALSE(isr->retrieve_basis(L"ã").own_extent());
    CHECK_FALSE(isr->retrieve_basis(L"ã").own_field());
    CHECK(isr->retrieve_basis(L"μ̃").own_extent() == 77);
    CHECK(isr->retrieve_basis(L"μ̃").own_field() == Field::Real);
    // â: registered with the space's current extent and field, which are
    // then its own and do not follow the space
    REQUIRE_NOTHROW(isr->add(L"â", IndexBasis{a, 6}, a.dimension(), a.field()));
    CHECK(isr->retrieve_basis(L"â").own_extent() == a.dimension());
    CHECK(isr->retrieve_basis(L"â").own_field() == a.field());
    // a negative extent is rejected, for a basis and for a space alike
    CHECK_THROWS_WITH(isr->add(L"ă", IndexBasis{a, 7}, -1),
                      Catch::Matchers::ContainsSubstring("negative"));
    CHECK_FALSE(isr->contains(L"ă"));
    CHECK_THROWS_WITH(isr->add(L"q", 0b100, a.qns(), -5),
                      Catch::Matchers::ContainsSubstring("negative"));
    CHECK_FALSE(isr->contains(L"q"));
    isr->extent(L"a", 33);
    CHECK(isr->retrieve_basis(L"â").extent() == a.dimension());
    CHECK(isr->retrieve(L"a").dimension() == 33);
    CHECK(std::ranges::find(isr->base_spaces(), a)->dimension() ==
          33);  // the memoized base spaces see the new size
    CHECK(isr->retrieve_basis(L"μ̃").extent() == 77);  // its own extent
    // a basis without an extent of its own has the space's, as set now
    CHECK(isr->retrieve_basis(L"ã").extent() == 33);
    CHECK(isr->retrieve_basis(L"ã").space().dimension() == 33);
    CHECK(isr->retrieve_basis(L"μ̃").space().dimension() == 33);
    const Field other_field =
        a.field() == Field::Real ? Field::Complex : Field::Real;
    isr->field(L"a", other_field);
    CHECK(isr->retrieve_basis(L"ã").field() == other_field);
    CHECK(isr->retrieve_basis(L"μ̃").field() == Field::Real);  // its own
    CHECK(isr->retrieve_basis(L"â").field() == a.field());    // its own
    isr->field(L"a", a.field());
    isr->extent(L"a", 34);
    CHECK(isr->retrieve(L"a").dimension() == 34);
    CHECK(isr->retrieve_basis(L"ã").extent() == 34);
    REQUIRE_NOTHROW(isr->remove(L"ã"));
    REQUIRE_NOTHROW(isr->remove(L"â"));
    CHECK_THROWS_AS(isr->extent(L"ζ", 1), IndexSpace::bad_key);
    // the deprecated spellings
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    isr->approximate_size(L"μ̃", 78);
    CHECK(isr->retrieve_basis(L"μ̃").extent() == 78);
    CHECK(isr->retrieve(L"a").approximate_size() ==
          isr->retrieve(L"a").dimension());
    IndexSpace s = isr->retrieve(L"a");
    s.approximate_size(44);
    CHECK(s.dimension() == 44);
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    isr->extent(L"μ̃", 77);

    // invariants
    CHECK_THROWS(isr->add(L"μ̃", IndexBasis{a, 5}));  // 1: label taken
    CHECK_THROWS(
        isr->add(IndexSpace{L"μ̃", 0b100}));  // 1: a space may not take a name
    CHECK_THROWS(
        isr->add(L"ν", IndexBasis{a, 2147483647}));  // 3: one name per basis
    CHECK_THROWS(isr->add(L"ν", IndexBasis{IndexSpace{L"q", 0b100},
                                           1}));  // 4: space not registered
    CHECK_THROWS(isr->add(L"ν", IndexBasis{a}));  // 5: needs an instance
    CHECK_THROWS(isr->add(
        L"a2", IndexBasis{a, 1}));  // 6: digit after the first character
    CHECK_THROWS(isr->add(L"μ̃_x", IndexBasis{a, 1}));  // 6: underscore
    CHECK_THROWS(isr->add(L"", IndexBasis{a, 1}));     // 6: empty
    // 6: no serializer syntax, for a name or a space alike
    for (auto bad : {L"ν<", L"ν}", L"ν;", L"ν,", L"ν ν", L"2ν"}) {
      CHECK_THROWS(isr->add(bad, IndexBasis{a, 1}));
      CHECK_THROWS(isr->add(IndexSpace{bad, 0b100000000}));
      CHECK_FALSE(isr->retrieve_basis_ptr(bad));
    }
    // 7: a named instance remains
    CHECK_THROWS_WITH(isr->remove(L"a"),
                      Catch::Matchers::ContainsSubstring("μ̃"));
    CHECK_THROWS(isr->remove(a));
    // a space under the key of a with another attr: removing or replacing it
    // would remove the space of μ̃
    const IndexSpace a_other(L"a", 0b100, a.qns());
    CHECK_THROWS_WITH(isr->remove(a_other),
                      Catch::Matchers::ContainsSubstring("μ̃"));
    CHECK_THROWS_WITH(isr->replace(a_other),
                      Catch::Matchers::ContainsSubstring("μ̃"));
    REQUIRE(std::as_const(*isr).retrieve_ptr(L"a") != nullptr);
    CHECK(isr->retrieve(L"a") == a);
    CHECK(same_instance(isr->retrieve_basis(L"μ̃"), pao));
    REQUIRE_NOTHROW(isr->remove(L"μ̃"));
    CHECK_FALSE(isr->contains(L"μ̃"));
    CHECK_FALSE(isr->basis_label(pao));  // early return, no scan
    REQUIRE_NOTHROW(isr->remove(L"a"));

    // the table is a value (#665): copies are deep and independent; equality
    // compares every entry, named ones and their metadata included; the table
    // ctor takes spaces and names together
    auto isr2 = sequant::mbpt::make_min_sr_spaces();
    isr2->add(L"μ̃", pao, 120ul);
    IndexBasisRegistry copy(*isr2);
    CHECK(copy == *isr2);
    // the name and its count travel
    CHECK(copy.basis_label(pao) == std::optional<std::wstring_view>{L"μ̃"});
    copy.extent(L"μ̃", 121);
    CHECK_FALSE(copy == *isr2);  // metadata of a named entry counts
    // ... and the copy is independent
    CHECK(isr2->retrieve_basis(L"μ̃").extent() == 120);
    copy.remove(L"μ̃");
    CHECK(isr2->contains(L"μ̃"));
    CHECK_FALSE(copy == *isr2);
    IndexBasisRegistry moved(std::move(copy));
    CHECK_FALSE(moved.contains(L"μ̃"));
    // the table ctor copies spaces and names together ...
    IndexBasisRegistry from_table(isr2->bases());
    CHECK(from_table.contains(L"μ̃"));
    CHECK(from_table.contains(L"a"));
    // ... and recounts the names
    CHECK(from_table.basis_label(pao) ==
          std::optional<std::wstring_view>{L"μ̃"});
    // the move members carry the count with the table and leave the source
    // naming nothing (its table is empty)
    IndexBasisRegistry moved_names(std::move(from_table));
    CHECK(moved_names.basis_label(pao) ==
          std::optional<std::wstring_view>{L"μ̃"});
    CHECK_FALSE(from_table.basis_label(pao));
    CHECK_FALSE(from_table.contains(L"μ̃"));
    from_table = std::move(moved_names);
    CHECK(from_table.basis_label(pao) ==
          std::optional<std::wstring_view>{L"μ̃"});
    CHECK_FALSE(moved_names.basis_label(pao));
    // self-move-assignment keeps the table and the count
    auto& from_table_alias = from_table;
    from_table = std::move(from_table_alias);
    CHECK(from_table.contains(L"μ̃"));
    CHECK(from_table.basis_label(pao) ==
          std::optional<std::wstring_view>{L"μ̃"});
    isr2->clear();
    CHECK(isr2->bases().size() == 1);     // the null space
    CHECK_FALSE(isr2->basis_label(pao));  // the count was reset with the table

    // remove(IndexSpace) never erases a named entry under the same key (R3)
    auto isr3 = sequant::mbpt::make_min_sr_spaces();
    isr3->add(L"μ̃", pao);
    // an unregistered space whose key is the name: no-op
    REQUIRE_NOTHROW(isr3->remove(IndexSpace{L"μ̃", 0b100}));
    CHECK(isr3->contains(L"μ̃"));
    // a registry inside a Context is immutable (#665): the context copies a
    // still-owned registry, equal by value
    const Context ctx({.index_basis_registry_shared_ptr = isr3});
    CHECK(ctx.index_basis_registry().get() != isr3.get());
    CHECK(*ctx.index_basis_registry() == *isr3);
    CHECK(same_instance(ctx.index_basis_registry()->retrieve_basis(L"μ̃"), pao));
    isr3->extent(L"μ̃", 5);
    CHECK_FALSE(*ctx.index_basis_registry() == *isr3);
  }

  SECTION("registry construction") {
    auto isr = std::make_shared<IndexBasisRegistry>();

    // similar make_sr_spaces, but no spin AND occupied and unoccupied bits are
    // NOT contiguous!
    REQUIRE_NOTHROW(
        isr->add(L"o", 0b0001, 3)       // dimension
            .add("i", 0b0100, is_hole)  // N.B. narrow string
            .add(L'a', 0b0010, Field::Real, is_particle, QuantumNumbersAttr{},
                 50)  // N.B. wchar_t + explicit quantum numbers + size
            .add('g', 0b1000, Field::Real)
            .add_union(L"m", {L"o", L"i"}, is_vacuum_occupied,
                       is_reference_occupied)
            .add_union(L"e", {L"a", L"g"})
            .add_unIon(L"p", {L"m", L"e"}, is_complete)  // N.B. unIon
    );
    REQUIRE(isr->retrieve(L"o").dimension() == 3);
    REQUIRE(isr->retrieve(L"o").field() == Field::Complex);
    REQUIRE(isr->retrieve(L"a").field() == Field::Real);
    REQUIRE(isr->retrieve(L"m").type().to_int32() == 0b0101);
    REQUIRE(isr->retrieve(L"m").field() == Field::Complex);
    REQUIRE(isr->retrieve(L"m").dimension() ==
            isr->retrieve(L"o").dimension() + isr->retrieve(L"i").dimension());
    REQUIRE(isr->retrieve(L"p").type().to_int32() == 0b1111);
    REQUIRE(isr->retrieve(L"p").dimension() ==
            isr->retrieve(L"m").dimension() + isr->retrieve(L"e").dimension());
    REQUIRE(isr->retrieve(L"p").dimension() == 73);
    REQUIRE_NOTHROW(
        isr->add_union(L"iag", {L"i", L"a", L"g"})
            .add_union(L"oia", {"o", "i", "a"})  // N.B. narrow strings
            .add_intersection(L"x", {L"oia", L"iag"}));
    REQUIRE(isr->retrieve(L"x").type().to_int32() == 0b0110);
    REQUIRE(isr->retrieve(L"x").dimension() ==
            isr->retrieve(L"i").dimension() + isr->retrieve(L"a").dimension());
    REQUIRE(isr->vacuum_occupied_space(IndexSpace::QuantumNumbers::null) ==
            isr->retrieve(L"m"));
    REQUIRE(isr->vacuum_unoccupied_space(IndexSpace::QuantumNumbers::null) ==
            isr->retrieve(L"e"));
    const auto& m = isr->retrieve(L"m");
    REQUIRE_NOTHROW(isr->vacuum_occupied_space(
        container::map<IndexSpace::QuantumNumbers, IndexSpace::Type>{
            {m.qns(), m.type()}}));
    REQUIRE(isr->vacuum_occupied_space(m.qns()) == m);
    REQUIRE(isr->retrieve("e").field() == Field::Real);
  }

  SECTION("equality") {
    auto sr_isr = sequant::mbpt::make_sr_spaces();
    REQUIRE(sr_isr->retrieve(L"i") == sr_isr->retrieve(L"i"));
    REQUIRE(sr_isr->retrieve(L"i") == IndexSpace(L"i"));
    REQUIRE(sr_isr->retrieve(L"i") != sr_isr->retrieve(L"a"));

    // registries compare their space specifications as well as their spaces
    REQUIRE(*sr_isr == IndexBasisRegistry(*sr_isr));
    IndexBasisRegistry other_hole = *sr_isr;
    other_hole.hole_space(L"o");
    REQUIRE(other_hole != *sr_isr);
    IndexBasisRegistry other_mask = *sr_isr;
    other_mask.physical_particle_attribute_mask(bitset::null);
    REQUIRE(other_mask != *sr_isr);

    // ... and the dimensions and fields of their spaces, which
    // IndexSpace equality ignores
    IndexBasisRegistry other_size = *sr_isr;
    other_size.extent(L"i", sr_isr->retrieve(L"i").dimension() + 1);
    REQUIRE(other_size.retrieve(L"i") == sr_isr->retrieve(L"i"));
    REQUIRE(other_size != *sr_isr);
    IndexBasisRegistry other_field = *sr_isr;
    other_field.field(L"i", Field::Real);
    REQUIRE(sr_isr->retrieve(L"i").field() == Field::Complex);
    REQUIRE(other_field != *sr_isr);
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

    // assignment, hence clear(), replaces the memoized base spaces
    const auto sr_base_space_types =
        sequant::mbpt::make_sr_spaces()->base_space_types();
    REQUIRE(sr_base_space_types.size() == 4);
    *isr = *sequant::mbpt::make_sr_spaces();
    REQUIRE(isr->base_space_types() == sr_base_space_types);
    isr->clear();
    REQUIRE(isr->base_space_types().empty());
    REQUIRE(isr->base_spaces().empty());

    // a registry that was moved from memoizes none of the spaces it gave up
    IndexBasisRegistry moved_from = *sequant::mbpt::make_sr_spaces();
    REQUIRE(moved_from.base_space_types() == sr_base_space_types);
    REQUIRE(!moved_from.base_spaces().empty());
    IndexBasisRegistry moved_to = std::move(moved_from);
    REQUIRE(moved_from.base_space_types().empty());
    REQUIRE(moved_from.base_spaces().empty());
    REQUIRE(moved_to.base_space_types() == sr_base_space_types);
    moved_from = std::move(moved_to);
    REQUIRE(moved_to.base_space_types().empty());
    REQUIRE(moved_to.base_spaces().empty());
    REQUIRE(moved_from.base_space_types() == sr_base_space_types);
  }

  SECTION("AO basis") {
    auto isr = sequant::mbpt::make_min_sr_spaces();
    REQUIRE_NOTHROW(mbpt::add_ao_basis(isr, mbpt::Spin::any));

    // the OBS AO basis spans the complete space, non-orthonormal
    REQUIRE(isr->contains(L"μ"));
    CHECK(std::as_const(*isr).retrieve_ptr(L"μ") ==
          nullptr);  // a name is not a space
    const IndexBasis ao = isr->retrieve_basis(L"μ");
    const auto p = isr->retrieve(L"p");
    CHECK(ao == IndexBasis(p, mbpt::default_ao_basis_instance, L"μ"));
    CHECK(ao.space() == p);
    CHECK(ao.metric() == IndexSpaceMetric::General);
    CHECK(ao.extent() == p.dimension());
    CHECK(mbpt::default_ao_basis_instance < mbpt::default_pao_basis_instance);

    // ordering: after every other basis of p
    CHECK((IndexBasis{p} < ao && IndexBasis{p, 0} < ao &&
           IndexBasis{p, 10} < ao));
    CHECK(Index(isr->retrieve(L"i"), 1) < Index(ao, 1));

    // the VBS and ABS variants need the union spaces registered: the F12
    // registry has p = m + e and κ = p + α' ...
    auto vbs = sequant::mbpt::make_F12_sr_spaces();
    REQUIRE_NOTHROW(mbpt::add_ao_basis(vbs, mbpt::Spin::any, /*vbs=*/true));
    for (const auto* label : {L"μ", L"Α", L"Γ"}) {
      CAPTURE(label);
      REQUIRE(vbs->contains(label));
      CHECK(vbs->retrieve_basis(label).metric() == IndexSpaceMetric::General);
    }
    CHECK(vbs->retrieve_basis(L"μ").space() == vbs->retrieve(L"m"));
    CHECK(vbs->retrieve_basis(L"Α").space() == vbs->retrieve(L"e"));
    CHECK(vbs->retrieve_basis(L"Γ").space() == vbs->retrieve(L"p"));
    auto abs = sequant::mbpt::make_F12_sr_spaces();
    REQUIRE_NOTHROW(mbpt::add_ao_basis(abs, mbpt::Spin::any, /*vbs=*/false,
                                       /*abs=*/true));
    CHECK(abs->retrieve_basis(L"μ").space() == abs->retrieve(L"p"));
    CHECK(abs->retrieve_basis(L"σ").space() == abs->retrieve(L"α'"));
    CHECK(abs->retrieve_basis(L"ρ").space() == abs->retrieve(L"κ"));
    // ... but no m + α', which ρ spans with both
    auto both = sequant::mbpt::make_F12_sr_spaces();
    CHECK_THROWS_WITH(mbpt::add_ao_basis(both, mbpt::Spin::any, true, true),
                      Catch::Matchers::ContainsSubstring("ρ"));
    // ... and a failure registers nothing, whether a space is missing or a
    // label (or the instance) is taken
    for (const auto* label : {L"μ", L"Α", L"Γ", L"σ", L"ρ", L"Ρ"}) {
      CAPTURE(label);
      CHECK_FALSE(both->contains(label));
    }
    auto taken = sequant::mbpt::make_F12_sr_spaces();
    taken->add(L"ρ", IndexBasis{taken->retrieve(L"κ"), 3});
    CHECK_THROWS_WITH(
        mbpt::add_ao_basis(taken, mbpt::Spin::any, /*vbs=*/false, /*abs=*/true),
        Catch::Matchers::ContainsSubstring("ρ"));
    CHECK_FALSE(taken->contains(L"μ"));
    CHECK_FALSE(taken->contains(L"σ"));
    auto named = sequant::mbpt::make_F12_sr_spaces();
    named->add(L"ν", IndexBasis{named->retrieve(L"α'"),
                                mbpt::default_ao_basis_instance});
    CHECK_THROWS_WITH(
        mbpt::add_ao_basis(named, mbpt::Spin::any, /*vbs=*/false, /*abs=*/true),
        Catch::Matchers::ContainsSubstring("ν"));
    CHECK_FALSE(named->contains(L"μ"));
    CHECK_FALSE(named->contains(L"σ"));
  }

  SECTION("PAO basis") {
    auto isr = sequant::mbpt::make_min_sr_spaces();
    REQUIRE_NOTHROW(mbpt::add_pao_basis(isr, mbpt::Spin::any));
    const auto uocc = isr->retrieve(isr->particle_space(), mbpt::Spin::any);
    REQUIRE(isr->contains(L"μ̃"));
    CHECK(std::as_const(*isr).retrieve_ptr(L"μ̃") == nullptr);
    const IndexBasis pao = isr->retrieve_basis(L"μ̃");
    CHECK(pao == IndexBasis(uocc, mbpt::default_pao_basis_instance, L"μ̃"));
    CHECK(
        same_instance(pao, IndexBasis{uocc, mbpt::default_pao_basis_instance}));
    CHECK(mbpt::default_pao_basis_instance ==
          std::numeric_limits<IndexBasis::instance_type>::max());
    CHECK(pao.space() == uocc);
    // the PAOs follow the OBS AO basis, registered on demand
    REQUIRE(isr->contains(L"μ"));
    const IndexBasis ao = isr->retrieve_basis(L"μ");
    CHECK(isr->follows(L"μ̃") == std::optional<std::wstring>{L"μ"});
    CHECK(pao.extent() == ao.extent());
    CHECK(pao.extent() == isr->retrieve(L"p").dimension());
    CHECK(pao.metric() == IndexSpaceMetric::General);
    CHECK(pao.field() == ao.field());
    // the spin-cased PAO bases: the same instance of a↑ and a↓, following μ
    for (const auto& [label, spin] :
         {std::pair{L"μ̃↑", mbpt::Spin::alpha}, {L"μ̃↓", mbpt::Spin::beta}}) {
      REQUIRE(isr->contains(label));
      const IndexBasis pao_spin = isr->retrieve_basis(label);
      CHECK(pao_spin.space() == isr->retrieve(isr->particle_space(), spin));
      CHECK(pao_spin.basis_instance() == pao.basis_instance());
      CHECK(pao_spin.metric() == IndexSpaceMetric::General);
      CHECK(pao_spin.extent() == ao.extent());
      CHECK(isr->follows(label) == std::optional<std::wstring>{L"μ"});
    }
    // an existing μ is used as is
    auto with_ao = sequant::mbpt::make_min_sr_spaces();
    mbpt::add_ao_basis(with_ao, mbpt::Spin::any);
    with_ao->extent(L"μ", 600);
    REQUIRE_NOTHROW(mbpt::add_pao_basis(with_ao, mbpt::Spin::any));
    CHECK(with_ao->retrieve_basis(L"μ̃").extent() == 600);
    // a registry that cannot take the AO basis is left untouched
    auto no_p = std::make_shared<IndexBasisRegistry>();
    const IndexSpace::QuantumNumbers any{mbpt::Spin::any};
    no_p->add(L"i", 0b01, any, is_hole).add(L"a", 0b10, any, is_particle);
    CHECK_THROWS_AS(mbpt::add_pao_basis(no_p, mbpt::Spin::any),
                    IndexSpace::bad_key);
    CHECK_FALSE(no_p->contains(L"μ̃"));
    CHECK_FALSE(no_p->contains(L"μ"));
    // ... as is one whose μ is not a named basis instance the PAOs can
    // follow: a space, or a follower itself
    auto mu_space = sequant::mbpt::make_min_sr_spaces();
    mu_space->add(L"μ", 0b100, any);
    CHECK_THROWS_WITH(mbpt::add_pao_basis(mu_space, mbpt::Spin::any),
                      Catch::Matchers::ContainsSubstring("space"));
    CHECK_FALSE(mu_space->contains(L"μ̃"));
    CHECK_FALSE(mu_space->contains(L"μ̃↑"));
    auto mu_follows = sequant::mbpt::make_min_sr_spaces();
    mbpt::add_ao_basis(mu_follows, mbpt::Spin::any);
    mu_follows->add(L"κ", IndexBasis{mu_follows->retrieve(L"p"), 3})
        .follow(L"μ", L"κ");
    CHECK_THROWS_WITH(mbpt::add_pao_basis(mu_follows, mbpt::Spin::any),
                      Catch::Matchers::ContainsSubstring("κ"));
    CHECK_FALSE(mu_follows->contains(L"μ̃"));
    CHECK_FALSE(mu_follows->contains(L"μ̃↑"));
    // ... and the PAO bases cannot take the label of the basis they follow
    auto mu_label = sequant::mbpt::make_min_sr_spaces();
    CHECK_THROWS_WITH(
        mbpt::add_pao_basis(mu_label, mbpt::Spin::any,
                            mbpt::default_pao_basis_instance, L"μ"),
        Catch::Matchers::ContainsSubstring("μ"));
    CHECK_FALSE(mu_label->contains(L"μ"));
    // ordering: after every basis of a (I2), between i and Κ as before
    CHECK((IndexBasis{uocc} < pao && IndexBasis{uocc, 0} < pao &&
           IndexBasis{uocc, 10} < pao));
    CHECK(Index(isr->retrieve(L"i"), 1) < Index(pao, 1));
    // the label is taken
    CHECK_THROWS(mbpt::add_pao_basis(isr, mbpt::Spin::any));
    // the deprecated spelling registers the same entry
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto isr2 = sequant::mbpt::make_min_sr_spaces();
    mbpt::add_pao_spaces(isr2, mbpt::Spin::any);
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    CHECK(isr2->retrieve_basis(L"μ̃") == pao);
    CHECK(*isr2 == *isr);
    // sized through μ, all three at once; not through their own labels
    isr->extent(L"μ", 500);
    for (const auto* label : {L"μ̃", L"μ̃↑", L"μ̃↓"}) {
      CAPTURE(label);
      CHECK(isr->retrieve_basis(label).extent() == 500);
      CHECK_THROWS_WITH(isr->extent(label, 1),
                        Catch::Matchers::ContainsSubstring("follows"));
    }
    CHECK_FALSE(*isr2 == *isr);  // equality sees the extent

    // a basis given to add() under another registry keeps its metadata (as
    // values of its own: the new registry has no μ to follow) ...
    isr->extent(L"μ", 1234).field(L"μ", Field::Real);
    auto isr3 = sequant::mbpt::make_min_sr_spaces();
    isr3->add(L"μ̃", isr->retrieve_basis(L"μ̃"));
    CHECK(isr3->retrieve_basis(L"μ̃") == pao);
    CHECK(isr3->retrieve_basis(L"μ̃").extent() == 1234);
    CHECK(isr3->retrieve_basis(L"μ̃").metric() == IndexSpaceMetric::General);
    CHECK(isr3->retrieve_basis(L"μ̃").field() == Field::Real);
    CHECK_FALSE(isr3->follows(L"μ̃"));
    // ... unless overridden by the arguments, and under the label given
    auto isr4 = sequant::mbpt::make_min_sr_spaces();
    isr4->add(L"ν̃", isr->retrieve_basis(L"μ̃"), 56ul, IndexSpaceMetric::Unit,
              Field::Complex);
    CHECK_FALSE(isr4->contains(L"μ̃"));
    CHECK(isr4->retrieve_basis(L"ν̃").name() == L"ν̃");
    CHECK(same_instance(isr4->retrieve_basis(L"ν̃"), pao));
    CHECK(isr4->retrieve_basis(L"ν̃").extent() == 56);
    CHECK(isr4->retrieve_basis(L"ν̃").metric() == IndexSpaceMetric::Unit);
    CHECK(isr4->retrieve_basis(L"ν̃").field() == Field::Complex);
    // a basis without metadata of its own takes that of this registry's space
    auto isr5 = sequant::mbpt::make_min_sr_spaces();
    isr5->extent(L"a", 99);
    isr5->add(L"μ̃", IndexBasis{uocc, mbpt::default_pao_basis_instance});
    CHECK(isr5->retrieve_basis(L"μ̃").extent() == 99);
    CHECK(isr5->retrieve_basis(L"μ̃").metric() == IndexSpaceMetric::Unit);
  }

  SECTION("basis follows basis") {
    auto isr = sequant::mbpt::make_min_sr_spaces();
    const IndexSpace p = isr->retrieve(L"p"), a = isr->retrieve(L"a");
    isr->add(L"μ", IndexBasis{p, 7}, 500ul, IndexSpaceMetric::General)
        .add(L"μ̃", IndexBasis{a, 7})
        .add(L"ν̃", IndexBasis{a, 8});
    REQUIRE(isr->retrieve_basis(L"μ̃").extent() == a.dimension());
    CHECK_FALSE(isr->follows(L"μ̃"));
    REQUIRE_NOTHROW(isr->follow(L"μ̃", L"μ"));
    CHECK(isr->follows(L"μ̃") == std::optional<std::wstring>{L"μ"});
    // the follower mirrors the source now and on every change of it
    CHECK(isr->retrieve_basis(L"μ̃").extent() == 500);
    CHECK(isr->retrieve_basis(L"μ̃").own_extent() == 500);
    CHECK(isr->retrieve_basis(L"μ̃").metric() == IndexSpaceMetric::General);
    CHECK(isr->retrieve_basis(L"μ̃").field() == p.field());
    isr->extent(L"μ", 600).metric(L"μ", IndexSpaceMetric::Unit);
    CHECK(isr->retrieve_basis(L"μ̃").extent() == 600);
    CHECK(isr->retrieve_basis(L"μ̃").metric() == IndexSpaceMetric::Unit);
    // ... also when the source follows its space
    isr->add(L"κ", IndexBasis{p, 9}).follow(L"ν̃", L"κ");
    CHECK(isr->retrieve_basis(L"ν̃").extent() == p.dimension());
    isr->extent(L"p", 700);
    CHECK(isr->retrieve_basis(L"κ").extent() == 700);
    CHECK(isr->retrieve_basis(L"ν̃").extent() == 700);
    const Field other = p.field() == Field::Real ? Field::Complex : Field::Real;
    isr->field(L"p", other);
    CHECK(isr->retrieve_basis(L"ν̃").field() == other);
    // a follower's metadata is set through the source
    CHECK_THROWS_WITH(isr->extent(L"μ̃", 1),
                      Catch::Matchers::ContainsSubstring("follows"));
    CHECK_THROWS_WITH(isr->field(L"μ̃", Field::Real),
                      Catch::Matchers::ContainsSubstring("follows"));
    CHECK_THROWS_WITH(isr->metric(L"μ̃", IndexSpaceMetric::General),
                      Catch::Matchers::ContainsSubstring("follows"));
    // no chains, no self, no spaces
    CHECK_THROWS(isr->follow(L"κ", L"ν̃"));  // ν̃ follows κ
    CHECK_THROWS(isr->follow(L"ν̃", L"μ̃"));  // μ̃ follows μ
    CHECK_THROWS(isr->follow(L"μ", L"μ"));
    CHECK_THROWS(isr->follow(L"μ̃", L"a"));
    CHECK_THROWS(isr->follow(L"a", L"μ"));
    CHECK_THROWS_AS(isr->follow(L"μ̃", L"ζ"), IndexSpace::bad_key);
    // re-pointing a follower
    REQUIRE_NOTHROW(isr->follow(L"μ̃", L"κ"));
    CHECK(isr->retrieve_basis(L"μ̃").extent() == 700);
    // copies follow alike, and equality sees who follows whom
    IndexBasisRegistry copy = *isr;
    CHECK(copy == *isr);
    CHECK(copy.follows(L"μ̃") == std::optional<std::wstring>{L"κ"});
    copy.extent(L"κ", 701);
    CHECK(copy.retrieve_basis(L"μ̃").extent() == 701);
    CHECK(isr->retrieve_basis(L"μ̃").extent() == 700);
    copy.extent(L"κ", 700);
    CHECK(copy == *isr);
    copy.remove(L"ν̃").add(L"ν̃", IndexBasis{a, 8}, 700ul);
    CHECK_FALSE(copy == *isr);  // same values, ν̃ no longer follows κ
    // removing the source detaches the followers, which keep the values
    isr->remove(L"κ");
    CHECK_FALSE(isr->follows(L"μ̃"));
    CHECK_FALSE(isr->follows(L"ν̃"));
    CHECK(isr->retrieve_basis(L"μ̃").extent() == 700);
    REQUIRE_NOTHROW(isr->extent(L"μ̃", 1));
    // removing a follower forgets it
    isr->follow(L"ν̃", L"μ").remove(L"ν̃");
    isr->add(L"ν̃", IndexBasis{a, 8});
    CHECK_FALSE(isr->follows(L"ν̃"));
  }

  SECTION("paper 1 example") {
    IndexBasisRegistry isr;
    using namespace sequant::mbpt;
    isr.add("i", 0b01,
            // fully occupied in vacuum state
            is_vacuum_occupied)
        .add("a", 0b10)
        .add_union("p", {"i", "a"});

    REQUIRE(isr.retrieve("i") == isr.intersection("p", "i"));
    REQUIRE(isr.retrieve("p") == isr.unIon("a", "i"));
    REQUIRE(isr.valid_intersection("i", "a") == false);
    REQUIRE(isr.is_base("a") == true);
    REQUIRE(isr.is_base("p") == false);
    REQUIRE(isr.is_pure_occupied("i") == true);
    REQUIRE(isr.is_pure_unoccupied("a") == true);
    REQUIRE(isr.contains_occupied("p") == true);

    // to use ISR load into default context
    auto _ = set_scoped_default_context(
        {.index_basis_registry = isr, .vacuum = Vacuum::SingleProduct});
    REQUIRE_NOTHROW(Index("a_1"));
    REQUIRE_NOTHROW(Index("p_1"));
    Index a1("a_1");
    Index p1("p_1");
  }

  SECTION("paper 2 example") {
    IndexBasisRegistry isr;
    using namespace sequant::mbpt;
    const QuantumNumbersAttr spin_any = 0b0011;
    isr.add("i", 0b001, spin_any,
            // fully occupied in vacuum state
            is_vacuum_occupied)
        .add("u", 0b010, spin_any)
        .add("a", 0b100, spin_any)
        // reference RDMs are nonzero in this space, annihilator = qp creator
        .add_union("I", {"i", "u"}, is_reference_occupied, is_hole)
        // creator = qp creator
        .add_union("A", {"u", "a"}, is_particle)
        // general particle operators (e.g. Hamiltonian) act on this space
        .add_union("p", {"I", "A"}, is_complete);

    // spaces whose QNs contain spin_any contain physical particles
    isr.physical_particle_attribute_mask(spin_any);

    // projected atomic orbitals = localized form of virtuals
    const QuantumNumbersAttr lcao = 0b0100;
    isr.add("μ̃", 0b100, lcao | spin_any);
    // nonparticle mode: density fitting, or Cholesky, or Hubbard-Stratonivich
    isr.add("x", 0b001, QuantumNumbersAttr{0b1000});

    REQUIRE(isr.retrieve("u") == isr.intersection("A", "I"));
    REQUIRE(isr.retrieve("p") == isr.unIon("A", "i"));
    REQUIRE(isr.valid_intersection("i", "x") == false);
    REQUIRE(isr.valid_intersection("a", "μ̃") == true);
    REQUIRE(isr.is_base("u") == true);
    REQUIRE(isr.is_base("A") == false);
    REQUIRE(isr.is_pure_occupied("i") == true);
    REQUIRE(isr.is_pure_unoccupied("μ̃") == true);
    REQUIRE(isr.contains_occupied("p") == true);
    REQUIRE(isr.contains_unoccupied("p") == true);

    auto _ = set_scoped_default_context(
        {.index_basis_registry = isr, .vacuum = Vacuum::SingleProduct});
    REQUIRE_NOTHROW(Index("A_1"));
    REQUIRE_NOTHROW(Index("μ̃_1"));
    REQUIRE_NOTHROW(Index("x_1"));
    Index A1("A_1");
    Index μ̃1("μ̃_1");
    Index x1("x_1");
  }
}
