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

namespace {
// whether the deprecated argument-list spellings of add_union() and
// add_intersection() accept Args
SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
template <typename... Args>
concept legacy_union_takes =
    requires(sequant::IndexBasisRegistry& r, Args... args) {
      r.add_union(L"R", {L"q", L"i"}, args...);
      r.add_intersection(L"R", {L"q", L"i"}, args...);
    };
SEQUANT_PRAGMA_IGNORE_DEPRECATED_END

// whether a brace list of a T initializes an IndexSpaceSpec / IndexBasisSpec
template <typename T>
concept space_spec_from = requires(T v) { sequant::IndexSpaceSpec{v}; };
template <typename T>
concept basis_spec_from = requires(T v) { sequant::IndexBasisSpec{v}; };
}  // namespace

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
    auto sr_ibr = sequant::mbpt::make_sr_spaces();
    REQUIRE_NOTHROW(sr_ibr->retrieve_space(L"i"));
    REQUIRE_NOTHROW(sr_ibr->remove(L"i"));
    IndexSpace active_occupied(L"i", 0b0010);
    REQUIRE_NOTHROW(sr_ibr->add(active_occupied));
    REQUIRE_THROWS(sr_ibr->add(
        active_occupied));  // cannot add a space that is already present
    IndexSpace jactive_occupied(L"j", 0b0010);
    REQUIRE_THROWS(sr_ibr->add(
        jactive_occupied));  // cannot add a space with duplicate attr

    // can use bytestrings too
    REQUIRE_NOTHROW(sr_ibr->retrieve_space("a"));  // N.B. string
    REQUIRE_NOTHROW(sr_ibr->remove('a'));          // N.B. char
  }

  // the table constructor takes the entries as given, but their labels must
  // still be base keys, or indices in them could not be parsed back
  SECTION("table constructor validates labels") {
    auto ibr = sequant::mbpt::make_min_sr_spaces();
    CHECK_NOTHROW(IndexBasisRegistry{ibr->bases()});  // the null space's
                                                      // empty key is fine
    auto table = ibr->bases();
    table.emplace(L"a1", IndexBasis{ibr->retrieve_space(L"a"), 1});
    CHECK_THROWS_AS(IndexBasisRegistry{table}, Exception);
    auto table2 = ibr->bases();
    table2.emplace(L"x_1", IndexBasis{IndexSpace{L"x_1", 0b100000}});
    CHECK_THROWS_AS(IndexBasisRegistry{table2}, Exception);
  }

  // the pre-rename spellings compile and reach the same registry
  SECTION("deprecated spellings") {
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto ibr = sequant::mbpt::make_min_sr_spaces();
    Context by_shared_ptr({.index_space_registry_shared_ptr = ibr});
    CHECK(*by_shared_ptr.index_space_registry() == *ibr);
    CHECK(by_shared_ptr.index_space_registry() ==
          by_shared_ptr.index_basis_registry());
    Context by_value({.index_space_registry = IndexBasisRegistry(*ibr)});
    CHECK(*by_value.index_space_registry() == *ibr);
    auto scoped = set_scoped_default_context(by_shared_ptr);
    CHECK(*get_default_index_space_registry() == *ibr);
    CHECK(get_default_index_space_registry() ==
          get_default_index_basis_registry());
    // the current spelling wins when both are given
    auto other = sequant::mbpt::make_sr_spaces();
    Context both({.index_basis_registry_shared_ptr = other,
                  .index_space_registry_shared_ptr = ibr});
    CHECK(*both.index_basis_registry() == *other);
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
  }

  SECTION("named basis instances") {
    auto ibr = sequant::mbpt::make_min_sr_spaces();
    const IndexSpace a = ibr->retrieve_space(L"a");
    const IndexBasis pao{a, 2147483647};
    const auto nbase = ibr->base_spaces().size();
    const auto nentries = ibr->bases().size();
    static_assert(std::is_same_v<decltype(ibr->bases()),
                                 const IndexBasisRegistry::table_type&>);

    // registration: label, space, instance, own size; the space keeps its own
    // size
    CHECK_FALSE(ibr->basis_label(pao));  // nothing named yet
    REQUIRE_NOTHROW(ibr->add(L"μ̃", a, 2147483647, {.extent = 120}));
    CHECK(ibr->bases().size() == nentries + 1);
    CHECK(ibr->contains(L"μ̃"));
    CHECK(ibr->contains(pao));
    CHECK(ibr->contains(IndexBasis{a}));
    CHECK_FALSE(ibr->contains(IndexBasis{a, 7}));
    CHECK(std::as_const(*ibr).retrieve_space_ptr(L"μ̃") ==
          nullptr);  // a name is not a space
    CHECK_THROWS_AS(ibr->retrieve_space(L"μ̃"), IndexBasisRegistry::not_a_space);
    CHECK_THROWS_AS(ibr->retrieve_space(L"μ̃"),
                    IndexSpace::bad_key);  // not_a_space is a bad_key
    REQUIRE(ibr->retrieve_ptr(L"μ̃") != nullptr);
    REQUIRE(ibr->retrieve_ptr(L"a") !=
            nullptr);  // both kinds live in the one table
    CHECK(*ibr->retrieve_ptr(L"a") == IndexBasis{a});
    CHECK(ibr->retrieve_ptr(L"ζ") == nullptr);
    const IndexBasis got =
        ibr->retrieve(L"μ̃_3");  // an Index label reduces to its base key
    CHECK(got == IndexBasis(a, 2147483647, L"μ̃"));  // the name is identity,
    CHECK(got != pao);                              // the unnamed pair is not
    CHECK(same_instance(got, pao));
    CHECK(got.extent() == 120);  // the entry carries the metadata
    CHECK(got.space().dimension() == a.dimension());
    CHECK(got.metric() == IndexBasisMetric::Unit);
    CHECK(got.field() == a.field());
    CHECK(ibr->retrieve_space(L"a").dimension() == a.dimension());
    CHECK(ibr->basis_label(pao) == std::optional<std::wstring_view>{L"μ̃"});
    CHECK_FALSE(ibr->basis_label(IndexBasis{a}));
    CHECK_FALSE(ibr->basis_label(IndexBasis{a, 7}));
    CHECK(ibr->resolve(IndexBasis{a, 2147483647}).extent() == 120);
    CHECK(ibr->resolve(IndexBasis{a, 7}) == IndexBasis{a, 7});
    CHECK(ibr->resolve(IndexBasis{a}) == IndexBasis{a});
    CHECK(ibr->bases().find(std::wstring_view(L"μ̃"))->second == got);

    // the space views do not see the name: the space a appears once, not also
    // as the named entry's copy of it
    const auto is_a = [](const IndexSpace& s) { return s.base_key() == L"a"; };
    CHECK(std::ranges::count_if(*ibr, is_a) == 1);
    CHECK(std::ranges::count_if(ibr->spaces(), is_a) == 1);
    CHECK(std::distance(ibr->spaces().begin(), ibr->spaces().end()) ==
          static_cast<std::ptrdiff_t>(nentries));
    CHECK(ibr->base_spaces().size() == nbase);  // unchanged by the registration
    // type/qns and attr lookups return the space entry, never the named
    // entry's copy
    REQUIRE(std::as_const(*ibr).retrieve_space_ptr(L"a") != nullptr);
    CHECK(ibr->retrieve_space_ptr(a.type(), a.qns()) ==
          std::as_const(*ibr).retrieve_space_ptr(L"a"));
    CHECK(ibr->retrieve_space_ptr(a.attr()) ==
          std::as_const(*ibr).retrieve_space_ptr(L"a"));
    // ... also when the name sorts before the space's key (Ĩ < i): localized
    // occupied orbitals as basis instance 1 of i, with a size of their own
    const IndexSpace i = ibr->retrieve_space(L"i");
    REQUIRE(i.dimension() != 50);
    REQUIRE_NOTHROW(ibr->add(L"Ĩ", i, 1, {.extent = 50}));
    REQUIRE(ibr->bases().find(std::wstring_view(L"Ĩ")) <
            ibr->bases().find(std::wstring_view(L"i")));
    REQUIRE(std::as_const(*ibr).retrieve_space_ptr(L"i") != nullptr);
    CHECK(ibr->retrieve_space_ptr(i.attr()) ==
          std::as_const(*ibr).retrieve_space_ptr(L"i"));
    CHECK(ibr->retrieve_space_ptr(i.type(), i.qns()) ==
          std::as_const(*ibr).retrieve_space_ptr(L"i"));
    CHECK(ibr->retrieve_space_ptr(i.attr())->dimension() == i.dimension());
    CHECK(ibr->retrieve_space_ptr(i.type(), i.qns())->dimension() ==
          i.dimension());
    REQUIRE_NOTHROW(ibr->remove(L"Ĩ"));

    // metadata by label, either kind
    ibr->extent(L"μ̃", 77)
        .field(L"μ̃", Field::Real)
        .metric(L"μ̃", IndexBasisMetric::General);
    CHECK(ibr->retrieve(L"μ̃").extent() == 77);
    CHECK(ibr->retrieve(L"μ̃").field() == Field::Real);
    CHECK(ibr->retrieve(L"μ̃").metric() == IndexBasisMetric::General);
    // the basis's metadata is its own: the space is untouched
    CHECK(ibr->retrieve(L"μ̃").space().dimension() == a.dimension());
    CHECK(ibr->retrieve(L"μ̃").space().field() == a.field());
    CHECK(ibr->retrieve_space(L"a").field() == a.field());
    // the own basis of a space is orthonormal
    CHECK_THROWS_AS(ibr->metric(L"a", IndexBasisMetric::General), Exception);
    REQUIRE(std::ranges::find(ibr->base_spaces(), a) !=
            ibr->base_spaces().end());  // memoizes the base spaces
    // ã: a named basis of a with no extent or field of its own
    REQUIRE_NOTHROW(ibr->add(L"ã", a, 5));
    REQUIRE(ibr->retrieve(L"ã").extent() == a.dimension());
    REQUIRE(ibr->retrieve(L"ã").field() == a.field());
    CHECK_FALSE(ibr->retrieve(L"ã").own_extent());
    CHECK_FALSE(ibr->retrieve(L"ã").own_field());
    CHECK(ibr->retrieve(L"μ̃").own_extent() == 77);
    CHECK(ibr->retrieve(L"μ̃").own_field() == Field::Real);
    // â: registered with the space's current extent and field, which are
    // then its own and do not follow the space
    REQUIRE_NOTHROW(
        ibr->add(L"â", a, 6, {.extent = a.dimension(), .field = a.field()}));
    CHECK(ibr->retrieve(L"â").own_extent() == a.dimension());
    CHECK(ibr->retrieve(L"â").own_field() == a.field());
    ibr->extent(L"a", 33);
    CHECK(ibr->retrieve(L"â").extent() == a.dimension());
    CHECK(ibr->retrieve_space(L"a").dimension() == 33);
    CHECK(std::ranges::find(ibr->base_spaces(), a)->dimension() ==
          33);  // the memoized base spaces see the new size
    CHECK(ibr->retrieve(L"μ̃").extent() == 77);  // its own extent
    // a basis without an extent of its own has the space's, as set now
    CHECK(ibr->retrieve(L"ã").extent() == 33);
    CHECK(ibr->retrieve(L"ã").space().dimension() == 33);
    CHECK(ibr->retrieve(L"μ̃").space().dimension() == 33);
    const Field other_field =
        a.field() == Field::Real ? Field::Complex : Field::Real;
    ibr->field(L"a", other_field);
    CHECK(ibr->retrieve(L"ã").field() == other_field);
    CHECK(ibr->retrieve(L"μ̃").field() == Field::Real);  // its own
    CHECK(ibr->retrieve(L"â").field() == a.field());    // its own
    ibr->field(L"a", a.field());
    ibr->extent(L"a", 34);
    CHECK(ibr->retrieve_space(L"a").dimension() == 34);
    CHECK(ibr->retrieve(L"ã").extent() == 34);
    REQUIRE_NOTHROW(ibr->remove(L"ã"));
    REQUIRE_NOTHROW(ibr->remove(L"â"));
    CHECK_THROWS_AS(ibr->extent(L"ζ", 1), IndexSpace::bad_key);
    // the deprecated spellings
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN {
      IndexBasisRegistry legacy = *ibr;
      // a negative dimension is rejected
      CHECK_THROWS_WITH(legacy.add(L"q", 0b100, a.qns(), -5),
                        Catch::Matchers::ContainsSubstring("negative"));
      CHECK_FALSE(legacy.contains(L"q"));
      legacy.add(L"q", 0b100, a.qns(), 5, is_particle);
      CHECK(legacy.retrieve_space(L"q").dimension() == 5);
      CHECK(legacy.particle_space() == legacy.retrieve_space(L"q").type());
      legacy.add_union(L"Q", {L"q", L"i"}, is_complete);
      CHECK(legacy.complete_space() == legacy.retrieve_space(L"Q").type());
      // ... and registers the same entries as the spec spelling: equality
      // sees the dimension and the roles
      IndexBasisRegistry by_spec = *ibr;
      by_spec
          .add(L"q", 0b100, {.dimension = 5, .qns = a.qns(), .particle = true})
          .add_union(L"Q", {L"q", L"i"}, {.complete = true});
      CHECK(by_spec == legacy);
      // the dimension is an integer, never a bool or a character
      static_assert(index_space_argument<int> && index_space_argument<long> &&
                    index_space_argument<unsigned long>);
      static_assert(!index_space_argument<bool> &&
                    !index_space_argument<char> &&
                    !index_space_argument<wchar_t>);
      // a union or an intersection takes role tags only
      static_assert(legacy_union_takes<IsHole, IsParticle>);
      static_assert(!legacy_union_takes<IndexSpace::QuantumNumbers>);
      static_assert(!legacy_union_takes<int>);
      static_assert(!legacy_union_takes<Field>);
    }
    ibr->approximate_size(L"μ̃", 78);
    CHECK(ibr->retrieve(L"μ̃").extent() == 78);
    CHECK(ibr->retrieve_space(L"a").approximate_size() ==
          ibr->retrieve_space(L"a").dimension());
    IndexSpace s = ibr->retrieve_space(L"a");
    s.approximate_size(44);
    CHECK(s.dimension() == 44);
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    ibr->extent(L"μ̃", 77);

    // the set algebra and the occupancy queries take the label of a named
    // basis as that of its space
    {
      const IndexSpace p = ibr->retrieve_space(L"p");
      CHECK(ibr->intersection(L"μ̃", L"p") == a);
      CHECK(ibr->intersection(L"μ̃", L"i") == IndexSpace::null);
      CHECK(ibr->valid_intersection(L"μ̃", L"p"));
      CHECK_FALSE(ibr->valid_intersection(L"μ̃", L"i"));
      CHECK(ibr->unIon(L"μ̃", L"i") == p);
      CHECK(ibr->valid_unIon(L"μ̃", L"i"));
      CHECK_FALSE(ibr->valid_unIon(L"μ̃", L"a"));  // the union of a with itself
      CHECK(ibr->is_base(L"μ̃"));
      CHECK_FALSE(ibr->is_pure_occupied(L"μ̃"));
      CHECK(ibr->is_pure_unoccupied(L"μ̃"));
      CHECK_FALSE(ibr->contains_occupied(L"μ̃"));
      CHECK(ibr->contains_unoccupied(L"μ̃"));
      // a label that is registered neither as a space nor as a name
      const auto unregistered =
          Catch::Matchers::ContainsSubstring("ζ is not registered");
      CHECK_THROWS_WITH(ibr->intersection(L"ζ", L"p"), unregistered);
      CHECK_THROWS_WITH(ibr->intersection(L"p", L"ζ"), unregistered);
      CHECK_THROWS_WITH(ibr->valid_intersection(L"ζ", L"p"), unregistered);
      CHECK_THROWS_WITH(ibr->unIon(L"ζ", L"i"), unregistered);
      CHECK_THROWS_WITH(ibr->valid_unIon(L"ζ", L"i"), unregistered);
      CHECK_THROWS_WITH(ibr->is_base(L"ζ"), unregistered);
      CHECK_THROWS_WITH(ibr->is_pure_occupied(L"ζ"), unregistered);
      CHECK_THROWS_WITH(ibr->is_pure_unoccupied(L"ζ"), unregistered);
      CHECK_THROWS_WITH(ibr->contains_occupied(L"ζ"), unregistered);
      CHECK_THROWS_WITH(ibr->contains_unoccupied(L"ζ"), unregistered);
      // the space-role setters designate a type: a name is not a space
      CHECK_THROWS_AS(ibr->vacuum_occupied_space(L"μ̃"),
                      IndexBasisRegistry::not_a_space);
      CHECK_THROWS_AS(ibr->particle_space(L"μ̃"),
                      IndexBasisRegistry::not_a_space);
    }

    // invariants
    CHECK_THROWS(ibr->add(L"μ̃", a, 5));  // 1: label taken
    CHECK_THROWS(
        ibr->add(IndexSpace{L"μ̃", 0b100}));  // 1: a space may not take a name
    CHECK_THROWS(ibr->add(L"ν", a, 2147483647));  // 3: one name per basis
    CHECK_THROWS(ibr->add(L"ν", IndexSpace{L"q", 0b100},
                          1));             // 4: space not registered
    CHECK_THROWS(ibr->add(L"a2", a, 1));   // 6: digit after the first character
    CHECK_THROWS(ibr->add(L"μ̃_x", a, 1));  // 6: underscore
    CHECK_THROWS(ibr->add(L"", a, 1));     // 6: empty
    // 6: no serializer syntax, for a name or a space alike
    for (auto bad : {L"ν<", L"ν}", L"ν;", L"ν,", L"ν ν", L"2ν"}) {
      CHECK_THROWS(ibr->add(bad, a, 1));
      CHECK_THROWS(ibr->add(IndexSpace{bad, 0b100000000}));
      CHECK_FALSE(ibr->retrieve_ptr(bad));
    }
    // 7: a named instance remains
    CHECK_THROWS_WITH(ibr->remove(L"a"),
                      Catch::Matchers::ContainsSubstring("μ̃"));
    CHECK_THROWS(ibr->remove(a));
    // a space under the key of a with another attr: removing or replacing it
    // would remove the space of μ̃
    const IndexSpace a_other(L"a", 0b100, a.qns());
    CHECK_THROWS_WITH(ibr->remove(a_other),
                      Catch::Matchers::ContainsSubstring("μ̃"));
    CHECK_THROWS_WITH(ibr->replace(a_other),
                      Catch::Matchers::ContainsSubstring("μ̃"));
    REQUIRE(std::as_const(*ibr).retrieve_space_ptr(L"a") != nullptr);
    CHECK(ibr->retrieve_space(L"a") == a);
    CHECK(same_instance(ibr->retrieve(L"μ̃"), pao));
    REQUIRE_NOTHROW(ibr->remove(L"μ̃"));
    CHECK_FALSE(ibr->contains(L"μ̃"));
    CHECK_FALSE(ibr->basis_label(pao));  // early return, no scan
    REQUIRE_NOTHROW(ibr->remove(L"a"));

    // the table is a value (#665): copies are deep and independent; equality
    // compares every entry, named ones and their metadata included; the table
    // ctor takes spaces and names together
    auto ibr2 = sequant::mbpt::make_min_sr_spaces();
    ibr2->add(L"μ̃", a, 2147483647, {.extent = 120});
    IndexBasisRegistry copy(*ibr2);
    CHECK(copy == *ibr2);
    // the name and its count travel
    CHECK(copy.basis_label(pao) == std::optional<std::wstring_view>{L"μ̃"});
    copy.extent(L"μ̃", 121);
    CHECK_FALSE(copy == *ibr2);  // metadata of a named entry counts
    // ... and the copy is independent
    CHECK(ibr2->retrieve(L"μ̃").extent() == 120);
    copy.remove(L"μ̃");
    CHECK(ibr2->contains(L"μ̃"));
    CHECK_FALSE(copy == *ibr2);
    IndexBasisRegistry moved(std::move(copy));
    CHECK_FALSE(moved.contains(L"μ̃"));
    // the table ctor copies spaces and names together ...
    IndexBasisRegistry from_table(ibr2->bases());
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
    ibr2->clear();
    CHECK(ibr2->bases().size() == 1);     // the null space
    CHECK_FALSE(ibr2->basis_label(pao));  // the count was reset with the table

    // remove(IndexSpace) never erases a named entry under the same key (R3)
    auto ibr3 = sequant::mbpt::make_min_sr_spaces();
    ibr3->add(L"μ̃", a, 2147483647);
    // an unregistered space whose key is the name: no-op
    REQUIRE_NOTHROW(ibr3->remove(IndexSpace{L"μ̃", 0b100}));
    CHECK(ibr3->contains(L"μ̃"));
    // a registry inside a Context is immutable (#665): the context copies a
    // still-owned registry, equal by value
    const Context ctx({.index_basis_registry_shared_ptr = ibr3});
    CHECK(ctx.index_basis_registry().get() != ibr3.get());
    CHECK(*ctx.index_basis_registry() == *ibr3);
    CHECK(same_instance(ctx.index_basis_registry()->retrieve(L"μ̃"), pao));
    ibr3->extent(L"μ̃", 5);
    CHECK_FALSE(*ctx.index_basis_registry() == *ibr3);
  }

  SECTION("registry construction") {
    auto ibr = std::make_shared<IndexBasisRegistry>();

    // similar make_sr_spaces, but no spin AND occupied and unoccupied bits are
    // NOT contiguous!
    REQUIRE_NOTHROW(
        ibr->add(L"o", 0b0001, {.dimension = 3})
            .add("i", 0b0100, {.hole = true})  // N.B. narrow string
            .add(L'a', 0b0010,
                 {.dimension = 50,
                  .qns = QuantumNumbersAttr{},
                  .field = Field::Real,
                  .particle = true})  // N.B. wchar_t + explicit quantum
                                      // numbers + size
            .add('g', 0b1000, {.field = Field::Real})
            .add_union(L"m", {L"o", L"i"},
                       {.vacuum_occupied = true, .reference_occupied = true})
            .add_union(L"e", {L"a", L"g"})
            .add_unIon(L"p", {L"m", L"e"}, {.complete = true})  // N.B. unIon
    );
    REQUIRE(ibr->retrieve_space(L"o").dimension() == 3);
    REQUIRE(ibr->retrieve_space(L"o").field() == Field::Complex);
    REQUIRE(ibr->retrieve_space(L"a").field() == Field::Real);
    REQUIRE(ibr->retrieve_space(L"m").type().to_int32() == 0b0101);
    REQUIRE(ibr->retrieve_space(L"m").field() == Field::Complex);
    REQUIRE(ibr->retrieve_space(L"m").dimension() ==
            ibr->retrieve_space(L"o").dimension() +
                ibr->retrieve_space(L"i").dimension());
    REQUIRE(ibr->retrieve_space(L"p").type().to_int32() == 0b1111);
    REQUIRE(ibr->retrieve_space(L"p").dimension() ==
            ibr->retrieve_space(L"m").dimension() +
                ibr->retrieve_space(L"e").dimension());
    REQUIRE(ibr->retrieve_space(L"p").dimension() == 73);
    REQUIRE_NOTHROW(
        ibr->add_union(L"iag", {L"i", L"a", L"g"})
            .add_union(L"oia", {"o", "i", "a"})  // N.B. narrow strings
            .add_intersection(L"x", {L"oia", L"iag"}));
    REQUIRE(ibr->retrieve_space(L"x").type().to_int32() == 0b0110);
    REQUIRE(ibr->retrieve_space(L"x").dimension() ==
            ibr->retrieve_space(L"i").dimension() +
                ibr->retrieve_space(L"a").dimension());
    REQUIRE(ibr->vacuum_occupied_space(IndexSpace::QuantumNumbers::null) ==
            ibr->retrieve_space(L"m"));
    REQUIRE(ibr->vacuum_unoccupied_space(IndexSpace::QuantumNumbers::null) ==
            ibr->retrieve_space(L"e"));
    const auto& m = ibr->retrieve_space(L"m");
    REQUIRE_NOTHROW(ibr->vacuum_occupied_space(
        container::map<IndexSpace::QuantumNumbers, IndexSpace::Type>{
            {m.qns(), m.type()}}));
    REQUIRE(ibr->vacuum_occupied_space(m.qns()) == m);
    REQUIRE(ibr->retrieve_space("e").field() == Field::Real);
  }

  SECTION("equality") {
    auto sr_ibr = sequant::mbpt::make_sr_spaces();
    REQUIRE(sr_ibr->retrieve_space(L"i") == sr_ibr->retrieve_space(L"i"));
    REQUIRE(sr_ibr->retrieve_space(L"i") == IndexSpace(L"i"));
    REQUIRE(sr_ibr->retrieve_space(L"i") != sr_ibr->retrieve_space(L"a"));

    // registries compare their space specifications as well as their spaces
    REQUIRE(*sr_ibr == IndexBasisRegistry(*sr_ibr));
    IndexBasisRegistry other_hole = *sr_ibr;
    other_hole.hole_space(L"o");
    REQUIRE(other_hole != *sr_ibr);
    IndexBasisRegistry other_mask = *sr_ibr;
    other_mask.physical_particle_attribute_mask(bitset::null);
    REQUIRE(other_mask != *sr_ibr);

    // ... and the dimensions and fields of their spaces, which
    // IndexSpace equality ignores
    IndexBasisRegistry other_size = *sr_ibr;
    other_size.extent(L"i", sr_ibr->retrieve_space(L"i").dimension() + 1);
    REQUIRE(other_size.retrieve_space(L"i") == sr_ibr->retrieve_space(L"i"));
    REQUIRE(other_size != *sr_ibr);
    IndexBasisRegistry other_field = *sr_ibr;
    other_field.field(L"i", Field::Real);
    REQUIRE(sr_ibr->retrieve_space(L"i").field() == Field::Complex);
    REQUIRE(other_field != *sr_ibr);
  }

  SECTION("ordering") {
    auto sr_ibr = sequant::mbpt::make_sr_spaces();
    REQUIRE(!(sr_ibr->retrieve_space(L"i") < sr_ibr->retrieve_space(L"i")));
    REQUIRE(sr_ibr->retrieve_space(L"i") < sr_ibr->retrieve_space(L"a"));
    REQUIRE(!(sr_ibr->retrieve_space(L"a") < sr_ibr->retrieve_space(L"i")));
    REQUIRE(sr_ibr->retrieve_space(L"i") < sr_ibr->retrieve_space(L"m"));
    REQUIRE(!(sr_ibr->retrieve_space(L"m") < sr_ibr->retrieve_space(L"m")));
    REQUIRE(sr_ibr->retrieve_space(L"m") < sr_ibr->retrieve_space(L"a"));
    REQUIRE(sr_ibr->retrieve_space(L"m") < sr_ibr->retrieve_space(L"p"));
    REQUIRE(sr_ibr->retrieve_space(L"a") < sr_ibr->retrieve_space(L"p"));

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
    auto ibr = sequant::mbpt::make_F12_sr_spaces();
    REQUIRE(ibr->retrieve_space(L"i") ==
            ibr->intersection(ibr->retrieve_space(L"i"),
                              ibr->retrieve_space(L"p")));
    REQUIRE(!ibr->intersection(ibr->retrieve_space(L"p↑"),
                               ibr->retrieve_space(L"p↓")));
    REQUIRE(ibr->intersection(ibr->retrieve_space(L"i↑"),
                              ibr->retrieve_space(L"p")) ==
            ibr->retrieve_space(L"i↑"));
    REQUIRE(!ibr->intersection(ibr->retrieve_space(L"a"),
                               ibr->retrieve_space(L"i")));
    REQUIRE(!ibr->intersection(ibr->retrieve_space(L"a"),
                               ibr->retrieve_space(L"α'")));

    REQUIRE(ibr->retrieve_space(L"κ") ==
            ibr->unIon(ibr->retrieve_space(L"p"), ibr->retrieve_space(L"α'")));

    REQUIRE(includes(ibr->retrieve_space(L"κ"), ibr->retrieve_space(L"m")));
    REQUIRE(!includes(ibr->retrieve_space(L"m"), ibr->retrieve_space(L"κ")));
    REQUIRE(includes(ibr->retrieve_space(L"α"), ibr->retrieve_space(L"a")));

    REQUIRE(ibr->valid_intersection(ibr->retrieve_space(L"i"),
                                    ibr->retrieve_space(L"p")));

    REQUIRE(
        ibr->valid_unIon(ibr->retrieve_space(L"i"), ibr->retrieve_space(L"a")));
    REQUIRE(ibr->valid_unIon(ibr->retrieve_space(L"i↑"),
                             ibr->retrieve_space(L"i↓")));
    REQUIRE(!ibr->valid_unIon(ibr->retrieve_space(L"i↑"),
                              ibr->retrieve_space(L"i↑")));
    REQUIRE(!ibr->valid_unIon(ibr->retrieve_space(L"p"),
                              ibr->retrieve_space(L"a")));
    REQUIRE(!ibr->valid_unIon(ibr->retrieve_space(L"p↑"),
                              ibr->retrieve_space(L"a↓")));
  }

  SECTION("occupancy_validation") {
    auto sr_ibr = sequant::mbpt::make_sr_spaces();
    auto mr_ibr = sequant::mbpt::make_mr_spaces();

    REQUIRE(sr_ibr->is_pure_occupied(sr_ibr->retrieve_space(L"i")));
    REQUIRE(sr_ibr->is_pure_unoccupied(sr_ibr->retrieve_space(L"a")));

    REQUIRE(mr_ibr->is_pure_occupied(mr_ibr->retrieve_space(L"i")));
    REQUIRE(mr_ibr->is_pure_unoccupied(mr_ibr->retrieve_space(L"a")));
    REQUIRE(!mr_ibr->is_pure_occupied(mr_ibr->retrieve_space(L"I")));
    // REQUIRE(!mr_ibr->is_pure_unoccupied(mr_ibr->retrieve_space(L"E")));

    REQUIRE(sr_ibr->contains_occupied(sr_ibr->retrieve_space(L"i")));
    REQUIRE(sr_ibr->contains_unoccupied(sr_ibr->retrieve_space(L"a")));

    REQUIRE(mr_ibr->contains_occupied(mr_ibr->retrieve_space(L"M")));
    REQUIRE(mr_ibr->contains_unoccupied(mr_ibr->retrieve_space(L"E")));
    REQUIRE(!mr_ibr->contains_occupied(mr_ibr->retrieve_space(L"a")));
    REQUIRE(!mr_ibr->contains_unoccupied(mr_ibr->retrieve_space(L"i")));
  }

  SECTION("base_space") {
    auto ibr = sequant::mbpt::make_F12_sr_spaces();
    const auto& f12_base_space_types = ibr->base_space_types();
    REQUIRE(f12_base_space_types.size() == 5);
    REQUIRE(f12_base_space_types[0] == 0b00001);
    REQUIRE(f12_base_space_types[1] == 0b00010);
    REQUIRE(f12_base_space_types[2] == 0b00100);
    REQUIRE(f12_base_space_types[3] == 0b01000);
    REQUIRE(f12_base_space_types[4] == 0b10000);
    const auto& f12_base_spaces = ibr->base_spaces();
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
    *ibr = *sequant::mbpt::make_sr_spaces();
    REQUIRE(ibr->base_space_types() == sr_base_space_types);
    ibr->clear();
    REQUIRE(ibr->base_space_types().empty());
    REQUIRE(ibr->base_spaces().empty());

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
    auto ibr = sequant::mbpt::make_min_sr_spaces();
    REQUIRE_NOTHROW(mbpt::add_ao_basis(ibr, mbpt::Spin::any));

    // the OBS AO basis spans the complete space, non-orthonormal
    REQUIRE(ibr->contains(L"μ"));
    CHECK(std::as_const(*ibr).retrieve_space_ptr(L"μ") ==
          nullptr);  // a name is not a space
    const IndexBasis ao = ibr->retrieve(L"μ");
    const auto p = ibr->retrieve_space(L"p");
    CHECK(ao == IndexBasis(p, mbpt::default_ao_basis_instance, L"μ"));
    CHECK(ao.space() == p);
    CHECK(ao.metric() == IndexBasisMetric::General);
    CHECK(ao.extent() == p.dimension());
    CHECK(mbpt::default_ao_basis_instance < mbpt::default_pao_basis_instance);

    // ordering: after every other basis of p
    CHECK((IndexBasis{p} < ao && IndexBasis{p, 0} < ao &&
           IndexBasis{p, 10} < ao));
    CHECK(Index(ibr->retrieve_space(L"i"), 1) < Index(ao, 1));

    // with a CABS the OBS AOs span the complete space less it: p, not κ
    auto f12 = sequant::mbpt::make_F12_sr_spaces();
    REQUIRE_NOTHROW(mbpt::add_ao_basis(f12, mbpt::Spin::any));
    CHECK(f12->retrieve(L"μ").space() == f12->retrieve_space(L"p"));
    // a failure registers nothing, whether that space is missing or the
    // label (or the instance) is taken
    auto no_p = sequant::mbpt::make_F12_sr_spaces();
    no_p->remove(L"p");
    CHECK_THROWS_WITH(mbpt::add_ao_basis(no_p, mbpt::Spin::any),
                      Catch::Matchers::ContainsSubstring("μ"));
    CHECK_FALSE(no_p->contains(L"μ"));
    auto taken = sequant::mbpt::make_F12_sr_spaces();
    taken->add(L"μ", L"κ", 3);
    CHECK_THROWS_WITH(mbpt::add_ao_basis(taken, mbpt::Spin::any),
                      Catch::Matchers::ContainsSubstring("μ"));
    CHECK(taken->retrieve(L"μ").space() == taken->retrieve_space(L"κ"));
    auto named = sequant::mbpt::make_F12_sr_spaces();
    named->add(L"ν", L"p", mbpt::default_ao_basis_instance);
    CHECK_THROWS_WITH(mbpt::add_ao_basis(named, mbpt::Spin::any),
                      Catch::Matchers::ContainsSubstring("ν"));
    CHECK_FALSE(named->contains(L"μ"));

    // the deprecated overload registers MPQC's VBS and ABS layout from the
    // labels of the F12 convention; without either it registers what
    // add_ao_basis(ibr, spin_any) does, and an integer third argument still
    // means vbs, not an instance
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto neither = sequant::mbpt::make_F12_sr_spaces();
    REQUIRE_NOTHROW(mbpt::add_ao_basis(neither, mbpt::Spin::any, false));
    CHECK(*neither == *f12);
    auto by_int = sequant::mbpt::make_F12_sr_spaces();
    REQUIRE_NOTHROW(mbpt::add_ao_basis(by_int, mbpt::Spin::any, 1));
    CHECK(by_int->contains(L"Α"));
    auto vbs = sequant::mbpt::make_F12_sr_spaces();
    REQUIRE_NOTHROW(mbpt::add_ao_basis(vbs, mbpt::Spin::any, /*vbs=*/true));
    for (const auto* label : {L"μ", L"Α", L"Γ"}) {
      CAPTURE(label);
      REQUIRE(vbs->contains(label));
      CHECK(vbs->retrieve(label).metric() == IndexBasisMetric::General);
    }
    CHECK(vbs->retrieve(L"μ").space() == vbs->retrieve_space(L"m"));
    CHECK(vbs->retrieve(L"Α").space() == vbs->retrieve_space(L"e"));
    CHECK(vbs->retrieve(L"Γ").space() == vbs->retrieve_space(L"p"));
    auto abs = sequant::mbpt::make_F12_sr_spaces();
    REQUIRE_NOTHROW(mbpt::add_ao_basis(abs, mbpt::Spin::any, /*vbs=*/false,
                                       /*abs=*/true));
    CHECK(abs->retrieve(L"μ").space() == abs->retrieve_space(L"p"));
    CHECK(abs->retrieve(L"σ").space() == abs->retrieve_space(L"α'"));
    CHECK(abs->retrieve(L"ρ").space() == abs->retrieve_space(L"κ"));
    // ... but no m + α', which ρ spans with both, and a failure registers
    // nothing
    auto both = sequant::mbpt::make_F12_sr_spaces();
    CHECK_THROWS_WITH(mbpt::add_ao_basis(both, mbpt::Spin::any, true, true),
                      Catch::Matchers::ContainsSubstring("ρ"));
    for (const auto* label : {L"μ", L"Α", L"Γ", L"σ", L"ρ", L"Ρ"}) {
      CAPTURE(label);
      CHECK_FALSE(both->contains(label));
    }
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
  }

  SECTION("PAO basis") {
    auto ibr = sequant::mbpt::make_min_sr_spaces();
    REQUIRE_NOTHROW(mbpt::add_pao_basis(ibr, mbpt::Spin::any));
    const auto uocc =
        ibr->retrieve_space(ibr->particle_space(), mbpt::Spin::any);
    REQUIRE(ibr->contains(L"μ̃"));
    CHECK(std::as_const(*ibr).retrieve_space_ptr(L"μ̃") == nullptr);
    const IndexBasis pao = ibr->retrieve(L"μ̃");
    CHECK(pao == IndexBasis(uocc, mbpt::default_pao_basis_instance, L"μ̃"));
    CHECK(
        same_instance(pao, IndexBasis{uocc, mbpt::default_pao_basis_instance}));
    CHECK(mbpt::default_pao_basis_instance ==
          std::numeric_limits<IndexBasis::instance_type>::max());
    CHECK(pao.space() == uocc);
    // the PAOs follow the OBS AO basis, registered on demand
    REQUIRE(ibr->contains(L"μ"));
    const IndexBasis ao = ibr->retrieve(L"μ");
    CHECK(ibr->follows(L"μ̃") == std::optional<std::wstring>{L"μ"});
    CHECK(pao.extent() == ao.extent());
    CHECK(pao.extent() == ibr->retrieve_space(L"p").dimension());
    CHECK(pao.metric() == IndexBasisMetric::General);
    CHECK(pao.field() == ao.field());
    // the spin-cased PAO bases: the same instance of a↑ and a↓, following μ
    for (const auto& [label, spin] :
         {std::pair{L"μ̃↑", mbpt::Spin::alpha}, {L"μ̃↓", mbpt::Spin::beta}}) {
      REQUIRE(ibr->contains(label));
      const IndexBasis pao_spin = ibr->retrieve(label);
      CHECK(pao_spin.space() ==
            ibr->retrieve_space(ibr->particle_space(), spin));
      CHECK(pao_spin.basis_instance() == pao.basis_instance());
      CHECK(pao_spin.metric() == IndexBasisMetric::General);
      CHECK(pao_spin.extent() == ao.extent());
      CHECK(ibr->follows(label) == std::optional<std::wstring>{L"μ"});
    }
    // an existing μ is used as is
    auto with_ao = sequant::mbpt::make_min_sr_spaces();
    mbpt::add_ao_basis(with_ao, mbpt::Spin::any);
    with_ao->extent(L"μ", 600);
    REQUIRE_NOTHROW(mbpt::add_pao_basis(with_ao, mbpt::Spin::any));
    CHECK(with_ao->retrieve(L"μ̃").extent() == 600);
    // the PAOs can be projected from another AO basis, e.g. a program's
    // OBS+VBS one, which must be registered
    auto from_vbs_plus = sequant::mbpt::make_min_sr_spaces();
    CHECK_THROWS_WITH(
        mbpt::add_pao_basis(from_vbs_plus, mbpt::Spin::any, {.ao_basis = L"Γ"}),
        Catch::Matchers::ContainsSubstring("Γ"));
    CHECK_FALSE(from_vbs_plus->contains(L"μ̃"));
    CHECK_FALSE(from_vbs_plus->contains(L"μ"));
    from_vbs_plus
        ->add(L"Γ", L"p", mbpt::default_ao_basis_instance,
              {.metric = IndexBasisMetric::General})
        .extent(L"Γ", 700);
    REQUIRE_NOTHROW(mbpt::add_pao_basis(from_vbs_plus, mbpt::Spin::any,
                                        {.ao_basis = L"Γ"}));
    CHECK_FALSE(from_vbs_plus->contains(L"μ"));
    for (const auto* l : {L"μ̃", L"μ̃↑", L"μ̃↓"}) {
      CAPTURE(l);
      CHECK(from_vbs_plus->follows(l) == std::optional<std::wstring>{L"Γ"});
      CHECK(from_vbs_plus->retrieve(l).extent() == 700);
    }
    // ... but not from one that does not span the PAOs
    auto from_occ = sequant::mbpt::make_min_sr_spaces();
    from_occ->add(L"ι", L"i", 4, {.metric = IndexBasisMetric::General});
    CHECK_THROWS_WITH(
        mbpt::add_pao_basis(from_occ, mbpt::Spin::any, {.ao_basis = L"ι"}),
        Catch::Matchers::ContainsSubstring("does not include"));
    CHECK_FALSE(from_occ->contains(L"μ̃"));
    CHECK_FALSE(from_occ->contains(L"μ̃↑"));
    // a registry that cannot take the AO basis (no complete space designated)
    // is left untouched
    auto no_p = std::make_shared<IndexBasisRegistry>();
    const IndexSpace::QuantumNumbers any{mbpt::Spin::any};
    no_p->add(L"i", 0b01, {.qns = any, .hole = true})
        .add(L"a", 0b10, {.qns = any, .particle = true});
    CHECK_THROWS_WITH(mbpt::add_pao_basis(no_p, mbpt::Spin::any),
                      Catch::Matchers::ContainsSubstring("complete space"));
    CHECK_FALSE(no_p->contains(L"μ̃"));
    CHECK_FALSE(no_p->contains(L"μ"));
    // ... as is one whose μ is not a named basis instance the PAOs can
    // follow: a space, or a follower itself
    auto mu_space = sequant::mbpt::make_min_sr_spaces();
    mu_space->add(L"μ", 0b100, {.qns = any});
    CHECK_THROWS_WITH(mbpt::add_pao_basis(mu_space, mbpt::Spin::any),
                      Catch::Matchers::ContainsSubstring("space"));
    CHECK_FALSE(mu_space->contains(L"μ̃"));
    CHECK_FALSE(mu_space->contains(L"μ̃↑"));
    auto mu_follows = sequant::mbpt::make_min_sr_spaces();
    mbpt::add_ao_basis(mu_follows, mbpt::Spin::any);
    mu_follows->add(L"κ", L"p", 3).follow(L"μ", L"κ");
    CHECK_THROWS_WITH(mbpt::add_pao_basis(mu_follows, mbpt::Spin::any),
                      Catch::Matchers::ContainsSubstring("κ"));
    CHECK_FALSE(mu_follows->contains(L"μ̃"));
    CHECK_FALSE(mu_follows->contains(L"μ̃↑"));
    // ... and the PAO bases cannot take the label of the basis they follow
    auto mu_label = sequant::mbpt::make_min_sr_spaces();
    CHECK_THROWS_WITH(
        mbpt::add_pao_basis(mu_label, mbpt::Spin::any, {.label = L"μ"}),
        Catch::Matchers::ContainsSubstring("μ"));
    CHECK_FALSE(mu_label->contains(L"μ"));
    // ordering: after every basis of a (I2), between i and Κ as before
    CHECK((IndexBasis{uocc} < pao && IndexBasis{uocc, 0} < pao &&
           IndexBasis{uocc, 10} < pao));
    CHECK(Index(ibr->retrieve_space(L"i"), 1) < Index(pao, 1));
    // the label is taken
    CHECK_THROWS(mbpt::add_pao_basis(ibr, mbpt::Spin::any));
    // the deprecated spelling registers the same entry
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto ibr2 = sequant::mbpt::make_min_sr_spaces();
    mbpt::add_pao_spaces(ibr2, mbpt::Spin::any);
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    CHECK(ibr2->retrieve(L"μ̃") == pao);
    CHECK(*ibr2 == *ibr);
    // sized through μ, all three at once; not through their own labels
    ibr->extent(L"μ", 500);
    for (const auto* label : {L"μ̃", L"μ̃↑", L"μ̃↓"}) {
      CAPTURE(label);
      CHECK(ibr->retrieve(label).extent() == 500);
      CHECK_THROWS_WITH(ibr->extent(label, 1),
                        Catch::Matchers::ContainsSubstring("follows"));
    }
    CHECK_FALSE(*ibr2 == *ibr);  // equality sees the extent

    // an entry of one registry is registered in another by its space,
    // instance and own metadata (as values of its own: the new registry has
    // no μ to follow) ...
    ibr->extent(L"μ", 1234).field(L"μ", Field::Real);
    const IndexBasis src = ibr->retrieve(L"μ̃");
    auto ibr3 = sequant::mbpt::make_min_sr_spaces();
    ibr3->add(L"μ̃", src.space(), *src.basis_instance(),
              {.extent = src.own_extent(),
               .metric = src.metric(),
               .field = src.own_field()});
    CHECK(ibr3->retrieve(L"μ̃") == pao);
    CHECK(ibr3->retrieve(L"μ̃").extent() == 1234);
    CHECK(ibr3->retrieve(L"μ̃").metric() == IndexBasisMetric::General);
    CHECK(ibr3->retrieve(L"μ̃").field() == Field::Real);
    CHECK_FALSE(ibr3->follows(L"μ̃"));
    // ... or with metadata and under a label of its own
    auto ibr4 = sequant::mbpt::make_min_sr_spaces();
    ibr4->add(L"ν̃", src.space(), *src.basis_instance(),
              {.extent = 56,
               .metric = IndexBasisMetric::Unit,
               .field = Field::Complex});
    CHECK_FALSE(ibr4->contains(L"μ̃"));
    CHECK(ibr4->retrieve(L"ν̃").name() == L"ν̃");
    CHECK(same_instance(ibr4->retrieve(L"ν̃"), pao));
    CHECK(ibr4->retrieve(L"ν̃").extent() == 56);
    CHECK(ibr4->retrieve(L"ν̃").metric() == IndexBasisMetric::Unit);
    CHECK(ibr4->retrieve(L"ν̃").field() == Field::Complex);
    // a basis without metadata of its own takes that of this registry's space
    auto ibr5 = sequant::mbpt::make_min_sr_spaces();
    ibr5->extent(L"a", 99);
    ibr5->add(L"μ̃", uocc, mbpt::default_pao_basis_instance);
    CHECK(ibr5->retrieve(L"μ̃").extent() == 99);
    CHECK(ibr5->retrieve(L"μ̃").metric() == IndexBasisMetric::Unit);
  }

  SECTION("basis follows basis") {
    auto ibr = sequant::mbpt::make_min_sr_spaces();
    const IndexSpace p = ibr->retrieve_space(L"p"),
                     a = ibr->retrieve_space(L"a");
    ibr->add(L"μ", p, 7, {.extent = 500, .metric = IndexBasisMetric::General})
        .add(L"μ̃", a, 7)
        .add(L"ν̃", a, 8);
    REQUIRE(ibr->retrieve(L"μ̃").extent() == a.dimension());
    CHECK_FALSE(ibr->follows(L"μ̃"));
    REQUIRE_NOTHROW(ibr->follow(L"μ̃", L"μ"));
    CHECK(ibr->follows(L"μ̃") == std::optional<std::wstring>{L"μ"});
    // the follower mirrors the source now and on every change of it
    CHECK(ibr->retrieve(L"μ̃").extent() == 500);
    CHECK(ibr->retrieve(L"μ̃").own_extent() == 500);
    CHECK(ibr->retrieve(L"μ̃").metric() == IndexBasisMetric::General);
    CHECK(ibr->retrieve(L"μ̃").field() == p.field());
    ibr->extent(L"μ", 600).metric(L"μ", IndexBasisMetric::Unit);
    CHECK(ibr->retrieve(L"μ̃").extent() == 600);
    CHECK(ibr->retrieve(L"μ̃").metric() == IndexBasisMetric::Unit);
    // ... also when the source follows its space
    ibr->add(L"κ", p, 9).follow(L"ν̃", L"κ");
    CHECK(ibr->retrieve(L"ν̃").extent() == p.dimension());
    ibr->extent(L"p", 700);
    CHECK(ibr->retrieve(L"κ").extent() == 700);
    CHECK(ibr->retrieve(L"ν̃").extent() == 700);
    const Field other = p.field() == Field::Real ? Field::Complex : Field::Real;
    ibr->field(L"p", other);
    CHECK(ibr->retrieve(L"ν̃").field() == other);
    // a follower's metadata is set through the source
    CHECK_THROWS_WITH(ibr->extent(L"μ̃", 1),
                      Catch::Matchers::ContainsSubstring("follows"));
    CHECK_THROWS_WITH(ibr->field(L"μ̃", Field::Real),
                      Catch::Matchers::ContainsSubstring("follows"));
    CHECK_THROWS_WITH(ibr->metric(L"μ̃", IndexBasisMetric::General),
                      Catch::Matchers::ContainsSubstring("follows"));
    // no chains, no self, no spaces
    CHECK_THROWS(ibr->follow(L"κ", L"ν̃"));  // ν̃ follows κ
    CHECK_THROWS(ibr->follow(L"ν̃", L"μ̃"));  // μ̃ follows μ
    CHECK_THROWS(ibr->follow(L"μ", L"μ"));
    CHECK_THROWS(ibr->follow(L"μ̃", L"a"));
    CHECK_THROWS(ibr->follow(L"a", L"μ"));
    CHECK_THROWS_AS(ibr->follow(L"μ̃", L"ζ"), IndexSpace::bad_key);
    // re-pointing a follower
    REQUIRE_NOTHROW(ibr->follow(L"μ̃", L"κ"));
    CHECK(ibr->retrieve(L"μ̃").extent() == 700);
    // copies follow alike, and equality sees who follows whom
    IndexBasisRegistry copy = *ibr;
    CHECK(copy == *ibr);
    CHECK(copy.follows(L"μ̃") == std::optional<std::wstring>{L"κ"});
    copy.extent(L"κ", 701);
    CHECK(copy.retrieve(L"μ̃").extent() == 701);
    CHECK(ibr->retrieve(L"μ̃").extent() == 700);
    copy.extent(L"κ", 700);
    CHECK(copy == *ibr);
    copy.remove(L"ν̃").add(L"ν̃", a, 8, {.extent = 700});
    CHECK_FALSE(copy == *ibr);  // same values, ν̃ no longer follows κ
    // removing the source detaches the followers, which keep the values
    ibr->remove(L"κ");
    CHECK_FALSE(ibr->follows(L"μ̃"));
    CHECK_FALSE(ibr->follows(L"ν̃"));
    CHECK(ibr->retrieve(L"μ̃").extent() == 700);
    REQUIRE_NOTHROW(ibr->extent(L"μ̃", 1));
    // removing a follower forgets it
    ibr->follow(L"ν̃", L"μ").remove(L"ν̃");
    ibr->add(L"ν̃", a, 8);
    CHECK_FALSE(ibr->follows(L"ν̃"));
  }

  SECTION("add with specs") {
    IndexBasisRegistry ibr;
    const IndexSpace::QuantumNumbers qn{0b11};
    // a space: its type, then anything else by name
    ibr.add(L"i", 0b01,
            {.dimension = 20,
             .qns = qn,
             .field = Field::Real,
             .vacuum_occupied = true,
             .reference_occupied = true,
             .hole = true})
        .add(L"a", 0b10, {.qns = qn, .particle = true})
        .add("x", 0b100);  // no spec: the defaults of an IndexSpace
    // the spec is written by name: a bare value in its brace list, a number
    // (which the deprecated argument-list spelling takes as the dimension), a
    // bool or a Spin, does not compile
    static_assert(IndexSpaceSpec{.dimension = 100}.dimension == 100);
    static_assert(!IndexSpaceSpec{.dimension = 100}.qns);
    static_assert(!space_spec_from<int> && !space_spec_from<bool> &&
                  !space_spec_from<mbpt::Spin>);
    static_assert(IndexBasisSpec{.extent = 100}.extent == 100);
    static_assert(!basis_spec_from<int> && !basis_spec_from<bool> &&
                  !basis_spec_from<IndexBasisMetric>);
    const IndexSpace i = ibr.retrieve_space(L"i"), a = ibr.retrieve_space(L"a"),
                     x = ibr.retrieve_space(L"x");
    CHECK(i.qns() == qn);
    CHECK(i.dimension() == 20);
    CHECK(i.field() == Field::Real);
    CHECK(a.dimension() == IndexSpace(L"a", 0b10).dimension());
    CHECK(a.field() == IndexSpace(L"a", 0b10).field());
    CHECK(x == IndexSpace(L"x", 0b100));
    CHECK(ibr.vacuum_occupied_space() == i.type());
    CHECK(ibr.reference_occupied_space() == i.type());
    CHECK(ibr.hole_space() == i.type());
    CHECK(ibr.particle_space() == a.type());
    CHECK(ibr.complete_space(/* nulltype_ok = */ true) ==
          IndexSpace::Type::null);
    // a union or an intersection computes its qns, dimension and field, so
    // its spec carries the roles only
    ibr.add_union(L"p", {L"i", L"a"}, {.complete = true});
    CHECK(ibr.complete_space() == ibr.retrieve_space(L"p").type());
    REQUIRE_NOTHROW(ibr.add_union("q", {"a", "x"}));
    CHECK_THROWS(ibr.add_union(L"r", {L"i", L"x"}, {.dimension = 5}));
    CHECK_THROWS(ibr.add_intersection(L"r", {L"q", L"p"}, {.qns = qn}));
    CHECK_THROWS(
        ibr.add_intersection(L"r", {L"q", L"p"}, {.field = Field::Real}));
    CHECK_FALSE(ibr.contains(L"r"));

    // a basis instance: its space, by label or by value, and its instance,
    // then the metadata by name
    ibr.add(L"μ̃", L"a", 1, {.metric = IndexBasisMetric::General});
    const IndexBasis pao = ibr.retrieve(L"μ̃");
    CHECK(pao == IndexBasis(a, 1, L"μ̃"));
    CHECK(pao.metric() == IndexBasisMetric::General);
    CHECK_FALSE(pao.own_extent());  // follows the space
    CHECK_FALSE(pao.own_field());
    ibr.add("ν̃", a, 2, {.extent = 120, .field = Field::Real});
    CHECK(ibr.retrieve(L"ν̃") == IndexBasis(a, 2, L"ν̃"));
    CHECK(ibr.retrieve(L"ν̃").own_extent() == 120);
    CHECK(ibr.retrieve(L"ν̃").metric() == IndexBasisMetric::Unit);
    CHECK(ibr.retrieve(L"ν̃").own_field() == Field::Real);
    ibr.add("ã", "a", 3);  // no spec: the space's own basis up to a rotation
    CHECK(ibr.retrieve(L"ã") == IndexBasis(a, 3, L"ã"));
    CHECK(ibr.retrieve(L"ã").extent() == a.dimension());
    CHECK(ibr.retrieve(L"ã").metric() == IndexBasisMetric::Unit);
    CHECK(ibr.retrieve(L"ã").field() == a.field());
    // a follower mirrors its source from the start
    ibr.add(L"κ", L"p", 7, {.extent = 500, .metric = IndexBasisMetric::General})
        .add(L"λ", L"a", 7, {.follows = L"κ"});
    CHECK(ibr.follows(L"λ") == std::optional<std::wstring>{L"κ"});
    CHECK(ibr.retrieve(L"λ").extent() == 500);
    CHECK(ibr.retrieve(L"λ").metric() == IndexBasisMetric::General);
    ibr.add("ξ", "a", 11, {.follows = "κ"});  // any string type names it
    CHECK(ibr.follows(L"ξ") == std::optional<std::wstring>{L"κ"});
    ibr.remove(L"ξ");
    // ... so its metadata cannot be given alongside
    CHECK_THROWS(ibr.add(L"ρ", L"a", 8, {.extent = 1, .follows = L"κ"}));
    CHECK_FALSE(ibr.contains(L"ρ"));
    // a negative dimension or extent is rejected, and registers nothing
    CHECK_THROWS_WITH(ibr.add(L"n", 0b1000, {.dimension = -1}),
                      Catch::Matchers::ContainsSubstring("negative"));
    CHECK_FALSE(ibr.contains(L"n"));
    CHECK_THROWS_WITH(ibr.add(L"ρ", L"a", 8, {.extent = -5}),
                      Catch::Matchers::ContainsSubstring("negative"));
    CHECK_FALSE(ibr.contains(L"ρ"));
    // a source that cannot be followed registers nothing
    CHECK_THROWS_AS(ibr.add(L"ρ", L"a", 8, {.follows = L"ζ"}),
                    IndexSpace::bad_key);
    CHECK_THROWS(ibr.add(L"ρ", L"a", 8, {.follows = L"λ"}));  // λ follows κ
    CHECK_THROWS(ibr.add(L"ρ", L"a", 8, {.follows = L"a"}));  // a space
    CHECK_FALSE(ibr.contains(L"ρ"));
    // the invariants of a named basis instance
    CHECK_THROWS(ibr.add(L"μ̃", L"a", 9));  // label taken
    CHECK_THROWS(ibr.add(L"ρ", L"a", 1));  // instance 1 of a is μ̃
    CHECK_THROWS_AS(ibr.add(L"ρ", L"ζ", 1), IndexSpace::bad_key);
    CHECK_THROWS_AS(ibr.add(L"ρ", L"μ̃", 1),
                    IndexBasisRegistry::not_a_space);  // a name is not a space
    CHECK_THROWS_WITH(ibr.add(L"ρ", IndexSpace{L"ζ", 0b1000}, 1),
                      Catch::Matchers::ContainsSubstring("not registered"));
    // a space of a registered type and qns under another label is not that
    // space
    CHECK_THROWS_WITH(ibr.add(L"ρ", IndexSpace{L"ζ", 0b10, qn}, 1),
                      Catch::Matchers::ContainsSubstring("differs from 'a'"));
    CHECK_THROWS(ibr.add(L"a2", L"a", 9));  // not a label
    CHECK_FALSE(ibr.contains(L"ρ"));
    CHECK(ibr.bases().size() == 1 + 5 + 5);  // null + spaces + bases
  }

  // the paper 1 example as published; it uses deprecated spellings
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
  SECTION("paper 1 example") {
    IndexSpaceRegistry isr;
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
        {.index_space_registry = isr, .vacuum = Vacuum::SingleProduct});
    REQUIRE_NOTHROW(Index("a_1"));
    REQUIRE_NOTHROW(Index("p_1"));
    Index a1("a_1");
    Index p1("p_1");
  }
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_END

  SECTION("paper 2 example") {
    IndexBasisRegistry ibr;
    using namespace sequant::mbpt;
    const QuantumNumbersAttr spin_any = 0b0011;
    ibr.add("i", 0b001,
            {.qns = spin_any,
             // fully occupied in vacuum state
             .vacuum_occupied = true})
        .add("u", 0b010, {.qns = spin_any})
        .add("a", 0b100, {.qns = spin_any})
        // reference RDMs are nonzero in this space, annihilator = qp creator
        .add_union("I", {"i", "u"}, {.reference_occupied = true, .hole = true})
        // creator = qp creator
        .add_union("A", {"u", "a"}, {.particle = true})
        // general particle operators (e.g. Hamiltonian) act on this space
        .add_union("p", {"I", "A"}, {.complete = true});

    // spaces whose QNs contain spin_any contain physical particles
    ibr.physical_particle_attribute_mask(spin_any);

    // projected atomic orbitals = localized, non-orthonormal basis of the
    // virtuals, registered as basis instance 1 of a
    ibr.add("μ̃", "a", 1, {.metric = IndexBasisMetric::General});
    // nonparticle mode: density fitting, or Cholesky, or Hubbard-Stratonovich
    ibr.add("x", 0b001, {.qns = 0b1000});

    REQUIRE(ibr.retrieve_space("u") == ibr.intersection("A", "I"));
    REQUIRE(ibr.retrieve_space("p") == ibr.unIon("A", "i"));
    REQUIRE(ibr.valid_intersection("i", "x") == false);
    REQUIRE(ibr.valid_intersection("a", "μ̃") == true);
    REQUIRE(ibr.is_base("u") == true);
    REQUIRE(ibr.is_base("A") == false);
    REQUIRE(ibr.is_pure_occupied("i") == true);
    REQUIRE(ibr.is_pure_unoccupied("μ̃") == true);
    REQUIRE(ibr.contains_occupied("p") == true);
    REQUIRE(ibr.contains_unoccupied("p") == true);
    REQUIRE(ibr.retrieve("μ̃").metric() == IndexBasisMetric::General);

    auto _ = set_scoped_default_context(
        {.index_basis_registry = ibr, .vacuum = Vacuum::SingleProduct});
    REQUIRE_NOTHROW(Index("A_1"));
    REQUIRE_NOTHROW(Index("μ̃_1"));
    REQUIRE_NOTHROW(Index("x_1"));
    Index A1("A_1");
    Index μ̃1("μ̃_1");
    Index x1("x_1");
  }
}
