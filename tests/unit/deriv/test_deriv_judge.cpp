#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/utility/expr.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <catch2/catch_test_macros.hpp>
#include "evaluator.hpp"
#include "judge.hpp"

#include <complex>
#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>

using namespace sequant;
namespace t = sequant::deriv::test;

namespace {

std::complex<double> at(std::wstring_view e, const t::Values& values,
                        const t::Data& data = {}) {
  return t::element_of(deserialize(e), values, data);
}

}  // namespace

TEST_CASE("deriv/judge", "[deriv][judge]") {
  // the basis: two values per occupied and three per virtual elementary
  // space, per spin where spins are read
  for (const auto& [label, spins, n] :
       {std::tuple{L"i_1", false, 2}, std::tuple{L"a_1", false, 3},
        std::tuple{L"m_1", false, 4}, std::tuple{L"i↑_1", true, 2},
        std::tuple{L"a↓_1", true, 3}, std::tuple{L"a_1", true, 6}}) {
    t::Domains domains;
    t::add_domain(Index(label), {}, domains, spins);
    INFO(toUtf8(label));
    REQUIRE(domains.front().second.size() == static_cast<std::size_t>(n));
  }
  // three occupied values where a rule asks for triples of distinct values
  t::Domains three;
  t::add_domain(Index(L"i_1"), {}, three, false, 3);
  REQUIRE(three.front().second.size() == 3);
  // each tensor carries exactly its declared symmetry
  const t::Values v{{L"i_1", 0}, {L"i_2", 1}, {L"a_1", 2}, {L"a_2", 3}};
  REQUIRE(at(L"t{a_1,a_2;i_1,i_2}:A-N-S", v) ==
          -at(L"t{a_2,a_1;i_1,i_2}:A-N-S", v));
  REQUIRE(at(L"t{a_1,a_2;i_1,i_2}:S-N-S", v) ==
          at(L"t{a_2,a_1;i_1,i_2}:S-N-S", v));
  REQUIRE(at(L"t{a_1,a_2;i_1,i_2}:N-N-S", v) ==
          at(L"t{a_2,a_1;i_2,i_1}:N-N-S", v));
  REQUIRE(at(L"t{a_1,a_2;i_1,i_2}:N-N-S", v) !=
          at(L"t{a_2,a_1;i_1,i_2}:N-N-S", v));
  const auto h = at(L"g{i_1,i_2;a_1,a_2}:A-C-S", v);
  REQUIRE(h.imag() != 0);
  REQUIRE(at(L"g{a_1,a_2;i_1,i_2}:A-C-S", v) == std::conj(h));
  REQUIRE(at(L"t⁺{i_1;a_1}:N-N-N", v) == std::conj(at(L"t{a_1;i_1}:N-N-N", v)));
  REQUIRE(at(L"x", {}, {.real = {L"x"}}).imag() == 0);
  REQUIRE(at(L"δ{i_1;i_2}", {{L"i_1", 1}, {L"i_2", 1}}) == 1.);
  REQUIRE(at(L"s{i_1;i_2}", {{L"i_1", 1}, {L"i_2", 0}}) == 0.);
  // protoindices are addressed by value
  const t::Values pairs{{L"a_1", 2}, {L"a_2", 2}, {L"i_1", 0}, {L"i_2", 1}};
  REQUIRE(at(L"s{a_1<i_1>;a_2<i_2>}:N-C-S", pairs) != 1.);
  REQUIRE(at(L"s{a_1<i_1>;a_2<i_1>}:N-C-S", pairs) == 1.);

  auto T = [](std::wstring_view s) { return deserialize(s); };
  // a hand-written derivative the judge accepts, and wrong ones it rejects
  const auto f = T(L"g{i_1;a_1}:N-N-N t{a_1;i_1}:N-N-N");
  const auto x = T(L"t{a_2;i_2}:N-N-N");
  REQUIRE(t::compare_tangents(f, T(L"g{i_2;a_2}:N-N-N"), x).difference.empty());
  REQUIRE_FALSE(
      t::compare_tangents(f, T(L"2 g{i_2;a_2}:N-N-N"), x).difference.empty());
  // a result whose free indices are not those of f and x is invalid
  REQUIRE(t::validity(f, x, T(L"g{i_2;a_2}:N-N-N")).empty());
  REQUIRE_FALSE(t::validity(f, x, T(L"g{i_3;a_2}:N-N-N")).empty());
  // a label carried by an index of a projector is free like one carried by
  // any other free index
  {
    const auto fp =
        T(L"u{;i_1}:N-N-N M{a_3;a_1<i_1>,a_2<i_1>}:N-N-N "
          L"U{a_1<i_1>,a_2<i_1>;a_3}:A-N-S v{i_1;}:N-N-N");
    const auto xp = T(L"U{a_5<i_5>,a_6<i_5>;a_7}:A-N-S");
    REQUIRE(t::validity(fp, xp,
                        T(L"2 Â{a_5<i_5>,a_6<i_5>;}:A-N-S u{;i_5}:N-N-N "
                          L"M{a_7;a_5<i_5>,a_6<i_5>}:N-N-N v{i_5;}:N-N-N"))
                .empty());
    // ... so one that is not kept is rejected, even where it stands in two
    // slots
    REQUIRE(t::validity(fp, xp,
                        T(L"2 Â{a_5<i_5>,a_6<i_5>;a_8<i_6>}:A-N-S "
                          L"u{;i_5,i_6}:N-N-N "
                          L"M{a_7,a_8<i_6>;a_5<i_5>,a_6<i_5>}:N-N-N "
                          L"v{i_6;}:N-N-N"))
                .starts_with("(d) i_6 is free and not kept"));
    // differentiated again: the labels f's own projector carries may stay
    // free, a foreign one may not
    const auto f2p =
        T(L"2 Â{a_5<i_5>,a_6<i_5>;}:A-N-S u{;i_5}:N-N-N "
          L"M{a_7;a_5<i_5>,a_6<i_5>}:N-N-N v{i_5;}:N-N-N");
    const auto x2p = T(L"u{;i_9}:N-N-N");
    const auto r2p =
        T(L"2 Â{a_5<i_5>,a_6<i_5>;}:A-N-S M{a_7;a_5<i_5>,a_6<i_5>}:N-N-N "
          L"v{i_5;}:N-N-N δ{i_9;i_5}:N-S-S");
    REQUIRE(t::compare_tangents(f2p, r2p, x2p).difference.empty());
    REQUIRE(t::validity(f2p, x2p, r2p).empty());
    REQUIRE(t::validity(f2p, x2p,
                        T(L"2 Â{a_5<i_5>,a_6<i_5>;a_8<i_8>}:A-N-S "
                          L"M{a_7,a_8<i_8>;a_5<i_5>,a_6<i_5>}:N-N-N "
                          L"v{i_5;}:N-N-N δ{i_9;i_5}:N-S-S"))
                .starts_with("(d) i_8 is free and not kept"));
  }
  // antisymmetric: the sum over the four operations
  const auto f2 = T(L"1/4 g{i_1,i_2;a_1,a_2}:A-C-S t{a_1,a_2;i_1,i_2}:A-N-S");
  const auto x2 = T(L"t{a_3,a_4;i_3,i_4}:A-N-S");
  REQUIRE(t::compare_tangents(
              f2, T(L"Â{a_3,a_4;i_3,i_4}:A-N-S g{i_3,i_4;a_3,a_4}:A-C-S"), x2)
              .difference.empty());
  REQUIRE_FALSE(t::compare_tangents(f2, T(L"1/4 g{i_3,i_4;a_3,a_4}:A-C-S"), x2)
                    .difference.empty());
  // a result without the symmetry of x is rejected
  REQUIRE_FALSE(t::compare_tangents(f2, T(L"g{i_3,i_4;a_3,a_4}:N-N-N"), x2)
                    .difference.empty());
  REQUIRE(t::group_order(x2) == 4);
  REQUIRE(t::block_order(x2) == 4);
  REQUIRE(t::group_order(T(L"g{i_1,a_1;a_2,i_2}:N-C-S")) == 2);
  REQUIRE(t::block_order(T(L"g{i_1,a_1;a_2,i_2}:N-C-S")) == 1);
  REQUIRE(t::group_order(T(L"g{i_1,a_1;i_2,i_3}:A-C-S")) == 4);
  REQUIRE(t::block_order(T(L"g{i_1,a_1;i_2,i_3}:A-C-S")) == 2);
  REQUIRE(t::group_order(x) == 1);
  // bra and ket of unequal ranks, or of rank 0, exchanged address the same
  // stored element: the exchange is no operation
  for (const auto& [fu, xu, accepted, rejected] :
       {std::tuple{L"c{i_1,i_2;}:N-N-N u{;i_1,i_2}:N-S-N", L"u{;i_3,i_4}:N-S-N",
                   L"c{i_3,i_4;}:N-N-N", L"2 c{i_3,i_4;}:N-N-N"},
        std::tuple{L"c{i_1;}:N-N-N u{;i_1}:N-S-N", L"u{;i_3}:N-S-N",
                   L"c{i_3;}:N-N-N", L"2 c{i_3;}:N-N-N"},
        std::tuple{L"c{;;x_1}:N-N-N u{;;x_1}:N-S-N", L"u{;;x_2}:N-S-N",
                   L"c{;;x_2}:N-N-N", L"2 c{;;x_2}:N-N-N"}}) {
    INFO(toUtf8(fu));
    REQUIRE(t::group_order(T(xu)) == 1);
    REQUIRE(t::compare_tangents(T(fu), T(accepted), T(xu)).difference.empty());
    REQUIRE_FALSE(
        t::compare_tangents(T(fu), T(rejected), T(xu)).difference.empty());
  }
  // a per-pair rotation: the two orders of its pair labels are operations;
  // the diagonal pair counts twice
  const auto fk =
      T(L"u{;i_1,i_2}:N-N-N M{a_3<i_1,i_2>;a_1<i_1,i_2>}:N-N-N "
        L"k{a_1<i_1,i_2>;a_3<i_1,i_2>}:N-N-N v{i_1,i_2;}:N-N-N");
  const auto xk = T(L"k{a_2<i_5,i_6>;a_4<i_5,i_6>}:N-N-N");
  REQUIRE(t::group_order(xk) == 2);
  REQUIRE(t::compare_tangents(
              fk,
              T(L"u{;i_5,i_6}:N-N-N M{a_4<i_5,i_6>;a_2<i_5,i_6>}:N-N-N "
                L"v{i_5,i_6;}:N-N-N + "
                L"u{;i_6,i_5}:N-N-N M{a_4<i_5,i_6>;a_2<i_5,i_6>}:N-N-N "
                L"v{i_6,i_5;}:N-N-N"),
              xk)
              .difference.empty());
  REQUIRE_FALSE(t::compare_tangents(
                    fk,
                    T(L"u{;i_5,i_6}:N-N-N M{a_4<i_5,i_6>;a_2<i_5,i_6>}:N-N-N "
                      L"v{i_5,i_6;}:N-N-N"),
                    xk)
                    .difference.empty());
  // a sample is a draw of the full comparison
  const auto r = T(L"g{i_2;a_2}:N-N-N");
  const auto all = t::compare_tangents(f, r, x);
  REQUIRE(all.elements == 6);
  REQUIRE(all.every_element);
  const auto drawn = t::compare_tangents(f, r, x, {}, {.count = 4, .seed = 7});
  REQUIRE(drawn.elements == 4);
  REQUIRE_FALSE(drawn.every_element);
  // the element a sample finds wrong is the same for one seed and depends on
  // it: a result wrong where i_1 = i_2 only
  const auto fs = T(L"t{a_1;i_1}:N-N-N");
  const auto wrong =
      T(L"δ{a_1;a_2}:N-S-S δ{i_2;i_1}:N-S-S + "
        L"δ{i_2;i_1}:N-S-S h{a_1;a_2}:N-N-N");
  const auto sample = [&](std::uint64_t seed) {
    return t::compare_tangents(fs, wrong, x, {}, {.count = 4, .seed = seed})
        .difference;
  };
  REQUIRE_FALSE(sample(7).empty());
  REQUIRE(sample(7) == sample(7));
  REQUIRE(sample(7) != sample(8));
  // the dummies of f named like x's indices are renamed apart from every
  // index of f, also where f holds the next temporary names
  const auto minted = [](const Index& like) {
    Index::reset_tmp_index();
    const auto first = Index::make_tmp_index(like.basis());
    return std::pair{first, Index::make_tmp_index(like.basis())};
  };
  const auto [i_tmp1, i_tmp2] = minted(Index(L"i_1"));
  const auto [a_tmp1, a_tmp2] = minted(Index(L"a_1"));
  const container::map<Index, Index> to_minted{{Index(L"i_7"), i_tmp1},
                                               {Index(L"i_8"), i_tmp2},
                                               {Index(L"a_7"), a_tmp1},
                                               {Index(L"a_8"), a_tmp2}};
  const auto uv = L" u{i_7,i_8;a_7,a_8}:N-N-N v{a_7,a_8;i_7,i_8}:N-N-N";
  const auto ft = transform_expr(
      T(std::wstring(L"g{i_1;a_1}:N-N-N t{a_1;i_1}:N-N-N") + uv), to_minted);
  const auto rt =
      transform_expr(T(std::wstring(L"g{i_1;a_1}:N-N-N") + uv), to_minted);
  const auto xt = T(L"t{a_1;i_1}:N-N-N");
  Index::reset_tmp_index();
  REQUIRE(t::compare_tangents(ft, rt, xt).difference.empty());
  Index::reset_tmp_index();
  REQUIRE_FALSE(
      t::compare_tangents(ft, ex<Constant>(2) * rt, xt).difference.empty());
}
