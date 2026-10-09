//
// The integral projection pass (mbpt/rules/pno.hpp) on terms SeQuant's CC
// derivations produce.
//

#include <catch2/catch_test_macros.hpp>
#include "catch2_sequant.hpp"
#include "csv_test_utils.hpp"

#include <SeQuant/core/basis.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/reserved.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/expr.hpp>
#include <SeQuant/core/utility/indices.hpp>
#include <SeQuant/core/utility/string.hpp>
#include <SeQuant/domain/mbpt/rules/df.hpp>
#include <SeQuant/domain/mbpt/rules/pno.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <functional>
#include <iterator>
#include <map>
#include <optional>
#include <regex>
#include <set>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

namespace sequant::tests::csv_rules {

using csv::tensors_labelled;

using mbpt::ProjectionBasis;
using mbpt::ProjectionDomain;
using mbpt::ProjectionOptions;
using mbpt::ProjectionTerms;

/// exchange cells -> @p x, Coulomb cells -> @p c
inline auto class_map(IndexBasis::instance_type x,
                      IndexBasis::instance_type c) {
  return [x, c](ProjectionTerms cell) {
    return mbpt::contains(ProjectionTerms::Exchange, cell) ? x : c;
  };
}

/// {terms, OwnPair, Integral, exchange -> 1, Coulomb -> 2}
inline ProjectionOptions own_pair(ProjectionTerms terms) {
  return {.terms = terms,
          .cell_instance = class_map(test_csv::exchange, test_csv::coulomb)};
}

/// project_integral_domains() of @p expr with tmp numbering reset: the
/// minted legs start at a_100
inline ExprPtr project(ExprPtr const& expr, ProjectionOptions const& opts) {
  Index::reset_tmp_index();
  return mbpt::project_integral_domains(expr, opts);
}

/// serialize() without the tensor symmetry annotations
inline std::wstring terse(ExprPtr const& expr) {
  return std::regex_replace(serialize(expr), std::wregex(L":\\w-\\w-\\w"), L"");
}

inline std::size_t overlaps(ExprPtr const& expr) {
  return tensors_labelled(expr, reserved::overlap_label()).size();
}

/// the number of `g` slots that differ between @p before and @p after
inline std::size_t moved_legs(ExprPtr const& before, ExprPtr const& after) {
  const auto b = tensors_labelled(before, L"g");
  const auto a = tensors_labelled(after, L"g");
  REQUIRE(a.size() == b.size());
  std::size_t n = 0;
  for (std::size_t k = 0; k != a.size(); ++k)
    for (std::size_t s = 0; s != a[k]->_num_slots(); ++s)
      if (!(a[k]->_slots()[s] == b[k]->_slots()[s])) ++n;
  return n;
}

/// checks that the overlap of each leg the pass minted on a `g` of @p expr
/// holds it on the other side than the `g` does
/// @return the number of such legs in a `g` bra
inline std::size_t minted_bra_legs(ExprPtr const& expr) {
  const auto metrics = tensors_labelled(expr, reserved::overlap_label());
  auto overlaps_with = [&metrics](Index const& idx, bool in_bra) {
    return std::ranges::count_if(metrics, [&](auto const* s) {
      return (in_bra ? s->_bra() : s->_ket())[0] == idx;
    });
  };
  std::size_t n = 0;
  for (AbstractTensor const* g : tensors_labelled(expr, L"g"))
    for (bool bra_leg : {true, false})
      for (Index const& idx : bra_leg ? g->_bra() : g->_ket()) {
        if (idx.ordinal() < Index::min_tmp_index()) continue;
        n += bra_leg;
        CHECK(overlaps_with(idx, !bra_leg) == 1);
        CHECK(overlaps_with(idx, bra_leg) == 0);
      }
  return n;
}

/// the ordinals of the terms of @p eq that @p opts changes
inline std::set<std::size_t> firing(ExprPtr const& eq,
                                    ProjectionOptions const& opts) {
  std::set<std::size_t> result;
  for (std::size_t t = 0; t != eq->size(); ++t)
    if (project(eq->at(t), opts).get() != eq->at(t).get()) result.insert(t);
  return result;
}

/// the serialized tensors of @p expr labelled @p label, sorted
inline std::vector<std::wstring> serialized(ExprPtr const& expr,
                                            std::wstring_view label) {
  std::vector<std::wstring> result;
  expr->visit(
      [&](ExprPtr const& x) {
        if (x->is<AbstractTensor>() &&
            std::wstring_view(x->as<AbstractTensor>()._label()) == label)
          result.push_back(serialize(x));
      },
      /* atoms_only = */ true);
  std::ranges::sort(result);
  return result;
}

/// @p expr with every index instance @p from[k] replaced by @p to[k]
inline ExprPtr relabel_instances(
    ExprPtr const& expr,
    std::map<IndexBasis::instance_type, IndexBasis::instance_type> const& to) {
  container::map<Index, Index> m;
  for (Index const& idx : get_used_indices(expr))
    if (auto const& inst = idx.basis().basis_instance();
        inst && to.contains(*inst))
      m.emplace(idx, idx.replace_basis_instance(to.at(*inst)));
  return transform_expr(expr, m);
}

/// checks that the terms of @p eq named by @p ordinals have the Mulliken class
/// and order (linear: <= 1 in the amplitudes) that @p scope fixes
inline void check_cell(ExprPtr const& eq, std::set<std::size_t> const& ordinals,
                       ProjectionTerms scope) {
  using enum ProjectionTerms;
  if (!mbpt::any(scope)) return;  // None fixes nothing
  const bool coulomb = mbpt::contains(Coulomb, scope);
  for (std::size_t ordinal : ordinals) {
    CAPTURE(ordinal);
    const auto gs = tensors_labelled(eq->at(ordinal), L"g");
    if (coulomb || mbpt::contains(Exchange, scope))
      CHECK(std::ranges::any_of(gs, [coulomb](auto const* g) {
        return coulomb ? mbpt::is_oovv_integral(*g)
                       : mbpt::is_ovov_integral(*g);
      }));
    const bool linear = tensors_labelled(eq->at(ordinal), L"t").size() <= 1;
    if (mbpt::contains(Linear, scope)) CHECK(linear);
    if (mbpt::contains(Nonlinear, scope)) CHECK_FALSE(linear);
  }
}

}  // namespace sequant::tests::csv_rules

TEST_CASE("csv-a-metric-between-two-amplitude-bases-stands",
          "[mbpt][csv][valgrind_skip]") {
  using namespace sequant;
  using namespace sequant::tests::csv;
  using namespace sequant::tests::csv_rules;

  // Λ-CCSD with t and λ in two families: its standing metrics join λ, t and
  // projector slots; doubles t/λ are (ov|ov)-shaped, but only g is projected
  auto ctx = set_scoped_default_context(csv_cc_context());
  const auto reg = granted_registry({{L"t", 1}, {L"λ", 2}});
  ScopedCsvContext scoped{reg};
  ExprPtr const r2 = derive_λ(reg).at(2);
  for (auto const& opts : {own_pair(ProjectionTerms::All),
                           ProjectionOptions{.terms = ProjectionTerms::All,
                           .basis = ProjectionBasis::Amplitude}}) {
    std::size_t changed = 0;
    for (auto const& term : *r2) {
      CAPTURE(toUtf8(serialize(term)));
      const auto out = project(term, opts);
      if (out.get() != term.get()) ++changed;
      for (std::wstring_view label : {L"t", L"λ"})
        CHECK(serialized(out, label) == serialized(term, label));
      const auto s_in = serialized(term, reserved::overlap_label());
      const auto s_out = serialized(out, reserved::overlap_label());
      CHECK(std::ranges::includes(s_out, s_in));
      CHECK(s_out.size() == s_in.size() + moved_legs(term, out));
    }
    CHECK(changed > 0);
  }
}

TEST_CASE("csv-shape-predicates", "[mbpt][csv]") {
  using namespace sequant;
  using namespace sequant::tests::csv;

  auto ctx_resetter = set_scoped_default_context(csv_cc_context());
  // {tensor, (ov|ov), (oo|vv)}; Dirac <p1 p2|p3 p4> is Mulliken (p1 p3|p2 p4)
  const std::vector<std::tuple<std::wstring, bool, bool>> cases{
      {L"g{i1,a1;i2,a2}", false, true},  // Coulomb
      {L"g{i1,a1;a2,i2}", true, false},  // exchange over the same spaces
      {L"g{a1,a2;i1,i2}", true, false},
      {L"g{i1,i2;a1,a2}", true, false},
      {L"g{i1,i2;i3,i4}", false, false},  // (oo|oo) ring
      {L"g{a1,a2;a3,a4}", false, false},  // (vv|vv) ladder
      {L"g{i1,a1;a2,a3}", false, false},
      {L"g{i1,i2;i3,a1}", false, false},
      {L"g{i1;a1}", false, false},
      {L"g{i1,i2,i3;a1,a2,a3}", false, false},
      {L"g{i1,i2;a1}", false, false},
      {L"g{i1;a1,a2}", false, false},
  };
  for (auto const& [tensor, ovov, oovv] : cases) {
    CAPTURE(toUtf8(tensor));
    ExprPtr owner = deserialize(tensor);
    CHECK(mbpt::is_ovov_integral(owner->as<Tensor>()) == ovov);
    CHECK(mbpt::is_oovv_integral(owner->as<Tensor>()) == oovv);
  }
}

TEST_CASE("csv-projection-partner-pair-integral-basis",
          "[mbpt][csv][valgrind_skip]") {
  using namespace sequant;
  using namespace sequant::tests::csv;
  using namespace sequant::tests::csv_rules;

  ScopedCsvContext scoped;
  const auto st = derive_t(granted_registry({}));
  REQUIRE(st.size() == 3);
  const ProjectionOptions opts{
      .terms = ProjectionTerms::All,
      .domain = ProjectionDomain::PartnerPair,
      .cell_instance = class_map(test_csv::exchange, test_csv::coulomb)};

  // neither the energy nor R1 is ever projected
  CHECK(project(st[0], opts).get() == st[0].get());
  CHECK(project(st[1], opts).get() == st[1].get());

  ExprPtr const& r2 = st[2];
  CHECK(firing(r2, opts).size() == 25 + 3);
  const auto out = project(r2, opts);
  std::map<IndexBasis::instance_type, std::size_t> stamped;
  for (std::size_t t = 0; t != out->size(); ++t) {
    ExprPtr const& term = out->at(t);
    CAPTURE(t);
    const auto externals = projector_slots(term);
    CHECK(projector_slots(r2->at(t)) == externals);
    const auto metrics = tensors_labelled(term, reserved::overlap_label());
    for (AbstractTensor const* g : tensors_labelled(term, L"g")) {
      const bool ovov = mbpt::is_ovov_integral(*g);
      if (!ovov && !mbpt::is_oovv_integral(*g)) continue;
      std::set<IndexBasis::instance_type> seen;
      for (Index const& idx : g->_slots()) {
        if (mbpt::is_occupied(idx)) {
          CHECK(!idx.basis().has_basis_instance());
          continue;
        }
        if (std::ranges::count(externals, idx)) {
          CHECK(!idx.basis().has_basis_instance());
          continue;
        }
        REQUIRE(idx.basis().has_basis_instance());
        CHECK(*idx.basis().basis_instance() ==
              (ovov ? test_csv::exchange : test_csv::coulomb));
        seen.insert(*idx.basis().basis_instance());
        // one overlap to its partner, which is on the same pair
        CHECK(std::ranges::count_if(metrics, [&idx](auto const* s) {
                return s->_bra()[0] == idx &&
                       s->_ket()[0].proto_indices() == idx.proto_indices() &&
                       !s->_ket()[0].basis().has_basis_instance();
              }) == 1);
      }
      for (auto inst : seen) ++stamped[inst];
    }
  }
  CHECK(stamped == std::map<IndexBasis::instance_type, std::size_t>{
                       {test_csv::exchange, 25}, {test_csv::coulomb, 3}});
  CHECK(overlaps(out) == overlaps(r2) + moved_legs(r2, out));
  for (auto const& term : *out) CHECK(minted_bra_legs(term) == 0);
  CHECK(mbpt::project_integral_domains(out, opts).get() == out.get());

  REQUIRE_THROWS_MATCHES(project(r2, {.terms = ProjectionTerms::All,
                                      .domain = ProjectionDomain::PartnerPair,
                                      .basis = ProjectionBasis::Amplitude}),
                         sequant::Exception,
                         message_contains("PartnerPair with"));
}

TEST_CASE("csv-pno-projection-fixtures", "[mbpt][csv]") {
  using namespace sequant;
  using namespace sequant::tests::csv;
  using namespace sequant::tests::csv_rules;
  using enum mbpt::ProjectionTerms;

  ScopedCsvContext scoped;
  // CSV-CCSD terms as derive_t() leaves them (term numbers of its R2)
  const std::wstring S2 = L"Ŝ{i_1,i_2;a_1<i_1,i_2>,a_2<i_1,i_2>} * ";
  const std::wstring tt_ss =
      L" * t{a_3<i_1,i_2>,a_4<i_1,i_2>;i_1,i_2}"
      L" * t{a_5<i_3,i_4>,a_6<i_3,i_4>;i_3,i_4}"
      L" * s{a_1<i_1,i_2>;a_5<i_3,i_4>} * s{a_2<i_1,i_2>;a_6<i_3,i_4>}";
  // #52: (ov|ov) on the foreign pair {i_1,i_2}, own pair {i_3,i_4}, order 2
  const std::wstring g_foreign = L"g{i_3,i_4;a_3<i_1,i_2>,a_4<i_1,i_2>}";
  const std::wstring foreign = S2 + g_foreign + tt_ss;
  // #15: linear, the external a_1 on the integral
  const std::wstring linear =
      L"-2 " + S2 +
      L"g{i_3,a_1<i_1,i_2>;a_3<i_2,i_3>,i_1}"
      L" * t{a_3<i_2,i_3>,a_4<i_2,i_3>;i_2,i_3} * s{a_2<i_1,i_2>;a_4<i_2,i_3>}";
  // #25: rank-1 partners
  const std::wstring rank1 =
      S2 +
      L"g{i_3,i_4;a_3<i_1>,a_4<i_2>} * t{a_3<i_1>;i_1} * t{a_4<i_2>;i_2}"
      L" * t{a_5<i_3,i_4>,a_6<i_3,i_4>;i_3,i_4}"
      L" * s{a_1<i_1,i_2>;a_5<i_3,i_4>} * s{a_2<i_1,i_2>;a_6<i_3,i_4>}";
  // #45: a_3 already on the own pair
  const std::wstring half =
      L"-4 " + S2 +
      L"g{i_3,i_4;a_3<i_3,i_4>,a_4<i_1,i_2>}"
      L" * t{a_1<i_1,i_2>,a_4<i_1,i_2>;i_1,i_2}"
      L" * t{a_3<i_3,i_4>,a_5<i_3,i_4>;i_3,i_4} * s{a_2<i_1,i_2>;a_5<i_3,i_4>}";
  // #18: (vv|vv) ladder
  const std::wstring ladder =
      S2 +
      L"g{a_1<i_1,i_2>,a_2<i_1,i_2>;a_3<i_1,i_2>,a_4<i_1,i_2>}"
      L" * t{a_3<i_1,i_2>,a_4<i_1,i_2>;i_1,i_2}";
  // E #0 and R1 #14: (ov|ov) integrals with legs off their own pairs
  const std::wstring energy =
      L"2 g{i_1,i_2;a_1<i_1>,a_2<i_2>} * t{a_1<i_1>;i_1} * t{a_2<i_2>;i_2}";
  const std::wstring singles =
      L"-2 Ŝ{i_1;a_1<i_1>} * g{i_2,i_3;a_2<i_2>,a_3<i_1>} * t{a_2<i_2>;i_2}"
      L" * t{a_3<i_1>;i_1} * t{a_4<i_3>;i_3} * s{a_1<i_1>;a_4<i_3>}";

  const auto amplitude = [](ProjectionTerms terms) {
    return ProjectionOptions{.terms = terms,
                             .basis = ProjectionBasis::Amplitude};
  };

  // {input, options, expected; `#` = the minted legs' instance suffix}
  const std::wstring moved_foreign =
      S2 +
      L"g{i_3,i_4;a_100<i_3,i_4#>,a_101<i_3,i_4#>}"
      L" * s{a_100<i_3,i_4#>;a_3<i_1,i_2>} * s{a_101<i_3,i_4#>;a_4<i_1,i_2>}" +
      tt_ss;
  const std::vector<std::tuple<char const*, std::wstring, std::wstring>> rows{
      {"foreign pair", foreign, moved_foreign},
      {"the external stays, the contracted leg moves", linear,
       L"-2 " + S2 +
           L"g{i_3,a_1<i_1,i_2>;a_100<i_1,i_3#>,i_1}"
           L" * s{a_100<i_1,i_3#>;a_3<i_2,i_3>}"
           L" * t{a_3<i_2,i_3>,a_4<i_2,i_3>;i_2,i_3}"
           L" * s{a_2<i_1,i_2>;a_4<i_2,i_3>}"},
      {"rank-1 partners move onto the rank-2 own pair", rank1,
       S2 + L"g{i_3,i_4;a_100<i_3,i_4#>,a_101<i_3,i_4#>}"
            L" * s{a_100<i_3,i_4#>;a_3<i_1>} * s{a_101<i_3,i_4#>;a_4<i_2>}"
            L" * t{a_3<i_1>;i_1} * t{a_4<i_2>;i_2}"
            L" * t{a_5<i_3,i_4>,a_6<i_3,i_4>;i_3,i_4}"
            L" * s{a_1<i_1,i_2>;a_5<i_3,i_4>} * s{a_2<i_1,i_2>;a_6<i_3,i_4>}"},
  };
  for (auto const& [what, in, expected] : rows) {
    CAPTURE(what);
    const ExprPtr input = deserialize(in);
    for (auto const& [opts, suffix] :
         {std::pair{own_pair(All), std::wstring(L";1")},
          std::pair{amplitude(All), std::wstring()}}) {
      const auto out = project(input, opts);
      CHECK(terse(out) ==
            std::regex_replace(expected, std::wregex(L"#"), suffix));
      CHECK(overlaps(out) == overlaps(input) + moved_legs(input, out));
      // idempotent; on the own pair, `Amplitude` keeps the Integral stamp
      CHECK(mbpt::project_integral_domains(out, opts).get() == out.get());
      CHECK(mbpt::project_integral_domains(out, amplitude(All)).get() ==
            out.get());
    }
  }

  // the cell decides: order and class
  for (auto const& [in, cell] : {std::pair{foreign, ExchangeNonlinear},
                                 std::pair{linear, ExchangeLinear},
                                 std::pair{rank1, ExchangeNonlinear}}) {
    const ExprPtr input = deserialize(in);
    for (auto scope :
         {ExchangeLinear, ExchangeNonlinear, CoulombLinear, CoulombNonlinear}) {
      CAPTURE(toUtf8(in), static_cast<int>(scope));
      CHECK((project(input, own_pair(scope)).get() != input.get()) ==
            (scope == cell));
    }
  }

  // a leg already on the own pair stays; only the other one moves
  {
    const ExprPtr input = deserialize(half);
    const auto out = project(input, amplitude(All));
    CHECK(terse(out) == L"-4 " + S2 +
                            L"g{i_3,i_4;a_3<i_3,i_4>,a_100<i_3,i_4>}"
                            L" * s{a_100<i_3,i_4>;a_4<i_1,i_2>}"
                            L" * t{a_1<i_1,i_2>,a_4<i_1,i_2>;i_1,i_2}"
                            L" * t{a_3<i_3,i_4>,a_5<i_3,i_4>;i_3,i_4}"
                            L" * s{a_2<i_1,i_2>;a_5<i_3,i_4>}");
    CHECK(overlaps(out) == overlaps(input) + 1);
  }

  // no option reaches the energy, a singles residual or a (vv|vv) ladder
  for (std::wstring const& text : {energy, singles, ladder}) {
    CAPTURE(toUtf8(text));
    const ExprPtr input = deserialize(text);
    for (auto const& opts : {own_pair(All), amplitude(All),
                             ProjectionOptions{.terms = All,
                             .domain = ProjectionDomain::PartnerPair,
                             .cell_instance = class_map(1, 2)}})
      CHECK(project(input, opts).get() == input.get());
  }

  // overlap reduction is contraction-time only, so simplify keeps them
  ExprPtr simplified = project(deserialize(foreign), own_pair(All));
  simplify(simplified);
  CHECK(standing_metrics(simplified) == 2);

  // a Sum nested in a Product is rewritten in place; the legs are shared
  {
    const ExprPtr input = ex<Product>(
        1,
        ExprPtrList{deserialize(S2.substr(0, S2.size() - 3)),
                    ex<Sum>(ExprPtrList{deserialize(g_foreign),
                                        deserialize(L"g{i_3,i_4;a_4<i_1,i_2>,"
                                                    L"a_3<i_1,i_2>}")}),
                    deserialize(tt_ss.substr(3))});
    CHECK(terse(project(input, own_pair(All))) ==
          L"Ŝ{i_1,i_2;a_1<i_1,i_2>,a_2<i_1,i_2>} * "
          L"(g{i_3,i_4;a_100<i_3,i_4;1>,a_101<i_3,i_4;1>}"
          L" * s{a_100<i_3,i_4;1>;a_3<i_1,i_2>}"
          L" * s{a_101<i_3,i_4;1>;a_4<i_1,i_2>}"
          L" + g{i_3,i_4;a_101<i_3,i_4;1>,a_100<i_3,i_4;1>}"
          L" * s{a_101<i_3,i_4;1>;a_4<i_1,i_2>}"
          L" * s{a_100<i_3,i_4;1>;a_3<i_1,i_2>})" +
              tt_ss);
  }

  // a leg shared by two integrals of one term moves once per integral, each
  // to its own fresh index; the overlaps chain through the leg
  {
    const std::wstring shared = S2 +
                                L"g{i_3,a_4<i_1,i_2>;a_3<i_1,i_2>,i_4} * "
                                L"g{i_5,i_6;a_4<i_1,i_2>,a_5<i_1,i_2>}"
                                L" * t{a_3<i_1,i_2>,a_5<i_1,i_2>;i_1,i_2}";
    const ExprPtr input = deserialize(shared);
    const auto out = project(input, own_pair(All));
    CHECK(terse(out) == S2 + L"g{i_3,a_100<i_3,i_4;1>;a_101<i_3,i_4;1>,i_4}"
                             L" * s{a_4<i_1,i_2>;a_100<i_3,i_4;1>}"
                             L" * s{a_101<i_3,i_4;1>;a_3<i_1,i_2>}"
                             L" * g{i_5,i_6;a_102<i_5,i_6;1>,a_103<i_5,i_6;1>}"
                             L" * s{a_102<i_5,i_6;1>;a_4<i_1,i_2>}"
                             L" * s{a_103<i_5,i_6;1>;a_5<i_1,i_2>}"
                             L" * t{a_3<i_1,i_2>,a_5<i_1,i_2>;i_1,i_2}");
    CHECK(overlaps(out) == 4);
    CHECK(minted_bra_legs(out) == 1);
    CHECK(mbpt::project_integral_domains(out, own_pair(All)).get() ==
          out.get());
  }

  // malformed input throws: density fitting already ran
  const IndexSpace aux =
      get_default_context().index_basis_registry()->retrieve(L"Κ");
  REQUIRE_THROWS_MATCHES(
      project(mbpt::density_fit(deserialize(foreign), aux, L"g", L"g"),
              own_pair(All)),
      sequant::Exception, message_contains("BEFORE density fitting"));
  // `Integral` without a cell -> instance map
  REQUIRE_THROWS_MATCHES(project(deserialize(foreign), {.terms = All}),
                         sequant::Exception,
                         message_contains("needs a cell_instance"));
}

TEST_CASE("csv-pno-projection-ccsd", "[mbpt][csv][valgrind_skip]") {
  using namespace sequant;
  using namespace sequant::tests::csv;
  using namespace sequant::tests::csv_rules;
  using enum mbpt::ProjectionTerms;

  ScopedCsvContext scoped;
  const auto st = derive_t(granted_registry({}));
  REQUIRE(st.size() == 3);

  // the number of R2 terms each scope fires on
  const std::map<ProjectionTerms, std::size_t> cells{
      {None, 0},          {ExchangeLinear, 2},   {ExchangeNonlinear, 23},
      {CoulombLinear, 2}, {CoulombNonlinear, 1}, {Exchange, 25},
      {Coulomb, 3},       {Linear, 4},           {Nonlinear, 24},
      {All, 28}};

  std::map<ProjectionTerms, std::set<std::size_t>> kept;
  for (auto const& [scope, size] : cells) {
    CAPTURE(static_cast<int>(scope));
    for (auto basis : {ProjectionBasis::Integral, ProjectionBasis::Amplitude}) {
      CAPTURE(static_cast<int>(basis));
      auto opts = own_pair(scope);
      opts.basis = basis;
      // neither the energy nor R1 is ever projected
      CHECK(project(st[0], opts).get() == st[0].get());
      CHECK(project(st[1], opts).get() == st[1].get());
      const auto fired = firing(st[2], opts);
      CHECK(fired.size() == size);
      if (basis == ProjectionBasis::Integral)
        kept[scope] = fired;
      else
        CHECK(fired == kept[scope]);
      // an insertion: one overlap per moved leg
      const auto out = project(st[2], opts);
      CHECK(overlaps(out) == overlaps(st[2]) + moved_legs(st[2], out));
      // t's legs face g's ket: every overlap is s{x';x}
      for (auto const& term : *out) CHECK(minted_bra_legs(term) == 0);
      CHECK(mbpt::project_integral_domains(out, opts).get() == out.get());
    }
    check_cell(st[2], kept[scope], scope);
  }

  // the partition: each pair is disjoint and unites to the third
  for (auto [a, b, both] :
       {std::tuple{ExchangeLinear, ExchangeNonlinear, Exchange},
        std::tuple{CoulombLinear, CoulombNonlinear, Coulomb},
        std::tuple{ExchangeLinear, CoulombLinear, Linear},
        std::tuple{ExchangeNonlinear, CoulombNonlinear, Nonlinear},
        std::tuple{Exchange, Coulomb, All},
        std::tuple{Linear, Nonlinear, All}}) {
    CAPTURE(static_cast<int>(a), static_cast<int>(b));
    std::set<std::size_t> uni, common;
    std::ranges::set_union(kept[a], kept[b], std::inserter(uni, uni.end()));
    std::ranges::set_intersection(kept[a], kept[b],
                                  std::inserter(common, common.end()));
    CHECK(uni == kept[both]);
    CHECK(common.empty());
  }
}

TEST_CASE("csv-pno-mixed-class-and-opaque-ordinals",
          "[mbpt][csv][valgrind_skip]") {
  using namespace sequant;
  using namespace sequant::tests::csv;
  using namespace sequant::tests::csv_rules;

  ScopedCsvContext scoped;
  ExprPtr const r2 = derive_t(granted_registry({})).at(2);
  const auto x_terms = firing(r2, own_pair(ProjectionTerms::ExchangeLinear));
  const auto c_terms = firing(r2, own_pair(ProjectionTerms::CoulombLinear));
  REQUIRE(!x_terms.empty());
  REQUIRE(!c_terms.empty());
  // {probe ordinal, its class's instance}
  const std::array<std::pair<std::size_t, IndexBasis::instance_type>, 2> probes{
      {{*x_terms.begin(), test_csv::exchange},
       {*c_terms.begin(), test_csv::coulomb}}};

  // a mixed selection gives each class's legs its own instance
  const ExprPtr mixed =
      ex<Sum>(ExprPtrList{r2->at(probes[0].first), r2->at(probes[1].first)});
  const ExprPtr projected = project(mixed, own_pair(ProjectionTerms::All));
  for (std::size_t t = 0; t != 2; ++t) {
    CAPTURE(t);
    ExprPtr const& before = mixed->at(t);
    ExprPtr const& after = projected->at(t);
    REQUIRE(moved_legs(before, after) > 0);
    const auto externals = projector_slots(after);
    for (AbstractTensor const* g : tensors_labelled(after, L"g"))
      for (Index const& idx : g->_slots()) {
        if (mbpt::is_occupied(idx) || std::ranges::count(externals, idx))
          continue;
        CHECK(idx.basis().basis_instance() == probes[t].second);
        CHECK(idx.proto_indices().size() == 2);
      }
  }

  // an opaque map is stamped verbatim; swapping it swaps the instances only
  for (auto scope :
       {ProjectionTerms::ExchangeLinear, ProjectionTerms::ExchangeNonlinear,
        ProjectionTerms::CoulombLinear, ProjectionTerms::CoulombNonlinear,
        ProjectionTerms::All}) {
    CAPTURE(static_cast<int>(scope));
    auto with_map = [scope](IndexBasis::instance_type x,
                            IndexBasis::instance_type c) {
      return ProjectionOptions{.terms = scope,
                               .cell_instance = class_map(x, c)};
    };
    const auto shape = project(r2, with_map(1, 2));
    const auto odd = project(r2, with_map(7, 9));
    const auto swapped = project(r2, with_map(2, 1));
    CHECK(*odd == *relabel_instances(shape, {{1, 7}, {2, 9}}));
    CHECK(*swapped == *relabel_instances(shape, {{1, 2}, {2, 1}}));
    CHECK(instance_histogram(odd).count(1) == 0);
    CHECK(instance_histogram(odd).count(2) == 0);
  }

  // nothing selected: the map is never asked
  CHECK(mbpt::project_integral_domains(r2, {}).get() == r2.get());
}

TEST_CASE("csv-pno-amplitude-basis-follows-the-partner",
          "[mbpt][csv][valgrind_skip]") {
  using namespace sequant;
  using namespace sequant::tests::csv;
  using namespace sequant::tests::csv_rules;

  auto ctx = set_scoped_default_context(csv_cc_context());
  const auto reg = granted_registry({{L"t", test_csv::pno_pert}});
  ScopedCsvContext scoped{reg};
  ExprPtr const r2 = derive_t(reg).at(2);
  const ProjectionOptions opts{.terms = ProjectionTerms::ExchangeNonlinear,
                               .basis = ProjectionBasis::Amplitude};
  CHECK(firing(r2, opts).size() == 23);
  const auto out = project(r2, opts);
  const auto minted = [](Index const& idx) {
    return idx.ordinal() >= Index::min_tmp_index();
  };
  std::size_t n = 0;
  for (AbstractTensor const* g : tensors_labelled(out, L"g"))
    for (Index const& idx : g->_slots())
      if (minted(idx)) {
        ++n;
        CHECK(idx.basis().basis_instance() == test_csv::pno_pert);
      }
  CHECK(n == moved_legs(r2, out));
  CHECK(n > 0);
}

TEST_CASE("csv-pno-order-counts-every-amplitude-family",
          "[mbpt][csv][valgrind_skip]") {
  using namespace sequant;
  using namespace sequant::tests::csv;
  using namespace sequant::tests::csv_rules;
  using enum mbpt::ProjectionTerms;

  ScopedCsvContext scoped;
  ExprPtr const r2 = derive_λ(granted_registry({})).at(2);
  // a λ g t term is non-linear
  std::size_t probed = 0;
  for (auto const& term : *r2) {
    if (tensors_labelled(term, L"λ").size() != 1 ||
        tensors_labelled(term, L"t").size() != 1 ||
        project(term, own_pair(All)).get() == term.get())
      continue;
    CAPTURE(toUtf8(serialize(term)));
    ++probed;
    CHECK(project(term, own_pair(Linear)).get() == term.get());
    CHECK(project(term, own_pair(Nonlinear)).get() != term.get());
  }
  CHECK(probed > 0);
}

TEST_CASE("csv-pno-projection-lambda-canonicalizes",
          "[mbpt][csv][valgrind_skip]") {
  using namespace sequant;
  using namespace sequant::tests::csv;
  using namespace sequant::tests::csv_rules;

  // λ's virtuals are kets, so the g legs facing them are bras: their overlaps
  // must be s{x;x'} for the term to stay bra-ket strict
  using Grants =
      std::vector<std::pair<std::wstring, IndexBasis::instance_type>>;
  for (auto const& grants : {Grants{}, Grants{{L"t", 1}, {L"λ", 2}}}) {
    CAPTURE(grants.size());
    auto ctx = set_scoped_default_context(csv_cc_context());
    const auto reg = granted_registry(grants);
    ScopedCsvContext scoped{reg};
    ExprPtr const r2 = derive_λ(reg).at(2);
    for (auto domain :
         {ProjectionDomain::PartnerPair, ProjectionDomain::OwnPair}) {
      CAPTURE(static_cast<int>(domain));
      auto opts = own_pair(ProjectionTerms::All);
      opts.domain = domain;
      const auto out = project(r2, opts);
      std::size_t bra_legs = 0;
      for (auto const& term : *out) {
        CAPTURE(toUtf8(serialize(term)));
        bra_legs += minted_bra_legs(term);
        ExprPtr canonical = term->clone();
        REQUIRE_NOTHROW(canonicalize(canonical));
      }
      CHECK(bra_legs > 0);
    }
  }
}

TEST_CASE("csv-pno-projection-nested-sum-then-sibling", "[mbpt][csv]") {
  using namespace sequant;
  using namespace sequant::tests::csv;
  using namespace sequant::tests::csv_rules;

  ScopedCsvContext scoped;
  // an integral after a nested sum of integrals that move the same leg: the
  // sum's replacement is not reused, so no minted index occurs more than twice
  // in a term
  const auto g = [](std::wstring const& s) { return deserialize<ExprPtr>(s); };
  const auto g_foreign = L"g{i_3,i_4;a_3<i_1,i_2>,a_4<i_1,i_2>}";
  const ExprPtr sum = g(g_foreign) + g(L"g{i_4,i_3;a_3<i_1,i_2>,a_4<i_1,i_2>}");
  const ExprPtr term = ex<Product>(
      ExprPtrList{g(L"Ŝ{i_1,i_2;a_1<i_1,i_2>,a_2<i_1,i_2>}"), sum, g(g_foreign),
                  g(L"t{a_3<i_1,i_2>,a_4<i_1,i_2>;i_1,i_2}"),
                  g(L"t{a_5<i_3,i_4>,a_6<i_3,i_4>;i_3,i_4}")});
  auto out = project(term, own_pair(ProjectionTerms::All));
  REQUIRE(out.get() != term.get());
  expand(out);
  REQUIRE(out->is<Sum>());
  for (auto const& summand : *out) {
    CAPTURE(toUtf8(serialize(summand)));
    container::map<Index, std::size_t> counts;
    for (AbstractTensor const* t : tensors_labelled(summand, L"g"))
      for (Index const& idx : t->_slots()) ++counts[idx];
    for (AbstractTensor const* t :
         tensors_labelled(summand, reserved::overlap_label()))
      for (Index const& idx : t->_slots()) ++counts[idx];
    for (auto const& [idx, n] : counts)
      if (idx.ordinal() >= Index::min_tmp_index()) CHECK(n <= 2);
  }
}
