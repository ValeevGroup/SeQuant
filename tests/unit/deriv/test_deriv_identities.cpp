#include <SeQuant/core/container.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/deriv/deriv.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/options.hpp>
#include <SeQuant/core/utility/string.hpp>
#include <SeQuant/domain/mbpt/models/cc.hpp>
#include <SeQuant/domain/mbpt/op_registry.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>
#include <catch2/catch_test_macros.hpp>
#include "csv_test_utils.hpp"
#include "evaluator.hpp"
#include "identities.hpp"
#include "judge.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <sstream>
#include <string>
#include <string_view>

using namespace sequant;
namespace t = sequant::deriv::test;
using V = container::svector<ExprPtr>;

namespace {

ExprPtr T(std::wstring_view text) { return deserialize(text); }

/// the seed of the pair (@p model , @p equation , @p x )
std::uint64_t seed_of(const t::Model& model, std::wstring_view equation,
                      const ExprPtr& x) {
  return t::stable_hash(t::name_of(model) + toUtf8(equation) + t::text_of(x));
}

bool is_g(const ExprPtr& x) {
  return x->is<Tensor>() && x->as<Tensor>().label() == L"g";
}

/// the size of a model's family, and how many pairs each identity checked
struct Counts {
  std::size_t equations = 0, variables = 0, amplitudes = 0;
  /// I1 and I6 pairs
  std::size_t pairs = 0;
  /// I3 and I4 pairs
  std::size_t taylor = 0, schwarz = 0;
  /// the I1, I3 and I4 pairs that read a nonzero element
  std::size_t pairs_read = 0, taylor_read = 0, schwarz_read = 0;
  /// the variables diff rejects, by the text of the rejection
  container::map<std::string, std::size_t> rejected;
  bool operator==(const Counts&) const = default;
};

}  // namespace

TEST_CASE("deriv/identities/over-the-family", "[deriv][identities]") {
  // I1 and I6 over every equation of every model by every checked variable,
  // I3 and I4 by every amplitude
  const auto diff = [](const ExprPtr& E, const ExprPtr& x) {
    return sequant::deriv::diff(E, {x});
  };
  container::map<std::string, Counts> counts;
  // the spellings of the sampled variables other than amplitudes
  container::map<std::string, container::svector<std::string>> sampled;
  // what each model rejected and checked, shown if the pinned counts differ
  std::ostringstream listing;
  for (const auto& m : t::family_models())
    t::with_family(m, [&](const t::Family& f) {
      auto& c = counts[t::name_of(m)];
      const auto checked = t::checked_variables(f);
      const auto& xs = checked.variables;
      V amplitudes;
      for (const auto& x : xs)
        if (t::is_amplitude(x)) amplitudes.push_back(x);
      c.equations = f.equations.size();
      c.variables = xs.size();
      c.amplitudes = amplitudes.size();
      for (const auto& x : xs)
        if (!t::is_amplitude(x))
          sampled[t::name_of(m)].push_back(t::text_of(x));
      {
        INFO(t::name_of(m));
        // open-shell samples an antisymmetric g over more than one space
        if (m.open_shell)
          REQUIRE(std::ranges::any_of(xs, [](const ExprPtr& x) {
            if (!is_g(x) || x->as<Tensor>().symmetry() != Symmetry::Antisymm)
              return false;
            const auto slots = x->as<Tensor>().const_braketaux();
            return !std::ranges::all_of(slots, [&](const Index& i) {
              return i.space() == slots.begin()->space();
            });
          }));
        // closed-shell samples a g of the symmetries N-C-S
        if (m.traced)
          REQUIRE(std::ranges::any_of(xs, [](const ExprPtr& x) {
            if (!is_g(x)) return false;
            const auto& g = x->as<Tensor>();
            return g.symmetry() == Symmetry::Nonsymm &&
                   g.braket_symmetry() == BraKetSymmetry::Conjugate &&
                   g.column_symmetry() == ColumnSymmetry::Symm;
          }));
      }
      for (const auto& [text, rejected] : checked.rejected) {
        c.rejected[text] = rejected.size();
        listing << t::name_of(m) << ": rejected " << rejected.size() << " by \""
                << text << "\":";
        for (const auto& x : rejected) listing << " " << t::text_of(x);
        listing << "\n";
        // the engine accepts every amplitude SeQuant generates
        for (const auto& x : rejected) {
          INFO(t::name_of(m) << ": " << t::text_of(x) << ": " << text);
          REQUIRE_FALSE(t::is_amplitude(x));
        }
      }
      for (const auto& [name, e] : f.equations) {
        INFO(t::name_of(m) << ", " << toUtf8(name));
        for (const auto& x : xs) {
          INFO("by " << t::text_of(x));
          const auto seed = seed_of(m, name, x);
          const auto r = diff(e, x);
          ++c.pairs;
          const auto i1 = t::definition(e, x, r, seed);
          INFO("I1: " << i1.difference);
          REQUIRE(i1.difference.empty());
          REQUIRE((i1.nonzero != 0 || r->is_zero()));
          c.pairs_read += i1.nonzero != 0;
          const auto i6 = t::validity(e, x, r);
          INFO("I6: " << i6);
          REQUIRE(i6.empty());
          if (!t::is_amplitude(x)) continue;
          const auto i3 = t::taylor(e, x, seed);
          ++c.taylor;
          INFO("I3: " << i3.difference);
          REQUIRE(i3.difference.empty());
          c.taylor_read += i3.nonzero != 0;
          const auto& y = t::schwarz_partner(r, x, amplitudes);
          INFO("and " << t::text_of(y));
          const auto i4 = t::schwarz(e, x, y, seed);
          ++c.schwarz;
          INFO("I4: " << i4.difference);
          REQUIRE(i4.difference.empty());
          c.schwarz_read += i4.nonzero != 0;
        }
      }
      INFO(t::name_of(m));
      // every equation by every checked variable, I3 and I4 by every amplitude
      REQUIRE(c.pairs == c.equations * c.variables);
      REQUIRE(c.taylor == c.equations * c.amplitudes);
      REQUIRE(c.schwarz == c.equations * c.amplitudes);
      // the models without CSV reject nothing; the CSV ones leave some out
      REQUIRE(c.rejected.empty() == !m.csv);
    });
  for (const auto& [model, c] : counts)
    listing << model << ": " << c.equations << " equations, " << c.variables
            << " variables, " << c.amplitudes << " amplitudes, " << c.pairs
            << " pairs, " << c.taylor << " taylor, " << c.schwarz
            << " schwarz; nonzero read by " << c.pairs_read << " pairs, "
            << c.taylor_read << " taylor, " << c.schwarz_read << " schwarz\n";
  for (const auto& [model, spellings] : sampled) {
    listing << model << ": sampled";
    for (const auto& x : spellings) listing << " " << x;
    listing << "\n";
  }
  INFO(listing.str());
  // the variables other than amplitudes cannot change unnoticed
  REQUIRE(
      sampled ==
      container::map<std::string, container::svector<std::string>>{
          {"CC{2} CSV spin-orbital",
           {"h¹{i_51;i_52}:A-C-S", "f{i_51;i_52}:A-C-S",
            "g{i_51,a_51<i_52,i_53>;i_52,i_53}:A-C-S",
            "f{a_51<i_51>;i_51}:A-C-S"}},
          {"CC{2} CSV traced",
           {"h¹{i_51;i_52}:N-C-S", "g{i_51,a_51<i_52,i_53>;i_53,i_52}:N-C-S",
            "g{i_51,a_51<i_52>;i_52,a_52<i_51>}:N-C-S",
            "g{i_51,a_51<i_52,i_53>;i_52,i_53}:N-C-S"}},
          {"CC{2} canonical spin-orbital",
           {"h¹{i_51;i_52}:A-C-S", "g{a_51,a_52;i_51,i_52}:A-C-S",
            "f{i_51;i_52}:A-C-S", "f{a_51;i_51}:A-C-S"}},
          {"CC{2} canonical traced",
           {"h¹{a_51;i_51}:N-C-S", "h¹{i_51;i_52}:N-C-S",
            "g{i_51,a_51;a_52,i_52}:N-C-S", "f{a_51;a_52}:N-C-S"}},
          {"CC{2} open-shell",
           {"h¹{i↓_51;i↓_52}:N-C-S", "g{i↑_51,a↑_51;i↑_52,i↑_53}:A-C-S",
            "g{i↑_51,a↑_51;i↑_52,a↑_52}:A-C-S",
            "g{i↓_51,i↓_52;i↓_53,a↓_51}:A-C-S"}}});
  // the family cannot shrink
  const std::string mixes =
      "an index of the variable mixes protoindices that are its slots with "
      "others";
  const std::string not_followed =
      "a leaf's protoindices do not follow its slots as the variable's do";
  const std::string partly_shared =
      "the protoindex lists of the variable's slots share some but not all "
      "labels";
  const std::string not_summed =
      "a protoindex label of a leaf of the variable must be summed over the "
      "space of the variable's label";
  const container::map<std::string, std::size_t> csv_rejected{
      {mixes, 20}, {not_followed, 25}, {partly_shared, 5}, {not_summed, 11}};
  const container::map<std::string, std::size_t> csv_traced_rejected{
      {mixes, 27}, {not_followed, 27}, {partly_shared, 7}, {not_summed, 11}};
  // equations, variables, amplitudes; pairs of I1, I3, I4; those that read a
  // nonzero element; rejections
  REQUIRE(counts ==
          container::map<std::string, Counts>{
              {"CC{2} CSV spin-orbital",
               {12, 14, 10, 168, 120, 120, 71, 120, 38, csv_rejected}},
              {"CC{2} CSV traced",
               {12, 17, 13, 204, 156, 156, 86, 156, 58, csv_traced_rejected}},
              {"CC{2} canonical spin-orbital",
               {12, 14, 10, 168, 120, 120, 68, 120, 37, {}}},
              {"CC{2} canonical traced",
               {12, 14, 10, 168, 120, 120, 73, 120, 44, {}}},
              {"CC{2} open-shell",
               {21, 24, 20, 504, 420, 420, 248, 420, 171, {}}}});
}

TEST_CASE("deriv/identities/own-equations", "[deriv][identities]") {
  // I2: the equations SeQuant generates are the derivatives that give them,
  // element by element, with constant 1
  container::map<std::string, std::size_t> compared;
  for (const auto& m : t::family_models())
    t::with_family(m, [&](const t::Family& f) {
      for (const auto& c : t::own_equations(f)) {
        INFO(t::name_of(m) << ", " << toUtf8(c.name));
        const auto same = t::compare_arrays(c.generated, c.expected, c.free,
                                            {.projectors = true});
        INFO(same.difference);
        CHECK(same.difference.empty());
        CHECK(same.nonzero != 0);
        ++compared[t::name_of(m)];
      }
    });
  // t[k], λ[k], tʼ[k], λʼ[k] for k = 1, 2; open-shell, each spin case
  REQUIRE(compared == container::map<std::string, std::size_t>{
                          {"CC{2} CSV spin-orbital", 8},
                          {"CC{2} CSV traced", 8},
                          {"CC{2} canonical spin-orbital", 8},
                          {"CC{2} canonical traced", 8},
                          {"CC{2} open-shell", 20}});
}

TEST_CASE("deriv/identities/each-identity-can-fail", "[deriv][identities]") {
  // one deliberately wrong result per identity, each caught, on ℒ by t₂ of
  // the canonical spin-orbital model
  const auto diff = [](const ExprPtr& E, const ExprPtr& x) {
    return sequant::deriv::diff(E, {x});
  };
  t::with_family({}, [&](const t::Family& f) {
    const auto L = t::equation_named(f, L"ℒ");
    const auto t2 = T(L"t{a_51,a_52;i_51,i_52}:A-N-S"),
               l2 = T(L"λ{i_51,i_52;a_51,a_52}:A-N-S");
    const auto r = diff(L, t2);
    // I1: twice the derivative; the derivative without its Â
    REQUIRE(t::definition(L, t2, r, 1).difference.empty());
    REQUIRE_FALSE(
        t::definition(L, t2, ex<Constant>(2) * r, 1).difference.empty());
    const auto unprojected =
        t::without_projectors_over(r, t::indices_of_leaf(t2));
    REQUIRE(*unprojected != *r);
    REQUIRE_FALSE(t::definition(L, t2, unprojected, 1).difference.empty());
    // I6: the derivative spelled with other indices than x's
    REQUIRE(t::validity(L, t2, r).empty());
    REQUIRE_FALSE(
        t::validity(L, t2, diff(L, t::renamed_apart(t2, 10))).empty());
    // I3: 1 for the number of operations that keep t₂'s block (4)
    REQUIRE(t::block_order(t2) == 4);
    REQUIRE(t::taylor(L, t2, 1).difference.empty());
    REQUIRE_FALSE(t::taylor(L, t2, 1, 1).difference.empty());
    // I4: one side doubled
    REQUIRE(t::schwarz(L, t2, l2, 1).difference.empty());
    REQUIRE_FALSE(t::schwarz(L, t2, l2, 1, 2).difference.empty());
    // I2: the constant 2
    const auto own = t::own_equations(f);
    const auto t2_it = std::ranges::find_if(
        own, [](const t::OwnEquation& c) { return c.name == L"t[2]"; });
    REQUIRE(t2_it != own.end());
    const auto& t2_eq = *t2_it;
    REQUIRE(t::compare_arrays(t2_eq.generated, t2_eq.expected, t2_eq.free,
                              {.projectors = true})
                .difference.empty());
    REQUIRE_FALSE(t::compare_arrays(t2_eq.generated,
                                    ex<Constant>(2) * t2_eq.expected,
                                    t2_eq.free, {.projectors = true})
                      .difference.empty());
  });
}

TEST_CASE("deriv/identities/exclusions-still-fail-in-core",
          "[deriv][identities]") {
  // the core defect that keeps open-shell CSV out of the family; when core is
  // fixed, this fails and the exclusion goes. Each slot of an amplitude that
  // carries protoindices carries the amplitude's occupied slots: so in the
  // CSV t()[1], not so in its open_shell_CC_spintrace cases
  REQUIRE(t::family_exclusions().size() == 1);
  const auto excluded = t::family_exclusions().front().model;
  REQUIRE((excluded.csv && excluded.open_shell && !excluded.traced));
  const auto pair_follows_slots = [](const Expr& x) {
    const auto slots = x.as<Tensor>().const_braketaux();
    container::set<Index> occupied;
    for (const auto& i : slots)
      if (!i.has_proto_indices()) occupied.insert(i);
    return std::ranges::all_of(slots, [&](const Index& i) {
      return !i.has_proto_indices() ||
             container::set<Index>(i.proto_indices().begin(),
                                   i.proto_indices().end()) == occupied;
    });
  };
  const auto amplitudes_follow = [&](const ExprPtr& e) {
    bool result = true;
    e->visit(
        [&](const ExprPtr& l) {
          result = result && (!t::is_amplitude(l) || pair_follows_slots(*l));
        },
        /*atoms_only=*/true);
    return result;
  };
  sequant::tests::csv::ScopedCsvContext scoped{mbpt::make_legacy_registry(),
                                               mbpt::CSV::Yes};
  // the Complete canonicalization open_shell_CC_spintrace requires
  auto context = get_default_context();
  context.set(CanonicalizeOptions{.method = CanonicalizationMethod::Complete});
  const auto complete = set_scoped_default_context(context);
  const auto r1 = mbpt::CC{2}.t()[1];
  REQUIRE(amplitudes_follow(r1));
  REQUIRE_FALSE(std::ranges::all_of(mbpt::open_shell_CC_spintrace(r1),
                                    amplitudes_follow));
}
