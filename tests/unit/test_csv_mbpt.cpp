//
// Basis grants: OpRegistry grants, OpMaker authoring, projectors and the
// amplitude validator, on operators and derivations SeQuant generates.
//

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include "catch2_sequant.hpp"
#include "csv_test_utils.hpp"

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/index_basis.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/reserved.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/expr.hpp>
#include <SeQuant/core/utility/string.hpp>
#include <SeQuant/domain/mbpt/basis_grants.hpp>
#include <SeQuant/domain/mbpt/context.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/models/cc.hpp>
#include <SeQuant/domain/mbpt/op.hpp>
#include <SeQuant/domain/mbpt/op_registry.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>
#include <SeQuant/domain/mbpt/utils.hpp>
#include <SeQuant/domain/mbpt/vac_av.hpp>

#include <algorithm>
#include <cstddef>
#include <functional>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace sequant::tests::csv_mbpt {

using Instances = std::set<IndexBasis::optional_instance>;

/// what a derived equation's terms hold, as the slots of its amplitude
/// (mbpt::is_amplitude_tensor) and projector tensors bound them
struct Census {
  Instances instances;                 ///< over every slot
  std::size_t max_per_term = 0;        ///< distinct non-null instances
  std::size_t metrics = 0;             ///< overlap tensors
  std::size_t metrics_off = 0;         ///< ... with an end on no such slot
  std::size_t integral_slots = 0;      ///< other tensors' slots with instance
  std::size_t integral_slots_off = 0;  ///< ... on no such slot
};

inline Census census(ExprPtr const& eq) {
  auto const& reg = *mbpt::get_default_mbpt_context().op_registry();
  auto is_projector = [](AbstractTensor const& t) {
    return t._label() == reserved::antisymm_label() ||
           t._label() == reserved::symm_label();
  };
  Census result;
  auto add_term = [&](ExprPtr const& term) {
    container::set<Index> bound;
    Instances term_instances;
    term->visit(
        [&](ExprPtr const& x) {
          if (!x->is<AbstractTensor>()) return;
          auto const& t = x->as<AbstractTensor>();
          for (Index const& idx : t._slots()) {
            result.instances.insert(idx.basis().basis_instance());
            if (idx.basis().has_basis_instance())
              term_instances.insert(idx.basis().basis_instance());
            if (mbpt::is_amplitude_tensor(t, reg) || is_projector(t))
              bound.insert(idx);
          }
        },
        /* atoms_only = */ true);
    result.max_per_term = std::max(result.max_per_term, term_instances.size());
    term->visit(
        [&](ExprPtr const& x) {
          if (!x->is<AbstractTensor>()) return;
          auto const& t = x->as<AbstractTensor>();
          if (mbpt::is_amplitude_tensor(t, reg) || is_projector(t)) return;
          if (t._label() == reserved::overlap_label()) {
            ++result.metrics;
            if (std::ranges::any_of(t._slots(), [&](Index const& idx) {
                  return !bound.contains(idx);
                }))
              ++result.metrics_off;
            return;
          }
          for (Index const& idx : t._slots())
            if (idx.basis().has_basis_instance()) {
              ++result.integral_slots;
              if (!bound.contains(idx)) ++result.integral_slots_off;
            }
        },
        /* atoms_only = */ true);
  };
  if (eq->is<Sum>())
    for (auto const& term : *eq) add_term(term);
  else
    add_term(eq);
  return result;
}

/// checks that every term of @p eq has its unoccupied externals (projector
/// slots) at @p expected and its occupied ones at none
inline void check_externals(ExprPtr const& eq,
                            IndexBasis::optional_instance expected) {
  auto const& isr = get_default_context().index_space_registry();
  REQUIRE(eq->is<Sum>());
  for (auto const& term : *eq) {
    Instances virt, occ;
    for (Index const& idx : csv::projector_slots(term))
      (isr->is_pure_occupied(idx.space()) ? occ : virt)
          .insert(idx.basis().basis_instance());
    REQUIRE(!virt.empty());
    CHECK(virt == Instances{expected});
    CHECK(occ == Instances{std::nullopt});
  }
}

/// @p tensor with its slot @p idx given @p instance
inline ExprPtr with_slot_instance(ExprPtr const& tensor, Index const& idx,
                                  IndexBasis::optional_instance instance) {
  return transform_expr(tensor, {{idx, idx.replace_basis_instance(instance)}});
}

/// the first tensor labelled @p label in @p eqs (with @p with_instance, the
/// first with a slot carrying a basis instance), as an expression
inline ExprPtr first_tensor(std::vector<ExprPtr> const& eqs,
                            std::wstring_view label,
                            bool with_instance = false) {
  ExprPtr result;
  for (auto const& eq : eqs)
    if (eq && !result)
      eq->visit(
          [&](ExprPtr const& x) {
            if (result || !x->is<AbstractTensor>()) return;
            auto const& t = x->as<AbstractTensor>();
            if (t._label() == label &&
                (!with_instance ||
                 std::ranges::any_of(t._slots(), [](Index const& idx) {
                   return idx.basis().has_basis_instance();
                 })))
              result = x->clone();
          },
          /* atoms_only = */ true);
  return result;
}

/// the first pure-unoccupied bra/ket slot of @p tensor
inline Index unoccupied_slot(ExprPtr const& tensor) {
  auto const& isr = get_default_context().index_space_registry();
  for (Index const& idx : tensor->as<AbstractTensor>()._braket())
    if (isr->is_pure_unoccupied(idx.space())) return idx;
  throw Exception("unoccupied_slot: none");
}

}  // namespace sequant::tests::csv_mbpt

TEST_CASE("basis-grants-registry", "[mbpt][csv]") {
  using namespace sequant;
  using sequant::tests::csv::csv_cc_context;
  using sequant::tests::csv::message_contains;

  auto ctx = set_scoped_default_context(csv_cc_context());
  auto const& isr = get_default_context().index_space_registry();
  const auto a = isr->retrieve(L"a");
  const auto i = isr->retrieve(L"i");

  auto reg = mbpt::make_minimal_registry();
  auto granted = std::make_shared<mbpt::OpRegistry>(reg->clone());
  CHECK(*granted == *reg);
  CHECK_FALSE(granted->has_basis_grants(L"t"));
  granted->grant_basis(L"t", a, 0);
  CHECK(granted->has_basis_grants(L"t"));
  CHECK(granted->basis_grant(L"t", a) == IndexBasis::optional_instance{0});
  CHECK_FALSE(granted->basis_grant(L"t", i).has_value());
  CHECK_FALSE(reg->basis_grant(L"t", a).has_value());
  CHECK_FALSE(*granted == *reg);
  granted->grant_basis(L"λ", a, -3);
  CHECK(granted->basis_grant(L"λ", a) == IndexBasis::optional_instance{-3});

  CHECK_THROWS_MATCHES(granted->grant_basis(L"nope", a, 1), Exception,
                       message_contains("does not exist in registry"));
  CHECK_THROWS_MATCHES(granted->grant_basis(L"g", a, 1), Exception,
                       message_contains("is a general operator"));
  CHECK_FALSE(granted->has_basis_grants(L"nope"));
  CHECK_FALSE(granted->basis_grant(L"nope", a).has_value());

  granted->remove(L"t").add(L"t", mbpt::OpClass::Ex);
  CHECK_FALSE(granted->has_basis_grants(L"t"));
  granted->purge();
  CHECK_FALSE(granted->has_basis_grants(L"λ"));
}

TEST_CASE("basis-grants-authoring", "[mbpt][csv]") {
  using namespace sequant;
  using namespace sequant::mbpt;
  using namespace sequant::tests::csv;
  namespace t = op::tensor;

  auto ctx = set_scoped_default_context(csv_cc_context());
  auto const& isr = get_default_context().index_space_registry();
  const auto a = isr->retrieve(L"a");
  const auto i = isr->retrieve(L"i");

  // every bra/ket slot of the tensors labelled `label` in `expr`
  auto slots_of = [](ExprPtr const& expr, std::wstring_view label) {
    std::vector<Index> result;
    for (auto const* tn : tensors_labelled(expr, label))
      for (Index const& idx : tn->_braket()) result.push_back(idx);
    return result;
  };
  auto all_null = [](ExprPtr const& expr) {
    const auto histogram = instance_histogram(expr);
    return histogram.size() == 1 && histogram.contains(std::nullopt);
  };

  SECTION("stamp = grant") {
    const auto csv = GENERATE(CSV::Yes, CSV::No);
    const bool yes = csv == CSV::Yes;
    using Grants =
        std::vector<std::pair<std::wstring, IndexBasis::instance_type>>;
    // the K1-K5 products of test_basis_wick.cpp each grant set authors
    for (auto const& [grants, products] :
         std::vector<std::pair<Grants, std::vector<int>>>{
             {{{L"t", 1}}, {1, 2, 5}},
             {{{L"t", 1}, {L"λ", 1}}, {3}},
             {{{L"t", 1}, {L"λ", 2}}, {4}}}) {
      INFO("CSV " << yes << ", " << grants.size() << " grants");
      const auto reg = granted_registry(grants);
      auto author = [csv](std::shared_ptr<OpRegistry> r,
                          std::wstring_view label) {
        auto mbpt_ctx = set_scoped_default_mbpt_context(
            mbpt::Context({.csv = csv, .op_registry_ptr = std::move(r)}));
        return label == L"t"   ? t::t(2)
               : label == L"λ" ? t::λ(2)
               : label == L"h" ? t::h(2)
               : label == L"f" ? t::h(1)
                               : t::P(nₚ(1), nₕ(1), {}, L"t");
      };
      // granted: authored under the grants; stamped: authored without, then
      // given each leg's grant
      auto granted = [&](std::wstring_view label) {
        return author(reg, label);
      };
      auto stamped = [&](std::wstring_view label) {
        const auto op = author(make_minimal_registry(), label);
        const std::wstring owner = label == L"λ"   ? L"λ"
                                   : label == L"P" ? L"t"
                                                   : std::wstring(label);
        return with_leg_instance(op, reg->basis_grant(owner, a));
      };
      for (auto label : {L"t", L"λ", L"P"}) {
        INFO(toUtf8(std::wstring(label)));
        Index::reset_tmp_index();
        const auto g = granted(label);
        Index::reset_tmp_index();
        const auto s = stamped(label);
        CHECK(*g == *s);
      }

      // factors authored left to right after one reset, so that both
      // products of a case mint the same indices
      using Maker = std::function<ExprPtr(std::wstring_view)>;
      auto product = [](int k, Maker const& make) {
        Index::reset_tmp_index();
        std::vector<std::wstring_view> labels =
            k == 1   ? std::vector<std::wstring_view>{L"h", L"t"}
            : k == 2 ? std::vector<std::wstring_view>{L"P", L"f", L"t"}
                     : std::vector<std::wstring_view>{L"λ", L"t"};
        ExprPtr result;
        for (auto label : labels) {
          auto factor = make(label);
          result = result ? result * factor : factor;
        }
        return result;
      };
      auto vac_av = [](ExprPtr const& p) {
        Index::reset_tmp_index();
        return t::vac_av(p);
      };
      for (int k : products) {
        INFO("K" << k);
        const auto g = product(k, granted);
        const auto s = product(k, stamped);
        if (k == 5 && yes) {
          // the guard through grants: λ ungranted next to a granted t
          CHECK_THROWS_MATCHES(
              vac_av(g), Exception,
              message_contains("carries proto indices and no basis instance"));
          CHECK_THROWS_AS(vac_av(s), Exception);
        } else {
          CHECK(serialize(vac_av(g)) == serialize(vac_av(s)));
        }
      }
    }
  }

  SECTION("integrals carry no instance") {
    for (auto const& grants : std::vector<
             std::vector<std::pair<std::wstring, IndexBasis::instance_type>>>{
             {}, {{L"t", 1}}, {{L"t", 1}, {L"λ", 2}}}) {
      ScopedCsvContext scoped{granted_registry(grants)};
      CHECK(all_null(t::H(2)));
    }
    CHECK_THROWS_AS(make_minimal_registry()->grant_basis(L"g", a, 1),
                    Exception);
  }

  SECTION("a grant stamps the legs of its leg space only") {
    for (bool occupied : {false, true}) {
      INFO("occupied " << occupied);
      auto reg = make_minimal_registry();
      reg->grant_basis(L"t", occupied ? i : a, 5);
      ScopedCsvContext scoped{reg};
      const auto slots = slots_of(t::t(2), L"t");
      REQUIRE(slots.size() == 4);
      for (Index const& idx : slots) {
        const bool occ = isr->is_pure_occupied(idx.space());
        CHECK(idx.basis().basis_instance() ==
              (occ == occupied ? IndexBasis::optional_instance{5}
                               : std::nullopt));
      }
    }
  }

  SECTION("a grant stamps only its decorated label") {
    const std::vector<std::wstring> amplitudes{L"t", L"λ", L"t¹", L"λ¹"};
    auto make = [](std::wstring_view label) {
      return label == L"t"    ? t::t(2)
             : label == L"λ"  ? t::λ(2)
             : label == L"t¹" ? t::tʼ(2)
                              : t::λʼ(2);
    };
    for (auto const& granted : amplitudes) {
      ScopedCsvContext scoped{granted_registry({{granted, 7}})};
      for (auto const& label : amplitudes) {
        INFO(toUtf8(granted) << " granted, " << toUtf8(label) << " built");
        std::size_t n = 0;
        for (Index const& idx : slots_of(make(label), label))
          n += idx.basis().basis_instance() == 7;
        CHECK(n == (label == granted ? 2 : 0));
      }
    }
  }

  SECTION("dependent indices are stamped as minted") {
    auto reg = make_minimal_registry();
    reg->grant_basis(L"t", a, 7);
    ScopedCsvContext scoped{reg};
    const auto amplitude = t::t(2);
    const auto ts = tensors_labelled(amplitude, L"t");
    REQUIRE(ts.size() == 1);
    for (Index const& idx : ts.front()->_bra()) {
      CHECK(idx.basis().basis_instance() == 7);
      REQUIRE(idx.proto_indices().size() == 2);
      for (Index const& p : idx.proto_indices())
        CHECK_FALSE(p.basis().has_basis_instance());
    }
    for (Index const& idx : ts.front()->_ket())
      CHECK_FALSE(idx.basis().has_basis_instance());

    // a granted independent index is the proto index of its dependents
    using Maker = OpMaker<Statistics::FermiDirac>;
    const container::svector<IndexSpace> uoccs{a};
    const auto op = Maker::make(
        uoccs, uoccs,
        [](auto const& creidxs, auto const& annidxs, Symmetry opsymm) {
          return ex<Tensor>(L"X", bra(creidxs), ket(annidxs), opsymm);
        },
        Maker::UseDepIdx::Bra, Normalization::Default,
        [](IndexSpace const&) { return IndexBasis::optional_instance{7}; });
    const auto xs = tensors_labelled(op, L"X");
    REQUIRE(xs.size() == 1);
    auto ket = xs.front()->_ket();
    const Index& ket0 = ket[0];
    auto bra = xs.front()->_bra();
    const Index& bra0 = bra[0];
    CHECK(ket0.basis().basis_instance() == 7);
    CHECK(bra0.proto_indices() == Index::index_vector{ket0});
  }

  SECTION("a grant with CSV::No stamps domainless legs") {
    auto reg = make_minimal_registry();
    reg->grant_basis(L"t", a, 7);
    ScopedCsvContext scoped{reg, CSV::No};
    for (auto const& [op, label] :
         {std::pair{t::t(2), std::wstring(L"t")},
          std::pair{t::P(nₚ(2), nₕ(2), {}, L"t"),
                    std::wstring(reserved::antisymm_label())}}) {
      std::size_t n = 0;
      for (Index const& idx : slots_of(op, label))
        if (isr->is_pure_unoccupied(idx.space())) {
          ++n;
          CHECK_FALSE(idx.has_proto_indices());
          CHECK(idx.full_label().ends_with(L"<;7>"));
        }
      CHECK(n == 2);
    }
    CHECK_NOTHROW(t::H(2));
  }

  SECTION("theta and EOM amplitudes are authored without an instance") {
    ScopedCsvContext scoped{granted_registry({{L"t", 1}, {L"λ", 1}})};
    for (auto const& expr :
         {t::θ(1), t::θ(2), t::r(nₚ(1), nₕ(1)), t::l(nₚ(1), nₕ(1))})
      CHECK(all_null(expr));
  }
}

TEST_CASE("basis-grants-validator", "[mbpt][csv][valgrind_skip]") {
  using namespace sequant;
  using namespace sequant::mbpt;
  using namespace sequant::tests::csv;
  using namespace sequant::tests::csv_mbpt;

  auto ctx = set_scoped_default_context(csv_cc_context());
  const auto ungranted = granted_registry({});
  const auto t_granted = granted_registry({{L"t", 1}});
  auto derive_r1 = [](std::shared_ptr<OpRegistry> reg, CSV csv) {
    ScopedCsvContext scoped{std::move(reg), csv};
    Index::reset_tmp_index();
    return CC{2}.t(1, 1);
  };
  auto validate = [](std::shared_ptr<OpRegistry> reg, ExprPtr const& expr) {
    auto mbpt_ctx = set_scoped_default_mbpt_context(
        mbpt::Context({.csv = CSV::Yes, .op_registry_ptr = std::move(reg)}));
    assert_amplitudes_carry_granted_basis(expr);
  };

  SECTION("slots with proto indices carry exactly their grant") {
    const auto t0 = first_tensor(derive_r1(ungranted, CSV::Yes), L"t");
    const auto t1 = first_tensor(derive_r1(t_granted, CSV::Yes), L"t");
    REQUIRE((t0 && t1));
    const Index a0 = unoccupied_slot(t0), a1 = unoccupied_slot(t1);
    REQUIRE(a0.has_proto_indices());
    REQUIRE(a1.basis().basis_instance() == 1);
    CHECK_NOTHROW(validate(ungranted, t0));
    CHECK_NOTHROW(validate(t_granted, t1));
    CHECK_THROWS_MATCHES(validate(ungranted, with_slot_instance(t0, a0, 1)),
                         Exception, message_contains("the grant of t"));
    // a copy path that lost the instance
    CHECK_THROWS_AS(validate(t_granted, with_slot_instance(t1, a1, {})),
                    Exception);
  }

  SECTION("domainless slots carry their grant when one exists") {
    const auto t0 = first_tensor(derive_r1(ungranted, CSV::No), L"t");
    const auto t1 = first_tensor(derive_r1(t_granted, CSV::No), L"t");
    REQUIRE((t0 && t1));
    const Index a0 = unoccupied_slot(t0), a1 = unoccupied_slot(t1);
    REQUIRE_FALSE(a0.has_proto_indices());
    // absorbed from a partner by Wick
    CHECK_NOTHROW(validate(ungranted, with_slot_instance(t0, a0, 7)));
    CHECK_NOTHROW(validate(t_granted, t1));
    CHECK_THROWS_AS(validate(t_granted, with_slot_instance(t1, a1, 7)),
                    Exception);
  }

  SECTION("integrals, reserved labels, perturbed and adjoint amplitudes") {
    const auto all =
        granted_registry({{L"t", 1}, {L"λ", 1}, {L"t¹", 10}, {L"λ¹", 10}});
    std::vector<ExprPtr> tp, l;
    {
      ScopedCsvContext scoped{all};
      Index::reset_tmp_index();
      tp = CC{2}.tʼ(1, 1);
      Index::reset_tmp_index();
      l = CC{2}.λ();
    }
    const auto t1 = first_tensor(tp, L"t¹");
    // integral and reserved-label tensors with the instances Wick gave them
    const auto g = first_tensor(tp, L"g", true);
    const auto s = first_tensor(tp, reserved::overlap_label(), true);
    const auto proj = first_tensor(tp, reserved::antisymm_label(), true);
    const auto ladj = first_tensor(l, std::wstring(L"λ") + adjoint_label);
    REQUIRE((t1 && g && s && proj && ladj));
    CHECK_NOTHROW(validate(all, t1));
    CHECK_THROWS_AS(
        validate(all, with_slot_instance(t1, unoccupied_slot(t1), 1)),
        Exception);
    for (auto const& x : {g, s, proj}) CHECK_NOTHROW(validate(all, x));
    // the adjoint is checked against the grant of its amplitude
    const Index la = unoccupied_slot(ladj);
    REQUIRE(la.basis().basis_instance() == 1);
    CHECK_NOTHROW(validate(all, ladj));
    CHECK_THROWS_AS(validate(all, with_slot_instance(ladj, la, 2)), Exception);
    CHECK_THROWS_AS(validate(ungranted, ladj), Exception);
  }

  SECTION("a partial grant set throws") {
    // R1 externals mix the granted and the ungranted family's basis
    auto partial = granted_registry({{L"t", 1}});
    partial->add(L"t¹", OpClass::Ex);
    const auto full = granted_registry({{L"t", 1}, {L"t¹", 10}});
    for (auto const& [reg, ok] :
         {std::pair{partial, false}, std::pair{full, true}}) {
      ScopedCsvContext scoped{reg, CSV::No};
      Index::reset_tmp_index();
      const auto eqs = CC{2}.tʼ(1, 1);
      for (std::size_t p : {1u, 2u}) {
        INFO("ok " << ok << ", R" << p);
        if (ok)
          CHECK_NOTHROW(assert_amplitudes_carry_granted_basis(eqs.at(p)));
        else
          CHECK_THROWS_MATCHES(
              assert_amplitudes_carry_granted_basis(eqs.at(p)), Exception,
              message_contains("amplitude t¹ has no basis grant"));
      }
    }
  }

  SECTION("CSV::No: ungranted domainless t legs absorb λ's instance") {
    const auto reg = granted_registry({{L"λ", 7}});
    ScopedCsvContext scoped{reg, CSV::No};
    Index::reset_tmp_index();
    const auto r1 = CC{2}.λ().at(1);
    std::size_t absorbed = 0;
    for (auto const* tn : tensors_labelled(r1, L"t"))
      for (Index const& idx : tn->_braket())
        absorbed += idx.basis().basis_instance() == 7;
    CHECK(absorbed == 38);
    // every slot rule passes; the partial-grant rule names t
    CHECK_THROWS_MATCHES(assert_amplitudes_carry_granted_basis(r1), Exception,
                         message_contains("amplitude t has no basis grant"));
  }
}

TEST_CASE("basis-grants-cc-t", "[mbpt][csv][valgrind_skip]") {
  using namespace sequant;
  using namespace sequant::mbpt;
  using namespace sequant::tests::csv;
  using namespace sequant::tests::csv_mbpt;

  auto ctx = set_scoped_default_context(csv_cc_context());
  auto const& isr = get_default_context().index_space_registry();

  SECTION("a grant changes no term or metric count") {
    std::vector<std::size_t> terms, metrics;
    for (auto const& reg :
         {granted_registry({}), granted_registry({{L"t", 1}})}) {
      const bool granted = reg->has_basis_grants(L"t");
      ScopedCsvContext scoped{reg};
      Index::reset_tmp_index();
      std::vector<std::size_t> n, m;
      for (auto const& eq : CC{2}.t(2, 0)) {
        n.push_back(term_count(eq));
        const auto c = census(eq);
        m.push_back(c.metrics);
        CHECK(c.metrics_off == 0);
        if (granted)
          for (auto const* g : tensors_labelled(eq, L"g"))
            for (Index const& idx : g->_braket())
              if (isr->is_pure_unoccupied(idx.space()))
                CHECK(idx.basis().basis_instance() == 1);
        CHECK_NOTHROW(assert_amplitudes_carry_granted_basis(eq));
      }
      CHECK(n == std::vector<std::size_t>{3, 14, 31});
      CHECK(m == std::vector<std::size_t>{0, 9, 41});
    }
  }

  SECTION("Bernoulli UCC projector legs carry the amplitude's grant") {
    ScopedCsvContext scoped{granted_registry({{L"t", 1}}), CSV::No};
    Index::reset_tmp_index();
    const CC ucc(2, {.ansatz = CC::Ansatz::U,
                     .hbar_comm_rank = 2,
                     .hbar_expansion = CC::HbarExpansion::Bernoulli});
    const auto eqs = ucc.t(2, 1);
    for (std::size_t p : {1u, 2u}) {
      INFO("R" << p);
      check_externals(eqs.at(p), 1);
    }
  }

  SECTION("projector legs carry the amplitude's grant") {
    for (auto const& inst :
         {IndexBasis::optional_instance{}, IndexBasis::optional_instance{10}}) {
      ScopedCsvContext scoped{inst ? granted_registry({{L"t", *inst}})
                                   : granted_registry({})};
      Index::reset_tmp_index();
      const auto eqs = CC{2}.t(2, 1);
      for (std::size_t p : {1u, 2u}) {
        INFO("R" << p);
        check_externals(eqs.at(p), inst);
      }
    }
  }

  SECTION("CSV without grants carries no instance") {
    ScopedCsvContext scoped;
    const auto eqs = derive_t(make_minimal_registry());
    const auto eqs_empty = derive_t(granted_registry({}));
    REQUIRE(eqs.size() == 3);
    for (std::size_t k = 0; k != eqs.size(); ++k) {
      CHECK(instance_histogram(eqs[k]).size() == 1);
      CHECK(instance_histogram(eqs[k]).contains(std::nullopt));
      CHECK(standing_metrics(eqs[k]) == 0);
      CHECK(*eqs[k] == *eqs_empty[k]);
    }
    // the t equations hold no ungranted family next to t
    CHECK_NOTHROW(derive_t(granted_registry({{L"t", 1}})));
  }

  SECTION("occupied legs granted") {
    const auto csv = GENERATE(CSV::Yes, CSV::No);
    ScopedCsvContext scoped{make_minimal_registry(), csv};
    auto reg = make_minimal_registry();
    reg->grant_basis(L"t", isr->retrieve(L"i"), 5);
    std::vector<ExprPtr> eqs;
    REQUIRE_NOTHROW(eqs = derive_t(reg));
    const auto plain = derive_t(granted_registry({}));
    for (std::size_t k = 0; k != eqs.size(); ++k) {
      INFO("CSV " << (csv == CSV::Yes) << ", equation " << k);
      CHECK(term_count(eqs[k]) == term_count(plain[k]));
      std::size_t n = 0;
      for (auto const* tn : tensors_labelled(eqs[k], L"t"))
        for (Index const& idx : tn->_braket()) {
          ++n;
          CHECK(idx.basis().basis_instance() ==
                (isr->is_pure_occupied(idx.space())
                     ? IndexBasis::optional_instance{5}
                     : std::nullopt));
        }
      CHECK(n > 0);
    }
  }
}

TEST_CASE("basis-grants-lambda", "[mbpt][csv][valgrind_skip]") {
  using namespace sequant;
  using namespace sequant::mbpt;
  using namespace sequant::tests::csv;
  using namespace sequant::tests::csv_mbpt;

  // Λ-CCSD with t and λ in two families
  auto ctx = set_scoped_default_context(csv_cc_context());
  const auto reg = granted_registry({{L"t", 1}, {L"λ", 2}});
  ScopedCsvContext scoped{reg};
  Index::reset_tmp_index();
  const auto raw = CC{2}.λ();
  std::map<IndexBasis::optional_instance, std::size_t> g_r2;
  for (auto const* g : tensors_labelled(raw.at(2), L"g"))
    for (Index const& idx : g->_braket()) ++g_r2[idx.basis().basis_instance()];
  CHECK(g_r2.at(1) == 19);
  CHECK(g_r2.at(2) == 35);

  const auto st = derive_λ(reg);
  std::vector<std::size_t> metrics;
  for (std::size_t p : {1u, 2u}) {
    const auto c = census(st.at(p));
    CHECK(c.instances == Instances{std::nullopt, 1, 2});
    CHECK(c.metrics_off == 0);
    metrics.push_back(c.metrics);
    CHECK_NOTHROW(assert_amplitudes_carry_granted_basis(st.at(p)));
  }
  CHECK(metrics == std::vector<std::size_t>{183, 88});
}

TEST_CASE("basis-grants-linear-response", "[mbpt][csv][valgrind_skip]") {
  using namespace sequant;
  using namespace sequant::mbpt;
  using namespace sequant::tests::csv;
  using namespace sequant::tests::csv_mbpt;

  auto ctx = set_scoped_default_context(csv_cc_context());
  const auto pno = test_csv::pno_pert;

  SECTION("perturbed amplitudes in their own family") {
    const auto reg =
        granted_registry({{L"t", 1}, {L"λ", 1}, {L"t¹", pno}, {L"λ¹", pno}});
    ScopedCsvContext scoped{reg};
    const CC cc{2};
    Index::reset_tmp_index();
    const auto tp = cc.tʼ(1, 1);
    Index::reset_tmp_index();
    const auto lp = cc.λʼ(1, 1);

    std::vector<std::size_t> int_slots, metrics;
    for (auto const* eqs : {&tp, &lp})
      for (std::size_t p : {1u, 2u}) {
        INFO("R" << p);
        ExprPtr const& eq = eqs->at(p);
        check_externals(eq, pno);
        const auto c = census(eq);
        CHECK(c.integral_slots_off == 0);
        CHECK(c.metrics_off == 0);
        int_slots.push_back(c.integral_slots);
        metrics.push_back(c.metrics);
        CHECK_NOTHROW(assert_amplitudes_carry_granted_basis(eq));
        // through the closed-shell spin trace
        CHECK(census(closed_shell_CC_spintrace_v2(eq)).instances ==
              Instances{std::nullopt, 1, pno});
      }
    CHECK(int_slots == std::vector<std::size_t>{42, 104, 191, 112});
    CHECK(metrics == std::vector<std::size_t>{19, 89, 138, 97});

    // the adjoint perturbed amplitude of the pseudoresponse keeps its grant
    Index::reset_tmp_index();
    const auto h1_bar =
        lst(op::Hʼ(1), op::T(2), 2, {.use_connected_form = true});
    const OpConnections<std::wstring> connect = {{L"h¹", L"t"}};
    const auto pseudo = closed_shell_CC_spintrace_v2(
        op::vac_av(adjoint(op::Tʼ(2)) * h1_bar, {.connect = connect}));
    std::size_t n = 0;
    auto const& isr = get_default_context().index_space_registry();
    for (auto const* tn :
         tensors_labelled(pseudo, std::wstring(L"t¹") + adjoint_label))
      for (Index const& idx : tn->_braket())
        if (isr->is_pure_unoccupied(idx.space())) {
          ++n;
          CHECK(idx.basis().basis_instance() == pno);
        }
    CHECK(n > 0);
  }

  SECTION("four amplitude families") {
    const auto reg =
        granted_registry({{L"t", 1}, {L"t¹", 10}, {L"λ", 11}, {L"λ¹", 12}});
    ScopedCsvContext scoped{reg};
    Index::reset_tmp_index();
    const auto lp = CC{2}.λʼ(1, 1);
    const Instances all{std::nullopt, 1, 10, 11, 12};
    std::size_t max_per_term = 0;
    std::vector<std::size_t> metrics;
    for (bool traced : {false, true})
      for (std::size_t p : {1u, 2u}) {
        INFO("traced " << traced << ", R" << p);
        const auto eq =
            traced ? closed_shell_CC_spintrace_v2(lp.at(p)) : lp.at(p);
        const auto c = census(eq);
        CHECK(c.instances == all);
        CHECK(c.metrics_off == 0);
        metrics.push_back(c.metrics);
        if (!traced) max_per_term = std::max(max_per_term, c.max_per_term);
        CHECK_NOTHROW(assert_amplitudes_carry_granted_basis(eq));
      }
    CHECK(metrics == std::vector<std::size_t>{141, 97, 486, 197});
    CHECK(max_per_term == 4);
  }
}
