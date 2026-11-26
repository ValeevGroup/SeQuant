#include <SeQuant/core/context.hpp>
#include <SeQuant/core/deriv/deriv.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/reserved.hpp>
#include <SeQuant/core/utility/string.hpp>
#include <SeQuant/domain/mbpt/context.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/models/hf.hpp>
#include <SeQuant/domain/mbpt/op.hpp>
#include <SeQuant/domain/mbpt/op_registry.hpp>
#include <SeQuant/domain/mbpt/series.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "contexts.hpp"
#include "evaluator.hpp"
#include "judge.hpp"

#include <string>
#include <string_view>

using namespace sequant;

namespace {

/// @p text with O and V spelled as the labels @p occupied and @p unoccupied
ExprPtr spelled(std::wstring text, std::wstring_view occupied,
                std::wstring_view unoccupied) {
  for (auto& c : text)
    if (c == L'O' || c == L'V') {
      const auto& label = c == L'O' ? occupied : unoccupied;
      REQUIRE(label.size() == 1);
      c = label.front();
    }
  return deserialize(text);
}

using deriv::test::require_same;

/// the checks over complex spaces, the rotation running from the space
/// labelled @p o to that labelled @p v
void check_complex(std::wstring_view o, std::wstring_view v) {
  using deriv::at_zero;
  using deriv::diff;
  const auto T = [&](std::wstring_view s) {
    return spelled(std::wstring(s), o, v);
  };
  const auto index = [](std::wstring_view label, int ordinal) {
    return Index(std::wstring(label) + L"_" + std::to_wstring(ordinal));
  };
  const auto o1 = index(o, 1), o2 = index(o, 2), v1 = index(v, 1),
             v2 = index(v, 2);
  const auto x = T(L"κ̂{V_1;O_1}:A-N-S"), y = T(L"κ̂{V_2;O_2}:A-N-S");

  REQUIRE(mbpt::hf::rotated_energy(0)->is_zero());
  require_same(
      mbpt::hf::rotated_energy(1),
      T(L"f{O_1;V_1}:A-C-S") * x + adjoint(x) * T(L"f{V_1;O_1}:A-C-S"));
  const auto E = mbpt::hf::rotated_energy(2);
  require_same(at_zero(diff(E, {x}), {x}), T(L"f{O_1;V_1}:A-C-S"), {o1, v1});
  require_same(diff(E, {x, adjoint(y)}),
               T(L"f{V_2;V_1}:A-C-S δ{O_1;O_2}:N-S-S - "
                 L"f{O_1;O_2}:A-C-S δ{V_2;V_1}:N-S-S + "
                 L"g{O_1,V_2;V_1,O_2}:A-C-S"),
               {o1, v1, o2, v2});
  require_same(diff(E, {x, y}), T(L"g{O_1,O_2;V_1,V_2}:A-C-S"),
               {o1, v1, o2, v2});
  require_same(diff(E, {adjoint(x), adjoint(y)}),
               T(L"g{V_1,V_2;O_1,O_2}:A-C-S"), {o1, v1, o2, v2});

  const auto closed_shell = mbpt::closed_shell_spintrace(E, {});
  const auto x_cs = T(L"κ̂{V_1;O_1}:N-N-S");
  require_same(at_zero(diff(closed_shell, {x_cs}), {x_cs}),
               T(L"2 f{O_1;V_1}:N-C-S"), {o1, v1});
}

}  // namespace

TEST_CASE("deriv/hf/orbital-rotation", "[deriv][hf]") {
  using Catch::Matchers::StartsWith;
  using deriv::at_zero;
  using deriv::diff;
  const auto T = [](std::wstring_view s) { return deserialize(s); };
  const Index m1(L"m_1"), m2(L"m_2"), e1(L"e_1"), e2(L"e_2");

  SECTION("registered") {
    for (const auto& registry :
         {mbpt::make_legacy_registry(), mbpt::make_minimal_registry()})
      REQUIRE(registry->to_class(L"κ̂") == mbpt::OpClass::Ex);
    REQUIRE(adjoint(T(L"κ̂{e_1;m_1}:A-N-S"))->as<Tensor>().label() == L"κ̂⁺");
  }

  SECTION("label") {
    REQUIRE(L"κ̂" != reserved::cumulant_label());
    REQUIRE(reserved::is_nonreserved(L"κ̂"));
    for (const auto* text : {L"κ̂{a_1;i_1}:A-N-S", L"κ̂¹{a_1;i_1}:A-N-S",
                             L"κ̂⁺{i_1;a_1}:A-N-S", L"κ̂¹⁺{i_1;a_1}:A-N-S"})
      REQUIRE(serialize(T(text)) == text);
    // the rotation, its adjoint and their first-order series coefficients
    const auto latex = [](const ExprPtr& e, std::wstring_view label) {
      std::wstring result;
      e->visit(
          [&](const ExprPtr& leaf) {
            if (leaf->is<Tensor>() && leaf->as<Tensor>().label() == label)
              result = leaf->to_latex();
          },
          /*atoms_only=*/true);
      REQUIRE_THAT(toUtf8(result),
                   Catch::Matchers::ContainsSubstring("\\hat{\\kappa}"));
      REQUIRE(result.find(L'κ') == std::wstring::npos);
      return result;
    };
    const auto E = mbpt::hf::rotated_energy(1);
    REQUIRE(latex(E, L"κ̂") == L"{\\hat{\\kappa}^{{m_1}}_{{e_1}}}");
    latex(E, L"κ̂⁺");
    const auto S = mbpt::series(E, L"κ̂", ex<Variable>(L"ε"), 1);
    latex(S, L"κ̂¹");
    latex(S, L"κ̂¹⁺");
  }

  SECTION("complex spaces") { check_complex(L"m", L"e"); }

  SECTION("minimal context") {
    auto context = get_default_context();
    context.set(mbpt::make_min_sr_spaces());
    const auto scope = set_scoped_default_context(context);
    const auto mbpt_scope = mbpt::set_scoped_default_mbpt_context(
        mbpt::Context({.op_registry_ptr = mbpt::make_minimal_registry()}));
    check_complex(L"i", L"a");
  }

  SECTION("real spaces") {
    const auto scope = deriv::test::real_spaces();
    const deriv::test::Data real{.real = {L"f", L"g", L"κ̂"}};
    const auto x = T(L"κ̂{e_1;m_1}:A-N-S"), y = T(L"κ̂{e_2;m_2}:A-N-S");
    const auto E = mbpt::hf::rotated_energy(2);
    require_same(at_zero(diff(E, {x}), {x}), T(L"2 f{e_1;m_1}:A-S-S"), {m1, e1},
                 real);
    const auto hessian =
        T(L"2 δ{m_1;m_2}:N-S-S f{e_1;e_2}:A-S-S - "
          L"2 δ{e_1;e_2}:N-S-S f{m_1;m_2}:A-S-S + "
          L"2 g{m_1,m_2;e_1,e_2}:A-S-S + 2 g{e_1,m_2;m_1,e_2}:A-S-S");
    require_same(diff(E, {x, y}), hessian, {m1, e1, m2, e2}, real);
    require_same(diff(E, {x, adjoint(y)}), hessian, {m1, e1, m2, e2}, real);
  }

  SECTION("closed shell") {
    // the variables are spelled over the real spaces
    const auto scope = deriv::test::real_spaces();
    const auto x = T(L"κ̂{e_1;m_1}:N-N-S"), y = T(L"κ̂{e_2;m_2}:N-N-S");
    const deriv::test::Data real{.real = {L"f", L"g", L"κ̂"}};
    const auto E =
        mbpt::closed_shell_spintrace(mbpt::hf::rotated_energy(2), {});
    require_same(at_zero(diff(E, {x}), {x}), T(L"4 f{e_1;m_1}:N-S-S"), {m1, e1},
                 real);
    require_same(diff(E, {x, y}),
                 T(L"8 g{m_1,m_2;e_1,e_2}:N-S-S - 4 g{m_2,m_1;e_1,e_2}:N-S-S "
                   L"- 4 g{m_2,e_1;m_1,e_2}:N-S-S + 8 g{e_1,m_2;m_1,e_2}:N-S-S "
                   L"+ 4 δ{m_1;m_2}:N-S-S f{e_1;e_2}:N-S-S - "
                   L"4 δ{e_1;e_2}:N-S-S f{m_1;m_2}:N-S-S"),
                 {m1, e1, m2, e2}, real);
  }

  SECTION("rejections") {
    const auto rejects = [] {
      REQUIRE_THROWS_WITH(mbpt::hf::rotated_energy(1),
                          StartsWith("mbpt::hf::rotated_energy: "));
    };
    {
      auto registry = mbpt::make_legacy_registry();
      registry->remove(L"κ̂");
      const auto scope = mbpt::set_scoped_default_mbpt_context(
          mbpt::Context({.op_registry_ptr = registry}));
      REQUIRE_THROWS_WITH(
          mbpt::hf::rotated_energy(1),
          StartsWith("mbpt::hf::rotated_energy: ") &&
              Catch::Matchers::ContainsSubstring("not registered"));
    }
    {
      // the reference is not the Wick vacuum
      auto context = get_default_context();
      context.set(mbpt::make_min_mr_spaces());
      const auto scope = set_scoped_default_context(context);
      rejects();
    }
    {
      auto context = get_default_context();
      context.set(Vacuum::Physical);
      const auto scope = set_scoped_default_context(context);
      rejects();
    }
    {
      auto context = get_default_context();
      context.set(mbpt::make_legacy_spaces());
      const auto scope = set_scoped_default_context(context);
      REQUIRE_THROWS_WITH(
          mbpt::hf::rotated_energy(1),
          StartsWith("mbpt::hf::rotated_energy: ") &&
              Catch::Matchers::ContainsSubstring("unoccupied space"));
    }
    {
      // pair-specific orbitals
      const auto scope = mbpt::set_scoped_default_mbpt_context(
          mbpt::Context({.csv = mbpt::CSV::Yes,
                         .op_registry_ptr = mbpt::make_legacy_registry()}));
      REQUIRE_THROWS_WITH(
          mbpt::hf::rotated_energy(1),
          StartsWith("mbpt::hf::rotated_energy: ") &&
              Catch::Matchers::ContainsSubstring("pair-specific"));
    }
  }
}

namespace {

/// the Fock matrix element f{P;Q} = h{P;Q} + Σ_m g{P,m;Q,m}, in h and g
std::wstring fock(std::wstring_view p, std::wstring_view q) {
  return L"(h{" + std::wstring(p) + L";" + std::wstring(q) + L"}:A-C-S + g{" +
         std::wstring(p) + L",m_9;" + std::wstring(q) + L",m_9}:A-C-S)";
}

/// the derivative of the energy along κ(ε), first-order coefficients carrying
/// the auxiliary indices @p z (";z_3" or empty)
ExprPtr nuclear_gradient(std::wstring_view z) {
  const std::wstring Z(z);
  const auto T = [](const std::wstring& s) { return deserialize(s); };
  return T(L"h¹{m_1;m_1" + Z + L"}:A-C-S + 1/2 g¹{m_1,m_2;m_1,m_2" + Z +
           L"}:A-C-S - 1/2 s¹{m_1;p_1" + Z + L"}:N-C-S " +
           fock(L"p_1", L"m_1") + L" - 1/2 " + fock(L"m_1", L"p_1") +
           L" s¹{p_1;m_1" + Z + L"}:N-C-S + " + fock(L"m_1", L"e_1") +
           L" κ̂¹{e_1;m_1" + Z + L"}:A-N-S") +
         adjoint(T(L"κ̂¹{e_1;m_1" + Z + L"}:A-N-S")) * T(fock(L"e_1", L"m_1"));
}

}  // namespace

TEST_CASE("deriv/hf/nuclear", "[deriv][hf]") {
  using Catch::Matchers::StartsWith;
  using deriv::at_zero;
  using deriv::diff;
  const auto T = [](std::wstring_view s) { return deserialize(s); };
  const auto eps = ex<Variable>(L"ε");
  const auto kappa = T(L"κ̂{e_9;m_9}:A-N-S");
  const Index m1(L"m_1"), e1(L"e_1");

  SECTION("energy") {
    require_same(mbpt::hf::energy(0),
                 T(L"h{m_1;m_1}:A-C-S + 1/2 g{m_1,m_2;m_1,m_2}:A-C-S"));
    require_same(
        at_zero(diff(mbpt::hf::energy(1), {T(L"κ̂{e_1;m_1}:A-N-S")}), {kappa}),
        T(L"h{m_1;e_1}:A-C-S + g{m_1,m_2;e_1,m_2}:A-C-S"), {m1, e1});
    // every h and every g is one variable
    const auto E = mbpt::hf::energy(2);
    REQUIRE_NOTHROW(diff(E, {T(L"h{p_1;p_2}:A-C-S")}));
    REQUIRE_NOTHROW(diff(E, {T(L"g{p_1,p_2;p_3,p_4}:A-C-S")}));
  }

  SECTION("one-electron") {
    const auto h = mbpt::hf::orthonormalized(T(L"h{m_1;m_1}:A-C-S"), eps, 1);
    REQUIRE_NOTHROW(diff(h, {T(L"s¹{p_1;p_2}:N-C-S")}));
    require_same(
        at_zero(diff(h, {eps}), {eps}),
        T(L"h¹{m_1;m_1}:A-C-S - 1/2 s¹{m_1;p_1}:N-C-S h{p_1;m_1}:A-C-S - "
          L"1/2 h{m_1;p_1}:A-C-S s¹{p_1;m_1}:N-C-S"));
    // a Variable labelled like the overlap, or an overlap under a power, is
    // not expanded
    const auto s = ex<Variable>(L"s");
    require_same(mbpt::hf::orthonormalized(s * T(L"h{m_1;m_1}:A-C-S"), eps, 1),
                 s * h);
    REQUIRE_NOTHROW(mbpt::hf::orthonormalized(
        ex<TensorPower>(T(L"s{m_2;m_3}:N-C-S"), -1) * T(L"h{m_1;m_1}:A-C-S"),
        eps, 1));
    {
      // there op::H()'s one-body tensor is labelled h
      auto context = get_default_context();
      context.set(Vacuum::Physical);
      const auto scope = set_scoped_default_context(context);
      require_same(
          at_zero(diff(mbpt::hf::orthonormalized(
                           T(L"h{p_1;p_2}:A-C-S D{p_2;p_1}:N-N-N"), eps, 1),
                       {eps}),
                  {eps}),
          T(L"(h¹{p_1;p_2}:A-C-S - 1/2 s¹{p_1;p_3}:N-C-S h{p_3;p_2}:A-C-S - "
            L"1/2 h{p_1;p_3}:A-C-S s¹{p_3;p_2}:N-C-S) D{p_2;p_1}:N-N-N"));
    }
    require_same(
        at_zero(diff(mbpt::hf::orthonormalized(
                         T(L"h{p_1;p_2}:A-C-S D{p_2;p_1}:N-N-N"), eps, 2),
                     {eps, eps}),
                {eps}),
        T(L"(h²{p_1;p_2}:A-C-S - s¹{p_1;p_3}:N-C-S h¹{p_3;p_2}:A-C-S - "
          L"h¹{p_1;p_3}:A-C-S s¹{p_3;p_2}:N-C-S - "
          L"1/2 s²{p_1;p_3}:N-C-S h{p_3;p_2}:A-C-S - "
          L"1/2 h{p_1;p_3}:A-C-S s²{p_3;p_2}:N-C-S + "
          L"3/4 s¹{p_1;p_3}:N-C-S s¹{p_3;p_4}:N-C-S h{p_4;p_2}:A-C-S + "
          L"3/4 h{p_1;p_3}:A-C-S s¹{p_3;p_4}:N-C-S s¹{p_4;p_2}:N-C-S + "
          L"1/2 s¹{p_1;p_3}:N-C-S h{p_3;p_4}:A-C-S s¹{p_4;p_2}:N-C-S) "
          L"D{p_2;p_1}:N-N-N"));
  }

  SECTION("an index named like a new one") {
    // t: the first temporary indices of the complete space, then the counter
    // back at them
    Index::reset_tmp_index();
    container::svector<Index> t;
    for (int k = 0; k != 4; ++k)
      t.push_back(Index::make_tmp_index(Index(L"p_1").space()));
    Index::reset_tmp_index();
    const Index p8(L"p_8");
    const auto tensor = [](std::wstring label, const Index& a, const Index& b,
                           BraKetSymmetry braket = BraKetSymmetry::Conjugate,
                           Symmetry symmetry = Symmetry::Antisymm) {
      return ex<Tensor>(std::move(label), bra{a}, ket{b}, symmetry, braket,
                        ColumnSymmetry::Symm);
    };
    const auto s1 = [&](const Index& a, const Index& b) {
      return tensor(L"s¹", a, b, BraKetSymmetry::Conjugate, Symmetry::Nonsymm);
    };
    const auto D = tensor(L"D", t[2], t[3], BraKetSymmetry::Nonsymm);
    const auto half = ex<Constant>(rational{1, 2});
    require_same(at_zero(diff(mbpt::hf::orthonormalized(
                                  tensor(L"h", t[0], t[1]) * D, eps, 1),
                              {eps}),
                         {eps}),
                 (tensor(L"h¹", t[0], t[1]) -
                  half * s1(t[0], p8) * tensor(L"h", p8, t[1]) -
                  half * tensor(L"h", t[0], p8) * s1(p8, t[1])) *
                     D,
                 {t.begin(), t.end()});
  }

  SECTION("gradient and response") {
    const auto E = mbpt::hf::orthonormalized(mbpt::hf::energy(2), eps, 2);
    const auto Et = at_zero(mbpt::series(E, L"κ̂", eps, 2), {kappa});
    require_same(at_zero(diff(Et, {eps}), {eps}), nuclear_gradient(L""));

    const auto G = diff(E, {adjoint(T(L"κ̂{e_1;m_1}:A-N-S"))});
    const auto R1 = at_zero(
        diff(at_zero(mbpt::series(G, L"κ̂", eps, 1), {kappa}), {eps}), {eps});
    const auto b =
        T(L"h¹{e_1;m_1}:A-C-S + g¹{e_1,m_2;m_1,m_2}:A-C-S - "
          L"1/2 s¹{e_1;p_1}:N-C-S " +
          fock(L"p_1", L"m_1") + L" - 1/2 " + fock(L"e_1", L"p_1") +
          L" s¹{p_1;m_1}:N-C-S - "
          L"1/2 s¹{m_2;p_1}:N-C-S g{e_1,p_1;m_1,m_2}:A-C-S - "
          L"1/2 g{e_1,m_2;m_1,p_1}:A-C-S s¹{p_1;m_2}:N-C-S");
    const auto expected =
        b +
        T(fock(L"e_1", L"e_2") + L" κ̂¹{e_2;m_1}:A-N-S - κ̂¹{e_1;m_2}:A-N-S " +
          fock(L"m_2", L"m_1") +
          L" + g{e_1,m_2;m_1,e_2}:A-C-S κ̂¹{e_2;m_2}:A-N-S") +
        T(L"g{e_1,e_2;m_1,m_2}:A-C-S") * adjoint(T(L"κ̂¹{e_2;m_2}:A-N-S"));
    require_same(R1, expected, {m1, e1});
  }

  SECTION("strength tensor") {
    const auto scope = deriv::test::with_batching_spaces();
    const auto strength = T(L"ε{;;z_1}:N-N-N");
    const Index z3(L"z_3");
    const auto E = mbpt::hf::orthonormalized(mbpt::hf::energy(2), strength, 2);
    const auto Et = at_zero(mbpt::series(E, L"κ̂", strength, 2), {kappa});
    require_same(at_zero(diff(Et, {T(L"ε{;;z_3}:N-N-N")}), {strength}),
                 nuclear_gradient(L";z_3"), {z3});
  }

  SECTION("rejections") {
    const auto rejects = [&](const ExprPtr& e) {
      REQUIRE_THROWS_WITH(mbpt::hf::orthonormalized(e, eps, 1),
                          StartsWith("mbpt::hf::orthonormalized: "));
    };
    rejects(mbpt::hf::rotated_energy(1));
    rejects(ex<TensorPower>(T(L"f{m_1;m_1}:A-C-S"), -1));
    rejects(mbpt::op::H());
    rejects(adjoint(T(L"h{e_1;m_1}:A-N-S")));
    rejects(adjoint(T(L"g{e_1,e_2;m_1,m_2}:A-N-S")));
    // pair-specific orbitals, spin orbitals of one spin
    rejects(T(L"h{e_1<m_2>;m_1}:A-C-S"));
    rejects(T(L"g{m_1,m_2;m_1,e_1<m_3>}:A-C-S"));
    rejects(T(L"h{m↑_1;m↑_1}:A-C-S"));
    for (const auto* leaf : {L"h{m_1;m_1}:A-C-S", L"g{m_1,m_2;m_1,m_2}:A-C-S"})
      REQUIRE_THROWS_WITH(
          mbpt::hf::orthonormalized(
              T(L"h{m_3;m_3}:A-C-S") * ex<TensorPower>(T(leaf), -1), eps, 1),
          StartsWith("mbpt::series: "));
    REQUIRE_THROWS_WITH(mbpt::hf::orthonormalized(nullptr, eps, 1),
                        StartsWith("mbpt::series: "));
    const auto h = T(L"h{p_1;p_2}:A-C-S D{p_2;p_1}:N-N-N");
    REQUIRE_THROWS_WITH(
        mbpt::hf::orthonormalized(h + T(L"h¹{p_1;p_2}:A-C-S"), eps, 1),
        StartsWith("mbpt::series: "));
    {
      const auto scope = deriv::test::with_batching_spaces();
      REQUIRE_THROWS_WITH(
          mbpt::hf::orthonormalized(h, T(L"ε{p_1;p_2}:N-N-N"), 1),
          StartsWith("mbpt::series: "));
    }
    {
      auto context = get_default_context();
      context.set(SPBasis::Spinfree);
      const auto scope = set_scoped_default_context(context);
      REQUIRE_THROWS_WITH(mbpt::hf::energy(1),
                          StartsWith("mbpt::hf::energy: "));
    }
    {
      // orthonormalized needs neither κ̂ nor the spaces of the rotation
      auto registry = mbpt::make_legacy_registry();
      registry->remove(L"κ̂");
      const auto scope = mbpt::set_scoped_default_mbpt_context(
          mbpt::Context({.csv = mbpt::CSV::Yes, .op_registry_ptr = registry}));
      REQUIRE_THROWS_WITH(mbpt::hf::energy(1),
                          StartsWith("mbpt::hf::energy: "));
      REQUIRE_NOTHROW(mbpt::hf::orthonormalized(h, eps, 1));
    }
    {
      auto context = get_default_context();
      context.set(mbpt::make_legacy_spaces());
      const auto scope = set_scoped_default_context(context);
      REQUIRE_NOTHROW(mbpt::hf::orthonormalized(
          T(L"h{p_1;p_2}:A-C-S D{p_2;p_1}:N-N-N"), eps, 1));
    }
  }
}
