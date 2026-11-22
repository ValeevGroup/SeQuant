#include <SeQuant/core/context.hpp>
#include <SeQuant/core/deriv/deriv.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/domain/mbpt/context.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/models/hf.hpp>
#include <SeQuant/domain/mbpt/op_registry.hpp>
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
  const auto x = T(L"κₒ{V_1;O_1}:A-N-S"), y = T(L"κₒ{V_2;O_2}:A-N-S");

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
  const auto x_cs = T(L"κₒ{V_1;O_1}:N-N-S");
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
      REQUIRE(registry->to_class(L"κₒ") == mbpt::OpClass::Ex);
    REQUIRE(adjoint(T(L"κₒ{e_1;m_1}:A-N-S"))->as<Tensor>().label() == L"κₒ⁺");
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
    const deriv::test::Data real{.real = {L"f", L"g", L"κₒ"}};
    const auto x = T(L"κₒ{e_1;m_1}:A-N-S"), y = T(L"κₒ{e_2;m_2}:A-N-S");
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
    const auto x = T(L"κₒ{e_1;m_1}:N-N-S"), y = T(L"κₒ{e_2;m_2}:N-N-S");
    const deriv::test::Data real{.real = {L"f", L"g", L"κₒ"}};
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
      registry->remove(L"κₒ");
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
