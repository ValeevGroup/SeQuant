#include <SeQuant/core/deriv/deriv.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/domain/mbpt/models/hf.hpp>
#include <SeQuant/domain/mbpt/series.hpp>

#include <catch2/catch_test_macros.hpp>
#include "contexts.hpp"
#include "evaluator.hpp"
#include "hf_model.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <string_view>
#include <utility>
#include <vector>

using namespace sequant;

namespace {

namespace model = deriv::test::hf_model;
using model::Number;

// The error of a five-point central difference of step d is truncation,
// O(d⁴), plus rounding, O(ε_mach |E| / d) for the first and
// O(ε_mach |E| / d²) for the second derivative; the rounding error is
// absolute, not relative to the derivative. The evaluator's data come from
// std::hash and differ between standard libraries, so the tolerances are
// scaled by max(10, |symbolic value|), an absolute floor, and set with room
// for data sets other than the one this library gives.
constexpr double step = 5e-4, first_order = 1e-8, second_order = 1e-6;
/// how many tolerances apart a wrong input must be
constexpr double margin = 30;

/// the scale of the tolerance for @p symbolic
double scale(Number symbolic) { return std::max(10., std::abs(symbolic)); }

/// f'(0) and f''(0) by five-point central differences of step @p d
std::pair<Number, Number> derivatives(const std::function<Number(double)>& f,
                                      double d) {
  const auto p2 = f(2 * d), p1 = f(d), z = f(0), m1 = f(-d), m2 = f(-2 * d);
  return {(-p2 + 8. * p1 - 8. * m1 + m2) / (12 * d),
          (-p2 + 16. * p1 - 30. * z + 16. * m1 - m2) / (12 * d * d)};
}

/// requires |symbolic − difference| ≤ tolerance · scale(symbolic)
void require_close(std::string_view what, Number symbolic, Number difference,
                   double tolerance) {
  const auto error = std::abs(symbolic - difference);
  INFO(what << ": symbolic " << symbolic << ", finite difference " << difference
            << ", |difference| " << error << ", tolerance " << tolerance
            << " · max(10, |symbolic|)");
  REQUIRE(error <= tolerance * scale(symbolic));
}

/// requires |symbolic − difference| > margin · tolerance · scale(symbolic)
void require_apart(std::string_view what, Number symbolic, Number difference,
                   double tolerance) {
  const auto error = std::abs(symbolic - difference);
  INFO(what << ": symbolic " << symbolic << ", finite difference " << difference
            << ", |difference| " << error << ", tolerance " << tolerance
            << " · max(10, |symbolic|), margin " << margin);
  REQUIRE(error > margin * tolerance * scale(symbolic));
}

template <typename T>
void constant(model::Series<T>& q) {
  q.first.clear();
  q.second.clear();
}

/// E(x) of @p m along a scalar strength
std::function<Number(double)> energy_of(
    const model::Model& m,
    model::Coefficients c = model::Coefficients::orthonormalized_then_rotated) {
  return [&m, c](double x) { return model::energy(m, {x}, c); };
}

}  // namespace

TEST_CASE("deriv/hf/finite-differences", "[deriv][hf]") {
  using deriv::at_zero;
  using deriv::diff;
  using deriv::test::element_of;
  const auto T = [](std::wstring_view s) { return deserialize(s); };
  const deriv::test::Data data{.real = {L"ε"}};
  const auto value = [&data](const ExprPtr& e) {
    return element_of(e, {}, data);
  };
  const auto kappa = T(L"κ̂{e_9;m_9}:A-N-S");

  SECTION("scalar strength") {
    const auto eps = ex<Variable>(L"ε");
    const auto m = model::read({}, data, /*quantum_numbers=*/false);

    // conventions, fixing transposition and conjugation: E(0); κ̂¹⁺ the
    // conjugate transpose of κ̂¹; the first-order change under a step in κ
    // alone, which needs both ∂E/∂κ̂ and ∂E/∂κ̂⁺
    require_close("E(0)", value(mbpt::hf::energy(0)), model::energy(m, {0.}),
                  1e-12);
    const auto k1 = T(L"κ̂¹{e_1;m_1}:A-N-S"), k2 = T(L"κ̂¹{e_2;m_2}:A-N-S");
    for (const auto a : m.orbitals.unoccupied)
      for (const auto i : m.orbitals.occupied)
        REQUIRE(element_of(adjoint(k1),
                           {{L"e_1", m.orbitals.values[a]},
                            {L"m_1", m.orbitals.values[i]}},
                           data) == std::conj(m.kappa.first.at(0)(a, i)));
    auto rotation = m;  // κ(x) = x κ̂¹, the integrals and S constant
    constant(rotation.h);
    constant(rotation.g);
    constant(rotation.s);
    rotation.kappa.second.clear();
    const auto x = T(L"κ̂{e_1;m_1}:A-N-S"), y = T(L"κ̂{e_2;m_2}:A-N-S");
    {
      const auto E = mbpt::hf::energy(1);
      require_close(
          "rotation, first order",
          value(at_zero(diff(E, {x}), {kappa}) * k1 +
                at_zero(diff(E, {adjoint(x)}), {kappa}) * adjoint(k1)),
          derivatives(energy_of(rotation), step).first, first_order);
    }

    // rotation only: the Hessian of energy(2) along κ̂¹
    {
      const auto E = mbpt::hf::energy(2);
      const auto hessian = at_zero(diff(E, {x, y}), {kappa}) * k1 * k2 +
                           ex<Constant>(2) *
                               at_zero(diff(E, {x, adjoint(y)}), {kappa}) * k1 *
                               adjoint(k2) +
                           at_zero(diff(E, {adjoint(x), adjoint(y)}), {kappa}) *
                               adjoint(k1) * adjoint(k2);
      require_close("rotation, second order", value(hessian),
                    derivatives(energy_of(rotation), step).second,
                    second_order);
    }

    // integrals only, no rotation
    {
      auto integrals = m;
      constant(integrals.kappa);
      const auto E = mbpt::hf::orthonormalized(mbpt::hf::energy(0), eps, 2);
      const auto [d1, d2] = derivatives(energy_of(integrals), step);
      require_close("integrals, first order",
                    value(at_zero(diff(E, {eps}), {eps})), d1, first_order);
      require_close("integrals, second order",
                    value(at_zero(diff(E, {eps, eps}), {eps})), d2,
                    second_order);
    }

    // the integrals, and the orbitals along κ(ε)
    const auto E = mbpt::hf::orthonormalized(mbpt::hf::energy(2), eps, 2);
    const auto Et = at_zero(mbpt::series(E, L"κ̂", eps, 2), {kappa});
    const auto E1 = value(at_zero(diff(Et, {eps}), {eps}));
    const auto hessian = at_zero(diff(Et, {eps, eps}), {eps});
    const auto E2 = value(hessian);
    const auto [d1, d2] = derivatives(energy_of(m), step);
    require_close("both, first order", E1, d1, first_order);
    require_close("both, second order", E2, d2, second_order);

    // wrong inputs: the κ̂² terms removed; the orbitals rotated, then
    // orthonormalized, which agrees at first order only
    require_apart("both, second order, without κ̂²",
                  value(at_zero(hessian, {T(L"κ̂²{e_9;m_9}:A-N-S")})), d2,
                  second_order);
    const auto [r1, r2] = derivatives(
        energy_of(m, model::Coefficients::rotated_then_orthonormalized), step);
    require_close("rotated, then orthonormalized, first order", E1, r1,
                  first_order);
    require_apart("rotated, then orthonormalized, second order", E2, r2,
                  second_order);
  }

  SECTION("strength tensor") {
    const auto scope = deriv::test::with_batching_spaces();
    // the evaluator distinguishes spin cases where the indices of its input
    // differ in quantum numbers, so the model has to as well
    const Index z1(L"z_1");
    const bool quantum_numbers =
        z1.space().qns() != Index(L"p_1").space().qns();
    deriv::test::Domains domains;
    deriv::test::add_domain(z1, {}, domains, quantum_numbers, data.occupied);
    const auto& z = domains.front().second;
    REQUIRE(z.size() == 2);
    const auto m = model::read(z, data, quantum_numbers);

    const auto strength = T(L"ε{;;z_1}:N-N-N");
    const auto E = mbpt::hf::orthonormalized(mbpt::hf::energy(2), strength, 2);
    const auto Et = at_zero(mbpt::series(E, L"κ̂", strength, 2), {kappa});
    const auto e3 = T(L"ε{;;z_3}:N-N-N"), e4 = T(L"ε{;;z_4}:N-N-N");
    const auto E1 = at_zero(diff(Et, {e3}), {strength});
    const auto E2 = at_zero(diff(Et, {e3, e4}), {strength});
    // E along x = t (e_a + sign e_b)
    const auto along = [&](std::size_t a, std::size_t b, double sign) {
      return [&m, &z, a, b, sign](double t) {
        std::vector<double> x(z.size(), 0.);
        x[a] += t;
        x[b] += sign * t;
        return model::energy(m, x);
      };
    };
#ifdef SEQUANT_SKIP_LONG_TESTS
    // the first value of z, with itself and with the other one
    const std::size_t first_values = 1;
#else
    const std::size_t first_values = z.size();
#endif
    for (std::size_t a = 0; a != first_values; ++a) {
      const auto [d1, daa] = derivatives(along(a, a, 0), step);
      require_close("E1{z}", element_of(E1, {{L"z_3", z[a]}}, data), d1,
                    first_order);
      for (std::size_t b = a; b != z.size(); ++b) {
        // ∂²E/∂x_a∂x_b = (D₊ − D₋)/4, D± the second derivative along
        // e_a ± e_b
        const auto d2 = a == b ? daa
                               : (derivatives(along(a, b, 1), step).second -
                                  derivatives(along(a, b, -1), step).second) /
                                     4.;
        for (const auto& [u, w] : {std::pair{a, b}, std::pair{b, a}})
          require_close("E2{z,w}",
                        element_of(E2, {{L"z_3", z[u]}, {L"z_4", z[w]}}, data),
                        d2, second_order);
      }
    }
  }
}
