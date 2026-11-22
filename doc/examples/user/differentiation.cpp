#include <SeQuant/core/context.hpp>
#include <SeQuant/core/deriv/deriv.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/rational.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/string.hpp>
#include <SeQuant/domain/mbpt/context.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/models/cc.hpp>
#include <SeQuant/domain/mbpt/models/hf.hpp>
#include <SeQuant/domain/mbpt/op.hpp>
#include <SeQuant/domain/mbpt/series.hpp>

#include <string>
#include <vector>

int main() {
  // start-snippet-0
  using namespace sequant;
  using deriv::at_zero;
  using deriv::diff;
  set_default_context(
      {.index_basis_registry_shared_ptr = mbpt::make_min_sr_spaces(),
       .vacuum = Vacuum::SingleProduct});
  mbpt::set_default_mbpt_context(
      {.op_registry_ptr = mbpt::make_legacy_registry()});

  // throws unless the two expressions have the same canonical form; the names
  // of external indices count, except those a symmetrizer carries, which
  // simplify renames
  const auto require_equivalent = [](ExprPtr result, ExprPtr expected) {
    auto context = get_default_context();
    context.set(CanonicalizeOptions::default_options().copy_and_set(
        CanonicalizeOptions::IgnoreNamedIndexLabel::No));
    const auto scope = set_scoped_default_context(context);
    simplify(result);
    simplify(expected);
    if (!(*result == *expected))
      throw Exception("unexpected result: " + toUtf8(serialize(result)));
  };
  // end-snippet-0

  {
    // start-snippet-1
    const auto x = ex<Variable>(L"x");
    const auto E = deserialize(L"3 x^(2) y + x");
    require_equivalent(diff(E, {x}), deserialize(L"6 x y + 1"));
    // a second derivative: the variables are taken left to right
    require_equivalent(diff(E, {x, x}), deserialize(L"6 y"));
    // x and its adjoint x^* are independent variables
    const auto x_star = adjoint(x);
    require_equivalent(diff(x * x_star, {x}), x_star);
    require_equivalent(diff(x * x_star, {x_star}), x);
    // end-snippet-1
  }

  {
    // start-snippet-2
    // the rest of the term, under the names of the variable's indices
    const auto t1 = deserialize(L"t{a_2;i_2}:N-N-N");
    require_equivalent(
        diff(deserialize(L"f{i_1;a_1}:N-N-N t{a_1;i_1}:N-N-N"), {t1}),
        deserialize(L"f{i_2;a_2}:N-N-N"));
    // the external indices of the function stay: Kronecker deltas
    require_equivalent(diff(deserialize(L"t{a_1;i_1}:N-N-N"), {t1}),
                       deserialize(L"δ{a_1;a_2}:N-S-S δ{i_2;i_1}:N-S-S"));
    // an antisymmetric variable: SeQuant's antisymmetrizer on the result
    const auto t2 = deserialize(L"t{a_3,a_4;i_3,i_4}:A-N-S");
    require_equivalent(
        diff(deserialize(L"t{a_1,a_2;i_1,i_2}:A-N-S"), {t2}),
        deserialize(L"4 Â{a_3,a_4;i_3,i_4}:A-N-S δ{a_1;a_3}:N-S-S "
                    L"δ{a_2;a_4}:N-S-S δ{i_3;i_1}:N-S-S δ{i_4;i_2}:N-S-S"));
    require_equivalent(
        diff(deserialize(
                 L"1/4 g{i_1,i_2;a_1,a_2}:A-C-S t{a_1,a_2;i_1,i_2}:A-N-S"),
             {t2}),
        deserialize(L"Â{a_3,a_4;i_3,i_4}:A-N-S g{i_3,i_4;a_3,a_4}:A-C-S"));
    // end-snippet-2
  }

  {
    // start-snippet-7
    // h over all orbitals holds the elements of h{i_3;a_4} in both slot
    // orders: one term per order
    require_equivalent(diff(deserialize(L"G{p_2;p_1}:N-N-N h{p_1;p_2}:N-S-N"),
                            {deserialize(L"h{i_3;a_4}:N-S-N")}),
                       deserialize(L"G{a_4;i_3}:N-N-N + G{i_3;a_4}:N-N-N"));
    // end-snippet-7
  }

  {
    // start-snippet-3
    // the coupled-cluster singles and doubles Lagrangian <0|(1 + Λ) H̄|0>
    mbpt::CC cc{2};
    auto L = mbpt::op::ref_av((1 + mbpt::op::Λ(2)) * cc.hbar(),
                              {.connect = cc.hbar_connections()});
    simplify(L);
    // the doubles multipliers and amplitudes
    const auto λ2 = deserialize(L"λ{i_1,i_2;a_1,a_2}:A-N-S");
    const auto t2 = deserialize(L"t{a_1,a_2;i_1,i_2}:A-N-S");
    require_equivalent(diff(L, {λ2}), cc.t()[2]);
    require_equivalent(diff(L, {t2}), cc.λ()[2]);
    // end-snippet-3
  }

  {
    // start-snippet-4
    // complex index spaces (the default): t and t⁺ are independent
    const auto E =
        L"f{i_1;a_1}:N-N-N t{a_1;i_1}:N-N-N + "
        L"w{a_1;i_1}:N-N-N t⁺{i_1;a_1}:N-N-N";
    const auto t = L"t{a_2;i_2}:N-N-N";
    require_equivalent(diff(deserialize(E), {deserialize(t)}),
                       deserialize(L"f{i_2;a_2}:N-N-N"));
    // real index spaces: t⁺ holds the elements of t
    auto registry = mbpt::make_min_sr_spaces();
    std::vector<std::wstring> keys;
    for (const auto& space : *registry) keys.push_back(space.base_key());
    for (const auto& key : keys)
      if (!key.empty()) registry->field(key, Field::Real);
    auto context = get_default_context();
    context.set(registry);
    const auto real = set_scoped_default_context(context);
    require_equivalent(diff(deserialize(E), {deserialize(t)}),
                       deserialize(L"f{i_2;a_2}:N-N-N + w{a_2;i_2}:N-N-N"));
    // a bra-ket symmetric tensor over one space, in any field: G + Gᵀ
    require_equivalent(diff(deserialize(L"G{p_2;p_1}:N-N-N h{p_1;p_2}:N-S-N"),
                            {deserialize(L"h{p_3;p_4}:N-S-N")}),
                       deserialize(L"G{p_4;p_3}:N-N-N + G{p_3;p_4}:N-N-N"));
    // end-snippet-4
  }

  {
    // start-snippet-5
    // a function of the strength ε of a perturbation u of the amplitudes t
    const auto ε = ex<Variable>(L"ε");
    const auto E = deserialize(
        L"f{i_1;a_1}:N-N-N (t{a_1;i_1}:N-N-N + ε u{a_1;i_1}:N-N-N) + "
        L"1/2 g{i_1,i_2;a_1,a_2}:A-C-S "
        L"(t{a_1;i_1}:N-N-N + ε u{a_1;i_1}:N-N-N) "
        L"(t{a_2;i_2}:N-N-N + ε u{a_2;i_2}:N-N-N)");
    // the first- and second-order responses: differentiate, then set ε = 0
    require_equivalent(at_zero(diff(E, {ε}), {ε}),
                       deserialize(L"f{i_1;a_1}:N-N-N u{a_1;i_1}:N-N-N + "
                                   L"g{i_1,i_2;a_1,a_2}:A-C-S t{a_1;i_1}:N-N-N "
                                   L"u{a_2;i_2}:N-N-N"));
    require_equivalent(at_zero(diff(E, {ε, ε}), {ε}),
                       deserialize(L"g{i_1,i_2;a_1,a_2}:A-C-S u{a_1;i_1}:N-N-N "
                                   L"u{a_2;i_2}:N-N-N"));
    // end-snippet-5
  }

  {
    // start-snippet-6
    // a variable under an inverse: D⁻¹ r by D
    const auto E =
        ex<TensorPower>(deserialize(L"D{i_1;a_1}:N-N-N"), rational(-1)) *
        deserialize(L"r{a_1;i_1}:N-N-N");
    std::string message;
    try {
      (void)diff(E, {deserialize(L"D{i_2;a_2}:N-N-N")});
    } catch (const Exception& error) {
      message = error.what();
    }
    if (message != "deriv::diff: the variable stands under a TensorPower")
      throw Exception("unexpected message: " + message);
    // at_zero: a leaf over all orbitals, with a variable over the virtual
    // ones, would be neither zero nor unchanged
    message.clear();
    try {
      (void)at_zero(deserialize(L"W{i_1;p_1}:N-N-N k{p_1;i_1}:N-N-N"),
                    {deserialize(L"k{a_2;i_2}:N-N-N")});
    } catch (const Exception& error) {
      message = error.what();
    }
    if (message !=
        "deriv::at_zero: a leaf holds elements of a variable and lies inside "
        "none of the variables")
      throw Exception("unexpected message: " + message);
    // end-snippet-6
  }

  {
    // start-snippet-8
    // the amplitudes t as a series in the strength ε of a perturbation
    const auto ε = ex<Variable>(L"ε");
    const auto E = deserialize(L"f{i_1;a_1}:N-C-S t{a_1;i_1}:N-N-N");
    const auto S = mbpt::series(E, L"t", ε, 2);
    require_equivalent(S, deserialize(L"f{i_1;a_1}:N-C-S (t{a_1;i_1}:N-N-N + "
                                      L"ε t¹{a_1;i_1}:N-N-N + "
                                      L"1/2 ε^(2) t²{a_1;i_1}:N-N-N)"));
    // the first derivative at ε = 0 has the first-order coefficient t¹
    require_equivalent(at_zero(diff(S, {ε}), {ε}),
                       deserialize(L"f{i_1;a_1}:N-C-S t¹{a_1;i_1}:N-N-N"));
    // end-snippet-8
  }

  {
    // start-snippet-9
    // the energy through second order in the rotation, and two elements of
    // the rotation, spelled with the symmetries of its leaves
    const auto E = mbpt::hf::rotated_energy(2);
    const auto x = deserialize(L"κ̂{a_1;i_1}:A-N-S");
    const auto y = deserialize(L"κ̂{a_2;i_2}:A-N-S");
    // the orbital gradient
    require_equivalent(at_zero(diff(E, {x}), {x}),
                       deserialize(L"f{i_1;a_1}:A-C-S"));
    // the blocks A and B of the orbital Hessian
    require_equivalent(diff(E, {x, adjoint(y)}),
                       deserialize(L"f{a_2;a_1}:A-C-S δ{i_1;i_2}:N-S-S - "
                                   L"f{i_1;i_2}:A-C-S δ{a_2;a_1}:N-S-S + "
                                   L"g{i_1,a_2;a_1,i_2}:A-C-S"));
    require_equivalent(diff(E, {x, y}),
                       deserialize(L"g{i_1,i_2;a_1,a_2}:A-C-S"));
    // end-snippet-9
  }

  {
    // start-snippet-10
    // the energy in h and g, in the integrals of orthonormalized orbitals,
    // through second order in the strength ε of a nuclear displacement
    const auto ε = ex<Variable>(L"ε");
    const auto E = mbpt::hf::orthonormalized(mbpt::hf::energy(2), ε, 2);
    // the orbitals follow the displacement: κ̂(ε) = ε κ̂¹ + 1/2 ε² κ̂²
    const auto κ = deserialize(L"κ̂{a_9;i_9}:A-N-S");
    const auto Et = at_zero(mbpt::series(E, L"κ̂", ε, 2), {κ});
    // the Fock matrix in h and g, F{P;Q} = h{P;Q} + Σ_j g{P,j;Q,j}
    const auto F = [](const std::wstring& p, const std::wstring& q) {
      return L"(h{" + p + L";" + q + L"}:A-C-S + g{" + p + L",i_9;" + q +
             L",i_9}:A-C-S)";
    };
    // the nuclear gradient dE/dε: the differentiated integrals, the
    // differentiated overlap with F, and the orbital gradient times κ̂¹ and
    // its adjoint, zero where F{i;a} = 0, kept
    require_equivalent(
        at_zero(diff(Et, {ε}), {ε}),
        deserialize(L"h¹{i_1;i_1}:A-C-S + 1/2 g¹{i_1,i_2;i_1,i_2}:A-C-S - "
                    L"1/2 s¹{i_1;p_1}:N-C-S " +
                    F(L"p_1", L"i_1") + L" - 1/2 " + F(L"i_1", L"p_1") +
                    L" s¹{p_1;i_1}:N-C-S + " + F(L"i_1", L"a_1") +
                    L" κ̂¹{a_1;i_1}:A-N-S + κ̂¹⁺{i_1;a_1}:A-N-S " +
                    F(L"a_1", L"i_1")));
    // the first-order response equation: the derivative of ∂E/∂κ̂⁺ along
    // κ̂(ε), A κ̂¹ + B κ̂¹⁺ + b = 0
    const auto G = diff(E, {adjoint(deserialize(L"κ̂{a_1;i_1}:A-N-S"))});
    const auto R1 =
        at_zero(diff(at_zero(mbpt::series(G, L"κ̂", ε, 1), {κ}), {ε}), {ε});
    require_equivalent(
        R1, deserialize(F(L"a_1", L"a_2") +
                        L" κ̂¹{a_2;i_1}:A-N-S - κ̂¹{a_1;i_2}:A-N-S " +
                        F(L"i_2", L"i_1") +
                        L" + g{a_1,i_2;i_1,a_2}:A-C-S κ̂¹{a_2;i_2}:A-N-S + "
                        L"g{a_1,a_2;i_1,i_2}:A-C-S κ̂¹⁺{i_2;a_2}:A-N-S + "
                        L"h¹{a_1;i_1}:A-C-S + g¹{a_1,i_2;i_1,i_2}:A-C-S - "
                        L"1/2 s¹{a_1;p_1}:N-C-S " +
                        F(L"p_1", L"i_1") + L" - 1/2 " + F(L"a_1", L"p_1") +
                        L" s¹{p_1;i_1}:N-C-S - "
                        L"1/2 s¹{i_2;p_1}:N-C-S g{a_1,p_1;i_1,i_2}:A-C-S - "
                        L"1/2 g{a_1,i_2;i_1,p_1}:A-C-S s¹{p_1;i_2}:N-C-S"));
    // end-snippet-10
  }

  {
    // start-snippet-11
    // a strength tensor over the batching space z, registered before the
    // context is set: one strength per nuclear coordinate
    auto registry = mbpt::make_min_sr_spaces();
    mbpt::add_batching_spaces(registry);
    auto context = get_default_context();
    context.set(registry);
    const auto scope = set_scoped_default_context(context);
    const auto ε = deserialize(L"ε{;;z_1}:N-N-N");
    const auto E = mbpt::hf::orthonormalized(mbpt::hf::energy(2), ε, 2);
    const auto Et = at_zero(mbpt::series(E, L"κ̂", ε, 2),
                            {deserialize(L"κ̂{a_9;i_9}:A-N-S")});
    // the gradient, over z_3, and the Hessian, over z_3 and z_4
    const auto ε3 = deserialize(L"ε{;;z_3}:N-N-N");
    const auto ε4 = deserialize(L"ε{;;z_4}:N-N-N");
    const auto gradient = at_zero(diff(Et, {ε3}), {ε});
    const auto hessian = at_zero(diff(Et, {ε3, ε4}), {ε});
    // the Hessian holds the orbital gradient times κ̂², kept
    bool κ2 = false;
    hessian->visit(
        [&κ2](const ExprPtr& leaf) {
          κ2 |= leaf->is<Tensor>() && leaf->as<Tensor>().label() == L"κ̂²";
        },
        /*atoms_only=*/true);
    if (gradient->is_zero() || !κ2) throw Exception("unexpected Hessian");
    // end-snippet-11
  }

  return 0;
}
