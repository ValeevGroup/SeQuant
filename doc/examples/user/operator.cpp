//
// Created by Ajay Melekamburath on 5/3/25.
//

#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/domain/mbpt/context.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/op.hpp>
#include <SeQuant/domain/mbpt/vac_av.hpp>

int main() {
  // start-snippet-0
  using namespace sequant;
  using namespace sequant::mbpt;
  set_default_context({.index_space_registry_shared_ptr = make_min_sr_spaces(),
                       .vacuum = Vacuum::SingleProduct});
  set_default_mbpt_context({.op_registry_ptr = make_legacy_registry()});
  // end-snippet-0

  // start-snippet-1
  using sequant::mbpt::op_t;
  using sequant::mbpt::qns_t;

  // constructor takes label, tensor form and QuantumNumberChange
  op_t fock_op([]() -> std::wstring_view { return L"f"; },
               []() -> ExprPtr {
                 return ex<Tensor>(L"f", bra{L"p_1"}, ket{L"p_2"}) *
                        ex<FNOperator>(cre({L"p_1"}), ann({L"p_2"}));
               },
               [](qns_t& qns) { qns += sequant::mbpt::general_type_qns(1); });
  // end-snippet-1

  // start-snippet-2
  // get the label
  auto label = fock_op.label();

  // apply operator to another state
  auto vac = qns_t{};  // vacuum state for example
  auto modified = fock_op(vac);

  // access the tensor form
  auto fock_tensor = fock_op.tensor_form();
  // end-snippet-2
  (void)label;

  // start-snippet-3
  using namespace sequant::mbpt::op;

  // Construct a double excitation operator (T2)
  auto T2 = T_(2);

  // Construct a two-body Hamiltonian
  auto H2 = H_(2);

  // Product: H2 * T2 * T2
  auto expr1 = H2 * T2 * T2;

  // Screening: can an expression raise the vacuum to double excitation?
  bool raises_to_double = raises_vacuum_to_rank(expr1, 2);  // true

  // Screening: can an expression raise the vacuum up to double excitation?
  bool raises_up_to_double = raises_vacuum_up_to_rank(expr1, 2);  // true

  // Screening: can an expression lower a double excitation to the vacuum?
  auto lambda2 = Λ_(2);
  auto expr2 = lambda2 * H_(1);
  bool lowers_to_vacuum = lowers_rank_to_vacuum(expr2, 2);  // true
  // end-snippet-3
  (void)raises_to_double;
  (void)raises_up_to_double;
  (void)lowers_to_vacuum;

  // start-snippet-4
  using namespace sequant::mbpt;

  auto expr = op::H(2) * op::T(2) * op::T(2);
  auto result = op::vac_av(op::P(2) * expr);
  // vac_av is equivalent to ref_av for single-determinant reference:
  // auto result = op::ref_av(op::P(2) * expr);

  std::wcout << "Result: " << to_latex(result) << "\n";
  // end-snippet-4

  return 0;
}
