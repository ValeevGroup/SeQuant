//
// Created by Ajay Melekamburath on 4/19/25.
//

#include <SeQuant/core/context.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/domain/mbpt/context.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/models/cc.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>

// SEQUANT_DOC_EXAMPLE: DO NOT RUN

int main() {
  // start-snippet-0
  using namespace sequant;
  using namespace sequant::mbpt;
  set_default_context({.index_basis_registry_shared_ptr = make_min_sr_spaces(),
                       .vacuum = Vacuum::SingleProduct});
  set_default_mbpt_context({.op_registry_ptr = make_legacy_registry()});
  // end-snippet-0

  // start-snippet-1
  // Traditional CCSD
  auto t_eqs = CC{2}.t();
  std::cout << "T1: " << to_latex(t_eqs[1]) << "\n"
            << "T2: " << to_latex(t_eqs[2]) << "\n";

  // Lambda equations
  auto l_eqs = CC{2}.λ();
  std::cout << "λ1: " << to_latex(l_eqs[1]) << "\n"
            << "λ2: " << to_latex(l_eqs[2]) << "\n";

  // Unitary CCSD
  auto Ut_eqs = CC(2, {.ansatz = CC::Ansatz::U, .hbar_comm_rank = 4})
                    .t();  // Use 4th-order commutator expansion
  std::cout << "T1 (UCC): " << to_latex(Ut_eqs[1]) << "\n"
            << "T2 (UCC): " << to_latex(Ut_eqs[2]) << "\n";

  // Orbital-optimized CCSD
  auto oT_eqs = CC(2, {.ansatz = CC::Ansatz::oT}).t();
  std::cout << "T2 (oT): " << to_latex(oT_eqs[2]) << "\n";
  // end-snippet-1

  // start-snippet-2
  // EE-EOM-CCSD (excitation energy)
  auto r_eqs = CC{2}.eom_r(nₚ(2), nₕ(2));
  std::cout << "R1: " << to_latex(r_eqs[1]) << "\n"
            << "R2: " << to_latex(r_eqs[2]) << "\n";

  //  EE-EOM-CCSD Left eigenvectors
  auto ee_l_eqs = CC{2}.eom_l(nₚ(2), nₕ(2));
  std::cout << "L1: " << to_latex(ee_l_eqs[1]) << "\n"
            << "L2: " << to_latex(ee_l_eqs[2]) << "\n";

  // IP-EOM-CCSD (ionization potential)
  auto ip_eqs = CC{2}.eom_r(nₚ(0), nₕ(1));
  std::cout << "IP-R1: " << to_latex(ip_eqs[0]) << "\n";

  // EA-EOM-CCSD (electron attachment)
  auto ea_eqs = CC{2}.eom_r(nₚ(1), nₕ(0));
  std::cout << "EA-R1: " << to_latex(ea_eqs[0]) << "\n";
  // end-snippet-2

  // start-snippet-3
  // First-order perturbed amplitude equations
  auto t_pt = CC{2}.tʼ(1, 1);
  std::cout << "t1 perturbed: " << to_latex(t_pt[1]) << "\n"
            << "t2 perturbed: " << to_latex(t_pt[2]) << "\n";

  // First-order perturbed Lambda amplitude equations
  auto l_pt = CC{2}.λʼ(1, 1);
  std::cout << "λ1 perturbed: " << to_latex(l_pt[1]) << "\n"
            << "λ2 perturbed: " << to_latex(l_pt[2]) << "\n";
  // end-snippet-3

  // start-snippet-4
  // Unitary CCSDT with 6th-order commutator expansion
  auto Ut_ccsdt = CC(3, {.ansatz = CC::Ansatz::U, .hbar_comm_rank = 6}).t();
  // end-snippet-4

  // start-snippet-5
  // Generate CCSD R2 equation
  auto ccsd_eqs = CC{2}.t();
  auto t2_eq = t_eqs[2];

  // Convert to the closed-shell spin-traced form
  auto t2_cs = closed_shell_CC_spintrace(t2_eq);
  std::cout << "Closed-shell spin-traced CCSD-R2: " << to_latex(t2_cs) << "\n";

  // Convert to the open-shell spin-traced form
  auto t2_os = open_shell_CC_spintrace(
      t2_eq);  // vector of expressions: {ɑɑ, ɑβ, ββ} blocks
  std::cout << "Open-shell spin-traced CCSD-R2:\n";
  for (const auto& eq : t2_os) {
    std::cout << to_latex(eq) << "\n";
  }
  // end-snippet-5

  // start-snippet-6
  // Evaluate H̄ outside the class with the connectivity that matches its form
  CC cc{2};
  auto hbar = cc.hbar();
  auto r2_eq =
      op::ref_av(op::P(nₚ(2)) * hbar, {.connect = cc.hbar_connections()});
  // end-snippet-6

  return 0;
}
