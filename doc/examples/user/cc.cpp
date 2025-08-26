//
// Created by Ajay Melekamburath on 4/19/25.
//

#include <SeQuant/core/context.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/models/cc.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>

// SEQUANT_DOC_EXAMPLE: DO NOT RUN

int main() {
  // start-snippet-0
  using namespace sequant;
  using namespace sequant::mbpt;
  set_default_context({.index_space_registry_shared_ptr = make_min_sr_spaces(),
                       .vacuum = Vacuum::SingleProduct});
  // end-snippet-0

  // start-snippet-1
  // Traditional CCSD
  auto t_eqs = CC{2}.t();
  std::wcout << "T1: " << to_latex(t_eqs[1]) << "\n"
             << "T2: " << to_latex(t_eqs[2]) << "\n";

  // Lambda equations
  auto l_eqs = CC{2}.λ();
  std::wcout << "λ1: " << to_latex(l_eqs[1]) << "\n"
             << "λ2: " << to_latex(l_eqs[2]) << "\n";

  // Unitary CCSD
  auto Ut_eqs =
      CC{2, CC::Ansatz::U}.t(4);  // Use 4th-order commutator expansion
  std::wcout << "T1 (UCC): " << to_latex(Ut_eqs[1]) << "\n"
             << "T2 (UCC): " << to_latex(Ut_eqs[2]) << "\n";

  // Orbital-optimized CCSD
  auto oT_eqs = CC{2, CC::Ansatz::oT}.t();
  std::wcout << "T2 (oT): " << to_latex(oT_eqs[2]) << "\n";
  // end-snippet-1

  // start-snippet-2
  // EE-EOM-CCSD (excitation energy)
  auto r_eqs = CC{2}.eom_r(nₚ(2), nₕ(2));
  std::wcout << "R1: " << to_latex(r_eqs[1]) << "\n"
             << "R2: " << to_latex(r_eqs[2]) << "\n";

  //  EE-EOM-CCSD Left eigenvectors
  auto ee_l_eqs = CC{2}.eom_l(nₚ(2), nₕ(2));
  std::wcout << "L1: " << to_latex(l_eqs[1]) << "\n"
             << "L2: " << to_latex(l_eqs[2]) << "\n";

  // IP-EOM-CCSD (ionization potential)
  auto ip_eqs = CC{2}.eom_r(nₚ(0), nₕ(1));
  std::wcout << "IP-R1: " << to_latex(ip_eqs[0]) << "\n";

  // EA-EOM-CCSD (electron attachment)
  auto ea_eqs = CC{2}.eom_r(nₚ(1), nₕ(0));
  std::wcout << "EA-R1: " << to_latex(ea_eqs[0]) << "\n";
  // end-snippet-2

  // start-snippet-3
  // First-order perturbed amplitude equations
  auto t_pt = CC{2}.t_pt(1, 1);
  std::wcout << "T1 perturbed: " << to_latex(t_pt[1]) << "\n"
             << "T2 perturbed: " << to_latex(t_pt[2]) << "\n";

  // First-order perturbed Lambda amplitude equations
  auto l_pt = CC{2}.λ_pt(1, 1);
  std::wcout << "λ1 perturbed: " << to_latex(l_pt[1]) << "\n"
             << "λ2 perturbed: " << to_latex(l_pt[2]) << "\n";
  // end-snippet-3

  // start-snippet-4
  // Unitary CCSDT with 6th-order commutator expansion
  auto Ut_ccsdt = CC{3, CC::Ansatz::U}.t(6);
  // end-snippet-4

  // start-snippet-5
  // Generate CCSD R2 equation
  auto ccsd_eqs = CC{2}.t();
  auto t2_eq = t_eqs[2];

  // Convert to the closed-shell spin-traced form
  auto t2_cs = closed_shell_CC_spintrace(t2_eq);
  std::wcout << "Closed-shell spin-traced CCSD-R2: " << to_latex(t2_cs) << "\n";

  // Convert to the open-shell spin-traced form
  auto t2_os = open_shell_CC_spintrace(
      t2_eq);  // vector of expressions: {ɑɑ, ɑβ, ββ} blocks
  std::wcout << "Open-shell spin-traced CCSD-R2:\n";
  for (const auto& eq : t2_os) {
    std::wcout << to_latex(eq) << "\n";
  }
  // end-snippet-5

  return 0;
}
