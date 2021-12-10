//
// Created by Eduard Valeyev on 2019-03-26.
//

#ifndef SEQUANT_DOMAIN_MBPT_OP_HPP
#define SEQUANT_DOMAIN_MBPT_OP_HPP

#include <string>
#include <vector>

namespace sequant {
namespace mbpt {

/// Operators
enum class OpType {
  h,   //!< 1-body Hamiltonian
  f,   //!< Fock operator
  g,   //!< 2-body Coulomb
  t,   //!< cluster amplitudes
  l,   //!< deexcitation cluster amplitudes
  A,   //!< antisymmetrizer
  S,   //!< particle symmetrizer
  L,   //!< left-hand eigenstate
  R,   //!< right-hand eigenstate
  R12, //!< geminal kernel
  GR   //!< GR kernel from f12 theory
};

/// Operator character relative to Fermi vacuum
enum class OpClass { ex, deex, gen };

/// @return the tensor labels in the cardinal order
std::vector<std::wstring> cardinal_tensor_labels();

std::wstring to_wstring(OpType op);
OpClass to_class(OpType op);

}  // namespace mbpt
}  // namespace sequant

#endif //SEQUANT_DOMAIN_MBPT_OP_HPP
