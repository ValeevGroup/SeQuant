//
// Created by Ajay Melekamburath on 7/12/23.
//

#ifndef SEQUANT_CCRESPONSE_HPP
#define SEQUANT_CCRESPONSE_HPP

#include <SeQuant/domain/mbpt/models/cc.hpp>

// clang-format off
// 1. implement the one-body operator to represent perturbations
// 2. implement the tensor/operator to represent perturbed operators
// 3. add a ccresponse class to derive the residual equations [x]
// 4. construct perturbed Hamiltonian (try first order for now)
// 5. solve perturbed t amplitudes (try first order for now)
// 6. solve perturbed lambda amplitudes (try first order for now)
// clang-format on

namespace sequant::mbpt::sr {

/// derives residual equations for coupled-cluster response theory
/// @param N coupled-cluster rank
/// @param R order of perturbation

class ccresponse {
  size_t N, R;

 public:
  // functions for deriving derive perturbed t and lambda amplitudes
  std::vector<sequant::ExprPtr> t();
  std::vector<sequant::ExprPtr> lambda();
};
}  // namespace sequant::mbpt::sr

#endif  // SEQUANT_CCRESPONSE_HPP
