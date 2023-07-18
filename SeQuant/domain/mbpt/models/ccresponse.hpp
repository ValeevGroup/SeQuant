//
// Created by Ajay Melekamburath on 7/12/23.
//

#ifndef SEQUANT_CCRESPONSE_HPP
#define SEQUANT_CCRESPONSE_HPP

#include <SeQuant/domain/mbpt/models/cc.hpp>
#include <SeQuant/domain/mbpt/op.hpp>
#include <SeQuant/domain/mbpt/sr/sr.hpp>

// clang-format off
// 1. implement the one-body operator to represent perturbations
// 2. implement the tensor/operator to represent perturbed operators
// 3. add a ccresponse class to derive the residual equations [x]
// 4. construct perturbed Hamiltonian (try first order for now)
// 5. solve perturbed t amplitudes (try first order for now)
// 6. solve perturbed lambda amplitudes (try first order for now)
// clang-format on

namespace sequant::mbpt::sr {
// Functions relating to perturbation and response

/// one-body perturbation operator of perturbation order \p r
ExprPtr V();

ExprPtr pT1_(std::size_t Nbra,
             std::size_t Nket = std::numeric_limits<std::size_t>::max());
ExprPtr pLambda1_(std::size_t Nbra,
                  std::size_t Nket = std::numeric_limits<std::size_t>::max());

namespace op {
// perturbation and response related operators
/// one-body perturbation operator of first order
ExprPtr V();

/// perturbed cluster amplitudes
ExprPtr pT1_(std::size_t K);
ExprPtr pT1(std::size_t K);

/// perturted de-excitation cluster amplitudes
ExprPtr pLambda1_(std::size_t K);
ExprPtr pLambda1(std::size_t K);

}  // namespace op

/// derives residual equations for coupled-cluster response theory
/// @param N coupled-cluster rank
/// @param R order of perturbation
class ccresponse {
  size_t N, R;

 public:
  ccresponse(size_t n, size_t r);

  // functions for deriving derive perturbed t and lambda amplitudes
  std::vector<sequant::ExprPtr> t();
  std::vector<sequant::ExprPtr> lambda();
};
}  // namespace sequant::mbpt::sr

#endif  // SEQUANT_CCRESPONSE_HPP
