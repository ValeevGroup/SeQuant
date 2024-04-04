#ifndef SEQUANT_DOMAIN_MBPT_MODELS_CC_HPP
#define SEQUANT_DOMAIN_MBPT_MODELS_CC_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr_fwd.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/timer.hpp>

namespace sequant::mbpt::sr {

/// CC is a derivation engine for the coupled-cluster method
class CC {
 public:
  enum Ansatz {
    /// traditional ansatz
    T,
    /// traditional orbital-optimized (singles-free) ansatz
    oT,
    /// unitary ansatz
    U,
    /// unitary orbital-optimized (singles-free) ansatz
    oU
  };

  /// @brief constructs CC engine
  /// @param N coupled cluster excitation rank
  /// @param ansatz the type of CC ansatz
  CC(size_t N, Ansatz ansatz = Ansatz::T);

  /// @return the type of ansatz
  Ansatz ansatz() const;

  /// @return true if the ansatz is unitary
  bool unitary() const;

  /// @brief derives similarity-transformed expressions of mpbt::Operators
  /// @param expr expression to be transformed
  /// @param r order of truncation
  /// @pre expr should be composed of mbpt::Operators
  /// @return transformed expression
  ExprPtr sim_tr(ExprPtr expr, size_t r);

  /// @brief derives t amplitude equations, \f$ \langle P|\bar{H}|0 \rangle = 0
  /// \f$
  /// @param commutator_rank rank of commutators included in \f$ \bar{H} \f$ ;
  /// must be specified for unitary ansatz and/or
  ///   Hamiltonians with particle rank != 2
  /// @param pmax highest particle rank of the projector manifold `\f \langle P
  /// | \f`; the default value is to use
  ///   the cluster operator rank of this engine
  /// @param pmin lowest particle rank of the projector manifold `\f \langle P |
  /// \f`; the default value is 0
  /// @return vector of t amplitude equations, with element `k` containing
  /// equation
  ///   \f$ \langle k |\bar{H}|0 \rangle = 0 \f$ for `k` in the [\p pmin,\p
  ///   pmax] range, and null value otherwise
  [[nodiscard]] std::vector<sequant::ExprPtr> t(
      size_t commutator_rank = 4,
      size_t pmax = std::numeric_limits<size_t>::max(), size_t pmin = 0);

  /// @brief derives λ amplitude equations,
  /// \f$ \langle 0| (1 + \hat{\Lambda}) \frac{d \bar{H}}{d \hat{T}_P} |0
  /// \rangle = 0 \f$
  /// @param commutator_rank rank of commutators included in \f$ \bar{H} \f$ ;
  /// must be specified for unitary ansatz and/or Hamiltonians with particle
  /// rank != 2
  /// @return vector of λ amplitude equations, with element `k` containing
  /// equation
  ///   \f$ \langle 0| (1 + \hat{\Lambda}) \frac{d \bar{H}}{d \hat{T}_k} |0
  ///   \rangle = 0 \f$ for `k` in
  /// the [1,N] range; element 0 is always null
  [[nodiscard]] std::vector<sequant::ExprPtr> λ(size_t commutator_rank = 4);

  // clang-format off
  /// @brief derives perturbed t amplitude equations
  /// @param order order of perturbation
  /// @param rank rank of perturbation operator. r = 1 means one-body perturbation operator
  /// @pre `rank==1 && order==1`, only first order perturbation and one-body perturbation operator is supported now
  /// @return std::vector of perturbed t amplitude equations
  // clang-format on
  [[nodiscard]] std::vector<sequant::ExprPtr> t_pt(size_t order = 1,
                                                   size_t rank = 1);

  // clang-format off
  /// @brief derives perturbed λ amplitude equations
  /// @param order order of perturbation
  /// @param rank rank of perturbation operator. r = 1 means one-body perturbation operator
  /// @pre `rank==1 && order==1`, only first order perturbation and one-body perturbation operator is supported now
  /// @return std::vector of perturbed λ amplitude equations
  // clang-format on
  [[nodiscard]] std::vector<sequant::ExprPtr> λ_pt(size_t order = 1,
                                                   size_t rank = 1);

  // clang-format off
  /// @brief derives right-side sigma equations for EOM-CC
  /// @param K_occ number of operators in the occupied space in R operator
  /// @param K_uocc number of operators in the unoccupied space in R operator
  /// @return vector of right side sigma equations, with projector corresponding to \p K_occ and \p K_uocc; element 0 is always null
  // clang-format on
  [[nodiscard]] std::vector<sequant::ExprPtr> R(size_t K_occ, size_t K_uocc);

  // clang-format off
  /// @brief derives left-side sigma equations for EOM-CC
  /// @param K_occ number of operators in the occupied space in L operator
  /// @param K_uocc number of operators in the unoccupied space in L operator
  /// @return vector of left side sigma equations, with projector corresponding to \p K_occ and \p K_uocc; element 0 is always null
  // clang-format on
  [[nodiscard]] std::vector<sequant::ExprPtr> L(size_t K_occ, size_t K_uocc);

 private:
  size_t N;
  Ansatz ansatz_ = Ansatz::T;
};  // class CC

}  // namespace sequant::mbpt::sr

#endif  // SEQUANT_DOMAIN_MBPT_MODELS_CC_HPP
