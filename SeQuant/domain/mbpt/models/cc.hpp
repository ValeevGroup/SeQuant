#ifndef SEQUANT_DOMAIN_MBPT_MODELS_CC_HPP
#define SEQUANT_DOMAIN_MBPT_MODELS_CC_HPP

#include <SeQuant/core/op.hpp>
#include <SeQuant/domain/mbpt/op.hpp>

#include <cstddef>
#include <limits>
#include <vector>

namespace sequant {
class ExprPtr;
}

namespace sequant::mbpt {

/// CC is a derivation engine for the coupled-cluster method
class CC {
 public:
  enum class Ansatz {
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
  /// @param screen if true, uses Operator level screening before applying
  /// WickTheorem
  /// @param use_topology if true, uses topological optimizations in WickTheorem
  /// @param use_connectivity if true, uses connectivity information in
  /// WickTheorem
  explicit CC(size_t N, Ansatz ansatz = Ansatz::T, bool screen = true,
              bool use_topology = true, bool use_connectivity = true);

  /// @return the type of ansatz
  Ansatz ansatz() const;

  /// @return true if the ansatz is unitary
  [[nodiscard]] bool unitary() const;

  /// @return whether screening is on or not
  [[nodiscard]] bool screen() const;

  /// @return whether topological optimization is used in WickTheorem
  [[nodiscard]] bool use_topology() const;

  /// @return whether connectivity information is used in WickTheorem
  [[nodiscard]] bool use_connectivity() const;

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
  [[nodiscard]] std::vector<ExprPtr> t(
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
  [[nodiscard]] std::vector<ExprPtr> λ(size_t commutator_rank = 4);

  // clang-format off
  /// @brief derives perturbed t amplitude equations
  /// @param rank rank of perturbation operator. r = 1 means one-body perturbation operator
  /// @param order order of perturbation
  /// @pre `rank==1 && order==1`, only first order perturbation and one-body perturbation operator is supported now
  /// @return std::vector of perturbed t amplitude equations
  // clang-format on
  [[nodiscard]] std::vector<ExprPtr> t_pt(size_t rank = 1,
                                          [[maybe_unused]] size_t order = 1);

  // clang-format off
  /// @brief derives perturbed λ amplitude equations
  /// @param rank rank of perturbation operator. r = 1 means one-body perturbation operator
  /// @param order order of perturbation
  /// @pre `rank==1 && order==1`, only first order perturbation and one-body perturbation operator is supported now
  /// @return std::vector of perturbed λ amplitude equations
  // clang-format on
  [[nodiscard]] std::vector<ExprPtr> λ_pt(size_t rank = 1,
                                          [[maybe_unused]] size_t order = 1);

  /// @brief derives right-side sigma equations for EOM-CC
  /// @param np number of particle creators in R operator
  /// @param nh number of hole creators in R operator
  /// @return vector of right side sigma equations, element 0 is always null
  [[nodiscard]] std::vector<ExprPtr> eom_r(nₚ np, nₕ nh);

  /// @brief derives left-side sigma equations for EOM-CC
  /// @param np number of particle annihilators in L operator
  /// @param nh number of hole annihilators in L operator
  /// @return vector of left side sigma equations, element 0 is always null
  [[nodiscard]] std::vector<ExprPtr> eom_l(nₚ np, nₕ nh);

 private:
  size_t N;
  Ansatz ansatz_ = Ansatz::T;
  bool screen_ = true;
  bool use_topology_ = true;
  bool use_connectivity_ = true;

  /// @brief computes reference expectation value of an expression. Dispatches
  /// to `mbpt::op::ref_av()`
  /// @param[in] expr input expression
  /// @param[in] op_connect list of pairs of operators to be connected. Default
  /// is given by `mbpt::default_op_connections()`.
  auto ref_av(const ExprPtr& expr,
              const OpConnections<mbpt::OpType>& op_connect =
                  default_op_connections()) const {
    return op::ref_av(expr, op_connect, this->use_topology(), this->screen(),
                      this->use_connectivity());
  }
};  // class CC

}  // namespace sequant::mbpt

#endif  // SEQUANT_DOMAIN_MBPT_MODELS_CC_HPP
