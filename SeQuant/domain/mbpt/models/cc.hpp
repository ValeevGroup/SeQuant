#ifndef SEQUANT_DOMAIN_MBPT_MODELS_CC_HPP
#define SEQUANT_DOMAIN_MBPT_MODELS_CC_HPP

#include <SeQuant/core/op.hpp>
#include <SeQuant/domain/mbpt/op.hpp>
#include <SeQuant/domain/mbpt/vac_av.hpp>

#include <cstddef>
#include <limits>
#include <optional>
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

  /// Configuration options for CC class
  struct Options {
    /// type of CC ansatz. see CC::Ansatz
    Ansatz ansatz = Ansatz::T;
    /// if true, singles amplitudes are excluded from \f$ \hat{T} \f$ and \f$
    /// \hat{\Lambda} \f$; if not specified, defaults to true for
    /// orbital-optimized ansätze (oT, oU) and false otherwise. Must be true
    /// for orbital-optimized ansätze.
    std::optional<bool> skip_singles = std::nullopt;
    /// if true, uses Operator level screening before applying WickTheorem.
    /// This propagates to all ref_av() calls
    bool screen = true;
    /// if true, uses topological optimizations in WickTheorem
    bool use_topology = true;
    /// maximum order of nested commutators in H̄; must be specified if unitary
    /// ansatz is used
    std::optional<size_t> hbar_comm_rank = std::nullopt;
    /// maximum order of nested commutators in the similarity transformed
    /// perturbation operator; must be specified if unitary ansatz is used in
    /// perturbed amplitude derivation
    std::optional<size_t> pertbar_comm_rank = std::nullopt;
  };

  /// @brief constructs CC engine with default options (traditional ansatz,
  /// screening enabled, topological optimization enabled)
  /// @param n coupled cluster excitation rank
  explicit CC(size_t n);

  /// @brief constructs CC engine with custom options
  /// @param n coupled cluster excitation rank
  /// @param opts configuration options @see CC::Options
  explicit CC(size_t n, const Options& opts);

  /// @return the type of ansatz
  [[nodiscard]] Ansatz ansatz() const;

  /// @return true if the ansatz is unitary
  [[nodiscard]] bool unitary() const;

  /// @return the maximum of nested commutators in H̄; returns std::nullopt if
  /// not set
  [[nodiscard]] std::optional<size_t> hbar_comm_rank() const;

  /// @return true if singles amplitudes are excluded from \f$ \hat{T} \f$ and
  /// \f$ \hat{\Lambda} \f$
  [[nodiscard]] bool skip_singles() const;

  /// @return whether screening is on or not
  [[nodiscard]] bool screen() const;

  /// @return whether topological optimization is used in WickTheorem
  [[nodiscard]] bool use_topology() const;

  /// @brief computes similarity transformed Hamiltonian, \f$ \bar{H} =
  /// e^{-\hat{\sigma}} \hat{H} e^{\hat{\sigma}} \f$. The form of \f$ \sigma \f$
  /// depends on the Ansatz choice.
  /// @param truncation_rank maximum order of nested commutators to include in
  /// the expansion; if not specified, will use the value of member
  /// `hbar_comm_rank`. If that is also not specified, will use 4 as the default
  /// value. If provided, will override all defaults.
  [[nodiscard]] ExprPtr hbar(
      std::optional<size_t> truncation_rank = std::nullopt) const;

  /// @brief derives the CC energy expression \f$ \langle 0|\bar{H}|0 \rangle
  /// \f$ at the requested commutator truncation, WITHOUT deriving the
  /// projected amplitude equations (avoids deriving the full `t()` manifold
  /// just to read `t().at(0)`).
  /// @param comm_rank optional H̄ commutator-truncation override, forwarded to
  ///   @ref hbar (defaults to the engine's `hbar_comm_rank`, else 4).
  /// @return the energy expression \f$ \langle 0|\bar{H}|0 \rangle \f$
  [[nodiscard]] ExprPtr energy(
      std::optional<size_t> comm_rank = std::nullopt) const;

  /// @brief derives t amplitude equations, \f$ \langle P|\bar{H}|0 \rangle = 0
  /// \f$
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
      size_t pmax = std::numeric_limits<size_t>::max(), size_t pmin = 0) const;

  /// @brief derives λ amplitude equations,
  /// \f$ \langle 0| (1 + \hat{\Lambda}) \frac{d \bar{H}}{d \hat{T}_P} |0
  /// \rangle = 0 \f$
  /// @return vector of λ amplitude equations, with element `k` containing
  /// equation
  ///   \f$ \langle 0| (1 + \hat{\Lambda}) \frac{d \bar{H}}{d \hat{T}_k} |0
  ///   \rangle = 0 \f$ for `k` in
  /// the [1,N] range; element 0 contains the λ pseudoenergy, computed as the
  /// CC energy with \f$ \hat{T} \f$ replaced by \f$ \hat{\Lambda}^{\dagger} \f$
  [[nodiscard]] std::vector<ExprPtr> λ() const;

  // clang-format off
  /// @brief derives perturbed t amplitude equations
  /// @param rank rank of perturbation operator. r = 1 means one-body perturbation operator
  /// @param order order of perturbation
  /// @param nbatch optional batching index rank for perturbation operators
  /// @pre `rank==1 && order==1`, only first order perturbation and one-body perturbation operator is supported now
  /// @return std::vector of perturbed t amplitude equations
  // clang-format on
  [[nodiscard]] std::vector<ExprPtr> tʼ(
      size_t rank = 1, size_t order = 1,
      std::optional<size_t> nbatch = std::nullopt) const;

  // clang-format off
  /// @brief derives perturbed λ amplitude equations
  /// @param rank rank of perturbation operator. r = 1 means one-body perturbation operator
  /// @param order order of perturbation
  /// @param nbatch optional batching index rank for perturbation operators
  /// @pre `rank==1 && order==1`, only first order perturbation and one-body perturbation operator is supported now
  /// @return std::vector of perturbed λ amplitude equations
  // clang-format on
  [[nodiscard]] std::vector<ExprPtr> λʼ(
      size_t rank = 1, size_t order = 1,
      std::optional<size_t> nbatch = std::nullopt) const;

  /// @brief derives right-side sigma equations for EOM-CC
  /// @param np number of particle creators in R operator
  /// @param nh number of hole creators in R operator
  /// @param nbatch if set, R carries this many batching (`z`) indices so the
  /// sigma equations can be evaluated for a batch of trial vectors at once;
  /// requires the batching space to be registered (see add_batching_spaces)
  /// @return vector of right side sigma equations
  [[nodiscard]] std::vector<ExprPtr> eom_r(
      nₚ np, nₕ nh, std::optional<size_t> nbatch = std::nullopt) const;

  /// @brief derives left-side sigma equations for EOM-CC
  /// @param np number of particle annihilators in L operator
  /// @param nh number of hole annihilators in L operator
  /// @param nbatch if set, L carries this many batching (`z`) indices so the
  /// sigma equations can be evaluated for a batch of trial vectors at once;
  /// requires the batching space to be registered (see add_batching_spaces)
  /// @return vector of left side sigma equations
  [[nodiscard]] std::vector<ExprPtr> eom_l(
      nₚ np, nₕ nh, std::optional<size_t> nbatch = std::nullopt) const;

 private:
  size_t N;
  Ansatz ansatz_ = Ansatz::T;
  bool skip_singles_ = false;
  bool screen_ = true;
  bool use_topology_ = true;
  std::optional<size_t> hbar_comm_rank_ = std::nullopt;
  std::optional<size_t> pertbar_comm_rank_ = std::nullopt;

  /// @brief computes reference expectation value of an expression. Dispatches
  /// to `mbpt::op::ref_av()`
  /// @param[in] expr input expression
  /// @param[in] connect list of operator label pairs to connect.
  /// @param[in] do_not_connect list of operator label pairs to never connect.
  /// @note Uses use_topology() and screen() from the CC instance to set other
  /// EVOptions
  auto ref_av(
      const ExprPtr& expr,
      const OpConnections<std::wstring>& connect = default_op_connections(),
      const OpConnections<std::wstring>& do_not_connect = {}) const {
    return op::ref_av(expr, {.connect = connect,
                             .do_not_connect = do_not_connect,
                             .screen = this->screen(),
                             .use_topology = this->use_topology()});
  }
};  // class CC

}  // namespace sequant::mbpt

#endif  // SEQUANT_DOMAIN_MBPT_MODELS_CC_HPP
