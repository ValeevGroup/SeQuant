#ifndef SEQUANT_DOMAIN_MBPT_MODELS_CC_HPP
#define SEQUANT_DOMAIN_MBPT_MODELS_CC_HPP

#include <SeQuant/core/op.hpp>
#include <SeQuant/core/utility/aggregate.hpp>
#include <SeQuant/domain/mbpt/op.hpp>
#include <SeQuant/domain/mbpt/utils.hpp>
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
    SEQUANT_DESIGNATED_INIT_ONLY;
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

  /// @brief returns a copy of this object whose H̄ commutator truncation rank
  /// is @p rank, leaving this engine unchanged; all other options are carried
  /// over
  [[nodiscard]] CC with_hbar_comm_rank(size_t rank) const;

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
  /// @note For a non-unitary ansatz each commutator is represented as a
  /// connected product, \f$ (\hat{A}\hat{B})_c \f$, equivalent to \f$
  /// [\hat{A},\hat{B}] \f$ only once operator connectivity is supplied
  /// downstream; this engine always supplies it, by handing `ref_av`
  /// `default_op_connections()` or a superset thereof. A unitary ansatz uses
  /// explicit commutators instead, and is handed empty connectivity.
  /// @warning Because of the above, the expression returned for a non-unitary
  /// ansatz is NOT self-contained: evaluating it with empty connectivity, e.g.
  /// `op::ref_av(P(nₚ(2)) * cc.hbar(), {.connect = {}})`, silently retains the
  /// disconnected terms that the commutator would have cancelled. Pass
  /// `default_op_connections()` (which `op::ref_av`/`op::vac_av` use when the
  /// options argument is omitted -- note that passing `{}` instead gives empty
  /// connections), or build H̄ yourself with `mbpt::lst(..., {})` if you need
  /// the explicit form. For a unitary ansatz the reverse holds: H̄ is already
  /// self-contained, so connectivity must be left empty. See the "Using H̄
  /// outside the CC class" section of the user guide.
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
  /// @return vector of right side sigma equations, element 0 is always null
  [[nodiscard]] std::vector<ExprPtr> eom_r(nₚ np, nₕ nh) const;

  /// @brief derives left-side sigma equations for EOM-CC
  /// @param np number of particle annihilators in L operator
  /// @param nh number of hole annihilators in L operator
  /// @return vector of left side sigma equations, element 0 is always null
  [[nodiscard]] std::vector<ExprPtr> eom_l(nₚ np, nₕ nh) const;

  /// @brief derives the @p rank -particle reduced density matrix (RDM) as the
  /// reference expectation value of a similarity-transformed number operator.
  /// Traditional ansatz:
  /// \f$ \gamma^{p_1 \dots p_r}_{q_1 \dots q_r} = \langle 0| (1 +
  /// \hat{\Lambda}) e^{-\hat{T}} \{ a^{p_1 \dots p_r}_{q_1 \dots q_r} \}
  /// e^{\hat{T}} |0 \rangle \f$; unitary ansatz drops \f$ \hat{\Lambda} \f$ and
  /// uses \f$ e^{-\hat{\sigma}} \dots e^{\hat{\sigma}} \f$ with \f$
  /// \hat{\sigma} = \hat{T} - \hat{T}^\dagger \f$. The number-operator indices
  /// are in the complete space, so they survive as the free indices of the
  /// result. Returns the *linked* density (the number
  /// operator is required to connect to the cluster amplitudes).
  /// @param rank particle rank of the RDM (1 = one-particle \f$ \gamma \f$,
  ///   2 = two-particle \f$ \Gamma \f$, ...)
  /// @param comm_rank maximum order of nested commutators in \f$ \bar{N} \f$;
  ///   if not specified, defaults to `min(2*rank, rank + N)` for the
  ///   traditional ansatz (where the expansion terminates exactly) and to
  ///   `hbar_comm_rank` for the unitary ansatz (where it does not). Pass an
  ///   explicit value to truncate earlier.
  /// @return the RDM expression (Fermi-vacuum normal-ordered / correlation
  /// part)
  [[nodiscard]] ExprPtr rdm(
      size_t rank = 1, std::optional<size_t> comm_rank = std::nullopt) const;

 private:
  size_t N;
  Ansatz ansatz_ = Ansatz::T;
  bool skip_singles_ = false;
  bool screen_ = true;
  bool use_topology_ = true;
  std::optional<size_t> hbar_comm_rank_ = std::nullopt;
  std::optional<size_t> pertbar_comm_rank_ = std::nullopt;

  /// @return the `LSTOptions` this engine uses for every `mbpt::lst()` call
  /// @note The choice of commutator representation is really a question of
  /// whether the caller supplies operator connectivity downstream; for this
  /// engine the two coincide. Every non-unitary path hands `ref_av` a
  /// connectivity map (`default_op_connections()`, or the λ⁺/perturbed
  /// supersets thereof), which is what makes the cheaper connected-product
  /// form equivalent to the explicit commutator. The unitary paths hand
  /// `ref_av` empty connectivity and so must use explicit commutators.
  [[nodiscard]] LSTOptions lst_options() const {
    return {.unitary = unitary(), .use_connected_form = !unitary()};
  }

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
