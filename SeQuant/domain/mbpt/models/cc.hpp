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

  enum class HbarExpansion {
    /// standard Baker-Campbell-Hausdorff commutator expansion
    BCH,
    /// Bernoulli expansion, 10.1063/1.5030344 (U ansatz only)
    Bernoulli
  };

  /// Assembly used for the right-hand UCC EOM equations.
  enum class UCCEOMAssembly {
    /// Project \f$ [\bar{H}, R] \f$.
    Commutator,
    /// Assemble the projected \f$ \bar{H} \f$ matrix.
    ProjectedHbar
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
    /// choice of H̄ expansion; Bernoulli requires the U ansatz, not merely a
    /// unitary one, and ignores `screen` and `use_topology`, see `hbar()`
    HbarExpansion hbar_expansion = HbarExpansion::BCH;
  };

  /// @brief constructs CC engine with default options (traditional ansatz,
  /// screening enabled, topological optimization enabled)
  /// @param n coupled cluster excitation rank
  explicit CC(size_t n);

  /// @brief constructs CC engine with custom options
  /// @param n coupled cluster excitation rank
  /// @param opts configuration options @see CC::Options
  /// @throw Exception if a unitary ansatz has no `hbar_comm_rank`, or if the
  /// Bernoulli expansion is requested with an ansatz other than `Ansatz::U`
  explicit CC(size_t n, const Options& opts);

  /// @return the type of ansatz
  [[nodiscard]] Ansatz ansatz() const;

  /// @return true if the ansatz is unitary
  [[nodiscard]] bool unitary() const;

  /// @return the maximum of nested commutators in H̄; returns std::nullopt if
  /// not set
  [[nodiscard]] std::optional<size_t> hbar_comm_rank() const;

  /// @return the choice of H̄ expansion
  [[nodiscard]] HbarExpansion hbar_expansion() const;

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
  /// @note The returned expression depends on the ansatz and expansion:
  ///   - A non-unitary ansatz represents each commutator as a connected
  ///     product,
  ///     \f$ (\hat{A}\hat{B})_c \f$. It is equivalent to \f$
  ///     [\hat{A},\hat{B}] \f$ only when operator connectivity is supplied
  ///     downstream.
  ///   - A unitary BCH ansatz uses explicit commutators and empty connectivity.
  ///   - `HbarExpansion::Bernoulli` returns a tensor-level expression for use
  ///     with `op::tensor` projectors and `op::tensor::ref_av`. It never
  ///     reaches `CC::ref_av`, so `screen` and `use_topology` have no effect.
  /// @warning A non-unitary H̄ is not self-contained. Evaluating it with empty
  ///   connectivity, e.g. `op::ref_av(P(nₚ(2)) * cc.hbar(), {.connect = {}})`,
  ///   retains disconnected terms. Pass `default_op_connections()` (the default
  ///   when the options argument is omitted), or build H̄ with
  ///   `mbpt::lst(..., {})` for an explicit form. A unitary H̄ is
  ///   self-contained, so its connectivity must be empty. See the "Using H̄
  ///   outside the CC class" section of the user guide.
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

  // clang-format off
  /// @brief derives right-side sigma equations for EOM-CC
  /// @param np number of particle creators in R operator
  /// @param nh number of hole creators in R operator
  /// @param block_ranks optional per-block H̄ truncation ranks:
  ///   - A `K`×`K` matrix is read row by row in ascending manifold rank. For
  ///     singles+doubles it is
  ///     | H_SS  H_SD |    e.g.  | 2  1 |
  ///     | H_DS  H_DD |          | 1  0 |
  ///     Thus `{2,1,1,0}` uses [[H,σ],σ] for H_SS, [H,σ] for H_SD and
  ///     H_DS, and bare Hamiltonian integrals for H_DD.
  ///   - The same ordering serves EE, IP, and EA: read S as 1h/1p and D as
  ///     2h1p/1h2p (10.1021/acs.jctc.5c01991, Fig. 1).
  ///   - Under `Bernoulli`, each rank is a Bernoulli order H̄^k whose
  ///     commutators are in V alone. These are the Bernoulli qUCCSD ranks
  ///     UCCSD[2|2,1,0] of 10.1063/5.0062090, Eqs. (29), (41), (44), (48). Its
  ///     BCH scheme needs F one rank above V, which one rank per block cannot
  ///     express.
  /// @param assembly optional explicit UCC assembly:
  ///   - If omitted, `Bernoulli` and non-empty `block_ranks` select
  ///     `ProjectedHbar`; otherwise `BCH` selects `Commutator`.
  ///   - The two BCH assemblies differ by ground-state amplitude-equation
  ///     residuals and agree when those vanish.
  ///   - Traditional CC accepts only `Commutator`, uses its connected H̄R
  ///     product, and does not support block ranks.
  /// @throw Exception if `Commutator` is used with the Bernoulli expansion, if
  ///   `ProjectedHbar` or non-empty `block_ranks` is used with a traditional
  ///   ansatz, or if `block_ranks` is not `K`×`K`
  /// @return projected sigma equations in a vector of size `min(np, nh) + 1`,
  ///   indexed by the smaller particle/hole count of each projection manifold:
  ///   - Element 0 is null iff `np == nh`.
  ///   - With projected-H̄ UCC assembly, each diagonal subtracts the scalar
  ///     carried by the same temporary \f$ \bar{H}^{(k)} \f$ used for that block,
  ///     leaving the per-block normal-ordered components of Eq. (10), not
  ///     distinct physical energy zeros.
  // clang-format on
  [[nodiscard]] std::vector<ExprPtr> eom_r(
      nₚ np, nₕ nh, const std::vector<std::size_t>& block_ranks = {},
      std::optional<UCCEOMAssembly> assembly = std::nullopt) const;

  /// @brief derives left-side sigma equations for EOM-CC
  /// @param np number of particle annihilators in L operator
  /// @param nh number of hole annihilators in L operator
  /// @return vector of left side sigma equations, element 0 is always null
  [[nodiscard]] std::vector<ExprPtr> eom_l(nₚ np, nₕ nh) const;

  /// @brief derives the reduced density matrix (RDM) of particle rank @p rank
  /// as the reference expectation value of a similarity-transformed
  /// replacement operator, whose free-index ordinals are fixed by op::ã.
  /// Traditional ansatz:
  /// \f$ \gamma^{p_1 \dots p_r}_{p_{r+1} \dots p_{2r}} = \langle 0| (1 +
  /// \hat{\Lambda}) e^{-\hat{T}} \{ \tilde{a}^{p_1 \dots p_r}_{p_{r+1} \dots
  /// p_{2r}} \} e^{\hat{T}} |0 \rangle \f$; unitary ansatz drops \f$
  /// \hat{\Lambda} \f$ and uses \f$ e^{-\hat{\sigma}} \dots e^{\hat{\sigma}}
  /// \f$ with \f$ \hat{\sigma} = \hat{T} - \hat{T}^\dagger \f$.
  /// @note Only the correlation contribution is returned: op::ã is
  ///   normal-ordered, so the reference contractions are absent.
  /// @note For @p rank >= 2 the result is not manifestly antisymmetric: the
  ///   replacement operator carries no antisymmetrizer.
  /// @note For the traditional ansatz this is the *linked* density; the
  ///   unitary ansatz uses no connectivity.
  /// @param rank particle rank of the RDM (1 = one-particle \f$ \gamma \f$,
  ///   2 = two-particle \f$ \Gamma \f$, ...)
  /// @param comm_rank maximum order of nested commutators in the
  ///   similarity transform of the replacement operator;
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
  HbarExpansion hbar_expansion_ = HbarExpansion::BCH;

  /// @brief assembles the right-hand UCC EOM equations
  /// @param block_ranks see `eom_r`; empty uses the configured H̄ rank, or 4
  /// @param assembly assembly to use
  /// @pre a unitary ansatz
  /// @note under the Bernoulli expansion each block's H̄ has its N part removed.
  ///   Where a block's rank equals `hbar_comm_rank` the removed terms vanish at
  ///   converged amplitudes, so its equations change but its values do not;
  ///   below that rank the values change too. `BCH` keeps them.
  [[nodiscard]] std::vector<ExprPtr> assemble_ucc_eom(
      nₚ np, nₕ nh, const std::vector<std::size_t>& block_ranks,
      UCCEOMAssembly assembly) const;

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
