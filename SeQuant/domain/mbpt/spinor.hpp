//
// Kramers tracing for closed-shell relativistic (2-/4-component) theories.
//
// Rewrites an all-spinor (Kramers-free) expression — e.g. the MP2 energy
// `1/4 g-bar t-bar` over full spinor indices — into a sum over the
// time-reversal-symmetry (TRS) canonical Kramers-block representatives.
//
// Design (round 2): a single function `closed_shell_kramers_trace` that
// mirrors the naive spin tracer `spintrace_impl` but, instead of the
// Ms-conservation filter (`can_expand`, which is invalid relativistically),
// folds whole Kramers configurations related by global time reversal
// (T = flip every Kramers label + complex-conjugate) into a single
// `2 Re[...]` representative. Particle-interchange (sigma) folding is left
// to `sequant::canonicalize`; external-index antisymmetry folding (CC,
// rank-general) and the per-column internal-T-reach are deferred.
//
// Kramers labels reuse the `Spin` quantum number (space_qns.hpp):
//   Spin::alpha (== Spin::up)   = Kramers-up   (label suffix U+2191)
//   Spin::beta  (== Spin::down) = Kramers-down (label suffix U+2193)
//
// `RealPart`/`ImagPart` are symbolic scalar Expr markers: they wrap a
// closed (scalar-valued) contraction so the evaluator takes its real/
// imaginary part. They are opaque to `simplify`/`canonicalize`, so the
// wrapped sub-expression is never re-permuted. SeQuant's own evaluation
// engine does not (yet) evaluate `Re()`/`Im()` of a tensor network — that
// is a TODO; for now the marker is consumed by the MPQC evaluator.
//

#ifndef SEQUANT_DOMAIN_MBPT_SPINOR_HPP
#define SEQUANT_DOMAIN_MBPT_SPINOR_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/hash.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <SeQuant/domain/mbpt/fwd.hpp>

#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>

namespace sequant::mbpt {

// =====================================================================
// Re/Im evaluation markers + complex-conjugation label convention.
//
// Convention: complex conjugation of a tensor is encoded as a `*` suffix
// on the tensor label (e.g. `g` -> `g*`). Evaluator-side dispatch
// translates the suffix into a `.conj()` call on the numeric tensor.
// =====================================================================

/// @brief Tests whether a tensor label encodes complex conjugation.
inline bool has_conj_suffix(std::wstring_view label) {
  return !label.empty() && label.back() == L'*';
}

/// @brief Symbolic real-part wrapper: `RealPart(E)` = `Re(E)` for a
///        scalar-valued Expr `E`.
class RealPart : public sequant::Expr {
 public:
  RealPart() = delete;
  RealPart(const RealPart&) = default;
  RealPart(RealPart&&) = default;
  ~RealPart() override = default;

  /// @pre @p inner is non-null and scalar-valued
  ///
  /// We do not enforce `Expr::is_scalar()` because Tensor (an atom)
  /// returns false unconditionally and a closed-contraction Product of
  /// two Tensors inherits the same answer — the typical inner here.
  explicit RealPart(ExprPtr inner) : inner_{std::move(inner)} {
    SEQUANT_ASSERT(inner_);
  }

  const ExprPtr& inner() const { return inner_; }
  bool is_scalar() const override { return true; }
  type_id_type type_id() const override { return get_type_id<RealPart>(); }
  ExprPtr clone() const override { return ex<RealPart>(inner_->clone()); }
  void adjoint() override {}  // Re(E) is real, self-adjoint
  std::wstring to_latex() const override {
    return L"\\Re\\left[" + inner_->to_latex() + L"\\right]";
  }

 private:
  ExprPtr inner_;

  hash_type memoizing_hash() const override {
    auto compute = [this]() {
      auto v = hash::value(*inner_);
      hash::combine(v, std::size_t{0xC0FFEE01ull});
      return v;
    };
    if (!hash_value_) hash_value_ = compute();
    return *hash_value_;
  }
  bool static_equal(const Expr& that) const override {
    return *inner_ == *static_cast<const RealPart&>(that).inner_;
  }
  bool static_less_than(const Expr& that) const override {
    return *inner_ < *static_cast<const RealPart&>(that).inner_;
  }
};

/// @brief Symbolic imaginary-part wrapper: `ImagPart(E)` = `Im(E)`.
class ImagPart : public sequant::Expr {
 public:
  ImagPart() = delete;
  ImagPart(const ImagPart&) = default;
  ImagPart(ImagPart&&) = default;
  ~ImagPart() override = default;

  /// @pre @p inner is non-null and scalar-valued (see RealPart for why
  ///      we don't use `Expr::is_scalar()` as the precondition)
  explicit ImagPart(ExprPtr inner) : inner_{std::move(inner)} {
    SEQUANT_ASSERT(inner_);
  }

  const ExprPtr& inner() const { return inner_; }
  bool is_scalar() const override { return true; }
  type_id_type type_id() const override { return get_type_id<ImagPart>(); }
  ExprPtr clone() const override { return ex<ImagPart>(inner_->clone()); }
  void adjoint() override {}  // Im(E) is real, self-adjoint
  std::wstring to_latex() const override {
    return L"\\Im\\left[" + inner_->to_latex() + L"\\right]";
  }

 private:
  ExprPtr inner_;

  hash_type memoizing_hash() const override {
    auto compute = [this]() {
      auto v = hash::value(*inner_);
      hash::combine(v, std::size_t{0xC0FFEE02ull});
      return v;
    };
    if (!hash_value_) hash_value_ = compute();
    return *hash_value_;
  }
  bool static_equal(const Expr& that) const override {
    return *inner_ == *static_cast<const ImagPart&>(that).inner_;
  }
  bool static_less_than(const Expr& that) const override {
    return *inner_ < *static_cast<const ImagPart&>(that).inner_;
  }
};

/// @brief Wraps @p expr as `Re(expr)`.
[[nodiscard]] inline ExprPtr real_part(ExprPtr expr) {
  return ex<RealPart>(std::move(expr));
}

/// @brief Wraps @p expr as `Im(expr)`.
[[nodiscard]] inline ExprPtr imaginary_part(ExprPtr expr) {
  return ex<ImagPart>(std::move(expr));
}

// =====================================================================
// Kramers tracer.
// =====================================================================

// clang-format off
/// @brief Traces an all-spinor (Kramers-free) closed-shell expression into a
///        sum over time-reversal-canonical Kramers-block representatives.
/// @details Mirrors the naive spin tracer `spintrace` in structure: it
/// enumerates the 2^n Kramers configurations of the term's indices and assigns
/// Kramers-up/down labels (reusing `Spin::alpha`/`Spin::beta`). Unlike
/// `spintrace`, it applies NO Ms-conservation filter (Kramers is not conserved
/// relativistically). Configurations related by global time reversal
/// (T: flip every Kramers label, complex-conjugate) are folded pairwise into a
/// single `RealPart`-wrapped representative (`A + A* = 2 Re A`). For the
/// integral `g` the bra antisymmetry is expanded (to the NonSymm form);
/// the amplitude `t` is kept antisymmetric (the "level-1" form). The caller is
/// expected to run `canonicalize` + `rapid_simplify` on the result to perform
/// the particle-interchange (sigma) merge.
/// @param expr an all-spinor expression (indices have no Kramers/Spin label)
/// @param ext_index_groups groups of external indices (empty for a fully
///        contracted energy; external-antisymmetry folding is not yet
///        implemented and these are treated like internal groups for now)
/// @param fold_T if true (default), fold each global-T (conjugate) pair into a
///        single `2 Re[...]` representative; if false, emit every configuration
///        (sigma-merged) verbatim with no `RealPart` wrapper — the sum is
///        complex and the caller takes its real part. The false form is for
///        evaluators that cannot evaluate `Re()` of a tensor network (e.g. the
///        CCk energy observable, which sums the blocks then takes the real part).
/// @param expand_g if true (default), expand the integral `g`'s antisymmetry
///        (the level-1 raw-g form); if false, keep `g` antisymmetric (ḡ) so the
///        evaluator fetches the factory [as] block and the cross-Kramers
///        antisymmetry is handled inside the integral. (g-expansion is a later
///        optimization stage.)
/// @return the Kramers-traced expression (unsimplified)
// clang-format on
ExprPtr closed_shell_kramers_trace(
    const ExprPtr& expr,
    const container::svector<container::svector<Index>>& ext_index_groups = {},
    bool fold_T = true, bool expand_g = true);

// clang-format off
/// @brief Orbits of the n-bit Kramers configurations under a group of bit
/// permutations and (optionally) global time reversal.
/// @details A configuration is an n-bit integer (bit p = Kramers label of index
/// slot p: 0 = up, 1 = down). The orbit group is generated by @p bit_perms
/// (each a permutation `q` of [0,n) acting as: bit p -> bit q[p]) together with,
/// if @p use_T, the global time-reversal T = full one's-complement. This is the
/// building block of the Kramers folds: the external-antisymmetry fold uses the
/// external transpositions P_ab / P_ij with T (doubles: 16 -> 5 blocks); the
/// internal/MP2 fold uses the particle interchange sigma with T. Mirrors
/// `~/code/KR-MP2/kramers_orbit_compress.py`.
/// @param n number of index slots (bits); must be <= 62
/// @param bit_perms generator permutations (each of length n)
/// @param use_T include global time reversal (complement) as a generator
/// @return one ascending-sorted vector of member configs per orbit; the
///         canonical representative is the front element (orbit-min)
// clang-format on
container::svector<container::svector<std::uint64_t>> kramers_config_orbits(
    std::size_t n,
    const container::svector<container::svector<std::size_t>>& bit_perms,
    bool use_T);

/// One member of a Kramers external block: its configuration plus the transform
/// that reconstructs its (residual) tensor from the orbit's canonical
/// representative,  block(config) = sign * [conj] * permute(block(canonical)).
struct KramersBlockMember {
  std::uint64_t config;  //!< this member's configuration
  int sign;              //!< +1 or -1
  bool conj;             //!< complex-conjugate the canonical block
  container::svector<std::size_t>
      perm;  //!< index permutation (slot p -> perm[p])
};

/// A symmetry-unique external block: the canonical (orbit-min) configuration
/// and every member with its reconstruction transform.
struct KramersBlock {
  std::uint64_t canonical;
  container::svector<KramersBlockMember> members;  //!< includes the canonical
};

// clang-format off
/// @brief External Kramers blocks with per-member reconstruction transforms.
/// @details Like kramers_config_orbits, but each generator carries a transform
/// so that every non-canonical member records how its tensor is obtained from
/// the canonical's: the transposition generators @p antisym_perms each act with
/// sign -1 (residual antisymmetry), the @p symm_perms each act with sign +1
/// (e.g. particle interchange sigma on a raw, non-antisymmetrized integral),
/// and (if @p use_T) global time reversal T acts as complex conjugation with
/// sign (-1)^(#down) and the identity slot permutation. This is the eval-time
/// external reconstruction (compute the canonical block, fill the rest by
/// sign/conj/perm): use {antisym_perms = external transpositions, T} for the
/// antisymmetric residual (doubles -> 5 blocks), or {symm_perms = sigma, T} for
/// a raw g leaf (doubles -> 6 blocks).
/// @param n number of external slots (bits), <= 62
/// @param antisym_perms transposition generators acting with sign -1
/// @param use_T include global time reversal
/// @param symm_perms permutation generators acting with sign +1 (default none)
/// @return one KramersBlock per symmetry-unique external representative
// clang-format on
container::svector<KramersBlock> kramers_external_blocks(
    std::size_t n,
    const container::svector<container::svector<std::size_t>>& antisym_perms,
    bool use_T,
    const container::svector<container::svector<std::size_t>>& symm_perms = {});

// clang-format off
/// @brief CC residual Kramers trace — external fold (stages 1-2).
/// @details Given a CC residual expression carrying a leading antisymmetrizer Â
/// (its bra = external virtuals, ket = external occupieds), folds the external
/// Kramers configurations under the residual's external antisymmetry (the S_k
/// adjacent transpositions of each external group, rank-general) together with
/// global time reversal T (kramers_config_orbits), and returns one block per
/// symmetry-unique external representative with the external indices
/// Kramers-labeled. Â is kept (not expanded) — this is stage 2 of the pipeline
/// (trace+fold under Â). Doubles -> 5 blocks.
/// @note Internal tracing/folding, A-expand, g-expansion and the further g TRS
///       folds (stages 3-5) are not yet implemented (TODO).
/// @param expr a CC residual term/expression with a leading Â
/// @return one Kramers-labeled block per external representative
// clang-format on
container::svector<ExprPtr> closed_shell_kramers_CC_trace(const ExprPtr& expr);

/// @brief True if @p expr contains an antisymmetrizer (Â) tensor anywhere.
/// @details Lets a caller (e.g. the CC spintrace dispatch) tell a residual
/// equation (carries a leading Â -> closed_shell_kramers_CC_trace) from a fully
/// contracted scalar like the energy (no Â -> closed_shell_kramers_trace).
bool has_antisymmetrizer(const ExprPtr& expr);

}  // namespace sequant::mbpt

#endif  // SEQUANT_DOMAIN_MBPT_SPINOR_HPP
