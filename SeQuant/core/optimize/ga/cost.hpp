#ifndef SEQUANT_CORE_OPTIMIZE_GA_COST_HPP
#define SEQUANT_CORE_OPTIMIZE_GA_COST_HPP

// Mask-based fast path for the per-merge flop count, reproducing
// opt::detail::flops_counter (with empty inner_pow) exactly: the extent
// product over the deduplicated union of the indices of lhs, rhs and result,
// including the proto-indices of any composite. Verified against
// flops_counter in the unit tests; keep the two in lockstep.

#include <SeQuant/core/optimize/ga/key_table.hpp>

namespace sequant::opt::ga {

/// Extent product over the union of the open indices of subsets \p a, \p b and
/// \p a|b of \p T. == flops_counter(F(a), F(b), F(a|b)) except for the
/// scalar-contraction rule, which the caller's CostModel applies.
inline double merge_volume(TermTable const& T, NodeMask a, NodeMask b) {
  const NodeMask ab = a | b;
  IndexMask m = T.outer_mask[a] | T.outer_mask[b] | T.outer_mask[ab] |
                T.inner_mask[a] | T.inner_mask[b] | T.inner_mask[ab];
  double v = 1.;
  for (; m; m &= m - 1) v *= T.extent[std::countr_zero(m)];
  return v;
}

/// Cost conventions of the two accounting modes. Native (default) is
/// SeQuant's flops_counter convention: a merge costs its extent-product,
/// scalar contractions are free, an L2 addition costs the face size.
/// Python-parity is the prototype's -- 2x the extent-product, no
/// scalar-contraction rule -- and exists only to validate against the
/// prototype's reference numbers. Both are pure FLOP models; nothing here
/// prices storage, and \ref volatile_weight is a REPLAY COUNT on the flops
/// already counted, not a second metric.
struct CostModel {
  double merge_factor = 1.0;
  bool zero_scalar_merges = true;
  /// Replay weight on volatile (amplitude-dependent) work: only that part of
  /// the DAG is rebuilt per solver iteration, so the honest objective is
  /// `persistent_flops + N * volatile_flops`. Mirrors
  /// CostParams::volatile_weight; 1.0 is the volatility-blind cost.
  double volatile_weight = 1.0;

  static CostModel native() { return {1.0, true}; }
  static CostModel python_parity() { return {2.0, false}; }

  /// Whether cluster \p S has to be rebuilt on every replay: true iff it
  /// contains an amplitude leaf. Monotone in S, so no propagation pass is
  /// needed for L1 clusters.
  static bool is_volatile(TermTable const& T, NodeMask S) {
    return (S & T.volatile_mask) != 0;
  }
  /// How many times work of this volatility is paid over the whole solve.
  double replays(bool volatile_work) const {
    return volatile_work ? volatile_weight : 1.0;
  }

  /// \p volatile_work: whether the merged value depends on the amplitudes.
  double merge(TermTable const& T, NodeMask a, NodeMask b,
               bool volatile_work) const {
    const double v = merge_volume(T, a, b);
    if (zero_scalar_merges && v == 1.) return 0.;
    return replays(volatile_work) * merge_factor * v;
  }
  /// L1 form: the merge produces cluster `a|b`, whose volatility is decided by
  /// its own tensors.
  double merge(TermTable const& T, NodeMask a, NodeMask b) const {
    return merge(T, a, b, is_volatile(T, a | b));
  }
  double addition(TermTable const& T, NodeMask S, bool volatile_work) const {
    return replays(volatile_work) * T.face_size(S);
  }
  double addition(TermTable const& T, NodeMask S) const {
    return addition(T, S, is_volatile(T, S));
  }
};

}  // namespace sequant::opt::ga

#endif  // SEQUANT_CORE_OPTIMIZE_GA_COST_HPP
