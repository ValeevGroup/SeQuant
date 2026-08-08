#ifndef SEQUANT_CORE_OPTIMIZE_GA_COST_HPP
#define SEQUANT_CORE_OPTIMIZE_GA_COST_HPP

// Mask-based fast path for the per-merge flop count, reproducing
// opt::detail::flops_counter (with empty inner_pow, i.e. flat composite
// extents) exactly: the extent product over the deduplicated union of the
// indices of lhs, rhs, and result -- including, via outer_mask, the
// proto-indices of any composite present. Verified against flops_counter in
// the unit tests; keep the two in lockstep.

#include <SeQuant/core/optimize/ga/key_table.hpp>

namespace sequant::opt::ga {

/// Extent product over the union of the open indices of subsets \p a, \p b,
/// and \p a|b of term \p T. Equals flops_counter(...)(F(a), F(b), F(a|b))
/// except for the scalar-contraction rule, applied by the caller's CostModel.
inline double merge_volume(TermTable const& T, NodeMask a, NodeMask b) {
  const NodeMask ab = a | b;
  IndexMask m = T.outer_mask[a] | T.outer_mask[b] | T.outer_mask[ab] |
                T.inner_mask[a] | T.inner_mask[b] | T.inner_mask[ab];
  double v = 1.;
  for (; m; m &= m - 1) v *= T.extent[std::countr_zero(m)];
  return v;
}

/// Cost conventions of the two supported accounting modes.
///   - Native (default): SeQuant's flops_counter convention -- a merge costs
///     its extent-product (1 flop per output element and contraction step),
///     scalar contractions are free; an L2 addition costs the face size.
///   - Python-parity: the prototype's convention -- a merge costs 2x the
///     extent-product (multiply + add), no scalar-contraction rule; an L2
///     addition costs the face size. Used only to validate the port against
///     the prototype's reference numbers.
///
/// Both are pure FLOP models: nothing here prices storage. What \ref
/// volatile_weight adds is not a second metric but a REPLAY COUNT on the flops
/// already counted -- the objective stays a flop-optimal DAG search, it just
/// stops pretending every merge is paid the same number of times.
struct CostModel {
  double merge_factor = 1.0;
  bool zero_scalar_merges = true;
  /// Replay weight on volatile (amplitude-dependent) work -- conceptually the
  /// expected number of times the schedule is replayed; mirrors
  /// CostParams::volatile_weight, which single-term optimization already
  /// applies. In an iterative solver the DAG is evaluated once per iteration,
  /// but only the part depending on the amplitudes is REBUILT each time: an
  /// amplitude-free array is produced once and reused, so the honest objective
  /// is `persistent_flops + N * volatile_flops`. Weighting each key by N iff
  /// its cluster contains a volatile leaf yields exactly that, because
  /// Fitness::resolve already visits every key once. 1.0 (the default) is the
  /// neutral, volatility-blind cost this optimizer used before, so an empty
  /// volatile-leaf predicate leaves every number unchanged.
  double volatile_weight = 1.0;

  static CostModel native() { return {1.0, true}; }
  static CostModel python_parity() { return {2.0, false}; }

  /// Whether cluster \p S has to be rebuilt on every replay: true iff it
  /// contains an amplitude leaf. Monotone in S -- any superset of a volatile
  /// cluster is volatile -- so no propagation pass is needed for L1 clusters.
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
